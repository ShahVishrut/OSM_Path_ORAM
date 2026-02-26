#include "client.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

Client::Client(size_t num_blocks, size_t blocks_per_bucket) {
  if (num_blocks > UINT32_MAX) {
    throw std::invalid_argument("Too many blocks!");
  }
  this->num_blocks = num_blocks;
  this->tree_height = std::ceil(std::log2(this->num_blocks));
  this->blocks_per_bucket = blocks_per_bucket;
  this->bucket_size_bytes = this->block_size_bytes * this->blocks_per_bucket;

  this->server = new Server(this->bucket_size_bytes, (1u << (tree_height + 1)));

  stash.resize(3 * (tree_height + 1));
  stash_full.resize(stash.size());

  ORAMBlock dummy;
  std::vector<uint8_t> data(bucket_size_bytes);
  for (size_t j = 0; j < blocks_per_bucket; j++) {
    std::memcpy(&data[j * block_size_bytes], &dummy, block_size_bytes);
  }
  for (size_t i = 1; i < (1u << (tree_height + 1)); i++) {
    server->write_buckets({i}, data);
  }
}

size_t Client::size(uint64_t key) {
  if (root.is_null) {
    return 0;
  }

  std::vector<ORAMBlock> history;
  ORAMBlock cur_read;
  cur_read.header = root;
  size_t result = 0;

  while (true) {
    cur_read = write_block(cur_read, false);
    history.push_back(cur_read);

    if (key > cur_read.data.key) {
      if (cur_read.data.r_child_ptr.is_null) break;
      cur_read.header = cur_read.data.r_child_ptr;
    } else if (key < cur_read.data.key) {
      if (cur_read.data.l_child_ptr.is_null) break;
      cur_read.header = cur_read.data.l_child_ptr;
    } else {
      result = 1 + cur_read.data.l_same_key_size + cur_read.data.r_same_key_size;
      break;
    }
  }

  // PathORAM Write-back: Update parent pointers with new leaf_labels
  for (int index = (int)history.size() - 1; index >= 0; index--) {
    uint32_t new_leaf_label = write_block(history[index], true).header.leaf_label;
    if (index > 0) {
      if (!history[index - 1].data.l_child_ptr.is_null &&
          history[index - 1].data.l_child_ptr.block_id == history[index].header.block_id) {
        history[index - 1].data.l_child_ptr.leaf_label = new_leaf_label;
      } else if (!history[index - 1].data.r_child_ptr.is_null &&
                 history[index - 1].data.r_child_ptr.block_id == history[index].header.block_id) {
        history[index - 1].data.r_child_ptr.leaf_label = new_leaf_label;
      }
    } else {
      root.leaf_label = new_leaf_label;
    }
  }

  return result;
}

std::vector<uint64_t> Client::find(uint64_t key, uint32_t i, uint32_t j) {
  if (root.is_null) return {};

  std::vector<ORAMBlock> common_path;
  ORAMBlock cur_read;
  cur_read.header = root;

  uint32_t recursive_i = i;
  uint32_t recursive_j = j;

  bool found_i = false;
  bool found_j = false;
  bool abort_search = false;

  ODSPointer next_header_i;
  ODSPointer next_header_j;

  // 1. Traverse together until i and j paths diverge
  while (true) {
    cur_read = write_block(cur_read, false);
    common_path.push_back(cur_read);

    uint32_t next_recursive_i = recursive_i;
    uint32_t next_recursive_j = recursive_j;

    // Determine next step for i
    if (key > cur_read.data.key) {
      if (cur_read.data.r_child_ptr.is_null) { abort_search = true; break; }
      next_header_i = cur_read.data.r_child_ptr;
    } else if (key < cur_read.data.key) {
      if (cur_read.data.l_child_ptr.is_null) { abort_search = true; break; }
      next_header_i = cur_read.data.l_child_ptr;
    } else {
      if (cur_read.data.l_same_key_size == recursive_i) {
        found_i = true;
      } else if (cur_read.data.l_same_key_size > recursive_i) {
        if (cur_read.data.l_child_ptr.is_null) { abort_search = true; break; }
        next_header_i = cur_read.data.l_child_ptr;
      } else {
        if (cur_read.data.r_child_ptr.is_null) { abort_search = true; break; }
        next_header_i = cur_read.data.r_child_ptr;
        next_recursive_i = recursive_i - (cur_read.data.l_same_key_size + 1);
      }
    }

    // Determine next step for j
    if (key > cur_read.data.key) {
      if (cur_read.data.r_child_ptr.is_null) { abort_search = true; break; }
      next_header_j = cur_read.data.r_child_ptr;
    } else if (key < cur_read.data.key) {
      if (cur_read.data.l_child_ptr.is_null) { abort_search = true; break; }
      next_header_j = cur_read.data.l_child_ptr;
    } else {
      if (cur_read.data.l_same_key_size == recursive_j) {
        found_j = true;
      } else if (cur_read.data.l_same_key_size > recursive_j) {
        if (cur_read.data.l_child_ptr.is_null) { abort_search = true; break; }
        next_header_j = cur_read.data.l_child_ptr;
      } else {
        if (cur_read.data.r_child_ptr.is_null) { abort_search = true; break; }
        next_header_j = cur_read.data.r_child_ptr;
        next_recursive_j = recursive_j - (cur_read.data.l_same_key_size + 1);
      }
    }

    // Continue shared path if both need the exact same next block
    if (!found_i && !found_j && next_header_i.block_id == next_header_j.block_id) {
      cur_read.header = next_header_i;
      recursive_i = next_recursive_i;
      recursive_j = next_recursive_j;
    } else {
      break; 
    }
  }

  std::vector<uint64_t> result;
  std::vector<ORAMBlock> left_bfs;
  std::vector<ORAMBlock> right_bfs;
  std::vector<ORAMBlock> full_bfs;

  // Only run BFS logic if the nodes were actually found
  if (!abort_search) {
    std::vector<uint32_t> i_values;
    std::vector<uint32_t> j_values;

    ORAMBlock divergence_node = common_path.back();

    if (divergence_node.data.key == key) {
      if (divergence_node.data.l_same_key_size >= recursive_i && 
          divergence_node.data.l_same_key_size <= recursive_j) {
        result.push_back(divergence_node.data.value);
      }
      
      if (!divergence_node.data.l_child_ptr.is_null) {
        cur_read.header = divergence_node.data.l_child_ptr;
        left_bfs.push_back(write_block(cur_read, false));
        i_values.push_back(recursive_i);
      }
      
      if (!divergence_node.data.r_child_ptr.is_null && divergence_node.data.l_same_key_size < recursive_j) {
        cur_read.header = divergence_node.data.r_child_ptr;
        right_bfs.push_back(write_block(cur_read, false));
        j_values.push_back(recursive_j - (divergence_node.data.l_same_key_size + 1));
      }
    }

    // 1. Process Left Path BFS
    size_t left_idx = 0;
    while (left_idx < left_bfs.size()) {
      ORAMBlock node = left_bfs[left_idx];
      uint32_t target_i = i_values[left_idx];
      left_idx++;

      if (node.data.key != key) {
        if (!node.data.r_child_ptr.is_null) {
          cur_read.header = node.data.r_child_ptr;
          left_bfs.push_back(write_block(cur_read, false));
          i_values.push_back(target_i);
        }
        continue;
      }

      uint32_t L = node.data.l_same_key_size;
      if (L < target_i) {
        if (!node.data.r_child_ptr.is_null) {
          cur_read.header = node.data.r_child_ptr;
          left_bfs.push_back(write_block(cur_read, false));
          i_values.push_back(target_i - (L + 1));
        }
      } else {
        result.push_back(node.data.value);
        if (!node.data.l_child_ptr.is_null) {
          cur_read.header = node.data.l_child_ptr;
          left_bfs.push_back(write_block(cur_read, false));
          i_values.push_back(target_i);
        }
        if (!node.data.r_child_ptr.is_null) {
          cur_read.header = node.data.r_child_ptr;
          full_bfs.push_back(write_block(cur_read, false));
        }
      }
    }

    // 2. Process Right Path BFS
    size_t right_idx = 0;
    while (right_idx < right_bfs.size()) {
      ORAMBlock node = right_bfs[right_idx];
      uint32_t target_j = j_values[right_idx];
      right_idx++;

      if (node.data.key != key) {
        if (!node.data.l_child_ptr.is_null) {
          cur_read.header = node.data.l_child_ptr;
          right_bfs.push_back(write_block(cur_read, false));
          j_values.push_back(target_j);
        }
        continue;
      }

      uint32_t L = node.data.l_same_key_size;
      if (L > target_j) {
        if (!node.data.l_child_ptr.is_null) {
          cur_read.header = node.data.l_child_ptr;
          right_bfs.push_back(write_block(cur_read, false));
          j_values.push_back(target_j);
        }
      } else {
        result.push_back(node.data.value);
        if (!node.data.l_child_ptr.is_null) {
          cur_read.header = node.data.l_child_ptr;
          full_bfs.push_back(write_block(cur_read, false));
        }
        if (L < target_j && !node.data.r_child_ptr.is_null) {
          cur_read.header = node.data.r_child_ptr;
          right_bfs.push_back(write_block(cur_read, false));
          j_values.push_back(target_j - (L + 1));
        }
      }
    }

    // 3. Process Full Subtree BFS
    size_t full_idx = 0;
    while (full_idx < full_bfs.size()) {
      ORAMBlock node = full_bfs[full_idx];
      full_idx++;

      if (node.data.key == key) {
        result.push_back(node.data.value);
        if (!node.data.l_child_ptr.is_null) {
          cur_read.header = node.data.l_child_ptr;
          full_bfs.push_back(write_block(cur_read, false));
        }
        if (!node.data.r_child_ptr.is_null) {
          cur_read.header = node.data.r_child_ptr;
          full_bfs.push_back(write_block(cur_read, false));
        }
      }
    }
  }

  // --- Write-back logic ---
  auto check_and_update = [](std::vector<ORAMBlock>& vec, uint32_t child_block_id, uint32_t new_leaf_label) {
    for (auto& p : vec) {
      if (!p.data.l_child_ptr.is_null && p.data.l_child_ptr.block_id == child_block_id) {
        p.data.l_child_ptr.leaf_label = new_leaf_label; 
        return true;
      }
      if (!p.data.r_child_ptr.is_null && p.data.r_child_ptr.block_id == child_block_id) {
        p.data.r_child_ptr.leaf_label = new_leaf_label; 
        return true;
      }
    }
    return false;
  };

  for (int idx = (int)full_bfs.size() - 1; idx >= 0; idx--) {
    uint32_t new_leaf = write_block(full_bfs[idx], true).header.leaf_label;
    uint32_t cid = full_bfs[idx].header.block_id;
    if (check_and_update(full_bfs, cid, new_leaf)) continue;
    if (check_and_update(left_bfs, cid, new_leaf)) continue;
    check_and_update(right_bfs, cid, new_leaf);
  }

  for (int idx = (int)right_bfs.size() - 1; idx >= 0; idx--) {
    uint32_t new_leaf = write_block(right_bfs[idx], true).header.leaf_label;
    uint32_t cid = right_bfs[idx].header.block_id;
    if (!check_and_update(right_bfs, cid, new_leaf)) {
      check_and_update(common_path, cid, new_leaf);
    }
  }

  for (int idx = (int)left_bfs.size() - 1; idx >= 0; idx--) {
    uint32_t new_leaf = write_block(left_bfs[idx], true).header.leaf_label;
    uint32_t cid = left_bfs[idx].header.block_id;
    if (!check_and_update(left_bfs, cid, new_leaf)) {
      check_and_update(common_path, cid, new_leaf);
    }
  }
  
  for (int idx = (int)common_path.size() - 1; idx >= 0; idx--) {
    uint32_t new_leaf = write_block(common_path[idx], true).header.leaf_label;
    if (idx > 0) {
      if (!common_path[idx - 1].data.l_child_ptr.is_null &&
          common_path[idx - 1].data.l_child_ptr.block_id == common_path[idx].header.block_id) {
        common_path[idx - 1].data.l_child_ptr.leaf_label = new_leaf;
      } else if (!common_path[idx - 1].data.r_child_ptr.is_null &&
                 common_path[idx - 1].data.r_child_ptr.block_id == common_path[idx].header.block_id) {
        common_path[idx - 1].data.r_child_ptr.leaf_label = new_leaf;
      }
    } else {
      root.leaf_label = new_leaf;
    }
  }

  return result;
}

void Client::insert(uint64_t key, uint64_t value) {
  if (root.is_null) {
    ORAMBlock to_write;
    to_write.header.block_id = next_available_block_id();
    to_write.header.is_null = false;
    to_write.data.key = key;
    to_write.data.value = value;
    root = write_block(to_write, true).header;
    return;
  }

  std::vector<ORAMBlock> avl_history;

  ORAMBlock cur_read;
  cur_read.header = root;
  bool duplicate = false;

  // Read AVL tree path from PathORAM and add new Node to end
  while (true) {
    cur_read = write_block(cur_read, false);
    avl_history.push_back(cur_read);

    if (key > cur_read.data.key ||
        key == cur_read.data.key && value > cur_read.data.value) {
      if (cur_read.data.r_child_ptr.is_null) {
        ORAMBlock to_write;
        to_write.header.block_id = next_available_block_id();
        to_write.header.is_null = false;
        to_write.data.key = key;
        to_write.data.value = value;

        avl_history.back().data.r_child_ptr = to_write.header;
        avl_history.push_back(to_write);
        break;
      } else {
        cur_read.header = cur_read.data.r_child_ptr;
      }
    } else if (key < cur_read.data.key ||
               key == cur_read.data.key && value < cur_read.data.value) {
      if (cur_read.data.l_child_ptr.is_null) {
        ORAMBlock to_write;
        to_write.header.block_id = next_available_block_id();
        to_write.header.is_null = false;
        to_write.data.key = key;
        to_write.data.value = value;

        avl_history.back().data.l_child_ptr = to_write.header;
        avl_history.push_back(to_write);
        break;
      } else {
        cur_read.header = cur_read.data.l_child_ptr;
      }
    } else {
      duplicate = true;
      if (!cur_read.data.l_child_ptr.is_null) {
        cur_read.header = cur_read.data.l_child_ptr;
      } else {
        break;
      }
    }
  }

  // Reverse traverse and update heights and rebalance the tree
  for (int height = 1; !duplicate && height < avl_history.size(); height++) {
    int cur_node_index = avl_history.size() - 1 - height;

    if (!avl_history[cur_node_index].data.r_child_ptr.is_null &&
        avl_history[cur_node_index].data.r_child_ptr.block_id ==
            avl_history[cur_node_index + 1].header.block_id) {
      avl_history[cur_node_index].data.r_height =
          1 + std::max(avl_history[cur_node_index + 1].data.l_height,
                       avl_history[cur_node_index + 1].data.r_height);
      if (key < avl_history[cur_node_index].data.r_min_key_subtree) {
        avl_history[cur_node_index].data.r_min_key_subtree = key;
        avl_history[cur_node_index].data.r_min_key_count = 1;
      } else if (key == avl_history[cur_node_index].data.r_min_key_subtree) {
        avl_history[cur_node_index].data.r_min_key_count++;
      }
      if (key > avl_history[cur_node_index].data.r_max_key_subtree) {
        avl_history[cur_node_index].data.r_max_key_subtree = key;
        avl_history[cur_node_index].data.r_max_key_count = 1;
      } else if (key == avl_history[cur_node_index].data.r_max_key_subtree) {
        avl_history[cur_node_index].data.r_max_key_count++;
      }
      if (key == avl_history[cur_node_index].data.key) {
        avl_history[cur_node_index].data.r_same_key_size++;
      }
    } else {
      avl_history[cur_node_index].data.l_height =
          1 + std::max(avl_history[cur_node_index + 1].data.l_height,
                       avl_history[cur_node_index + 1].data.r_height);
      if (key < avl_history[cur_node_index].data.l_min_key_subtree) {
        avl_history[cur_node_index].data.l_min_key_subtree = key;
        avl_history[cur_node_index].data.l_min_key_count = 1;
      } else if (key == avl_history[cur_node_index].data.l_min_key_subtree) {
        avl_history[cur_node_index].data.l_min_key_count++;
      }
      if (key > avl_history[cur_node_index].data.l_max_key_subtree) {
        avl_history[cur_node_index].data.l_max_key_subtree = key;
        avl_history[cur_node_index].data.l_max_key_count = 1;
      } else if (key == avl_history[cur_node_index].data.l_max_key_subtree) {
        avl_history[cur_node_index].data.l_max_key_count++;
      }
      if (key == avl_history[cur_node_index].data.key) {
        avl_history[cur_node_index].data.l_same_key_size++;
      }
    }

    int balance_factor = avl_history[cur_node_index].data.r_height -
                         avl_history[cur_node_index].data.l_height;

    if (std::abs(balance_factor) > 1) {
      bool child_is_right =
          (!avl_history[cur_node_index].data.r_child_ptr.is_null &&
           avl_history[cur_node_index].data.r_child_ptr.block_id ==
               avl_history[cur_node_index + 1].header.block_id);

      bool grandchild_is_right =
          (!avl_history[cur_node_index + 1].data.r_child_ptr.is_null &&
           avl_history[cur_node_index + 1].data.r_child_ptr.block_id ==
               avl_history[cur_node_index + 2].header.block_id);

      if (child_is_right && grandchild_is_right) {
        rotate_right_right(avl_history, cur_node_index);
      } else if (!child_is_right && !grandchild_is_right) {
        rotate_left_left(avl_history, cur_node_index);
      } else if (child_is_right && !grandchild_is_right) {
        rotate_right_left(avl_history, cur_node_index);
      } else {
        rotate_left_right(avl_history, cur_node_index);
      }
    }
  }

  // Write back to PathORAM in reverse order, updating new leaf references
  for (int index = avl_history.size() - 1; index >= 0; index--) {
    uint32_t new_leaf_label =
        write_block(avl_history[index], true).header.leaf_label;

    if (index > 0 && !avl_history[index - 1].data.l_child_ptr.is_null &&
        avl_history[index - 1].data.l_child_ptr.block_id ==
            avl_history[index].header.block_id) {
      avl_history[index - 1].data.l_child_ptr.leaf_label = new_leaf_label;
    } else if (index > 0 && !avl_history[index - 1].data.r_child_ptr.is_null &&
               avl_history[index - 1].data.r_child_ptr.block_id ==
                   avl_history[index].header.block_id) {
      avl_history[index - 1].data.r_child_ptr.leaf_label = new_leaf_label;
    } else if (index > 1 && !avl_history[index - 2].data.l_child_ptr.is_null &&
               avl_history[index - 2].data.l_child_ptr.block_id ==
                   avl_history[index].header.block_id) {
      avl_history[index - 2].data.l_child_ptr.leaf_label = new_leaf_label;
    } else if (index > 1 && !avl_history[index - 2].data.r_child_ptr.is_null &&
               avl_history[index - 2].data.r_child_ptr.block_id ==
                   avl_history[index].header.block_id) {
      avl_history[index - 2].data.r_child_ptr.leaf_label = new_leaf_label;
    } else if (index == 0) {
      root.leaf_label = new_leaf_label;
    }
  }
}

void Client::rotate_right_right(std::vector<ORAMBlock>& avl_history,
                                int cur_node_index) {
  avl_history[cur_node_index].data.r_child_ptr =
      avl_history[cur_node_index + 1].data.l_child_ptr;
  avl_history[cur_node_index].data.r_height =
      avl_history[cur_node_index + 1].data.l_height;
  avl_history[cur_node_index].data.r_max_key_count =
      avl_history[cur_node_index + 1].data.l_max_key_count;
  avl_history[cur_node_index].data.r_max_key_subtree =
      avl_history[cur_node_index + 1].data.l_max_key_subtree;
  avl_history[cur_node_index].data.r_min_key_count =
      avl_history[cur_node_index + 1].data.l_min_key_count;
  avl_history[cur_node_index].data.r_min_key_subtree =
      avl_history[cur_node_index + 1].data.l_min_key_subtree;
  if (avl_history[cur_node_index].data.r_min_key_subtree ==
      avl_history[cur_node_index].data.key) {
    avl_history[cur_node_index].data.r_same_key_size =
        avl_history[cur_node_index].data.r_min_key_count;
  } else {
    avl_history[cur_node_index].data.r_same_key_size = 0;
  }

  avl_history[cur_node_index + 1].data.l_child_ptr =
      avl_history[cur_node_index].header;
  avl_history[cur_node_index + 1].data.l_height =
      1 + std::max(avl_history[cur_node_index].data.l_height,
                   avl_history[cur_node_index].data.r_height);
  if (avl_history[cur_node_index].data.r_max_key_count != 0) {
    if (avl_history[cur_node_index].data.r_max_key_subtree ==
        avl_history[cur_node_index + 1].data.key) {
      if (avl_history[cur_node_index].data.key ==
          avl_history[cur_node_index].data.r_max_key_subtree) {
        avl_history[cur_node_index + 1].data.l_same_key_size =
            1 + avl_history[cur_node_index].data.l_same_key_size +
            avl_history[cur_node_index].data.r_same_key_size;
      } else {
        avl_history[cur_node_index + 1].data.l_same_key_size =
            avl_history[cur_node_index].data.r_max_key_count;
      }
    } else {
      avl_history[cur_node_index + 1].data.l_same_key_size = 0;
    }
  } else {
    if (avl_history[cur_node_index].data.key ==
        avl_history[cur_node_index + 1].data.key) {
      avl_history[cur_node_index + 1].data.l_same_key_size =
          1 + avl_history[cur_node_index].data.l_same_key_size +
          avl_history[cur_node_index].data.r_same_key_size;
    } else {
      avl_history[cur_node_index + 1].data.l_same_key_size = 0;
    }
  }

  // 1. Calculate MAX Key for B's new Left Child (A)
  // Logic: Check A's Right Subtree (which was just updated).
  // If it exists, that's the Max. Otherwise, A itself is the Max.
  avl_history[cur_node_index + 1].data.l_max_key_subtree =
      (avl_history[cur_node_index].data.r_max_key_count > 0)
          ? avl_history[cur_node_index].data.r_max_key_subtree
          : avl_history[cur_node_index].data.key;

  // Calculate Count for that Max Key
  if (avl_history[cur_node_index].data.r_max_key_count > 0) {
    // STRADDLE CHECK: Is the Right Subtree's Max Key equal to A's Key?
    if (avl_history[cur_node_index].data.r_max_key_subtree ==
        avl_history[cur_node_index].data.key) {
      // Yes: The 'Max' group includes A. Sum A's self(1) + Left Matches + Right
      // Matches. Note: Ensure A.r_same_key_size was updated to
      // B.l_same_key_size prior to this if relying on it.
      avl_history[cur_node_index + 1].data.l_max_key_count =
          1 + avl_history[cur_node_index].data.l_same_key_size +
          avl_history[cur_node_index].data.r_same_key_size;
    } else {
      // No: The Max is strictly in the right subtree.
      avl_history[cur_node_index + 1].data.l_max_key_count =
          avl_history[cur_node_index].data.r_max_key_count;
    }
  } else {
    // No right child: A is the Max. Count is Self(1) + Left Matches.
    avl_history[cur_node_index + 1].data.l_max_key_count =
        1 + avl_history[cur_node_index].data.l_same_key_size;
  }

  // 2. Calculate MIN Key for B's new Left Child (A)
  // Logic: Check A's Left Subtree (untouched by rotation).
  avl_history[cur_node_index + 1].data.l_min_key_subtree =
      (avl_history[cur_node_index].data.l_min_key_count > 0)
          ? avl_history[cur_node_index].data.l_min_key_subtree
          : avl_history[cur_node_index].data.key;

  // Calculate Count for that Min Key
  if (avl_history[cur_node_index].data.l_min_key_count > 0) {
    // STRADDLE CHECK: Is the Left Subtree's Min Key equal to A's Key?
    if (avl_history[cur_node_index].data.l_min_key_subtree ==
        avl_history[cur_node_index].data.key) {
      // Yes: The 'Min' group includes A. Sum Total Matches.
      avl_history[cur_node_index + 1].data.l_min_key_count =
          1 + avl_history[cur_node_index].data.l_same_key_size +
          avl_history[cur_node_index].data.r_same_key_size;
    } else {
      // No: The Min is strictly in the left subtree.
      avl_history[cur_node_index + 1].data.l_min_key_count =
          avl_history[cur_node_index].data.l_min_key_count;
    }
  } else {
    // No left child: A is the Min. Count is Self(1) + Right Matches.
    avl_history[cur_node_index + 1].data.l_min_key_count =
        1 + avl_history[cur_node_index].data.r_same_key_size;
  }

  std::swap(avl_history[cur_node_index + 1], avl_history[cur_node_index]);

  if (cur_node_index == 0) {
    root = avl_history[cur_node_index].header;
  } else {
    if (!avl_history[cur_node_index - 1].data.r_child_ptr.is_null &&
        avl_history[cur_node_index - 1].data.r_child_ptr.block_id ==
            avl_history[cur_node_index + 1].header.block_id) {
      avl_history[cur_node_index - 1].data.r_child_ptr =
          avl_history[cur_node_index].header;
    } else {
      avl_history[cur_node_index - 1].data.l_child_ptr =
          avl_history[cur_node_index].header;
    }
  }
}

void Client::rotate_left_left(std::vector<ORAMBlock>& avl_history,
                              int cur_node_index) {
  // 1. A (old root) adopts B's right child as its new left child
  avl_history[cur_node_index].data.l_child_ptr =
      avl_history[cur_node_index + 1].data.r_child_ptr;
  avl_history[cur_node_index].data.l_height =
      avl_history[cur_node_index + 1].data.r_height;
  avl_history[cur_node_index].data.l_max_key_count =
      avl_history[cur_node_index + 1].data.r_max_key_count;
  avl_history[cur_node_index].data.l_max_key_subtree =
      avl_history[cur_node_index + 1].data.r_max_key_subtree;
  avl_history[cur_node_index].data.l_min_key_count =
      avl_history[cur_node_index + 1].data.r_min_key_count;
  avl_history[cur_node_index].data.l_min_key_subtree =
      avl_history[cur_node_index + 1].data.r_min_key_subtree;

  // Check if the max of the new left subtree straddles A's key
  if (avl_history[cur_node_index].data.l_max_key_subtree ==
      avl_history[cur_node_index].data.key) {
    avl_history[cur_node_index].data.l_same_key_size =
        avl_history[cur_node_index].data.l_max_key_count;
  } else {
    avl_history[cur_node_index].data.l_same_key_size = 0;
  }

  // 2. B (new root) adopts A as its new right child
  avl_history[cur_node_index + 1].data.r_child_ptr =
      avl_history[cur_node_index].header;
  avl_history[cur_node_index + 1].data.r_height =
      1 + std::max(avl_history[cur_node_index].data.l_height,
                   avl_history[cur_node_index].data.r_height);

  // Calculate B's new r_same_key_size by looking at A's min key
  if (avl_history[cur_node_index].data.l_min_key_count != 0) {
    if (avl_history[cur_node_index].data.l_min_key_subtree ==
        avl_history[cur_node_index + 1].data.key) {
      if (avl_history[cur_node_index].data.key ==
          avl_history[cur_node_index].data.l_min_key_subtree) {
        avl_history[cur_node_index + 1].data.r_same_key_size =
            1 + avl_history[cur_node_index].data.l_same_key_size +
            avl_history[cur_node_index].data.r_same_key_size;
      } else {
        avl_history[cur_node_index + 1].data.r_same_key_size =
            avl_history[cur_node_index].data.l_min_key_count;
      }
    } else {
      avl_history[cur_node_index + 1].data.r_same_key_size = 0;
    }
  } else {
    if (avl_history[cur_node_index].data.key ==
        avl_history[cur_node_index + 1].data.key) {
      avl_history[cur_node_index + 1].data.r_same_key_size =
          1 + avl_history[cur_node_index].data.l_same_key_size +
          avl_history[cur_node_index].data.r_same_key_size;
    } else {
      avl_history[cur_node_index + 1].data.r_same_key_size = 0;
    }
  }

  // 3. Calculate MIN Key for B's new Right Child (A)
  // Logic: Check A's Left Subtree (which was just updated).
  avl_history[cur_node_index + 1].data.r_min_key_subtree =
      (avl_history[cur_node_index].data.l_min_key_count > 0)
          ? avl_history[cur_node_index].data.l_min_key_subtree
          : avl_history[cur_node_index].data.key;

  if (avl_history[cur_node_index].data.l_min_key_count > 0) {
    // STRADDLE CHECK: Is the Left Subtree's Min Key equal to A's Key?
    if (avl_history[cur_node_index].data.l_min_key_subtree ==
        avl_history[cur_node_index].data.key) {
      avl_history[cur_node_index + 1].data.r_min_key_count =
          1 + avl_history[cur_node_index].data.l_same_key_size +
          avl_history[cur_node_index].data.r_same_key_size;
    } else {
      avl_history[cur_node_index + 1].data.r_min_key_count =
          avl_history[cur_node_index].data.l_min_key_count;
    }
  } else {
    // No left child: A is the Min. Count is Self(1) + Right Matches.
    avl_history[cur_node_index + 1].data.r_min_key_count =
        1 + avl_history[cur_node_index].data.r_same_key_size;
  }

  // 4. Calculate MAX Key for B's new Right Child (A)
  // Logic: Check A's Right Subtree (untouched by rotation).
  avl_history[cur_node_index + 1].data.r_max_key_subtree =
      (avl_history[cur_node_index].data.r_max_key_count > 0)
          ? avl_history[cur_node_index].data.r_max_key_subtree
          : avl_history[cur_node_index].data.key;

  if (avl_history[cur_node_index].data.r_max_key_count > 0) {
    // STRADDLE CHECK: Is the Right Subtree's Max Key equal to A's Key?
    if (avl_history[cur_node_index].data.r_max_key_subtree ==
        avl_history[cur_node_index].data.key) {
      avl_history[cur_node_index + 1].data.r_max_key_count =
          1 + avl_history[cur_node_index].data.l_same_key_size +
          avl_history[cur_node_index].data.r_same_key_size;
    } else {
      avl_history[cur_node_index + 1].data.r_max_key_count =
          avl_history[cur_node_index].data.r_max_key_count;
    }
  } else {
    // No right child: A is the Max. Count is Self(1) + Left Matches.
    avl_history[cur_node_index + 1].data.r_max_key_count =
        1 + avl_history[cur_node_index].data.l_same_key_size;
  }

  std::swap(avl_history[cur_node_index + 1], avl_history[cur_node_index]);

  if (cur_node_index == 0) {
    root = avl_history[cur_node_index].header;
  } else {
    if (!avl_history[cur_node_index - 1].data.r_child_ptr.is_null &&
        avl_history[cur_node_index - 1].data.r_child_ptr.block_id ==
            avl_history[cur_node_index + 1].header.block_id) {
      avl_history[cur_node_index - 1].data.r_child_ptr =
          avl_history[cur_node_index].header;
    } else {
      avl_history[cur_node_index - 1].data.l_child_ptr =
          avl_history[cur_node_index].header;
    }
  }
}

void Client::rotate_right_left(std::vector<ORAMBlock>& avl_history,
                               int cur_node_index) {
  avl_history[cur_node_index].data.r_child_ptr =
      avl_history[cur_node_index + 2].data.l_child_ptr;
  avl_history[cur_node_index].data.r_height =
      avl_history[cur_node_index + 2].data.l_height;
  avl_history[cur_node_index].data.r_min_key_subtree =
      avl_history[cur_node_index + 2].data.l_min_key_subtree;
  avl_history[cur_node_index].data.r_min_key_count =
      avl_history[cur_node_index + 2].data.l_min_key_count;
  avl_history[cur_node_index].data.r_max_key_subtree =
      avl_history[cur_node_index + 2].data.l_max_key_subtree;
  avl_history[cur_node_index].data.r_max_key_count =
      avl_history[cur_node_index + 2].data.l_max_key_count;
  if (avl_history[cur_node_index].data.r_min_key_subtree ==
      avl_history[cur_node_index].data.key) {
    avl_history[cur_node_index].data.r_same_key_size =
        avl_history[cur_node_index].data.r_min_key_count;
  } else {
    avl_history[cur_node_index].data.r_same_key_size = 0;
  }

  avl_history[cur_node_index + 1].data.l_child_ptr =
      avl_history[cur_node_index + 2].data.r_child_ptr;
  avl_history[cur_node_index + 1].data.l_height =
      avl_history[cur_node_index + 2].data.r_height;
  avl_history[cur_node_index + 1].data.l_min_key_subtree =
      avl_history[cur_node_index + 2].data.r_min_key_subtree;
  avl_history[cur_node_index + 1].data.l_min_key_count =
      avl_history[cur_node_index + 2].data.r_min_key_count;
  avl_history[cur_node_index + 1].data.l_max_key_subtree =
      avl_history[cur_node_index + 2].data.r_max_key_subtree;
  avl_history[cur_node_index + 1].data.l_max_key_count =
      avl_history[cur_node_index + 2].data.r_max_key_count;
  if (avl_history[cur_node_index + 1].data.l_max_key_subtree ==
      avl_history[cur_node_index + 1].data.key) {
    avl_history[cur_node_index + 1].data.l_same_key_size =
        avl_history[cur_node_index + 1].data.l_max_key_count;
  } else {
    avl_history[cur_node_index + 1].data.l_same_key_size = 0;
  }

  // 3. C (cur_node_index + 2) adopts A as its left child
  avl_history[cur_node_index + 2].data.l_child_ptr =
      avl_history[cur_node_index].header;
  avl_history[cur_node_index + 2].data.l_height =
      1 + std::max(avl_history[cur_node_index].data.l_height,
                   avl_history[cur_node_index].data.r_height);

  // Calculate C's l_same_key_size
  if (avl_history[cur_node_index].data.r_max_key_count != 0) {
    if (avl_history[cur_node_index].data.r_max_key_subtree ==
        avl_history[cur_node_index + 2].data.key) {
      if (avl_history[cur_node_index].data.key ==
          avl_history[cur_node_index].data.r_max_key_subtree) {
        avl_history[cur_node_index + 2].data.l_same_key_size =
            1 + avl_history[cur_node_index].data.l_same_key_size +
            avl_history[cur_node_index].data.r_same_key_size;
      } else {
        avl_history[cur_node_index + 2].data.l_same_key_size =
            avl_history[cur_node_index].data.r_max_key_count;
      }
    } else {
      avl_history[cur_node_index + 2].data.l_same_key_size = 0;
    }
  } else {
    if (avl_history[cur_node_index].data.key ==
        avl_history[cur_node_index + 2].data.key) {
      avl_history[cur_node_index + 2].data.l_same_key_size =
          1 + avl_history[cur_node_index].data.l_same_key_size +
          avl_history[cur_node_index].data.r_same_key_size;
    } else {
      avl_history[cur_node_index + 2].data.l_same_key_size = 0;
    }
  }

  // Calculate C's l_max
  avl_history[cur_node_index + 2].data.l_max_key_subtree =
      (avl_history[cur_node_index].data.r_max_key_count > 0)
          ? avl_history[cur_node_index].data.r_max_key_subtree
          : avl_history[cur_node_index].data.key;
  if (avl_history[cur_node_index].data.r_max_key_count > 0) {
    if (avl_history[cur_node_index].data.r_max_key_subtree ==
        avl_history[cur_node_index].data.key) {
      avl_history[cur_node_index + 2].data.l_max_key_count =
          1 + avl_history[cur_node_index].data.l_same_key_size +
          avl_history[cur_node_index].data.r_same_key_size;
    } else {
      avl_history[cur_node_index + 2].data.l_max_key_count =
          avl_history[cur_node_index].data.r_max_key_count;
    }
  } else {
    avl_history[cur_node_index + 2].data.l_max_key_count =
        1 + avl_history[cur_node_index].data.l_same_key_size;
  }

  // Calculate C's l_min
  avl_history[cur_node_index + 2].data.l_min_key_subtree =
      (avl_history[cur_node_index].data.l_min_key_count > 0)
          ? avl_history[cur_node_index].data.l_min_key_subtree
          : avl_history[cur_node_index].data.key;
  if (avl_history[cur_node_index].data.l_min_key_count > 0) {
    if (avl_history[cur_node_index].data.l_min_key_subtree ==
        avl_history[cur_node_index].data.key) {
      avl_history[cur_node_index + 2].data.l_min_key_count =
          1 + avl_history[cur_node_index].data.l_same_key_size +
          avl_history[cur_node_index].data.r_same_key_size;
    } else {
      avl_history[cur_node_index + 2].data.l_min_key_count =
          avl_history[cur_node_index].data.l_min_key_count;
    }
  } else {
    avl_history[cur_node_index + 2].data.l_min_key_count =
        1 + avl_history[cur_node_index].data.r_same_key_size;
  }

  // 4. C (cur_node_index + 2) adopts B as its right child
  avl_history[cur_node_index + 2].data.r_child_ptr =
      avl_history[cur_node_index + 1].header;
  avl_history[cur_node_index + 2].data.r_height =
      1 + std::max(avl_history[cur_node_index + 1].data.l_height,
                   avl_history[cur_node_index + 1].data.r_height);

  // Calculate C's r_same_key_size
  if (avl_history[cur_node_index + 1].data.l_min_key_count != 0) {
    if (avl_history[cur_node_index + 1].data.l_min_key_subtree ==
        avl_history[cur_node_index + 2].data.key) {
      if (avl_history[cur_node_index + 1].data.key ==
          avl_history[cur_node_index + 1].data.l_min_key_subtree) {
        avl_history[cur_node_index + 2].data.r_same_key_size =
            1 + avl_history[cur_node_index + 1].data.l_same_key_size +
            avl_history[cur_node_index + 1].data.r_same_key_size;
      } else {
        avl_history[cur_node_index + 2].data.r_same_key_size =
            avl_history[cur_node_index + 1].data.l_min_key_count;
      }
    } else {
      avl_history[cur_node_index + 2].data.r_same_key_size = 0;
    }
  } else {
    if (avl_history[cur_node_index + 1].data.key ==
        avl_history[cur_node_index + 2].data.key) {
      avl_history[cur_node_index + 2].data.r_same_key_size =
          1 + avl_history[cur_node_index + 1].data.l_same_key_size +
          avl_history[cur_node_index + 1].data.r_same_key_size;
    } else {
      avl_history[cur_node_index + 2].data.r_same_key_size = 0;
    }
  }

  // Calculate C's r_min
  avl_history[cur_node_index + 2].data.r_min_key_subtree =
      (avl_history[cur_node_index + 1].data.l_min_key_count > 0)
          ? avl_history[cur_node_index + 1].data.l_min_key_subtree
          : avl_history[cur_node_index + 1].data.key;
  if (avl_history[cur_node_index + 1].data.l_min_key_count > 0) {
    if (avl_history[cur_node_index + 1].data.l_min_key_subtree ==
        avl_history[cur_node_index + 1].data.key) {
      avl_history[cur_node_index + 2].data.r_min_key_count =
          1 + avl_history[cur_node_index + 1].data.l_same_key_size +
          avl_history[cur_node_index + 1].data.r_same_key_size;
    } else {
      avl_history[cur_node_index + 2].data.r_min_key_count =
          avl_history[cur_node_index + 1].data.l_min_key_count;
    }
  } else {
    avl_history[cur_node_index + 2].data.r_min_key_count =
        1 + avl_history[cur_node_index + 1].data.r_same_key_size;
  }

  // Calculate C's r_max
  avl_history[cur_node_index + 2].data.r_max_key_subtree =
      (avl_history[cur_node_index + 1].data.r_max_key_count > 0)
          ? avl_history[cur_node_index + 1].data.r_max_key_subtree
          : avl_history[cur_node_index + 1].data.key;
  if (avl_history[cur_node_index + 1].data.r_max_key_count > 0) {
    if (avl_history[cur_node_index + 1].data.r_max_key_subtree ==
        avl_history[cur_node_index + 1].data.key) {
      avl_history[cur_node_index + 2].data.r_max_key_count =
          1 + avl_history[cur_node_index + 1].data.l_same_key_size +
          avl_history[cur_node_index + 1].data.r_same_key_size;
    } else {
      avl_history[cur_node_index + 2].data.r_max_key_count =
          avl_history[cur_node_index + 1].data.r_max_key_count;
    }
  } else {
    avl_history[cur_node_index + 2].data.r_max_key_count =
        1 + avl_history[cur_node_index + 1].data.l_same_key_size;
  }

  std::swap(avl_history[cur_node_index + 2], avl_history[cur_node_index]);

  if (cur_node_index == 0) {
    root = avl_history[cur_node_index].header;
  } else {
    if (!avl_history[cur_node_index - 1].data.r_child_ptr.is_null &&
        avl_history[cur_node_index - 1].data.r_child_ptr.block_id ==
            avl_history[cur_node_index + 2].header.block_id) {
      avl_history[cur_node_index - 1].data.r_child_ptr =
          avl_history[cur_node_index].header;
    } else {
      avl_history[cur_node_index - 1].data.l_child_ptr =
          avl_history[cur_node_index].header;
    }
  }
}

void Client::rotate_left_right(std::vector<ORAMBlock>& avl_history,
                               int cur_node_index) {
  // 1. A (cur_node_index) adopts C's right child as its new left child
  avl_history[cur_node_index].data.l_child_ptr =
      avl_history[cur_node_index + 2].data.r_child_ptr;
  avl_history[cur_node_index].data.l_height =
      avl_history[cur_node_index + 2].data.r_height;
  avl_history[cur_node_index].data.l_min_key_subtree =
      avl_history[cur_node_index + 2].data.r_min_key_subtree;
  avl_history[cur_node_index].data.l_min_key_count =
      avl_history[cur_node_index + 2].data.r_min_key_count;
  avl_history[cur_node_index].data.l_max_key_subtree =
      avl_history[cur_node_index + 2].data.r_max_key_subtree;
  avl_history[cur_node_index].data.l_max_key_count =
      avl_history[cur_node_index + 2].data.r_max_key_count;

  // Check if the max of the new left subtree straddles A's key
  if (avl_history[cur_node_index].data.l_max_key_subtree ==
      avl_history[cur_node_index].data.key) {
    avl_history[cur_node_index].data.l_same_key_size =
        avl_history[cur_node_index].data.l_max_key_count;
  } else {
    avl_history[cur_node_index].data.l_same_key_size = 0;
  }

  // 2. B (cur_node_index + 1) adopts C's left child as its new right child
  avl_history[cur_node_index + 1].data.r_child_ptr =
      avl_history[cur_node_index + 2].data.l_child_ptr;
  avl_history[cur_node_index + 1].data.r_height =
      avl_history[cur_node_index + 2].data.l_height;
  avl_history[cur_node_index + 1].data.r_min_key_subtree =
      avl_history[cur_node_index + 2].data.l_min_key_subtree;
  avl_history[cur_node_index + 1].data.r_min_key_count =
      avl_history[cur_node_index + 2].data.l_min_key_count;
  avl_history[cur_node_index + 1].data.r_max_key_subtree =
      avl_history[cur_node_index + 2].data.l_max_key_subtree;
  avl_history[cur_node_index + 1].data.r_max_key_count =
      avl_history[cur_node_index + 2].data.l_max_key_count;

  // Check if the min of the new right subtree straddles B's key
  if (avl_history[cur_node_index + 1].data.r_min_key_subtree ==
      avl_history[cur_node_index + 1].data.key) {
    avl_history[cur_node_index + 1].data.r_same_key_size =
        avl_history[cur_node_index + 1].data.r_min_key_count;
  } else {
    avl_history[cur_node_index + 1].data.r_same_key_size = 0;
  }

  // 3. C (cur_node_index + 2) adopts B as its left child
  avl_history[cur_node_index + 2].data.l_child_ptr =
      avl_history[cur_node_index + 1].header;
  avl_history[cur_node_index + 2].data.l_height =
      1 + std::max(avl_history[cur_node_index + 1].data.l_height,
                   avl_history[cur_node_index + 1].data.r_height);

  // Calculate C's l_same_key_size (looking at B's right side)
  if (avl_history[cur_node_index + 1].data.r_max_key_count != 0) {
    if (avl_history[cur_node_index + 1].data.r_max_key_subtree ==
        avl_history[cur_node_index + 2].data.key) {
      if (avl_history[cur_node_index + 1].data.key ==
          avl_history[cur_node_index + 1].data.r_max_key_subtree) {
        avl_history[cur_node_index + 2].data.l_same_key_size =
            1 + avl_history[cur_node_index + 1].data.l_same_key_size +
            avl_history[cur_node_index + 1].data.r_same_key_size;
      } else {
        avl_history[cur_node_index + 2].data.l_same_key_size =
            avl_history[cur_node_index + 1].data.r_max_key_count;
      }
    } else {
      avl_history[cur_node_index + 2].data.l_same_key_size = 0;
    }
  } else {
    if (avl_history[cur_node_index + 1].data.key ==
        avl_history[cur_node_index + 2].data.key) {
      avl_history[cur_node_index + 2].data.l_same_key_size =
          1 + avl_history[cur_node_index + 1].data.l_same_key_size +
          avl_history[cur_node_index + 1].data.r_same_key_size;
    } else {
      avl_history[cur_node_index + 2].data.l_same_key_size = 0;
    }
  }

  // Calculate C's l_max
  avl_history[cur_node_index + 2].data.l_max_key_subtree =
      (avl_history[cur_node_index + 1].data.r_max_key_count > 0)
          ? avl_history[cur_node_index + 1].data.r_max_key_subtree
          : avl_history[cur_node_index + 1].data.key;
  if (avl_history[cur_node_index + 1].data.r_max_key_count > 0) {
    if (avl_history[cur_node_index + 1].data.r_max_key_subtree ==
        avl_history[cur_node_index + 1].data.key) {
      avl_history[cur_node_index + 2].data.l_max_key_count =
          1 + avl_history[cur_node_index + 1].data.l_same_key_size +
          avl_history[cur_node_index + 1].data.r_same_key_size;
    } else {
      avl_history[cur_node_index + 2].data.l_max_key_count =
          avl_history[cur_node_index + 1].data.r_max_key_count;
    }
  } else {
    avl_history[cur_node_index + 2].data.l_max_key_count =
        1 + avl_history[cur_node_index + 1].data.l_same_key_size;
  }

  // Calculate C's l_min
  avl_history[cur_node_index + 2].data.l_min_key_subtree =
      (avl_history[cur_node_index + 1].data.l_min_key_count > 0)
          ? avl_history[cur_node_index + 1].data.l_min_key_subtree
          : avl_history[cur_node_index + 1].data.key;
  if (avl_history[cur_node_index + 1].data.l_min_key_count > 0) {
    if (avl_history[cur_node_index + 1].data.l_min_key_subtree ==
        avl_history[cur_node_index + 1].data.key) {
      avl_history[cur_node_index + 2].data.l_min_key_count =
          1 + avl_history[cur_node_index + 1].data.l_same_key_size +
          avl_history[cur_node_index + 1].data.r_same_key_size;
    } else {
      avl_history[cur_node_index + 2].data.l_min_key_count =
          avl_history[cur_node_index + 1].data.l_min_key_count;
    }
  } else {
    avl_history[cur_node_index + 2].data.l_min_key_count =
        1 + avl_history[cur_node_index + 1].data.r_same_key_size;
  }

  // 4. C (cur_node_index + 2) adopts A as its right child
  avl_history[cur_node_index + 2].data.r_child_ptr =
      avl_history[cur_node_index].header;
  avl_history[cur_node_index + 2].data.r_height =
      1 + std::max(avl_history[cur_node_index].data.l_height,
                   avl_history[cur_node_index].data.r_height);

  // Calculate C's r_same_key_size (looking at A's left side)
  if (avl_history[cur_node_index].data.l_min_key_count != 0) {
    if (avl_history[cur_node_index].data.l_min_key_subtree ==
        avl_history[cur_node_index + 2].data.key) {
      if (avl_history[cur_node_index].data.key ==
          avl_history[cur_node_index].data.l_min_key_subtree) {
        avl_history[cur_node_index + 2].data.r_same_key_size =
            1 + avl_history[cur_node_index].data.l_same_key_size +
            avl_history[cur_node_index].data.r_same_key_size;
      } else {
        avl_history[cur_node_index + 2].data.r_same_key_size =
            avl_history[cur_node_index].data.l_min_key_count;
      }
    } else {
      avl_history[cur_node_index + 2].data.r_same_key_size = 0;
    }
  } else {
    if (avl_history[cur_node_index].data.key ==
        avl_history[cur_node_index + 2].data.key) {
      avl_history[cur_node_index + 2].data.r_same_key_size =
          1 + avl_history[cur_node_index].data.l_same_key_size +
          avl_history[cur_node_index].data.r_same_key_size;
    } else {
      avl_history[cur_node_index + 2].data.r_same_key_size = 0;
    }
  }

  // Calculate C's r_min
  avl_history[cur_node_index + 2].data.r_min_key_subtree =
      (avl_history[cur_node_index].data.l_min_key_count > 0)
          ? avl_history[cur_node_index].data.l_min_key_subtree
          : avl_history[cur_node_index].data.key;
  if (avl_history[cur_node_index].data.l_min_key_count > 0) {
    if (avl_history[cur_node_index].data.l_min_key_subtree ==
        avl_history[cur_node_index].data.key) {
      avl_history[cur_node_index + 2].data.r_min_key_count =
          1 + avl_history[cur_node_index].data.l_same_key_size +
          avl_history[cur_node_index].data.r_same_key_size;
    } else {
      avl_history[cur_node_index + 2].data.r_min_key_count =
          avl_history[cur_node_index].data.l_min_key_count;
    }
  } else {
    avl_history[cur_node_index + 2].data.r_min_key_count =
        1 + avl_history[cur_node_index].data.r_same_key_size;
  }

  // Calculate C's r_max
  avl_history[cur_node_index + 2].data.r_max_key_subtree =
      (avl_history[cur_node_index].data.r_max_key_count > 0)
          ? avl_history[cur_node_index].data.r_max_key_subtree
          : avl_history[cur_node_index].data.key;
  if (avl_history[cur_node_index].data.r_max_key_count > 0) {
    if (avl_history[cur_node_index].data.r_max_key_subtree ==
        avl_history[cur_node_index].data.key) {
      avl_history[cur_node_index + 2].data.r_max_key_count =
          1 + avl_history[cur_node_index].data.l_same_key_size +
          avl_history[cur_node_index].data.r_same_key_size;
    } else {
      avl_history[cur_node_index + 2].data.r_max_key_count =
          avl_history[cur_node_index].data.r_max_key_count;
    }
  } else {
    avl_history[cur_node_index + 2].data.r_max_key_count =
        1 + avl_history[cur_node_index].data.l_same_key_size;
  }

  std::swap(avl_history[cur_node_index + 2], avl_history[cur_node_index]);

  if (cur_node_index == 0) {
    root = avl_history[cur_node_index].header;
  } else {
    if (!avl_history[cur_node_index - 1].data.r_child_ptr.is_null &&
        avl_history[cur_node_index - 1].data.r_child_ptr.block_id ==
            avl_history[cur_node_index + 2].header.block_id) {
      avl_history[cur_node_index - 1].data.r_child_ptr =
          avl_history[cur_node_index].header;
    } else {
      avl_history[cur_node_index - 1].data.l_child_ptr =
          avl_history[cur_node_index].header;
    }
  }
}

ORAMBlock Client::write_block(ORAMBlock to_write, bool write) {
  uint32_t initial_leaf_num = to_write.header.leaf_label;
  read_block_path(initial_leaf_num);

  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_int_distribution<> distr(0, (1u << tree_height) - 1);

  bool block_found = false;
  for (size_t i = 0; i < stash_full.size(); i++) {
    if (stash_full[i]) {
      if (stash[i].header.block_id == to_write.header.block_id) {
        block_found = true;
        if (write) {
          stash[i] = to_write;
        }
        stash[i].header.leaf_label = distr(gen);
        to_write = stash[i];
      }
    }
  }
  if (!block_found) {
    if (!write) {
      to_write.header.is_null = true;
    } else {
      for (size_t i = 0; i < stash_full.size(); i++) {
        if (!stash_full[i]) {
          stash[i] = to_write;
          stash[i].header.leaf_label = distr(gen);
          stash_full[i] = true;
          to_write = stash[i];
          break;
        }
      }
    }
  }
  evict(initial_leaf_num);
  return to_write;
}

void Client::read_block_path(uint32_t leaf_num) {
  std::vector<size_t> index_list;
  size_t leaf_path = (1u << tree_height) + leaf_num;
  while (leaf_path > 0) {
    index_list.push_back(leaf_path);
    leaf_path /= 2;
  }
  std::vector<uint8_t> temp = server->read_buckets(index_list);
  size_t stash_entry = 0;
  for (int cur_bucket = 0; cur_bucket < index_list.size(); cur_bucket++) {
    for (int cur_block = 0; cur_block < blocks_per_bucket; cur_block++) {
      ORAMBlock temp_block;
      std::memcpy(&temp_block,
                  &temp[(cur_bucket * blocks_per_bucket + cur_block) *
                        block_size_bytes],
                  sizeof(ORAMBlock));
      if (temp_block.header.is_null) {
        continue;
      }
      while (stash_entry < stash.size() && stash_full[stash_entry]) {
        stash_entry++;
      }
      if (stash_entry < stash.size()) {
        stash[stash_entry] = temp_block;
        stash_full[stash_entry] = true;
      } else {
        throw std::runtime_error("Stash overflow");
      }
    }
  }
}

bool Client::is_on_path(size_t index, size_t leaf_num) {
  size_t leaf_path = (1u << tree_height) + leaf_num;
  while (leaf_path >= index && leaf_path > 0) {
    if (leaf_path == index) {
      return true;
    }
    leaf_path /= 2;
  }
  return false;
}

void Client::evict(uint32_t leaf_num) {
  std::vector<size_t> index_list;
  std::vector<uint8_t> data((tree_height + 1) * bucket_size_bytes);
  size_t leaf_path = (1u << tree_height) + leaf_num;

  while (leaf_path > 0) {
    index_list.push_back(leaf_path);
    leaf_path /= 2;
  }
  for (int i = 0; i < index_list.size(); i++) {
    size_t cur_block_in_bucket = 0;
    for (int j = 0;
         j < stash_full.size() && cur_block_in_bucket < blocks_per_bucket;
         j++) {
      if (stash_full[j]) {
        uint32_t mapped_to_leaf = stash[j].header.leaf_label;
        if (is_on_path(index_list[i], mapped_to_leaf)) {
          std::memcpy(&data[(i * blocks_per_bucket + cur_block_in_bucket) *
                            block_size_bytes],
                      &stash[j], block_size_bytes);
          stash_full[j] = false;
          cur_block_in_bucket++;
        }
      }
    }
    while (cur_block_in_bucket < blocks_per_bucket) {
      ORAMBlock dummy;
      std::memcpy(&data[(i * blocks_per_bucket + cur_block_in_bucket) *
                        block_size_bytes],
                  &dummy, block_size_bytes);
      cur_block_in_bucket++;
    }
  }
  server->write_buckets(index_list, data);
}

uint32_t Client::next_available_block_id() { return cur_block_id++; }