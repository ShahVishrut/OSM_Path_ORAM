#include "client.h"
#include <cmath>
#include <random>
#include <stdexcept>

#include <vector>
#include <cstring>
#include <cstdint>
#include <array>


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
}


void Client::write_data(uint64_t key, uint64_t value) {
    if (root.is_null) {
        ORAMBlock to_write;
        to_write.header.block_id = next_available_block_id();
        to_write.header.is_null = false;
        to_write.data.key = key;
        to_write.data.value = value;
        root = write_block(to_write, true).header;
        return;
    }

    ORAMBlock cur_read;
    cur_read.header = root;
    std::vector<char> seq;
    
    // 1. Read Path and Append New Node
    while (true) {
        cur_read = write_block(cur_read, false);
        avl_history.push_back(cur_read);
        
        if (key > cur_read.data.key) {
            seq.push_back('R');
            if (cur_read.data.r_child_ptr.is_null) {
                ORAMBlock to_write;
                to_write.header.block_id = next_available_block_id();
                to_write.header.is_null = false;
                to_write.data.key = key;
                to_write.data.value = value;
                
                // Update parent pointer in history
                avl_history.back().data.r_child_ptr = to_write.header;
                avl_history.push_back(to_write);
                break;
            } else {
                cur_read.header = cur_read.data.r_child_ptr;
            }
        } else {
            seq.push_back('L');
            if (cur_read.data.l_child_ptr.is_null) {
                ORAMBlock to_write;
                to_write.header.block_id = next_available_block_id();
                to_write.header.is_null = false;
                to_write.data.key = key;
                to_write.data.value = value;
                
                // Update parent pointer in history
                avl_history.back().data.l_child_ptr = to_write.header;
                avl_history.push_back(to_write);
                break;
            } else {
                cur_read.header = cur_read.data.l_child_ptr;
            }
        }
    }

    // 2. Traverse Upwards: Update Heights and Rebalance
    for (int height = 1; height < avl_history.size(); height++) {
        int index_A = avl_history.size() - 1 - height; // Parent (Unbalanced Candidate)
        int index_B = avl_history.size() - height;     // Child
        
        // Update Height for Node A
        if (seq[seq.size() - height] == 'R') {
            avl_history[index_A].data.r_height = 1 + std::max(avl_history[index_B].data.l_height, avl_history[index_B].data.r_height);
        } else {
            avl_history[index_A].data.l_height = 1 + std::max(avl_history[index_B].data.l_height, avl_history[index_B].data.r_height);
        }

        // Check Balance Factor
        int balance_factor = avl_history[index_A].data.r_height - avl_history[index_A].data.l_height;

        if (std::abs(balance_factor) > 1) {
            // Check for RR Case: Right heavy parent, Right heavy (or neutral) child
            // Note: We check 'seq' to see if the insertion happened on the Right side
            bool child_is_right = (seq[seq.size() - height] == 'R');
            bool grandchild_is_right = (seq.size() - height + 1 < seq.size()) && (seq[seq.size() - height + 1] == 'R');

            if (child_is_right && grandchild_is_right) {
                
                // --- Step 1: Pointer & Height Updates (Standard AVL) ---
                
                // A adopts B's left child
                avl_history[index_A].data.r_child_ptr = avl_history[index_B].data.l_child_ptr;
                // A's right height becomes B's left height
                avl_history[index_A].data.r_height = avl_history[index_B].data.l_height;

                // B adopts A (A becomes Left child of B)
                avl_history[index_B].data.l_child_ptr = avl_history[index_A].header;
                // Recalculate B's left height (based on A's new height)
                avl_history[index_B].data.l_height = 1 + std::max(
                    avl_history[index_A].data.l_height, 
                    avl_history[index_A].data.r_height
                );

                // --- Step 2: Swap in History Vector (For Write-Back Order) ---
                // We swap so the new parent (B) is 'above' the new child (A) in the vector
                ORAMBlock swap = avl_history[index_B];
                avl_history[index_B] = avl_history[index_A];
                avl_history[index_A] = swap;
                
                // NOTE: After swap:
                // index_A contains Node B (The new Subtree Root)
                // index_B contains Node A (The new Child)

                // --- Step 3: Update Grandparent Pointer ---
                if (index_A == 0) {
                    // A was the root. Now B (at index_A) is the root.
                    this->root = avl_history[index_A].header;
                } else {
                    int index_GP = index_A - 1;
                    // Check if Grandparent pointed Right or Left to get to A
                    if (seq[index_GP] == 'R') { // Use index_GP for the seq check
                        avl_history[index_GP].data.r_child_ptr = avl_history[index_A].header;
                    } else {
                        avl_history[index_GP].data.l_child_ptr = avl_history[index_A].header;
                    }
                }
            }
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
            std::memcpy(&temp_block, &temp[(cur_bucket * blocks_per_bucket + cur_block) * block_size_bytes], sizeof(ORAMBlock));
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
        for (int j = 0; j < stash_full.size() && cur_block_in_bucket < blocks_per_bucket; j++) {
            if (stash_full[j]) {
                uint32_t mapped_to_leaf = stash[j].header.leaf_label;
                if (is_on_path(index_list[i], mapped_to_leaf)) {
                    std::memcpy(&data[(i * blocks_per_bucket + cur_block_in_bucket) * block_size_bytes], &stash[j], block_size_bytes);
                    stash_full[j] = false;
                    cur_block_in_bucket++;
                }
            }
        }
        while (cur_block_in_bucket < blocks_per_bucket) {
            ORAMBlock dummy;
            std::memcpy(&data[(i * blocks_per_bucket + cur_block_in_bucket) * block_size_bytes], &dummy, block_size_bytes);
            cur_block_in_bucket++;
        }
    }
    server->write_buckets(index_list, data);
}

uint32_t Client::next_available_block_id() {
    return cur_block_id++;
}