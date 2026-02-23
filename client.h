#ifndef CLIENT_H
#define CLIENT_H

#include "server.h"
#include <vector>
#include <cstdint>
#include <stack>

struct ODSPointer {
    uint32_t block_id; 
    uint32_t leaf_label;
    bool is_null;

    ODSPointer() : block_id(0), leaf_label(0), is_null(true) {}
    ODSPointer(uint32_t b, uint32_t l) : block_id(b), leaf_label(l), is_null(false) {}
};

struct OSMNode {
    uint64_t key;
    uint64_t value;

    ODSPointer l_child_ptr;
    ODSPointer r_child_ptr;

    uint8_t l_height;
    uint8_t r_height;
    uint32_t l_same_key_size;
    uint32_t r_same_key_size;

    uint64_t l_max_key_subtree;
    uint32_t l_max_key_count;

    uint64_t l_min_key_subtree;
    uint32_t l_min_key_count;

    uint64_t r_max_key_subtree;
    uint32_t r_max_key_count;

    uint64_t r_min_key_subtree;
    uint32_t r_min_key_count;

    OSMNode() : key(0), value(0), l_height(0), r_height(0), l_same_key_size(0), r_same_key_size(0), l_max_key_subtree(0), l_max_key_count(0), l_min_key_subtree(0), l_min_key_count(0), r_max_key_subtree(0), r_max_key_count(0), r_min_key_count(0), r_min_key_subtree(0) {}
};

struct ORAMBlock {
    ODSPointer header;
    OSMNode data;
};

class Client {
public:
    Client(size_t num_blocks, size_t blocks_per_bucket);
    void insert(uint64_t key, uint64_t value);
    size_t size(uint64_t key);

private:
    ORAMBlock write_block(ORAMBlock to_write, bool write);
    size_t block_size_bytes = sizeof(ORAMBlock);
    ODSPointer root;
    size_t num_blocks;
    size_t blocks_per_bucket;
    size_t bucket_size_bytes;
    size_t tree_height;
    std::vector<ORAMBlock> stash;
    std::vector<bool> stash_full;
    Server* server;
    void evict(uint32_t leaf_num);
    bool is_on_path(size_t index, size_t leaf_num);
    void read_block_path(uint32_t leaf_num);
    uint32_t next_available_block_id();
    uint32_t cur_block_id = 0;
    void rotate_right_right(std::vector<ORAMBlock>& avl_history, int cur_node_index);
    void rotate_left_left(std::vector<ORAMBlock>& avl_history, int cur_node_index);
    void rotate_right_left(std::vector<ORAMBlock>& avl_history, int cur_node_index);
    void rotate_left_right(std::vector<ORAMBlock>& avl_history, int cur_node_index);
};

#endif