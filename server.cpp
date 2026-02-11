#include "server.h"
#include <cstring>

Server::Server(size_t bucket_size_bytes, size_t num_buckets) {
    this->bucket_size_bytes = bucket_size_bytes;
    this->num_buckets = num_buckets;
    this->database.resize(bucket_size_bytes * num_buckets);
}

bool Server::write_buckets(std::vector<size_t> index_list, std::vector<uint8_t> data) {
    if (bucket_size_bytes * index_list.size() != data.size()) {
        return false;
    }
    
    for (size_t idx : index_list) {
        if (idx >= num_buckets) return false; 
    }

    for (size_t i = 0; i < index_list.size(); i++) {
        size_t db_offset = index_list[i] * bucket_size_bytes;
        size_t data_offset = i * bucket_size_bytes;

        std::memcpy(&database[db_offset], &data[data_offset], bucket_size_bytes);
    }
    return true;
}

std::vector<uint8_t> Server::read_buckets(std::vector<size_t> index_list) {
    for (size_t idx : index_list) {
        if (idx >= num_buckets) return {};
    }

    std::vector<uint8_t> data(index_list.size() * bucket_size_bytes);

    for (size_t i = 0; i < index_list.size(); i++) {
        size_t db_offset = index_list[i] * bucket_size_bytes;
        size_t data_offset = i * bucket_size_bytes;
        
        std::memcpy(&data[data_offset], &database[db_offset], bucket_size_bytes);
    }
    return data;
}