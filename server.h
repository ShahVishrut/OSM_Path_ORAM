#ifndef SERVER_H
#define SERVER_H

#include <vector>
#include <cstdint>

class Server {
public:
    Server(size_t bucket_size_bytes, size_t num_buckets);
    bool write_buckets(std::vector<size_t> index_list, std::vector<uint8_t> data);
    std::vector<uint8_t> read_buckets(std::vector<size_t> index_list);

private:
    std::vector<uint8_t> database;
    size_t bucket_size_bytes;
    size_t num_buckets;
};


#endif