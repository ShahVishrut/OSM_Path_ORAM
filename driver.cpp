#include "client.h"
#include <iostream>

int main() {
    Client *my_client = new Client(200, 10, 4);
    my_client->init_position_map();
    std::vector<uint8_t> data(2);
    data[0] = 23;
    data[1] = 255;
    my_client->write_block(20, data);
    std::cout << (my_client->read_block(20))[0] + 0;
    std::cout << (my_client->read_block(20))[1] + 0;
}