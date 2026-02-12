#include "client.h"
#include <iostream>

int main() {
    Client *my_client = new Client(200, 4);
    my_client->write_data(30, 67);
    my_client->write_data(10, 78);
    my_client->write_data(20, 88);
    my_client->write_data(40, 98);
}