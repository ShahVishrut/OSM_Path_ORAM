#include "client.h"
#include <iostream>

int main() {
    Client *my_client = new Client(200, 4);
    my_client->write_data(50, 67);
    my_client->write_data(10, 78);
    my_client->write_data(30, 88);
    my_client->write_data(75, 987);
    my_client->write_data(40, 984);
    my_client->write_data(45, 982);
    my_client->write_data(5, 981);
    my_client->write_data(42, 9889);
}