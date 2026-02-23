#include "client.h"
#include <iostream>

int main() {
    Client *my_client = new Client(200, 4);
    my_client->insert(50, 66);
    my_client->insert(10, 78);
    my_client->insert(30, 88);
    my_client->insert(75, 987);
    my_client->insert(40, 984);
    my_client->insert(45, 982);
    my_client->insert(5, 981);
    my_client->insert(50, 67);

    std::cout << my_client->size(10);
}