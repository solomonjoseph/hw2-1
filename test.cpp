// #include <stdio.h>
// #include <ostream>
#include <iostream>

int main(int argc, char* argv[]) {
    std::cout << "program" << std::endl;
    for (int i=0; i<1000000; i++) {
        if (999983 % i == 0) {
            std::cout << i << std::endl;
        }
    }
    return 0;
}
