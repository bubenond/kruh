// Tento soubor je určen pro C++ kód.
// Můžete začít psát svůj C++ program zde.

#include <iostream>

int main() {
    const int staticSize = 5;
    int staticArray[staticSize] = {1, 2, 3, 4, 5};

std::cout << "Static Array Example" << std::endl;
for (int i = 0; i < staticSize; ++i) {
    std::cout << staticArray[i] << " ";
}
std::cout << std::endl;
return 0;
}
