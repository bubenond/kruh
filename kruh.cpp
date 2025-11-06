// objekty velkym pismenem a typy malym pismenem.


#include <iostream>


class kruh 
{
private:
    float polomer;
    float pi = 3.14159f;

public: 
void nastavPolomer(float r) 
{
    polomer = r;
}

float spocitejObsah() 
{    
    return pi * polomer * polomer;
}

float spocitejObvod() 
{
    return 2 * pi * polomer;
}
};
int main() 
{
    kruh mujprvniKruh;
    mujprvniKruh.nastavPolomer(5.0f);

    std::cout << "Obsah kruhu: " << mujprvniKruh.spocitejObsah() << std::endl;
    std::cout << "Obvod kruhu: " << mujprvniKruh.spocitejObvod() << std::endl;
    
kruh mujDruhyKruh;
mujDruhyKruh.nastavPolomer(10.0f);

std::cout << "Obsah druheho kruhu: " << mujDruhyKruh.spocitejObsah() << std::endl;
std::cout << "Obvod druheho kruhu: " << mujDruhyKruh.spocitejObvod() << std::endl;


    return 0;
}