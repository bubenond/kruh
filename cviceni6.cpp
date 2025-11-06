#include <iostream>
#include <vector>
#include <cmath>
#include <string>

const float pi=3.14f;

class Kruh
{
    public: 
 void setRadius ( float novyRadius)
 {
   radius = novyRadius;
 }

 float getRadius()
 {
    return radius; 
  }

  void print()
{
std::cout << "radius je: " << radius << std::endl;    
}

float spocitejPlochu() const{
return pi*radius*radius ;
}


float spocitejObvod() const{
return 2*pi*radius ;
}

    private:

      float radius;
};


int main()
{
   float prvniRadius = 2;
   Kruh prvniKruh;
   prvniKruh.setRadius( prvniRadius );

 //  float druhyRaidus = 4;
 //  Kruh druhyKruh( druhyRaidus );
   std::cout << "Plocha kruhu:" << prvniKruh.spocitejPlochu() << std::endl;
   std::cout << "Obvod kruhu:" << prvniKruh.spocitejObvod() << std::endl;

  // druhyKruh.vypisUdaje();

   // BONUS: radius kruhu bude sablonovy parametr, a udela pole kruhu a pro kazdy kruh spocita vse

   return 0;
}