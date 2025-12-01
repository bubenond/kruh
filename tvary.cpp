#include <vector>
#include <string>
#include <iostream>

class Tvar {
public:
Tvar() {};
  virtual ~Tvar() {};
  virtual float obvod() = 0;
  virtual float obsah() = 0;
  virtual std::string jmeno() = 0;

};

// class Kruh odvozena z tvaru
class Kruh : public Tvar { 
private:
    float polomer;
    const float pi = 3.14159f;

public:
    Kruh(float r) : polomer(r) {}
    float obvod() override {
        return 2 * pi * polomer;
    }

    float obsah() override {
        return pi * polomer * polomer;
    }

    std::string jmeno() override {
        return "Kruh";
    }
};
// class Ctverec odvozena z tvaru
class Ctverec : public Tvar {
private:
    float strana;
public:
    Ctverec(float a) : strana(a) {}
    float obvod() override {
        return 4 * strana;
    }
    float obsah() override {
        return strana * strana;
    }
    std::string jmeno() override {
        return "Ctverec";
    }
};


int main()
{
    std::vector<Tvar*> tvary;
    tvary.push_back(new Kruh(5.0f));
    tvary.push_back(new Ctverec(4.0f));
    tvary.push_back(new Kruh(10.0f));

    for (Tvar* tvar : tvary) {
        std::cout << tvar->jmeno() << " - Obvod: " << tvar->obvod() << ", Obsah: " << tvar->obsah() << std::endl;
    }

    for (Tvar* tvar : tvary) {
        delete tvar;
    }

   return 0;
};
