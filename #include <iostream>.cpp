#include <iostream>
#include <vector>
#include <cmath>
#include <string>

// Konstantní hodnota PI
const float PI = 3.14159f;

// --- Třída s šablonovým (template) poloměrem pro BONUS ---
template <typename T>
class Kruh {
private:
    T radius;

public:
    // 1. Parametrický Konstruktor (pro Kruh druhyKruh(druhyRaidus);)
    Kruh(T r) : radius(r) {
        // Inicializační seznam ( : radius(r) ) je moderní a doporučený způsob
    }

    // 2. Výchozí Konstruktor (pro Kruh prvniKruh;)
    Kruh() : radius(T()) {
        // Inicializuje radius na výchozí hodnotu (0 pro numerické typy)
    }

    // Setter
    void setRadius(T r) {
        if (r >= 0) {
            radius = r;
        } else {
            // Jednoduchá kontrola, aby poloměr nebyl záporný
            std::cerr << "Chyba: Poloměr musí být nezáporné číslo." << std::endl;
        }
    }

    // Metoda pro výpočet plochy
    float spocitejPlochu() const {
        return PI * radius * radius;
    }

    // Metoda pro výpočet obvodu
    float spocitejObvod() const {
        return 2 * PI * radius;
    }

    // Metoda pro výpis všech údajů
    void vypisUdaje() const {
        std::cout << "\n--- Údaje o Kruhu ---" << std::endl;
        std::cout << "Poloměr (Radius): " << radius << std::endl;
        std::cout << "Plocha: " << spocitejPlochu() << std::endl;
        std::cout << "Obvod: " << spocitejObvod() << std::endl;
        std::cout << "----------------------" << std::endl;
    }
};

// =======================================================
// HLAVNÍ FUNKCE MAIN
// =======================================================

int main() {
    // Použijeme Kruh s typem poloměru 'float'
    using KruhFloat = Kruh<float>; 
    
    // --- Práce s prvním objektem (podle původní kostry) ---

    float prvniRadius = 2.0f; // Přidáno .0f pro float literál
    KruhFloat prvniKruh; // Volá se Výchozí Konstruktor
    prvniKruh.setRadius(prvniRadius); // Volá se Setter

    float druhyRaidus = 4.0f; 
    KruhFloat druhyKruh(druhyRaidus); // Volá se Parametrický Konstruktor
    
    // Původní výstupy (opravena chyba v zápisu std::stdl na std::endl)
    std::cout << "Plocha kruhu (první): " << prvniKruh.spocitejPlochu() << std::endl;
    std::cout << "Obvod kruhu (první): " << prvniKruh.spocitejObvod() << std::endl;

    druhyKruh.vypisUdaje(); // Vypíše všechny údaje o druhém kruhu

    // =======================================================
    // BONUS: Pole kruhů se šablonovým poloměrem
    // =======================================================
    
    std::cout << "\n========== BONUS - Pole Kruhů ==========" << std::endl;

    // Vytvoření vektoru (pole) kruhů s poloměrem typu 'double'
    using KruhDouble = Kruh<double>;
    std::vector<KruhDouble> poleKruhu;
    
    // Přidání kruhů s různými poloměry
    poleKruhu.emplace_back(1.5); // Kruh s r=1.5
    poleKruhu.emplace_back(3.0); // Kruh s r=3.0
    poleKruhu.emplace_back(5.2); // Kruh s r=5.2

    int pocitadlo = 1;
    for (const auto& k : poleKruhu) {
        std::cout << "\n--- Kruh č. " << pocitadlo++ << " ---" << std::endl;
        // Zde voláme výpočet pro každý kruh a vypisujeme ho
        std::cout << "  Plocha: " << k.spocitejPlochu() << std::endl;
        std::cout << "  Obvod:  " << k.spocitejObvod() << std::endl;
    }
    
    std::cout << "\n========================================" << std::endl;

    return 0;
}