#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm> // Pro std::max

// Třída pro řešení 3D STACIONÁRNÍ rovnice vedení tepla (Poissonova rovnice)
class ThermalSolver3D {
private:
    // **FYZIKÁLNÍ PARAMETRY**
    double k = 50.0;             // Tepelná vodivost [W/(m*K)]
    double q_dot = 200000.0;     // Hustota objemového zdroje tepla [W/m^3] (Jouleovo teplo)
    double h = 10.0;             // Součinitel přestupu tepla konvekcí [W/(m^2*K)]
    
    // **OKRAJOVÉ & POČÁTEČNÍ PODMÍNKY**
    double T_cooler = 300.0;     // Teplota chladiče (Dirichletova na x=0) [K]
    double T_infinity = 293.0;   // Teplota okolí (Robinova) [K]
    
    // **GEOMETRIE A DISKRETIZACE PROSTORU**
    int Nx, Ny, Nz;              // Počet diskretizačních bodů (uzlů)
    double Lx, Ly, Lz;           // Rozměry domény [m]
    double dx;                   // Velikost kroku sítě (předpokládáme dx = dy = dz) [m]

    // **NUMERICKÉ PARAMETRY**
    double tolerance = 1e-5;     // Tolerance konvergence (max. rozdíl T_new - T_old) [K]
    int max_iter = 50000;        // Maximální počet iterací

    // 3D pole pro teploty T(i, j, k)
    std::vector<std::vector<std::vector<double>>> T;

public:
    ThermalSolver3D(double Lx_in, double Ly_in, double Lz_in, int N_in) :
        Lx(Lx_in), Ly(Ly_in), Lz(Lz_in), Nx(N_in), Ny(N_in), Nz(N_in) 
    {
        // Nastavení kroku sítě
        dx = Lx / (Nx - 1); 

        // Inicializace 3D pole teplot s počáteční hodnotou T_infinity (první odhad)
        T.resize(Nx, std::vector<std::vector<double>>(Ny, std::vector<double>(Nz, T_infinity)));
        
        // Aplikace Dirichletovy BC na x=0 hned na začátku
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
                T[0][j][k] = T_cooler;
            }
        }
    }

    /**
     * @brief Aktualizuje teploty na okrajích na základě Robinovy BC (Konvekce).
     */
    void update_boundaries() {
        double C_rob = (h * dx / k); // Biotovo číslo na okrajích (pro zjednodušenou FDM)
        double q_term_edge = q_dot * dx * dx / k; // Člen zdroje na hranici
        
        // Okraj i = Nx-1 (x=Lx) - Robinova BC
        for (int j = 1; j < Ny - 1; ++j) {
            for (int k = 1; k < Nz - 1; ++k) {
                // FDM schéma pro okraj s Robinovou BC, kde je vliv sousedů a zdroje
                T[Nx-1][j][k] = (
                    2.0 * T[Nx-2][j][k] + // Vliv vnitřního souseda (z Robinovy aproximace)
                    T[Nx-1][j+1][k] + T[Nx-1][j-1][k] + // Vliv sousedů v y-směru
                    T[Nx-1][j][k+1] + T[Nx-1][j][k-1] + // Vliv sousedů v z-směru
                    2.0 * C_rob * T_infinity + // Vliv konvekce do okolí
                    q_term_edge // Vliv Jouleova tepla
                ) / (6.0 + 2.0 * C_rob);
            }
        }
        
        // Pro kompletní řešení by se přidala analogická schémata pro okraje j=0, j=Ny-1, k=0, k=Nz-1.
    }

    /**
     * @brief Řeší soustavu rovnic Gauss-Seidelovou metodou, dokud nekonverguje.
     */
    void solve() {
        double error = tolerance + 1.0;
        int iter_count = 0;

        double h_sq = dx * dx; 
        double q_term = q_dot * h_sq / k; // Člen zdroje pro vnitřní uzel
        
        std::cout << "Zahajuji řešení STACIONÁRNÍHO stavu (FDM + Gauss-Seidel)..." << std::endl;
        std::cout << "Síť: " << Nx << "x" << Ny << "x" << Nz << " uzlů." << std::endl;

        // Vlastní iterační proces
        while (error > tolerance && iter_count < max_iter) {
            error = 0.0;
            iter_count++;

            // 1. Iterace přes Vnitřní uzly
            for (int i = 1; i < Nx - 1; ++i) {
                for (int j = 1; j < Ny - 1; ++j) {
                    for (int k = 1; k < Nz - 1; ++k) {
                        double T_old = T[i][j][k];

                        // Diskrétní rovnice pro vnitřní bod: T_new = 1/6 * (součet 6 sousedů) + q_term/6
                        double T_new = (
                            T[i + 1][j][k] + T[i - 1][j][k] +
                            T[i][j + 1][k] + T[i][j - 1][k] +
                            T[i][j][k + 1] + T[i][j][k - 1] + q_term
                        ) / 6.0;

                        error = std::max(error, std::abs(T_new - T_old));

                        // Aktualizace teploty (Gauss-Seidel)
                        T[i][j][k] = T_new;
                    }
                }
            }
            
            // 2. Aktualizace okrajů (BCs)
            update_boundaries();

            if (iter_count % 5000 == 0) {
                std::cout << "Iterace: " << iter_count << ", Max. chyba: " << error << std::endl;
            }
        }

        std::cout << "\n--------------------------------------------------" << std::endl;
        if (error <= tolerance) {
            std::cout << "Řešení konvergovalo po " << iter_count << " iteracích k ustálenému stavu." << std::endl;
        } else {
            std::cout << "Řešení NEKONVERGOVALO po " << max_iter << " iteracích." << std::endl;
        }
        std::cout << "--------------------------------------------------" << std::endl;
    }

    /**
     * @brief Vypíše výslednou teplotu v centrálním průřezu.
     */
    void print_results() const {
        std::cout << "\nVýsledné STACIONÁRNÍ teploty (průřez j = Ny/2 a k = Nz/2):" << std::endl;
        std::cout << std::fixed << std::setprecision(3);
        
        int j_center = Ny / 2;
        int k_center = Nz / 2;

        std::cout << "x/dx | T(x) [K]" << std::endl;
        std::cout << "------|------------" << std::endl;

        double T_max = T_cooler;
        for (int i = 0; i < Nx; ++i) {
            double T_val = T[i][j_center][k_center];
            std::cout << std::setw(4) << i << " | " << T_val << std::endl;
            T_max = std::max(T_max, T_val);
        }
        
        std::cout << "\nMaximální dosažená teplota v ustáleném stavu: " << T_max << " K" << std::endl;
    }
};

int main() {
    // Rozměry: 0.1m x 0.1m x 0.1m (10cm krychle)
    double Lx = 0.1, Ly = 0.1, Lz = 0.1; 
    // Síť 21x21x21 uzlů (Nx=Ny=Nz=21)
    int N = 21; 

    // Inicializace a spuštění řešiče
    ThermalSolver3D solver(Lx, Ly, Lz, N);
    solver.solve();

    // Vypis ustáleného stavu
    solver.print_results();

    return 0;
}