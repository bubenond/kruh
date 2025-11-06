#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

int main() {
    // -----------------------------
    // Geometrie a síť
    // -----------------------------
    const double Lx = 0.1, Ly = 0.06, Lz = 0.01; 
    const int Nx = 200, Ny = 100, Nz = 50;
    const double dx = Lx / (Nx - 1);
    const double dy = Ly / (Ny - 1);
    const double dz = Lz / (Nz - 1);

    // -----------------------------
    // Fyzikální parametry
    // -----------------------------
    const double lambda = 0.4;     // [W/mK]
    const double rho = 1200.0;     // [kg/m3]
    const double cp = 1000.0;      // [J/kgK]
    const double Tinf = 293.15;    // [K]
    const double Q0 = 5e4;         // [W/m3]

    // Rozdílné konvekční koeficienty podle směru
    const double h_x = 2.0;   // např. chlazení na bocích
    const double h_y = 1.0;   // slabší chlazení
    const double h_z = 3.0;   // silnější chlazení (např. kontakt s deskou)

    // -----------------------------
    // Pomocná funkce pro index
    // -----------------------------
    auto idx = [Nx, Ny](int i, int j, int k) { return i + Nx * (j + Ny * k); };

    // -----------------------------
    // Inicializace teploty a zdroje
    // -----------------------------
    std::vector<double> T(Nx * Ny * Nz, Tinf);
    std::vector<double> Q(Nx * Ny * Nz, 0.0);

    // Zdroj uprostřed tělesa (1/3 objemu)
    for (int k = Nz / 3; k < 2 * Nz / 3; k++)
        for (int j = Ny / 3; j < 2 * Ny / 3; j++)
            for (int i = Nx / 3; i < 2 * Nx / 3; i++)
                Q[idx(i, j, k)] = Q0;

    // -----------------------------
    // Iterační výpočet (Gauss–Seidel)
    // -----------------------------
    const int maxIter = 10000;
    const double tol = 1e-6;

    for (int iter = 0; iter < maxIter; iter++) {
        double maxChange = 0.0;

        for (int k = 1; k < Nz - 1; k++)
            for (int j = 1; j < Ny - 1; j++)
                for (int i = 1; i < Nx - 1; i++) {
                    int p = idx(i, j, k);

                    double Tnew =
                        (T[idx(i + 1, j, k)] + T[idx(i - 1, j, k)]) / (dx * dx) +
                        (T[idx(i, j + 1, k)] + T[idx(i, j - 1, k)]) / (dy * dy) +
                        (T[idx(i, j, k + 1)] + T[idx(i, j, k - 1)]) / (dz * dz) +
                        Q[p] / lambda;

                    double denom = 2.0 / (dx * dx) + 2.0 / (dy * dy) + 2.0 / (dz * dz);
                    Tnew /= denom;

                    double diff = std::abs(Tnew - T[p]);
                    if (diff > maxChange) maxChange = diff;

                    T[p] = Tnew;
                }

        // -----------------------------
        // Robinovy okraje – různé h
        // -----------------------------
        // Směr X
        for (int j = 0; j < Ny; j++)
            for (int k = 0; k < Nz; k++) {
                T[idx(0, j, k)] =
                    (lambda / dx * T[idx(1, j, k)] + h_x * Tinf) / (lambda / dx + h_x);
                T[idx(Nx - 1, j, k)] =
                    (lambda / dx * T[idx(Nx - 2, j, k)] + h_x * Tinf) / (lambda / dx + h_x);
            }

        // Směr Y
        for (int i = 0; i < Nx; i++)
            for (int k = 0; k < Nz; k++) {
                T[idx(i, 0, k)] =
                    (lambda / dy * T[idx(i, 1, k)] + h_y * Tinf) / (lambda / dy + h_y);
                T[idx(i, Ny - 1, k)] =
                    (lambda / dy * T[idx(i, Ny - 2, k)] + h_y * Tinf) / (lambda / dy + h_y);
            }

        // Směr Z
        for (int i = 0; i < Nx; i++)
            for (int j = 0; j < Ny; j++) {
                T[idx(i, j, 0)] =
                    (lambda / dz * T[idx(i, j, 1)] + h_z * Tinf) / (lambda / dz + h_z);
                T[idx(i, j, Nz - 1)] =
                    (lambda / dz * T[idx(i, j, Nz - 2)] + h_z * Tinf) / (lambda / dz + h_z);
            }

        if (iter % 200 == 0)
            std::cout << "Iterace " << iter << ", maxChange=" << maxChange << std::endl;

        if (maxChange < tol) {
            std::cout << "Konvergence po " << iter << " iteracích" << std::endl;
            break;
        }
    }

    // -----------------------------
    // Min a max teplota
    // -----------------------------
    double Tmin = T[0], Tmax = T[0];
    for (double val : T) {
        if (val < Tmin) Tmin = val;
        if (val > Tmax) Tmax = val;
    }
    std::cout << "Min T = " << Tmin << " K, Max T = " << Tmax << " K" << std::endl;

    // -----------------------------
    // Uložení CSV
    // -----------------------------
    std::ofstream fxy("plane_xy.csv"), fxz("plane_xz.csv"), fyz("plane_yz.csv");

    fxy << "x,y,T\n";
    for (int j = 0; j < Ny; j++)
        for (int i = 0; i < Nx; i++)
            fxy << i * dx << "," << j * dy << "," << T[idx(i, j, Nz / 2)] << "\n";

    fxz << "x,z,T\n";
    for (int k = 0; k < Nz; k++)
        for (int i = 0; i < Nx; i++)
            fxz << i * dx << "," << k * dz << "," << T[idx(i, Ny / 2, k)] << "\n";

    fyz << "y,z,T\n";
    for (int k = 0; k < Nz; k++)
        for (int j = 0; j < Ny; j++)
            fyz << j * dy << "," << k * dz << "," << T[idx(Nx / 2, j, k)] << "\n";

    std::cout << "Výsledky uloženy do CSV." << std::endl;
    return 0;
}
