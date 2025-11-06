#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

int main()
{
    // --- Geometrie a fyzika
    const double x1 = 0.0;
    const double x2 = 1.0;
    const int    nx = 201;                       // počet uzlů (=> dx = (x2-x1)/(nx-1))
    const double c = 1.0;                       // rychlost advekce (>0 znamená proudění zleva doprava)

    // --- Diskretizace
    const double dx = (x2 - x1) / (nx - 1);
    const double CFL_target = 0.9;               // cílové CFL (<=1 pro upwind stabilní)
    double dt = CFL_target * dx / c;             // dt z CFL
    const double t_end = 0.5;                    // konec simulace (s)
    int nt = static_cast<int>(std::ceil(t_end / dt));
    dt = t_end / nt;                             // dorovnáme dt, aby poslední krok dopadl přesně v t_end
    const double lambda = c * dt / dx;           // skutečné CFL

    // --- Počáteční podmínka: Gauss kolem x=0.3
    const double x_center = 0.30;
    const double sigma = 0.05;

    // --- Pole řešení (staré a nové časové patro)
    std::vector<double> u(nx, 0.0), u_new(nx, 0.0), u0(nx, 0.0);

    // Naplnění počáteční podmínky
    for (int i = 0; i < nx; ++i)
    {
        const double x = x1 + i * dx;
        u[i] = std::exp(-(x - x_center) * (x - x_center) / (2.0 * sigma * sigma));
        u0[i] = u[i];
    }

    // --- Časová smyčka (explicitní upwind pro c > 0):
    // u_j^{n+1} = u_j^n - lambda * (u_j^n - u_{j-1}^n),  j = 1..nx-1
    for (int n = 0; n < nt; ++n)
    {
        // Levá okrajová podmínka: u(t, x=0) = 0
        u_new[0] = 0.0;

        for (int j = 1; j < nx; ++j)
        {
            u_new[j] = u[j] - lambda * (u[j] - u[j - 1]);
        }

        // Prohození vrstev
        std::swap(u, u_new);

        // (volitelné) odkomentuj, pokud chceš vidět čísla kroků
        // std::cout << "n = " << (n+1) << " / " << nt << "\r" << std::flush;
    }

    // --- Uložení výsledku (x, u0, u_final) do CSV
    std::ofstream fout("advekce_solution.csv");
    fout << std::setprecision(10);
    fout << "x,u0,u_final\n";
    for (int i = 0; i < nx; ++i)
    {
        const double x = x1 + i * dx;
        fout << x << "," << u0[i] << "," << u[i] << "\n";
    }
    fout.close();

    // Krátká rekapitulace do konzole
    std::cout << "Hotovo.\n"
        << "  dx = " << dx << ", dt = " << dt << ", CFL = " << lambda << "\n"
        << "  nt = " << nt << " kroku, vysledek ulozen do advekce_solution.csv\n";

    return 0;
}