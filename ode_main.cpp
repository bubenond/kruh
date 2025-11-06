#include <iostream>
#include <vector>
#include <cmath>

// Solve dy/dt = -k*y using forward Euler method

int main() {
    double k = 1.0;          // decay rate
    double y0 = 1.0;         // initial condition
    double t0 = 0.0;
    double tmax = 5.0;
    double dt = 0.1;

    int nsteps = static_cast<int>((tmax - t0) / dt);
    double y = y0;
    double t = t0;

    std::cout << "t y_euler y_exact abs_error\n";
    for (int i = 0; i <= nsteps; ++i) {
        double y_exact = y0 * std::exp(-k * t);
        double err = std::abs(y - y_exact);
        std::cout << t << " " << y << " " << y_exact << " " << err << "\n";

        // Euler step
        y = y + dt * (-k * y);
        t += dt;
    }

    return 0;
}
