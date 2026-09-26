
#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

constexpr double PI_GLOBAL = 3.14159265358979323846;

struct Candidate {
    std::string label;
    double target_phase_deg;
    int width_cells;
    int height_cells;
};

struct CPML1D {
    std::vector<double> kappa_e, b_e, c_e;
    std::vector<double> kappa_h, b_h, c_h;
};

struct LinePhasor {
    std::vector<std::complex<double>> E;
    std::vector<std::complex<double>> H;
};

struct Checkpoint {
    LinePhasor ref;
    LinePhasor trans_far;
};

struct RunOutput {
    std::vector<Checkpoint> checkpoints;
};

inline int id(int x, int y, int Nx) {
    return y * Nx + x;
}

double cpml_rho_e(int y, int Ny, int npml) {
    if (y < npml) return static_cast<double>(npml - y) / npml;
    if (y >= Ny - npml) return static_cast<double>(y - (Ny - npml - 1)) / npml;
    return 0.0;
}

double cpml_rho_h(int y, int Ny, int npml) {
    const double yh = y + 0.5;
    if (yh < npml) return (npml - yh) / npml;
    if (yh >= Ny - npml) return (yh - (Ny - npml)) / npml;
    return 0.0;
}

CPML1D build_cpml(
    int Ny,
    int npml,
    double dx,
    double dt,
    int m = 3,
    double kappa_max = 6.0,
    double target_R = 1.0e-8
) {
    constexpr double c0   = 299792458.0;
    constexpr double mu0  = 4.0e-7 * PI_GLOBAL;
    constexpr double eps0 = 1.0 / (mu0 * c0 * c0);
    const double eta0 = std::sqrt(mu0 / eps0);

    CPML1D p;
    p.kappa_e.assign(Ny, 1.0);
    p.b_e.assign(Ny, 1.0);
    p.c_e.assign(Ny, 0.0);
    p.kappa_h.assign(Ny, 1.0);
    p.b_h.assign(Ny, 1.0);
    p.c_h.assign(Ny, 0.0);

    const double L = npml * dx;
    const double sigma_max =
        -(m + 1.0) * std::log(target_R) / (2.0 * eta0 * L);
    const double alpha_max = 0.05 * sigma_max;

    auto coeff = [&](double rho, bool magnetic, double& kappa, double& b, double& c) {
        if (rho <= 0.0) {
            kappa = 1.0;
            b = 1.0;
            c = 0.0;
            return;
        }

        const double rpow = std::pow(rho, m);
        kappa = 1.0 + (kappa_max - 1.0) * rpow;

        double sigma = sigma_max * rpow;
        double alpha = alpha_max * (1.0 - rho);
        double medium = eps0;

        if (magnetic) {
            const double scale = mu0 / eps0;
            sigma *= scale;
            alpha *= scale;
            medium = mu0;
        }

        b = std::exp(-(sigma / kappa + alpha) * dt / medium);

        const double denom = sigma * kappa + kappa * kappa * alpha;
        c = (std::abs(denom) < 1.0e-30)
            ? 0.0
            : sigma * (b - 1.0) / denom;
    };

    for (int y = 0; y < Ny; ++y) {
        coeff(cpml_rho_e(y, Ny, npml), false, p.kappa_e[y], p.b_e[y], p.c_e[y]);
        coeff(cpml_rho_h(y, Ny, npml), true,  p.kappa_h[y], p.b_h[y], p.c_h[y]);
    }

    return p;
}

LinePhasor zero_line(int Nx) {
    return {
        std::vector<std::complex<double>>(Nx, {0.0, 0.0}),
        std::vector<std::complex<double>>(Nx, {0.0, 0.0})
    };
}

Checkpoint zero_checkpoint(int Nx) {
    return {zero_line(Nx), zero_line(Nx)};
}

void accumulate_line(
    LinePhasor& line,
    const std::vector<double>& Ez,
    const std::vector<double>& Hx,
    int y,
    int Nx,
    const std::complex<double>& kernel
) {
    for (int x = 0; x < Nx; ++x) {
        const double E = Ez[id(x, y, Nx)];
        const double H = 0.5 * (
            Hx[id(x, y - 1, Nx)] + Hx[id(x, y, Nx)]
        );

        line.E[x] += E * kernel;
        line.H[x] += H * kernel;
    }
}

RunOutput run_cell(
    bool with_pillar,
    int Nx,
    int Ny,
    int Nt,
    double dx,
    double dt,
    int source_y,
    int reflection_y,
    int trans_far_y,
    int pillar_x0,
    int pillar_x1,
    int pillar_y0,
    int pillar_y1,
    int npml,
    double lambda0,
    double pillar_n,
    const std::vector<double>& checkpoint_fs
) {
    constexpr double c0   = 299792458.0;
    constexpr double mu0  = 4.0e-7 * PI_GLOBAL;
    constexpr double eps0 = 1.0 / (mu0 * c0 * c0);

    const int N = Nx * Ny;

    std::vector<double> Ez(N, 0.0);
    std::vector<double> Hx(N, 0.0);
    std::vector<double> Hy(N, 0.0);
    std::vector<double> eps_r(N, 1.0);
    std::vector<double> psi_hx_y(N, 0.0);
    std::vector<double> psi_ez_y(N, 0.0);

    if (with_pillar) {
        if (pillar_x0 < 0 || pillar_x1 > Nx ||
            pillar_y0 < 0 || pillar_y1 >= trans_far_y) {
            throw std::runtime_error("Pillar geometry exceeds valid region.");
        }

        for (int y = pillar_y0; y < pillar_y1; ++y) {
            for (int x = pillar_x0; x < pillar_x1; ++x) {
                eps_r[id(x, y, Nx)] = pillar_n * pillar_n;
            }
        }
    }

    const CPML1D pml = build_cpml(Ny, npml, dx, dt);

    const double f0 = c0 / lambda0;
    const double omega0 = 2.0 * PI_GLOBAL * f0;
    const double sigma_t = 4.5e-15;
    const double t0 = 6.0 * sigma_t;

    std::vector<int> checkpoint_steps;
    for (double t_fs : checkpoint_fs) {
        checkpoint_steps.push_back(
            static_cast<int>(std::round(t_fs * 1e-15 / dt))
        );
    }

    Checkpoint accum = zero_checkpoint(Nx);
    RunOutput out;
    out.checkpoints.assign(checkpoint_fs.size(), zero_checkpoint(Nx));
    std::size_t next_cp = 0;

    for (int n = 0; n < Nt; ++n) {
        for (int y = 0; y < Ny - 1; ++y) {
            const double inv_kappa_y = 1.0 / pml.kappa_h[y];

            for (int x = 0; x < Nx; ++x) {
                const int xp = (x + 1) % Nx;
                const int k = id(x, y, Nx);

                const double dEz_dy =
                    (Ez[id(x, y + 1, Nx)] - Ez[k]) / dx;

                psi_hx_y[k] =
                    pml.b_h[y] * psi_hx_y[k]
                    + pml.c_h[y] * dEz_dy;

                Hx[k] -= (dt / mu0) * (
                    inv_kappa_y * dEz_dy + psi_hx_y[k]
                );

                const double dEz_dx =
                    (Ez[id(xp, y, Nx)] - Ez[k]) / dx;

                Hy[k] += (dt / mu0) * dEz_dx;
            }
        }

        for (int y = 1; y < Ny - 1; ++y) {
            const double inv_kappa_y = 1.0 / pml.kappa_e[y];

            for (int x = 0; x < Nx; ++x) {
                const int xm = (x - 1 + Nx) % Nx;
                const int k = id(x, y, Nx);

                const double dHy_dx =
                    (Hy[k] - Hy[id(xm, y, Nx)]) / dx;

                const double dHx_dy =
                    (Hx[k] - Hx[id(x, y - 1, Nx)]) / dx;

                psi_ez_y[k] =
                    pml.b_e[y] * psi_ez_y[k]
                    + pml.c_e[y] * dHx_dy;

                Ez[k] += (dt / (eps0 * eps_r[k])) * (
                    dHy_dx
                    - (inv_kappa_y * dHx_dy + psi_ez_y[k])
                );
            }
        }

        const double time = n * dt;
        const double tau = time - t0;

        const double envelope =
            std::exp(-0.5 * tau * tau / (sigma_t * sigma_t));

        const double source =
            envelope * std::sin(omega0 * tau);

        for (int x = 0; x < Nx; ++x) {
            Ez[id(x, source_y, Nx)] += source;
        }

        for (int x = 0; x < Nx; ++x) {
            Ez[id(x, 0, Nx)] = 0.0;
            Ez[id(x, Ny - 1, Nx)] = 0.0;
        }

        const double phase = -omega0 * time;
        const std::complex<double> kernel(
            std::cos(phase),
            std::sin(phase)
        );

        accumulate_line(
            accum.ref,
            Ez, Hx,
            reflection_y,
            Nx,
            kernel
        );

        accumulate_line(
            accum.trans_far,
            Ez, Hx,
            trans_far_y,
            Nx,
            kernel
        );

        while (
            next_cp < checkpoint_steps.size()
            && n >= checkpoint_steps[next_cp]
        ) {
            out.checkpoints[next_cp] = accum;
            ++next_cp;
        }
    }

    while (next_cp < out.checkpoints.size()) {
        out.checkpoints[next_cp] = accum;
        ++next_cp;
    }

    const std::complex<double> half_step_phase(
        std::cos(omega0 * dt / 2.0),
        std::sin(omega0 * dt / 2.0)
    );

    for (auto& cp : out.checkpoints) {
        for (int x = 0; x < Nx; ++x) {
            cp.ref.H[x] *= half_step_phase;
            cp.trans_far.H[x] *= half_step_phase;
        }
    }

    return out;
}

LinePhasor subtract_lines(
    const LinePhasor& a,
    const LinePhasor& b
) {
    LinePhasor out = zero_line(static_cast<int>(a.E.size()));

    for (std::size_t x = 0; x < a.E.size(); ++x) {
        out.E[x] = a.E[x] - b.E[x];
        out.H[x] = a.H[x] - b.H[x];
    }

    return out;
}

double integrated_power(
    const LinePhasor& line,
    double dx
) {
    double sum = 0.0;

    for (std::size_t x = 0; x < line.E.size(); ++x) {
        sum += 0.5 * std::real(
            line.E[x] * std::conj(line.H[x])
        );
    }

    return sum * dx;
}

std::complex<double> zeroth_order_E(
    const LinePhasor& line
) {
    std::complex<double> sum(0.0, 0.0);

    for (const auto& v : line.E) {
        sum += v;
    }

    return sum / static_cast<double>(line.E.size());
}

double wrapped_phase_deg(
    const std::complex<double>& z
) {
    double p = std::arg(z) * 180.0 / PI_GLOBAL;

    while (p < 0.0) p += 360.0;
    while (p >= 360.0) p -= 360.0;

    return p;
}

double circular_error_deg(
    double a,
    double b
) {
    double d = std::fmod(a - b + 540.0, 360.0) - 180.0;
    return std::abs(d);
}

int main() {
    constexpr double c0 = 299792458.0;

    const double period = 1.0e-6;
    const int Nx = 60;
    const double dx = period / Nx;

    const double S = 0.70;
    const double dt = S * dx / c0;

    const double lambda0 = 1.55e-6;
    const double pillar_n = 2.0;

    const double domain_height_um = 11.0;
    const double cpml_um = 1.2;

    const int Ny = static_cast<int>(
        std::round(domain_height_um * 1e-6 / dx)
    );

    const int npml = static_cast<int>(
        std::round(cpml_um * 1e-6 / dx)
    );

    const int source_y =
        static_cast<int>(std::round(2.2e-6 / dx));

    const int reflection_y =
        static_cast<int>(std::round(2.8e-6 / dx));

    const int pillar_y0 =
        static_cast<int>(std::round(3.6e-6 / dx));

    const int trans_far_y =
        static_cast<int>(std::round(7.5e-6 / dx));

    const std::vector<double> checkpoint_fs = {
        420.0,
        600.0
    };

    const double total_time_fs = 600.0;

    const int Nt = static_cast<int>(
        std::ceil(total_time_fs * 1e-15 / dt)
    ) + 2;

    const std::vector<Candidate> candidates = {
        {"S0", 0.0, 36, 108},
        {"S45", 45.0, 27, 110},
        {"S90", 90.0, 42, 78},
        {"S135", 135.0, 18, 96},
        {"S180", 180.0, 30, 60},
        {"S225", 225.0, 18, 60},
        {"S270", 270.0, 18, 37},
        {"S315", 315.0, 11, 32},
    };

    std::cout << "Step 5L: uniform final verification of all 8 phase states\n";
    std::cout << "--------------------------------------------------------\n";
    std::cout << "candidate count = " << candidates.size() << "\n";
    std::cout << "checkpoints = 420 and 600 fs\n";
    std::cout << "dx = " << dx * 1e9 << " nm\n";

#ifdef _OPENMP
    std::cout << "OpenMP enabled, max threads = "
              << omp_get_max_threads() << "\n";
#endif

    std::cout << "\nRunning 600-fs vacuum reference...\n";

    const auto vac = run_cell(
        false,
        Nx, Ny, Nt,
        dx, dt,
        source_y,
        reflection_y,
        trans_far_y,
        0, 0,
        pillar_y0,
        pillar_y0 + 1,
        npml,
        lambda0,
        pillar_n,
        checkpoint_fs
    );

    std::vector<RunOutput> sims(candidates.size());

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < static_cast<int>(candidates.size()); ++i) {
        const Candidate c = candidates[i];

        const int x0 = (Nx - c.width_cells) / 2;
        const int x1 = x0 + c.width_cells;
        const int y1 = pillar_y0 + c.height_cells;

        sims[i] = run_cell(
            true,
            Nx, Ny, Nt,
            dx, dt,
            source_y,
            reflection_y,
            trans_far_y,
            x0, x1,
            pillar_y0, y1,
            npml,
            lambda0,
            pillar_n,
            checkpoint_fs
        );

#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            std::cout << "Completed "
                      << c.label
                      << " target=" << c.target_phase_deg
                      << " deg\n";
        }
    }

    std::ofstream out("plots/final_8state_verification_step5l.csv");

    out << "checkpoint_fs,label,target_phase_deg,width_cells,height_cells,"
           "width_nm,height_nm,phase_deg,phase_error_deg,R,T,R_plus_T,"
           "energy_error,passes_T,passes_phase,passes_energy,passes_all\n";

    out << std::setprecision(15);

    for (int i = 0; i < static_cast<int>(candidates.size()); ++i) {
        const Candidate c = candidates[i];

        for (int j = 0; j < static_cast<int>(checkpoint_fs.size()); ++j) {
            const auto& v = vac.checkpoints[j];
            const auto& s = sims[i].checkpoints[j];

            const LinePhasor refl =
                subtract_lines(s.ref, v.ref);

            const double S_inc =
                integrated_power(v.ref, dx);

            const double S_ref =
                integrated_power(refl, dx);

            const double S_trans =
                integrated_power(s.trans_far, dx);

            const double S_vac_trans =
                integrated_power(v.trans_far, dx);

            const double eps = 1e-30;

            const double R =
                -S_ref / (S_inc + eps);

            const double T =
                S_trans / (S_vac_trans + eps);

            const auto E0_struct =
                zeroth_order_E(s.trans_far);

            const auto E0_vac =
                zeroth_order_E(v.trans_far);

            const double phase =
                wrapped_phase_deg(
                    E0_struct
                    / (E0_vac + std::complex<double>(eps, 0.0))
                );

            const double phase_error =
                circular_error_deg(
                    phase,
                    c.target_phase_deg
                );

            const double energy_error =
                std::abs(R + T - 1.0);

            const bool passT = T >= 0.80;
            const bool passP = phase_error <= 10.0;
            const bool passE = energy_error <= 0.02;
            const bool passAll = passT && passP && passE;

            out << checkpoint_fs[j] << ","
                << c.label << ","
                << c.target_phase_deg << ","
                << c.width_cells << ","
                << c.height_cells << ","
                << c.width_cells * dx * 1e9 << ","
                << c.height_cells * dx * 1e9 << ","
                << phase << ","
                << phase_error << ","
                << R << ","
                << T << ","
                << R + T << ","
                << energy_error << ","
                << passT << ","
                << passP << ","
                << passE << ","
                << passAll << "\n";
        }
    }

    std::cout << "\nStep 5L completed successfully.\n";
    std::cout << "Output: plots/final_8state_verification_step5l.csv\n";

    return 0;
}
