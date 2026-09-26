
#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

constexpr double PI_GLOBAL = 3.14159265358979323846;

struct State {
    std::string name;
    double verified_phase_deg;
    double width_nm;
    double height_nm;
    double T;
};

struct Cell {
    double x_um;
    double ideal_phase_deg;
    State state;
};

struct AxisCPML {
    std::vector<double> kappa_e, b_e, c_e;
    std::vector<double> kappa_h, b_h, c_h;
};

struct MonitorLine {
    std::vector<std::complex<double>> E;
    std::vector<std::complex<double>> H;
};

struct SimulationOutput {
    MonitorLine incident;
    MonitorLine postlens;
    MonitorLine focus;
};

inline int idx(int x, int y, int Nx) {
    return y * Nx + x;
}

double wrap360(double deg) {
    double p = std::fmod(deg, 360.0);
    if (p < 0.0) p += 360.0;
    return p;
}

double circular_error(double a, double b) {
    double d = std::fmod(a - b + 540.0, 360.0) - 180.0;
    return std::abs(d);
}

AxisCPML make_cpml_axis(
    int N, int npml, double dx, double dt,
    int m = 3, double kappa_max = 6.0, double target_R = 1.0e-8
) {
    constexpr double c0   = 299792458.0;
    constexpr double mu0  = 4.0e-7 * PI_GLOBAL;
    constexpr double eps0 = 1.0 / (mu0 * c0 * c0);
    const double eta0 = std::sqrt(mu0 / eps0);

    AxisCPML p;
    p.kappa_e.assign(N, 1.0); p.b_e.assign(N, 1.0); p.c_e.assign(N, 0.0);
    p.kappa_h.assign(N, 1.0); p.b_h.assign(N, 1.0); p.c_h.assign(N, 0.0);

    const double L = npml * dx;
    const double sigma_max =
        -(m + 1.0) * std::log(target_R) / (2.0 * eta0 * L);
    const double alpha_max = 0.05 * sigma_max;

    auto rho_e = [&](int i) {
        if (i < npml) return static_cast<double>(npml - i) / npml;
        if (i >= N - npml) return static_cast<double>(i - (N - npml - 1)) / npml;
        return 0.0;
    };

    auto rho_h = [&](int i) {
        const double ih = i + 0.5;
        if (ih < npml) return (npml - ih) / npml;
        if (ih >= N - npml) return (ih - (N - npml)) / npml;
        return 0.0;
    };

    auto coeff = [&](double rho, bool magnetic, double& kappa, double& b, double& c) {
        if (rho <= 0.0) {
            kappa = 1.0; b = 1.0; c = 0.0; return;
        }

        const double rp = std::pow(rho, m);
        kappa = 1.0 + (kappa_max - 1.0) * rp;

        double sigma = sigma_max * rp;
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
        c = (std::abs(denom) < 1e-30) ? 0.0 : sigma * (b - 1.0) / denom;
    };

    for (int i = 0; i < N; ++i) {
        coeff(rho_e(i), false, p.kappa_e[i], p.b_e[i], p.c_e[i]);
        coeff(rho_h(i), true,  p.kappa_h[i], p.b_h[i], p.c_h[i]);
    }

    return p;
}

MonitorLine zero_monitor(int Nx) {
    return {
        std::vector<std::complex<double>>(Nx, {0.0, 0.0}),
        std::vector<std::complex<double>>(Nx, {0.0, 0.0})
    };
}

void accumulate_monitor(
    MonitorLine& line,
    const std::vector<double>& Ez,
    const std::vector<double>& Hx,
    int y,
    int Nx,
    const std::complex<double>& kernel
) {
    for (int ix = 0; ix < Nx; ++ix) {
        const double E = Ez[idx(ix, y, Nx)];
        const double H = 0.5 * (
            Hx[idx(ix, y - 1, Nx)] + Hx[idx(ix, y, Nx)]
        );

        line.E[ix] += E * kernel;
        line.H[ix] += H * kernel;
    }
}

SimulationOutput run_simulation(
    bool with_lens,
    const std::vector<Cell>& cells,
    int Nx, int Ny, int Nt, int npml,
    double dx, double dx_um, double dt,
    int source_y, int incident_y, int postlens_y, int focus_y,
    double omega0,
    double lens_base_y_um,
    double domain_x_um,
    double pillar_n
) {
    constexpr double mu0  = 4.0e-7 * PI_GLOBAL;
    constexpr double c0   = 299792458.0;
    constexpr double eps0 = 1.0 / (mu0 * c0 * c0);

    const int N = Nx * Ny;

    std::vector<double> Ez(N, 0.0);
    std::vector<double> Hx(N, 0.0);
    std::vector<double> Hy(N, 0.0);
    std::vector<double> eps_r(N, 1.0);

    std::vector<double> psi_hx_y(N, 0.0);
    std::vector<double> psi_hy_x(N, 0.0);
    std::vector<double> psi_ez_x(N, 0.0);
    std::vector<double> psi_ez_y(N, 0.0);

    const AxisCPML px = make_cpml_axis(Nx, npml, dx, dt);
    const AxisCPML py = make_cpml_axis(Ny, npml, dx, dt);

    auto x_coord_um = [&](int ix) {
        return (ix + 0.5) * dx_um - 0.5 * domain_x_um;
    };

    auto y_coord_um = [&](int iy) {
        return (iy + 0.5) * dx_um;
    };

    if (with_lens) {
        // Exact rectangle/cell overlap for subpixel Ez-cell permittivity.
        // For TMz here, Ez is tangential to the pillar interfaces, so an
        // area-weighted arithmetic epsilon is a reasonable first-order
        // subpixel representation.
        const double eps_pillar = pillar_n * pillar_n;

        auto overlap_1d = [](double a0, double a1, double b0, double b1) {
            return std::max(0.0, std::min(a1, b1) - std::max(a0, b0));
        };

        for (const auto& cell : cells) {
            const double w_um = cell.state.width_nm * 1e-3;
            const double h_um = cell.state.height_nm * 1e-3;

            const double rx0 = cell.x_um - 0.5 * w_um;
            const double rx1 = cell.x_um + 0.5 * w_um;
            const double ry0 = lens_base_y_um;
            const double ry1 = lens_base_y_um + h_um;

            int ix0 = std::max(0, static_cast<int>(
                std::floor((rx0 + 0.5 * domain_x_um) / dx_um)) - 2);
            int ix1 = std::min(Nx - 1, static_cast<int>(
                std::ceil((rx1 + 0.5 * domain_x_um) / dx_um)) + 2);
            int iy0 = std::max(0, static_cast<int>(std::floor(ry0 / dx_um)) - 2);
            int iy1 = std::min(Ny - 1, static_cast<int>(std::ceil(ry1 / dx_um)) + 2);

            for (int iy = iy0; iy <= iy1; ++iy) {
                const double yc = y_coord_um(iy);
                const double cy0 = yc - 0.5 * dx_um;
                const double cy1 = yc + 0.5 * dx_um;
                const double oy = overlap_1d(cy0, cy1, ry0, ry1);
                if (oy <= 0.0) continue;

                for (int ix = ix0; ix <= ix1; ++ix) {
                    const double xc = x_coord_um(ix);
                    const double cx0 = xc - 0.5 * dx_um;
                    const double cx1 = xc + 0.5 * dx_um;
                    const double ox = overlap_1d(cx0, cx1, rx0, rx1);
                    if (ox <= 0.0) continue;

                    const double fill = (ox * oy) / (dx_um * dx_um);
                    const double eps_eff = 1.0 + fill * (eps_pillar - 1.0);
                    const int k = idx(ix, iy, Nx);
                    eps_r[k] = std::max(eps_r[k], eps_eff);
                }
            }
        }
    }

    MonitorLine incident = zero_monitor(Nx);
    MonitorLine postlens = zero_monitor(Nx);
    MonitorLine focus = zero_monitor(Nx);

    const double sigma_t = 10.0e-15;
    const double t0 = 6.0 * sigma_t;

    std::cout << (with_lens ? "Lens" : "Vacuum")
              << " simulation started...\n";

    for (int n = 0; n < Nt; ++n) {
        #pragma omp parallel for schedule(static)
        for (int iy = 0; iy < Ny - 1; ++iy) {
            const double inv_ky = 1.0 / py.kappa_h[iy];

            for (int ix = 0; ix < Nx - 1; ++ix) {
                const int k = idx(ix, iy, Nx);

                const double dEz_dy =
                    (Ez[idx(ix, iy + 1, Nx)] - Ez[k]) / dx;

                psi_hx_y[k] =
                    py.b_h[iy] * psi_hx_y[k]
                    + py.c_h[iy] * dEz_dy;

                Hx[k] -= (dt / mu0) * (
                    inv_ky * dEz_dy + psi_hx_y[k]
                );

                const double dEz_dx =
                    (Ez[idx(ix + 1, iy, Nx)] - Ez[k]) / dx;

                psi_hy_x[k] =
                    px.b_h[ix] * psi_hy_x[k]
                    + px.c_h[ix] * dEz_dx;

                Hy[k] += (dt / mu0) * (
                    (1.0 / px.kappa_h[ix]) * dEz_dx
                    + psi_hy_x[k]
                );
            }
        }

        #pragma omp parallel for schedule(static)
        for (int iy = 1; iy < Ny - 1; ++iy) {
            const double inv_ky = 1.0 / py.kappa_e[iy];

            for (int ix = 1; ix < Nx - 1; ++ix) {
                const int k = idx(ix, iy, Nx);

                const double dHy_dx =
                    (Hy[k] - Hy[idx(ix - 1, iy, Nx)]) / dx;

                psi_ez_x[k] =
                    px.b_e[ix] * psi_ez_x[k]
                    + px.c_e[ix] * dHy_dx;

                const double term_x =
                    (1.0 / px.kappa_e[ix]) * dHy_dx
                    + psi_ez_x[k];

                const double dHx_dy =
                    (Hx[k] - Hx[idx(ix, iy - 1, Nx)]) / dx;

                psi_ez_y[k] =
                    py.b_e[iy] * psi_ez_y[k]
                    + py.c_e[iy] * dHx_dy;

                const double term_y =
                    inv_ky * dHx_dy
                    + psi_ez_y[k];

                Ez[k] +=
                    (dt / (eps0 * eps_r[k]))
                    * (term_x - term_y);
            }
        }

        const double time = n * dt;
        const double tau = time - t0;
        const double envelope =
            std::exp(-0.5 * tau * tau / (sigma_t * sigma_t));
        const double src =
            envelope * std::sin(omega0 * tau);

        for (int ix = npml + 2; ix < Nx - npml - 2; ++ix) {
            const double x = std::abs(x_coord_um(ix));
            double w = 0.0;

            if (x <= 11.0) {
                w = 1.0;
            } else if (x < 12.5) {
                const double u = (x - 11.0) / 1.5;
                w = 0.5 * (1.0 + std::cos(PI_GLOBAL * u));
            }

            Ez[idx(ix, source_y, Nx)] += w * src;
        }

        for (int ix = 0; ix < Nx; ++ix) {
            Ez[idx(ix, 0, Nx)] = 0.0;
            Ez[idx(ix, Ny - 1, Nx)] = 0.0;
        }

        for (int iy = 0; iy < Ny; ++iy) {
            Ez[idx(0, iy, Nx)] = 0.0;
            Ez[idx(Nx - 1, iy, Nx)] = 0.0;
        }

        const double phase = -omega0 * time;
        const std::complex<double> kernel(
            std::cos(phase), std::sin(phase)
        );

        accumulate_monitor(incident, Ez, Hx, incident_y, Nx, kernel);
        accumulate_monitor(postlens, Ez, Hx, postlens_y, Nx, kernel);
        accumulate_monitor(focus, Ez, Hx, focus_y, Nx, kernel);

        if (n % std::max(1, Nt / 20) == 0) {
            std::cout << "  progress "
                      << std::fixed << std::setprecision(1)
                      << 100.0 * n / Nt << "%\n";
        }
    }

    // Yee half-step H correction.
    const std::complex<double> half_step_phase(
        std::cos(omega0 * dt / 2.0),
        std::sin(omega0 * dt / 2.0)
    );

    for (int ix = 0; ix < Nx; ++ix) {
        incident.H[ix] *= half_step_phase;
        postlens.H[ix] *= half_step_phase;
        focus.H[ix] *= half_step_phase;
    }

    return {incident, postlens, focus};
}

double sy_at(const MonitorLine& line, int ix) {
    return 0.5 * std::real(
        line.E[ix] * std::conj(line.H[ix])
    );
}

double integrate_sy(
    const MonitorLine& line,
    double dx,
    double dx_um,
    double domain_x_um,
    double x_abs_max_um
) {
    double P = 0.0;

    for (int ix = 0; ix < static_cast<int>(line.E.size()); ++ix) {
        const double x =
            (ix + 0.5) * dx_um - 0.5 * domain_x_um;

        if (std::abs(x) <= x_abs_max_um) {
            P += sy_at(line, ix) * dx;
        }
    }

    return P;
}

int main() {
    constexpr double c0 = 299792458.0;

    // ------------------------------------------------------------
    // 33.333-nm comparison point for efficiency convergence.
    // ------------------------------------------------------------
    const double grid_nm = 100.0 / 3.0;
    const double dx_um = grid_nm * 1e-3;
    const double dx = dx_um * 1e-6;

    const double S = 0.70;
    const double dt = S * dx / c0;

    const double lambda_um = 1.55;
    const double lambda0 = lambda_um * 1e-6;
    const double f0 = c0 / lambda0;
    const double omega0 = 2.0 * PI_GLOBAL * f0;

    const double domain_x_um = 30.0;
    const double domain_y_um = 34.0;
    const double cpml_um = 2.0;

    const int Nx = static_cast<int>(std::round(domain_x_um / dx_um));
    const int Ny = static_cast<int>(std::round(domain_y_um / dx_um));
    const int npml = static_cast<int>(std::round(cpml_um / dx_um));

    const double source_y_um = 3.0;
    const double incident_monitor_y_um = 4.5;
    const double lens_base_y_um = 6.0;
    const double postlens_monitor_y_um = 9.0;

    // Step 6C result at this grid.
    const double requested_focus_y_um = 25.6500;
    const double stabilized_fwhm_um = 1.3752;

    const double pillar_n = 2.0;
    const double focal_distance_um = 20.0;
    const double period_um = 1.0;

    auto nearest_y_index = [&](double y_um) {
        int iy = static_cast<int>(
            std::round(y_um / dx_um - 0.5)
        );
        return std::max(1, std::min(Ny - 2, iy));
    };

    auto y_coord_um = [&](int iy) {
        return (iy + 0.5) * dx_um;
    };

    const int source_y = nearest_y_index(source_y_um);
    const int incident_y = nearest_y_index(incident_monitor_y_um);
    const int postlens_y = nearest_y_index(postlens_monitor_y_um);
    const int focus_y = nearest_y_index(requested_focus_y_um);

    const double actual_focus_monitor_y_um = y_coord_um(focus_y);

    const double total_time_fs = 420.0;
    const int Nt = static_cast<int>(
        std::ceil(total_time_fs * 1e-15 / dt)
    ) + 2;

    const std::vector<State> library = {
        {"S0", 7.654000, 600.00000, 1800.00000, 0.923399000},
        {"S45", 44.568000, 450.00000, 1833.33000, 0.998266000},
        {"S90", 87.311000, 666.67000, 1333.33000, 0.940155000},
        {"S135", 140.156000, 300.00000, 1600.00000, 0.974974000},
        {"S180", 174.155000, 500.00000, 1000.00000, 0.909153000},
        {"S225", 224.540000, 300.00000, 1000.00000, 0.976975000},
        {"S270", 271.134000, 300.00000, 616.67000, 0.981414000},
        {"S315", 313.983000, 183.33000, 533.33000, 0.960412000},
    };

    const double k0_um = 2.0 * PI_GLOBAL / lambda_um;
    std::vector<Cell> cells;

    for (int m = -10; m <= 10; ++m) {
        const double x_um = m * period_um;

        const double opd =
            std::sqrt(x_um*x_um + focal_distance_um*focal_distance_um)
            - focal_distance_um;

        const double ideal_phase_deg =
            wrap360(+k0_um * opd * 180.0 / PI_GLOBAL);

        int best = 0;
        double best_err = std::numeric_limits<double>::infinity();

        for (int i = 0; i < static_cast<int>(library.size()); ++i) {
            const double e = circular_error(
                library[i].verified_phase_deg,
                ideal_phase_deg
            );
            if (e < best_err) {
                best_err = e;
                best = i;
            }
        }

        cells.push_back({x_um, ideal_phase_deg, library[best]});
    }

    std::cout << "Step 6G 33p333nm: subpixel efficiency run\n";
    std::cout << "============================================\n";
    std::cout << "Grid case: 33p333nm.\n";
    std::cout << "Using exact rectangle/cell overlap subpixel permittivity.\n";
    std::cout << "Nx x Ny = " << Nx << " x " << Ny
              << ", Nt = " << Nt << "\n";

#ifdef _OPENMP
    std::cout << "OpenMP enabled, max threads = "
              << omp_get_max_threads() << "\n";
#endif

    const auto vac = run_simulation(
        false, cells,
        Nx, Ny, Nt, npml,
        dx, dx_um, dt,
        source_y, incident_y, postlens_y, focus_y,
        omega0,
        lens_base_y_um,
        domain_x_um,
        pillar_n
    );

    const auto lens = run_simulation(
        true, cells,
        Nx, Ny, Nt, npml,
        dx, dx_um, dt,
        source_y, incident_y, postlens_y, focus_y,
        omega0,
        lens_base_y_um,
        domain_x_um,
        pillar_n
    );

    const double aperture_half_width_um = 10.5;
    const double source_half_width_um = 12.5;

    const double P_inc_aperture =
        integrate_sy(vac.incident, dx, dx_um, domain_x_um, aperture_half_width_um);

    const double P_inc_source =
        integrate_sy(vac.incident, dx, dx_um, domain_x_um, source_half_width_um);

    const double P_postlens =
        integrate_sy(lens.postlens, dx, dx_um, domain_x_um, source_half_width_um);

    const double P_focus_total =
        integrate_sy(lens.focus, dx, dx_um, domain_x_um, source_half_width_um);

    const double half1 = 0.5 * stabilized_fwhm_um;
    const double half2 = 1.0 * stabilized_fwhm_um;
    const double half3 = 1.5 * stabilized_fwhm_um;

    const double P1 =
        integrate_sy(lens.focus, dx, dx_um, domain_x_um, half1);

    const double P2 =
        integrate_sy(lens.focus, dx, dx_um, domain_x_um, half2);

    const double P3 =
        integrate_sy(lens.focus, dx, dx_um, domain_x_um, half3);

    const double eps = 1e-30;

    const double eta1 = P1 / (P_inc_aperture + eps);
    const double eta2 = P2 / (P_inc_aperture + eps);
    const double eta3 = P3 / (P_inc_aperture + eps);

    const double postT =
        P_postlens / (P_inc_source + eps);

    const double focusFrac =
        P_focus_total / (P_inc_source + eps);

    const double capture3 =
        P3 / (P_focus_total + eps);

    {
        std::ofstream out("plots/subpixel_summary_33p333nm_step6g.csv");

        out << "grid_nm,focus_monitor_y_um,fwhm_um,"
               "incident_aperture_power,incident_source_power,"
               "postlens_power,focus_total_power,"
               "postlens_transmission,focus_plane_forward_fraction,"
               "eta_1xFWHM,eta_2xFWHM,eta_3xFWHM,capture_3xFWHM\n";

        out << std::setprecision(15)
            << grid_nm << ","
            << actual_focus_monitor_y_um << ","
            << stabilized_fwhm_um << ","
            << P_inc_aperture << ","
            << P_inc_source << ","
            << P_postlens << ","
            << P_focus_total << ","
            << postT << ","
            << focusFrac << ","
            << eta1 << ","
            << eta2 << ","
            << eta3 << ","
            << capture3 << "\n";
    }

    {
        std::ofstream profile("plots/subpixel_poynting_33p333nm_step6g.csv");
        profile << "x_um,Sy_focus,Sy_postlens,Sy_incident_vacuum\n";
        profile << std::setprecision(15);

        for (int ix = 0; ix < Nx; ++ix) {
            const double x =
                (ix + 0.5) * dx_um - 0.5 * domain_x_um;

            profile << x << ","
                    << sy_at(lens.focus, ix) << ","
                    << sy_at(lens.postlens, ix) << ","
                    << sy_at(vac.incident, ix) << "\n";
        }
    }

    std::cout << "\nStep 6G 33p333nm simulation completed.\n";
    std::cout << "33p333nm eta_3xFWHM = "
              << eta3 << "\n";
    std::cout << "Run the second grid case, then python/plot_results.py.\n";

    return 0;
}
