
#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
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
    State state;
};

struct AxisCPML {
    std::vector<double> kappa_e, b_e, c_e;
    std::vector<double> kappa_h, b_h, c_h;
};

struct LineDFT {
    std::vector<std::complex<double>> E;
    std::vector<std::complex<double>> H;
};

struct MultiFreqOutput {
    std::vector<LineDFT> incident;
    std::vector<LineDFT> postlens;
    std::vector<std::vector<LineDFT>> focus;
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
    const double d = std::fmod(a - b + 540.0, 360.0) - 180.0;
    return std::abs(d);
}

double overlap_1d(double a0, double a1, double b0, double b1) {
    return std::max(0.0, std::min(a1, b1) - std::max(a0, b0));
}

AxisCPML make_cpml_axis(
    int N, int npml, double dx, double dt,
    int m = 3, double kappa_max = 6.0, double target_R = 1e-8
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
        if (i >= N - npml)
            return static_cast<double>(i - (N - npml - 1)) / npml;
        return 0.0;
    };

    auto rho_h = [&](int i) {
        const double ih = i + 0.5;
        if (ih < npml) return (npml - ih) / npml;
        if (ih >= N - npml) return (ih - (N - npml)) / npml;
        return 0.0;
    };

    auto coeff = [&](double rho, bool magnetic,
                     double& kappa, double& b, double& c) {
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

        c = (std::abs(denom) < 1e-30)
            ? 0.0
            : sigma * (b - 1.0) / denom;
    };

    for (int i = 0; i < N; ++i) {
        coeff(rho_e(i), false, p.kappa_e[i], p.b_e[i], p.c_e[i]);
        coeff(rho_h(i), true,  p.kappa_h[i], p.b_h[i], p.c_h[i]);
    }

    return p;
}

LineDFT zero_line(int Nx) {
    return {
        std::vector<std::complex<double>>(Nx, {0.0, 0.0}),
        std::vector<std::complex<double>>(Nx, {0.0, 0.0})
    };
}

void accumulate_line(
    LineDFT& line,
    const std::vector<double>& Ez,
    const std::vector<double>& Hx,
    int y, int Nx,
    const std::complex<double>& kernel
) {
    for (int ix = 0; ix < Nx; ++ix) {
        const double E = Ez[idx(ix, y, Nx)];
        const double H = 0.5 * (
            Hx[idx(ix, y - 1, Nx)] +
            Hx[idx(ix, y, Nx)]
        );
        line.E[ix] += E * kernel;
        line.H[ix] += H * kernel;
    }
}

double sy_at(const LineDFT& line, int ix) {
    return 0.5 * std::real(line.E[ix] * std::conj(line.H[ix]));
}

double integrate_sy(
    const LineDFT& line,
    double dx, double dx_um,
    double domain_x_um,
    double half_width_um
) {
    double P = 0.0;

    for (int ix = 0; ix < static_cast<int>(line.E.size()); ++ix) {
        const double x =
            (ix + 0.5) * dx_um - 0.5 * domain_x_um;

        if (std::abs(x) <= half_width_um) {
            P += sy_at(line, ix) * dx;
        }
    }

    return P;
}

double center_pair_intensity(const LineDFT& line) {
    const int Nx = static_cast<int>(line.E.size());
    const int r = Nx / 2;
    const int l = r - 1;

    return 0.5 * (
        std::norm(line.E[l]) +
        std::norm(line.E[r])
    );
}

double parabolic_peak(
    double ym, double y0, double yp,
    double Im, double I0, double Ip
) {
    const double denom = Im - 2.0 * I0 + Ip;
    if (std::abs(denom) < 1e-30) return y0;

    double delta = 0.5 * (Im - Ip) / denom;
    delta = std::max(-1.0, std::min(1.0, delta));

    return y0 + delta * (yp - y0);
}

double estimate_fwhm(
    const LineDFT& line,
    double dx_um,
    double domain_x_um
) {
    const int Nx = static_cast<int>(line.E.size());
    std::vector<double> I(Nx, 0.0);

    for (int ix = 0; ix < Nx; ++ix) {
        I[ix] = std::norm(line.E[ix]);
    }

    int center = Nx / 2;

    for (int ix = std::max(1, center - 7);
         ix <= std::min(Nx - 2, center + 7);
         ++ix) {
        if (I[ix] > I[center]) center = ix;
    }

    const double peak = I[center];
    const double half = 0.5 * peak;

    int il = center;
    while (il > 1 && I[il] >= half) --il;

    int ir = center;
    while (ir < Nx - 2 && I[ir] >= half) ++ir;

    auto xcoord = [&](int ix) {
        return (ix + 0.5) * dx_um - 0.5 * domain_x_um;
    };

    auto cross = [&](int i1, int i2) {
        const double y1 = I[i1];
        const double y2 = I[i2];

        if (std::abs(y2 - y1) < 1e-30) {
            return 0.5 * (xcoord(i1) + xcoord(i2));
        }

        double t = (half - y1) / (y2 - y1);
        t = std::max(0.0, std::min(1.0, t));

        return xcoord(i1) + t * (xcoord(i2) - xcoord(i1));
    };

    return std::max(0.0, cross(ir - 1, ir) - cross(il, il + 1));
}

MultiFreqOutput run_sim(
    bool with_lens,
    const std::vector<Cell>& cells,
    const std::vector<double>& wavelengths_nm,
    int Nx, int Ny, int Nt, int npml,
    double dx, double dx_um, double dt,
    int source_y, int incident_y, int postlens_y,
    const std::vector<int>& focus_y_indices,
    double lens_base_y_um,
    double domain_x_um,
    double pillar_n
) {
    constexpr double c0   = 299792458.0;
    constexpr double mu0  = 4.0e-7 * PI_GLOBAL;
    constexpr double eps0 = 1.0 / (mu0 * c0 * c0);

    const int N = Nx * Ny;
    const int NF = static_cast<int>(wavelengths_nm.size());

    std::vector<double> omega(NF);

    for (int f = 0; f < NF; ++f) {
        omega[f] =
            2.0 * PI_GLOBAL * c0 / (wavelengths_nm[f] * 1e-9);
    }

    std::vector<double> Ez(N, 0.0), Hx(N, 0.0), Hy(N, 0.0);
    std::vector<double> eps_r(N, 1.0);

    std::vector<double> psi_hx_y(N, 0.0), psi_hy_x(N, 0.0);
    std::vector<double> psi_ez_x(N, 0.0), psi_ez_y(N, 0.0);

    const AxisCPML px = make_cpml_axis(Nx, npml, dx, dt);
    const AxisCPML py = make_cpml_axis(Ny, npml, dx, dt);

    auto xcenter = [&](int ix) {
        return (ix + 0.5) * dx_um - 0.5 * domain_x_um;
    };

    auto ycenter = [&](int iy) {
        return (iy + 0.5) * dx_um;
    };

    if (with_lens) {
        const double eps_p = pillar_n * pillar_n;

        for (const auto& cell : cells) {
            const double w_um = cell.state.width_nm * 1e-3;
            const double h_um = cell.state.height_nm * 1e-3;

            const double rx0 = cell.x_um - 0.5 * w_um;
            const double rx1 = cell.x_um + 0.5 * w_um;
            const double ry0 = lens_base_y_um;
            const double ry1 = lens_base_y_um + h_um;

            int ix0 = std::max(
                0, static_cast<int>(
                    std::floor((rx0 + 0.5 * domain_x_um) / dx_um)
                ) - 2
            );

            int ix1 = std::min(
                Nx - 1, static_cast<int>(
                    std::ceil((rx1 + 0.5 * domain_x_um) / dx_um)
                ) + 2
            );

            int iy0 = std::max(
                0, static_cast<int>(std::floor(ry0 / dx_um)) - 2
            );

            int iy1 = std::min(
                Ny - 1, static_cast<int>(std::ceil(ry1 / dx_um)) + 2
            );

            for (int iy = iy0; iy <= iy1; ++iy) {
                const double yc = ycenter(iy);

                const double oy = overlap_1d(
                    yc - 0.5 * dx_um,
                    yc + 0.5 * dx_um,
                    ry0, ry1
                );

                if (oy <= 0.0) continue;

                for (int ix = ix0; ix <= ix1; ++ix) {
                    const double xc = xcenter(ix);

                    const double ox = overlap_1d(
                        xc - 0.5 * dx_um,
                        xc + 0.5 * dx_um,
                        rx0, rx1
                    );

                    if (ox <= 0.0) continue;

                    const double fill =
                        (ox * oy) / (dx_um * dx_um);

                    const double eps_eff =
                        1.0 + fill * (eps_p - 1.0);

                    eps_r[idx(ix, iy, Nx)] =
                        std::max(
                            eps_r[idx(ix, iy, Nx)],
                            eps_eff
                        );
                }
            }
        }
    }

    MultiFreqOutput out;
    out.incident.assign(NF, zero_line(Nx));
    out.postlens.assign(NF, zero_line(Nx));
    out.focus.resize(NF);

    for (int f = 0; f < NF; ++f) {
        out.focus[f].assign(
            focus_y_indices.size(),
            zero_line(Nx)
        );
    }

    // Center on the notch/re-entry region for strong support across 2350-2450.
    const double lambda_c = 2400.0e-9;
    const double omega_c = 2.0 * PI_GLOBAL * c0 / lambda_c;
    const double sigma_t = 5.0e-15;
    const double t0 = 6.0 * sigma_t;

    std::cout
        << (with_lens ? "Subpixel metalens" : "Vacuum")
        << " dense upper-edge run...\n";

    for (int n = 0; n < Nt; ++n) {
        #pragma omp parallel for schedule(static)
        for (int iy = 0; iy < Ny - 1; ++iy) {
            const double inv_ky = 1.0 / py.kappa_h[iy];

            for (int ix = 0; ix < Nx - 1; ++ix) {
                const int k = idx(ix, iy, Nx);

                const double dEz_dy =
                    (Ez[idx(ix, iy + 1, Nx)] - Ez[k]) / dx;

                psi_hx_y[k] =
                    py.b_h[iy] * psi_hx_y[k] +
                    py.c_h[iy] * dEz_dy;

                Hx[k] -= (dt / mu0) * (
                    inv_ky * dEz_dy +
                    psi_hx_y[k]
                );

                const double dEz_dx =
                    (Ez[idx(ix + 1, iy, Nx)] - Ez[k]) / dx;

                psi_hy_x[k] =
                    px.b_h[ix] * psi_hy_x[k] +
                    px.c_h[ix] * dEz_dx;

                Hy[k] += (dt / mu0) * (
                    (1.0 / px.kappa_h[ix]) * dEz_dx +
                    psi_hy_x[k]
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
                    px.b_e[ix] * psi_ez_x[k] +
                    px.c_e[ix] * dHy_dx;

                const double dHx_dy =
                    (Hx[k] - Hx[idx(ix, iy - 1, Nx)]) / dx;

                psi_ez_y[k] =
                    py.b_e[iy] * psi_ez_y[k] +
                    py.c_e[iy] * dHx_dy;

                Ez[k] +=
                    (dt / (eps0 * eps_r[k])) *
                    (
                        (1.0 / px.kappa_e[ix]) * dHy_dx +
                        psi_ez_x[k] -
                        inv_ky * dHx_dy -
                        psi_ez_y[k]
                    );
            }
        }

        const double time = n * dt;
        const double tau = time - t0;

        const double envelope =
            std::exp(-0.5 * tau * tau / (sigma_t * sigma_t));

        const double src =
            envelope * std::sin(omega_c * tau);

        for (int ix = npml + 2; ix < Nx - npml - 2; ++ix) {
            const double ax = std::abs(xcenter(ix));
            double w = 0.0;

            if (ax <= 11.0) {
                w = 1.0;
            } else if (ax < 12.5) {
                const double u = (ax - 11.0) / 1.5;
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

        for (int f = 0; f < NF; ++f) {
            const double phase = -omega[f] * time;

            const std::complex<double> kernel(
                std::cos(phase),
                std::sin(phase)
            );

            accumulate_line(
                out.incident[f],
                Ez, Hx,
                incident_y, Nx,
                kernel
            );

            accumulate_line(
                out.postlens[f],
                Ez, Hx,
                postlens_y, Nx,
                kernel
            );

            for (std::size_t j = 0;
                 j < focus_y_indices.size();
                 ++j) {
                accumulate_line(
                    out.focus[f][j],
                    Ez, Hx,
                    focus_y_indices[j],
                    Nx,
                    kernel
                );
            }
        }

        if (n % std::max(1, Nt / 20) == 0) {
            std::cout
                << "  progress "
                << std::fixed
                << std::setprecision(1)
                << 100.0 * n / Nt
                << "%\n";
        }
    }

    for (int f = 0; f < NF; ++f) {
        const std::complex<double> corr(
            std::cos(omega[f] * dt / 2.0),
            std::sin(omega[f] * dt / 2.0)
        );

        for (int ix = 0; ix < Nx; ++ix) {
            out.incident[f].H[ix] *= corr;
            out.postlens[f].H[ix] *= corr;
        }

        for (auto& line : out.focus[f]) {
            for (int ix = 0; ix < Nx; ++ix) {
                line.H[ix] *= corr;
            }
        }
    }

    return out;
}

int main() {
    constexpr double c0 = 299792458.0;

    const double grid_nm = 100.0 / 3.0;
    const double dx_um = grid_nm * 1e-3;
    const double dx = dx_um * 1e-6;

    const double S = 0.70;
    const double dt = S * dx / c0;

    const double domain_x_um = 30.0;
    const double domain_y_um = 34.0;
    const double cpml_um = 2.0;

    const int Nx =
        static_cast<int>(std::round(domain_x_um / dx_um));

    const int Ny =
        static_cast<int>(std::round(domain_y_um / dx_um));

    const int npml =
        static_cast<int>(std::round(cpml_um / dx_um));

    const double source_y_um = 3.0;
    const double incident_y_um = 4.5;
    const double lens_base_y_um = 6.0;
    const double postlens_y_um = 9.0;
    const double pillar_n = 2.0;

    const double total_time_fs = 420.0;

    const int Nt =
        static_cast<int>(
            std::ceil(total_time_fs * 1e-15 / dt)
        ) + 2;

    auto nearest_y_index = [&](double y_um) {
        int iy = static_cast<int>(
            std::round(y_um / dx_um - 0.5)
        );

        return std::max(
            1,
            std::min(Ny - 2, iy)
        );
    };

    auto ycoord = [&](int iy) {
        return (iy + 0.5) * dx_um;
    };

    const int source_y =
        nearest_y_index(source_y_um);

    const int incident_y =
        nearest_y_index(incident_y_um);

    const int postlens_y =
        nearest_y_index(postlens_y_um);

    // Dense 10-nm sampling through the dip and recovery seen in Step 6K.
    std::vector<double> wavelengths_nm;
    for (double wl = 2350.0; wl <= 2450.0001; wl += 10.0) {
        wavelengths_nm.push_back(wl);
    }

    // Focus is around 20 um in this wavelength range.
    std::set<int> focus_set;

    for (double y = 18.5; y <= 21.5001; y += 0.25) {
        focus_set.insert(
            nearest_y_index(y)
        );
    }

    std::vector<int> focus_y_indices(
        focus_set.begin(),
        focus_set.end()
    );

    std::vector<double> focus_y_um;

    for (int iy : focus_y_indices) {
        focus_y_um.push_back(
            ycoord(iy)
        );
    }

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

    const double design_lambda_um = 1.55;
    const double design_focal_um = 20.0;
    const double period_um = 1.0;

    const double k0_design =
        2.0 * PI_GLOBAL / design_lambda_um;

    std::vector<Cell> cells;

    for (int m = -10; m <= 10; ++m) {
        const double x_um = m * period_um;

        const double opd =
            std::sqrt(
                x_um*x_um +
                design_focal_um*design_focal_um
            ) - design_focal_um;

        const double ideal_phase =
            wrap360(
                +k0_design * opd
                * 180.0 / PI_GLOBAL
            );

        int best = 0;
        double best_err =
            std::numeric_limits<double>::infinity();

        for (int i = 0;
             i < static_cast<int>(library.size());
             ++i) {
            const double e =
                circular_error(
                    library[i].verified_phase_deg,
                    ideal_phase
                );

            if (e < best_err) {
                best_err = e;
                best = i;
            }
        }

        cells.push_back({
            x_um,
            library[best]
        });
    }

    std::cout
        << "Step 6L: dense upper-threshold notch localization\n"
        << "================================================\n"
        << "Subpixel grid = " << grid_nm << " nm\n"
        << "Source center = 2400 nm, sigma_t = 5 fs\n"
        << "Wavelengths = 2350...2450 nm in 10-nm steps\n"
        << "Purpose: locate BOTH the downward and recovery crossings\n"
        << "Only two FDTD runs: vacuum + metalens\n";

#ifdef _OPENMP
    std::cout
        << "OpenMP enabled, max threads = "
        << omp_get_max_threads() << "\n";
#endif

    const auto vac = run_sim(
        false, cells, wavelengths_nm,
        Nx, Ny, Nt, npml,
        dx, dx_um, dt,
        source_y, incident_y, postlens_y,
        focus_y_indices,
        lens_base_y_um,
        domain_x_um,
        pillar_n
    );

    const auto lens = run_sim(
        true, cells, wavelengths_nm,
        Nx, Ny, Nt, npml,
        dx, dx_um, dt,
        source_y, incident_y, postlens_y,
        focus_y_indices,
        lens_base_y_um,
        domain_x_um,
        pillar_n
    );

    const double aperture_half_um = 10.5;
    const double common_half_um = 12.5;

    // 2400-nm source support reference for this local run.
    int i2400 = 0;

    for (int f = 0;
         f < static_cast<int>(wavelengths_nm.size());
         ++f) {
        if (std::abs(wavelengths_nm[f] - 2400.0) < 1e-9) {
            i2400 = f;
        }
    }

    const double Psrc2400 =
        integrate_sy(
            vac.incident[i2400],
            dx, dx_um,
            domain_x_um,
            common_half_um
        );

    std::ofstream out(
        "plots/dense_upper_notch_step6l.csv"
    );

    out
        << "wavelength_nm,selected_focus_y_um,interpolated_focus_y_um,"
           "fwhm_um,postlens_transmission,focus_plane_forward_fraction,"
           "eta_1xFWHM,eta_2xFWHM,eta_3xFWHM,capture_3xFWHM,"
           "source_power_relative_to_2400\n";

    out << std::setprecision(15);

    for (int f = 0;
         f < static_cast<int>(wavelengths_nm.size());
         ++f) {
        int bestj = 0;
        double bestI = -1.0;

        for (int j = 0;
             j < static_cast<int>(focus_y_indices.size());
             ++j) {
            const double I =
                center_pair_intensity(
                    lens.focus[f][j]
                );

            if (I > bestI) {
                bestI = I;
                bestj = j;
            }
        }

        const LineDFT& focus =
            lens.focus[f][bestj];

        const double ysel =
            focus_y_um[bestj];

        double yinterp =
            ysel;

        if (bestj > 0 &&
            bestj + 1 < static_cast<int>(focus_y_um.size())) {
            yinterp = parabolic_peak(
                focus_y_um[bestj - 1],
                focus_y_um[bestj],
                focus_y_um[bestj + 1],
                center_pair_intensity(lens.focus[f][bestj - 1]),
                center_pair_intensity(lens.focus[f][bestj]),
                center_pair_intensity(lens.focus[f][bestj + 1])
            );
        }

        const double fwhm =
            estimate_fwhm(
                focus,
                dx_um,
                domain_x_um
            );

        const double PincA =
            integrate_sy(
                vac.incident[f],
                dx, dx_um,
                domain_x_um,
                aperture_half_um
            );

        const double PincS =
            integrate_sy(
                vac.incident[f],
                dx, dx_um,
                domain_x_um,
                common_half_um
            );

        const double Ppost =
            integrate_sy(
                lens.postlens[f],
                dx, dx_um,
                domain_x_um,
                common_half_um
            );

        const double Pfocus =
            integrate_sy(
                focus,
                dx, dx_um,
                domain_x_um,
                common_half_um
            );

        const double P1 =
            integrate_sy(
                focus,
                dx, dx_um,
                domain_x_um,
                0.5 * fwhm
            );

        const double P2 =
            integrate_sy(
                focus,
                dx, dx_um,
                domain_x_um,
                1.0 * fwhm
            );

        const double P3 =
            integrate_sy(
                focus,
                dx, dx_um,
                domain_x_um,
                1.5 * fwhm
            );

        const double eps = 1e-30;

        out
            << wavelengths_nm[f] << ","
            << ysel << ","
            << yinterp << ","
            << fwhm << ","
            << Ppost / (PincS + eps) << ","
            << Pfocus / (PincS + eps) << ","
            << P1 / (PincA + eps) << ","
            << P2 / (PincA + eps) << ","
            << P3 / (PincA + eps) << ","
            << P3 / (Pfocus + eps) << ","
            << PincS / (Psrc2400 + eps)
            << "\n";
    }

    std::cout
        << "\nStep 6L simulations completed.\n"
        << "Run python/plot_results.py.\n";

    return 0;
}
