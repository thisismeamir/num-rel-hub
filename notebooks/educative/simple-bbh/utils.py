import numpy as np
import numba 
@numba.jit(nopython=True)
def compute_christoffel_symbols(gamma_xx, gamma_xy, gamma_yy, dx, dy):
    """Compute Christoffel symbols using Numba"""
    nx, ny = gamma_xx.shape

    # Preallocate arrays for derivatives
    gamma_xx_dx = np.zeros_like(gamma_xx)
    gamma_xy_dx = np.zeros_like(gamma_xx)
    gamma_yy_dx = np.zeros_like(gamma_xx)
    gamma_xx_dy = np.zeros_like(gamma_xx)
    gamma_xy_dy = np.zeros_like(gamma_xx)
    gamma_yy_dy = np.zeros_like(gamma_xx)

    # Compute derivatives using finite differences
    for i in range(1, nx - 1):
        for j in range(ny):
            gamma_xx_dx[i, j] = (gamma_xx[i + 1, j] - gamma_xx[i - 1, j]) / (2 * dx)
            gamma_xy_dx[i, j] = (gamma_xy[i + 1, j] - gamma_xy[i - 1, j]) / (2 * dx)
            gamma_yy_dx[i, j] = (gamma_yy[i + 1, j] - gamma_yy[i - 1, j]) / (2 * dx)

    for i in range(nx):
        for j in range(1, ny - 1):
            gamma_xx_dy[i, j] = (gamma_xx[i, j + 1] - gamma_xx[i, j - 1]) / (2 * dy)
            gamma_xy_dy[i, j] = (gamma_xy[i, j + 1] - gamma_xy[i, j - 1]) / (2 * dy)
            gamma_yy_dy[i, j] = (gamma_yy[i, j + 1] - gamma_yy[i, j - 1]) / (2 * dy)

    return gamma_xx_dx, gamma_xy_dx, gamma_yy_dx, gamma_xx_dy, gamma_xy_dy, gamma_yy_dy


@numba.jit(nopython=True)
def evolve_fields(chi, gamma_xx, gamma_xy, gamma_yy, K, A_xx, A_xy, A_yy,
                  alpha, beta_x, beta_y, dx, dy, dt):
    """Evolve BSSN fields one timestep using Numba"""
    nx, ny = chi.shape

    # Preallocate arrays for new values
    chi_new = np.copy(chi)
    gamma_xx_new = np.copy(gamma_xx)
    gamma_xy_new = np.copy(gamma_xy)
    gamma_yy_new = np.copy(gamma_yy)
    K_new = np.copy(K)
    A_xx_new = np.copy(A_xx)
    A_xy_new = np.copy(A_xy)
    A_yy_new = np.copy(A_yy)

    # Compute derivatives
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            # Compute derivatives
            chi_dx = (chi[i + 1, j] - chi[i - 1, j]) / (2 * dx)
            chi_dy = (chi[i, j + 1] - chi[i, j - 1]) / (2 * dy)

            beta_x_dx = (beta_x[i + 1, j] - beta_x[i - 1, j]) / (2 * dx)
            beta_y_dy = (beta_y[i, j + 1] - beta_y[i, j - 1]) / (2 * dy)

            # Evolution equations
            # Modified chi evolution with momentum terms
            chi_dt = (2.0 / 3.0) * chi[i, j] * (
                    alpha[i, j] * K[i, j] -
                    beta_x_dx - beta_y_dy +
                    beta_x[i, j] * chi_dx + beta_y[i, j] * chi_dy
            )

            # Modified conformal metric evolution with momentum terms
            gamma_xx_dt = (-2 * alpha[i, j] * A_xx[i, j] +
                           beta_x[i, j] * (gamma_xx[i + 1, j] - gamma_xx[i - 1, j]) / (2 * dx) +
                           beta_y[i, j] * (gamma_xx[i, j + 1] - gamma_xx[i, j - 1]) / (2 * dy))

            gamma_xy_dt = (-2 * alpha[i, j] * A_xy[i, j] +
                           beta_x[i, j] * (gamma_xy[i + 1, j] - gamma_xy[i - 1, j]) / (2 * dx) +
                           beta_y[i, j] * (gamma_xy[i, j + 1] - gamma_xy[i, j - 1]) / (2 * dy))

            gamma_yy_dt = (-2 * alpha[i, j] * A_yy[i, j] +
                           beta_x[i, j] * (gamma_yy[i + 1, j] - gamma_yy[i - 1, j]) / (2 * dx) +
                           beta_y[i, j] * (gamma_yy[i, j + 1] - gamma_yy[i, j - 1]) / (2 * dy))

            # Evolution equation for K with momentum terms
            K_dt = (alpha[i, j] * (A_xx[i, j] ** 2 + 2 * A_xy[i, j] ** 2 + A_yy[i, j] ** 2) -
                    alpha[i, j] * K[i, j] ** 2 / 3 +
                    beta_x[i, j] * (K[i + 1, j] - K[i - 1, j]) / (2 * dx) +
                    beta_y[i, j] * (K[i, j + 1] - K[i, j - 1]) / (2 * dy))

            # Update fields
            chi_new[i, j] = chi[i, j] + dt * chi_dt
            gamma_xx_new[i, j] = gamma_xx[i, j] + dt * gamma_xx_dt
            gamma_xy_new[i, j] = gamma_xy[i, j] + dt * gamma_xy_dt
            gamma_yy_new[i, j] = gamma_yy[i, j] + dt * gamma_yy_dt
            K_new[i, j] = K[i, j] + dt * K_dt

    return chi_new, gamma_xx_new, gamma_xy_new, gamma_yy_new, K_new
