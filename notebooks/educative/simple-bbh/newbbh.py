import numpy as np
import numba
from scipy.ndimage import maximum_filter
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from dataclasses import dataclass
from typing import Tuple, List


@dataclass
class SimulationParameters:
    """Parameters for the binary black hole simulation"""
    M1: float = 1.0  # Mass of first black hole
    M2: float = 1.0  # Mass of second black hole
    d: float = 3.0  # Initial separation
    nx: int = 200  # Grid points in x
    ny: int = 200  # Grid points in y
    dx: float = 0.1  # Grid spacing
    dy: float = 0.1  # Grid spacing
    dt: float = 0.01  # Time step
    x_min: float = -10.0
    x_max: float = 10.0
    y_min: float = -10.0
    y_max: float = 10.0



def compute_derivatives(field: np.ndarray, dx: float, dy: float) -> Tuple[np.ndarray, np.ndarray]:
    """Compute spatial derivatives using 4th order finite differences"""
    dx_field = np.zeros_like(field)
    dy_field = np.zeros_like(field)

    # 4th order central differences for bulk
    for i in range(2, field.shape[0] - 2):
        for j in range(field.shape[1]):
            dx_field[i, j] = (-field[i + 2, j] + 8 * field[i + 1, j] - 8 * field[i - 1, j] + field[i - 2, j]) / (
                        12 * dx)

    for i in range(field.shape[0]):
        for j in range(2, field.shape[1] - 2):
            dy_field[i, j] = (-field[i, j + 2] + 8 * field[i, j + 1] - 8 * field[i, j - 1] + field[i, j - 2]) / (
                        12 * dy)

    return dx_field, dy_field


class ImprovedBBH:
    def __init__(self, params: SimulationParameters):
        self.params = params
        self.t = 0.0
        self.setup_grid()
        self.initialize_fields()

    def setup_grid(self):
        """Initialize the computational grid"""
        self.x = np.linspace(self.params.x_min, self.params.x_max, self.params.nx)
        self.y = np.linspace(self.params.y_min, self.params.y_max, self.params.ny)
        self.X, self.Y = np.meshgrid(self.x, self.y)

    def initialize_fields(self):
        """Initialize BSSN fields with improved initial data"""
        self.chi = np.ones((self.params.nx, self.params.ny))
        self.gamma_xx = np.ones((self.params.nx, self.params.ny))
        self.gamma_xy = np.zeros((self.params.nx, self.params.ny))
        self.gamma_yy = np.ones((self.params.nx, self.params.ny))
        self.K = np.zeros((self.params.nx, self.params.ny))
        self.A_xx = np.zeros((self.params.nx, self.params.ny))
        self.A_xy = np.zeros((self.params.nx, self.params.ny))
        self.A_yy = np.zeros((self.params.nx, self.params.ny))
        self.alpha = np.ones((self.params.nx, self.params.ny))
        self.beta_x = np.zeros((self.params.nx, self.params.ny))
        self.beta_y = np.zeros((self.params.nx, self.params.ny))

        self.setup_binary_black_holes()

    def setup_binary_black_holes(self):
        """Set up initial data for binary black holes with improved orbital parameters"""
        M_total = self.params.M1 + self.params.M2
        eta = (self.params.M1 * self.params.M2) / (M_total * M_total)

        # Improved orbital frequency with 2PN corrections
        omega = np.sqrt(M_total / self.params.d ** 3) * (
                1 + (-0.75 - eta / 12) * (M_total / self.params.d) +
                (-(3.375 + (19 / 8) * eta - eta * eta / 24)) * (M_total / self.params.d) ** 2
        )

        # Calculate positions with time evolution
        phi = omega * self.t
        x1 = -self.params.d / 2 * np.cos(phi)
        y1 = -self.params.d / 2 * np.sin(phi)
        x2 = self.params.d / 2 * np.cos(phi)
        y2 = self.params.d / 2 * np.sin(phi)

        # Improved momentum parameters
        p_tang = eta * np.sqrt(M_total * self.params.d) * (
                1 + (-0.25 - eta / 12) * (M_total / self.params.d) +
                (-(3.375 + (19 / 8) * eta - eta * eta / 24)) * (M_total / self.params.d) ** 2
        )

        # Update fields with improved initial data
        for i in range(self.params.nx):
            for j in range(self.params.ny):
                x, y = self.X[i, j], self.Y[i, j]

                # Calculate distances with smoothing
                r1 = np.sqrt((x - x1) ** 2 + (y - y1) ** 2 + 1e-6)
                r2 = np.sqrt((x - x2) ** 2 + (y - y2) ** 2 + 1e-6)

                # Improved conformal factor
                self.chi[i, j] = 1.0 + self.params.M1 / (2 * r1) + self.params.M2 / (2 * r2)

                # Improved lapse function with better behavior near horizons
                self.alpha[i, j] = (1.0 / self.chi[i, j] ** 2) * (
                        1.0 - 0.5 * self.params.M1 / r1 - 0.5 * self.params.M2 / r2
                )

                # Improved shift vector
                self.beta_x[i, j] = -omega * y
                self.beta_y[i, j] = omega * x

                # Improved extrinsic curvature
                r = np.sqrt(x * x + y * y + 1e-6)
                if r > 0.5:
                    self.K[i, j] = -p_tang * (x * self.beta_y[i, j] - y * self.beta_x[i, j]) / r ** 3

                    # Improved traceless part of extrinsic curvature
                    self.A_xx[i, j] = p_tang * (x * y) / r ** 3
                    self.A_xy[i, j] = -0.5 * p_tang * (x * x - y * y) / r ** 3
                    self.A_yy[i, j] = -self.A_xx[i, j]

    
    def evolve_BSSN(self):
        """Evolve the BSSN equations one timestep with improved numerical methods"""
        # Update time
        self.t += self.params.dt

        # Compute derivatives
        chi_dx, chi_dy = compute_derivatives(self.chi, self.params.dx, self.params.dy)
        K_dx, K_dy = compute_derivatives(self.K, self.params.dx, self.params.dy)

        # Evolution equations with improved handling of constraints
        self.chi = self.evolve_field(self.chi,
                                     2.0 / 3.0 * self.chi * (
                                             self.alpha * self.K -
                                             np.gradient(self.beta_x, self.params.dx, axis=0) -
                                             np.gradient(self.beta_y, self.params.dy, axis=1)
                                     ))

        self.K = self.evolve_field(self.K,
                                   self.alpha * (self.A_xx ** 2 + 2 * self.A_xy ** 2 + self.A_yy ** 2) -
                                   self.alpha * self.K ** 2 / 3)

        # Update gauge conditions
        self.update_gauge()

        # Apply improved boundary conditions
        self.apply_boundary_conditions()

        # Update orbital motion
        self.setup_binary_black_holes()

    def evolve_field(self, field: np.ndarray, rhs: np.ndarray) -> np.ndarray:
        """Evolve a field using RK4 time integration"""
        k1 = self.params.dt * rhs
        k2 = self.params.dt * (rhs + 0.5 * k1)
        k3 = self.params.dt * (rhs + 0.5 * k2)
        k4 = self.params.dt * (rhs + k3)

        return field + (k1 + 2 * k2 + 2 * k3 + k4) / 6

    def update_gauge(self):
        """Update gauge conditions with improved damping terms"""
        # Improved 1+log slicing condition
        self.alpha = self.evolve_field(
            self.alpha,
            -2 * self.alpha * self.K + self.beta_x * np.gradient(self.alpha, self.params.dx, axis=0) +
            self.beta_y * np.gradient(self.alpha, self.params.dy, axis=1)
        )

        # Improved Gamma-driver shift condition
        self.beta_x = self.evolve_field(
            self.beta_x,
            0.75 * self.alpha * np.gradient(self.gamma_xx, self.params.dx, axis=0)
        )

        self.beta_y = self.evolve_field(
            self.beta_y,
            0.75 * self.alpha * np.gradient(self.gamma_yy, self.params.dy, axis=1)
        )

    def apply_boundary_conditions(self):
        """Apply improved boundary conditions"""
        fields = [self.chi, self.gamma_xx, self.gamma_xy, self.gamma_yy,
                  self.K, self.alpha, self.beta_x, self.beta_y]

        for field in fields:
            # Improved extrapolation at boundaries
            field[0, :] = 3 * field[1, :] / 2 - field[2, :] / 2
            field[-1, :] = 3 * field[-2, :] / 2 - field[-3, :] / 2
            field[:, 0] = 3 * field[:, 1] / 2 - field[:, 2] / 2
            field[:, -1] = 3 * field[:, -2] / 2 - field[:, -3] / 2

    def get_horizon_locations(self) -> List[Tuple[float, float]]:
        """Detect apparent horizons with improved algorithm"""
        # Use chi field to detect horizons
        local_max = maximum_filter(self.chi, size=5) == self.chi
        maxima = np.where(local_max & (self.chi > 1.5))
        return list(zip(self.X[maxima], self.Y[maxima]))


def create_improved_animation(solver: ImprovedBBH, num_frames: int = 200) -> FuncAnimation:
    """Create smoother animation with improved visualization"""
    fig, ax = plt.subplots(figsize=(12, 12))

    def update(frame):
        ax.clear()
        solver.evolve_BSSN()

        # Improved visualization with log scale for better contrast
        chi_plot = np.log10(solver.chi)
        contour = ax.contourf(solver.X, solver.Y, chi_plot,
                              levels=np.linspace(chi_plot.min(), chi_plot.max(), 30),
                              cmap='magma')

        # Plot horizons
        horizons = solver.get_horizon_locations()
        for x, y in horizons:
            ax.plot(x, y, 'wo', markersize=8, markeredgecolor='black')

        # Improved aesthetics
        ax.set_xlim(solver.params.x_min, solver.params.x_max)
        ax.set_ylim(solver.params.y_min, solver.params.y_max)
        ax.set_title(f'Binary Black Hole Evolution (t = {solver.t:.2f})')
        ax.set_xlabel('x')
        ax.set_ylabel('y')

        return contour,

    anim = FuncAnimation(fig, update, frames=num_frames, interval=20, blit=True)
    plt.close()
    return anim