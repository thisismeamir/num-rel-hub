from scipy.ndimage import maximum_filter
import numba
import numpy as np

class bbh:
    t = 0

    def __init__(self, nx=100, ny=100, dx=0.1, dy=0.1, dt=0.01):
        self.nx = nx
        self.ny = ny
        self.dx = dx
        self.dy = dy
        self.dt = dt

        # Initialize fields and setup grid
        self.x = np.linspace(-5, 5, nx)
        self.y = np.linspace(-5, 5, ny)
        self.X, self.Y = np.meshgrid(self.x, self.y)
        self.initialize_fields()
        self.t = 0

    def initialize_fields(self):
        """Initialize the BSSN variables with orbital motion"""
        self.chi = np.ones((self.nx, self.ny))
        self.gamma_xx = np.ones((self.nx, self.ny))
        self.gamma_xy = np.zeros((self.nx, self.ny))
        self.gamma_yy = np.ones((self.nx, self.ny))
        self.K = np.zeros((self.nx, self.ny))
        self.A_xx = np.zeros((self.nx, self.ny))
        self.A_xy = np.zeros((self.nx, self.ny))
        self.A_yy = np.zeros((self.nx, self.ny))
        self.alpha = np.ones((self.nx, self.ny))
        self.beta_x = np.zeros((self.nx, self.ny))
        self.beta_y = np.zeros((self.nx, self.ny))

        self.setup_binary_black_holes()

    def setup_binary_black_holes(self):
        """Set up initial data for binary black holes with orbital motion"""
        M1 = 1.0
        M2 = 1.0
        d = 2.0  # Initial separation

        # Initial orbital angular velocity (Kepler approximation)
        omega = np.sqrt((M1 + M2) / d ** 3)

        # Set up initial positions and momenta
        for i in range(self.nx):
            for j in range(self.ny):
                # Calculate positions relative to origin
                x = self.X[i, j]
                y = self.Y[i, j]

                # Rotating coordinate positions of the black holes
                x1 = -d / 2 * np.cos(self.t * omega)
                y1 = -d / 2 * np.sin(self.t * omega)
                x2 = d / 2 * np.cos(self.t * omega)
                y2 = d / 2 * np.sin(self.t * omega)

                # Calculate distances to both black holes
                r1 = np.sqrt((x - x1) ** 2 + (y - y1) ** 2)
                r2 = np.sqrt((x - x2) ** 2 + (y - y2) ** 2)

                # Avoid singularities
                r1 = max(r1, 0.1)
                r2 = max(r2, 0.1)

                # Conformal factor (Brill-Lindquist with momentum)
                self.chi[i, j] = 1.0 + M1 / (2 * r1) + M2 / (2 * r2)

                # Initial lapse
                self.alpha[i, j] = 1.0 / self.chi[i, j] ** 2

                # Initial shift (orbital motion)
                self.beta_x[i, j] = -omega * y
                self.beta_y[i, j] = omega * x

                # Initial extrinsic curvature for orbital motion
                r = np.sqrt(x ** 2 + y ** 2)
                if r > 0.2:
                    self.K[i, j] = -omega * (x * self.beta_y[i, j] - y * self.beta_x[i, j]) / r ** 2

                    # Add momentum terms to A_ij
                    self.A_xx[i, j] = omega * (x * y) / r ** 2
                    self.A_xy[i, j] = -omega * (x ** 2 - y ** 2) / (2 * r ** 2)
                    self.A_yy[i, j] = -self.A_xx[i, j]

    def update_gauge(self):
        """Update gauge conditions"""
        # Modified 1+log slicing with orbital terms
        alpha_dx = np.gradient(self.alpha, self.dx, axis=0)
        alpha_dy = np.gradient(self.alpha, self.dy, axis=1)

        self.alpha += self.dt * (
                -2 * self.alpha * self.K +
                self.beta_x * alpha_dx +
                self.beta_y * alpha_dy
        )

        # Gamma-driver shift with orbital terms
        self.beta_x += self.dt * (
                0.75 * self.alpha * np.gradient(self.gamma_xx, self.dx, axis=0) +
                0.25 * np.gradient(self.beta_x, self.dx, axis=0)
        )
        self.beta_y += self.dt * (
                0.75 * self.alpha * np.gradient(self.gamma_yy, self.dy, axis=1) +
                0.25 * np.gradient(self.beta_y, self.dy, axis=1)
        )

    def evolve_BSSN(self):
        # Update time and evolve fields
        self.t += self.dt
        self.chi, self.gamma_xx, self.gamma_xy, self.gamma_yy, self.K = evolve_fields(
            self.chi, self.gamma_xx, self.gamma_xy, self.gamma_yy,
            self.K, self.A_xx, self.A_xy, self.A_yy,
            self.alpha, self.beta_x, self.beta_y,
            self.dx, self.dy, self.dt
        )
        self.update_gauge()
        self.apply_boundary_conditions()
        self.setup_binary_black_holes()

    def apply_boundary_conditions(self):
        # Implementing boundary conditions
        fields = [self.chi, self.gamma_xx, self.gamma_xy, self.gamma_yy,
                  self.K, self.alpha, self.beta_x, self.beta_y]
        for field in fields:
            field[:, 0] = field[:, 1]
            field[:, -1] = field[:, -2]
            field[0, :] = field[1, :]
            field[-1, :] = field[-2, :]

    def get_horizon_locations(self):
        # Detecting apparent horizons
        local_max = maximum_filter(self.chi, size=3) == self.chi
        maxima = np.where(local_max & (self.chi > 1.5))
        return list(zip(self.X[maxima], self.Y[maxima]))

