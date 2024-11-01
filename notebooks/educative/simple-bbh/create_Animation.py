import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
def create_animation(solver, num_frames=200):
    fig, ax = plt.subplots(figsize=(10, 10))

    def update(frame):
        ax.clear()
        solver.evolve_BSSN()

        # Plot conformal factor and horizons
        contour = ax.contourf(solver.X, solver.Y, solver.chi,
                              levels=np.linspace(1, 2, 20), cmap='viridis')

        horizons = solver.get_horizon_locations()
        for x, y in horizons:
            ax.plot(x, y, 'ko', markersize=10)

        ax.set_xlim(-5, 5)
        ax.set_ylim(-5, 5)
        ax.set_title(f'Time = {solver.t:.2f}')

        return contour,

    anim = FuncAnimation(fig, update, frames=num_frames, interval=50, blit=True)
    plt.close()
    return anim
