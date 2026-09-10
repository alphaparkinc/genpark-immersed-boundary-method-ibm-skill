import math

class ImmersedBoundaryMethod:
    """
    Immersed Boundary Method (IBM) discrete Dirac delta kernel.
    Interpolates Eulerian fluid grid velocities to Lagrangian solid points,
    and spreads Lagrangian forces back to Eulerian grid.
    """
    def dirac_delta_1d(self, r):
        abs_r = abs(r)
        if abs_r <= 1.0:
            return 0.125 * (3.0 - 2.0*abs_r + math.sqrt(1.0 + 4.0*abs_r - 4.0*(abs_r**2)))
        elif abs_r <= 2.0:
            return 0.125 * (5.0 - 2.0*abs_r - math.sqrt(-7.0 + 12.0*abs_r - 4.0*(abs_r**2)))
        return 0.0

    def dirac_delta_2d(self, rx, ry):
        return self.dirac_delta_1d(rx) * self.dirac_delta_1d(ry)

    def interpolate_velocity(self, grid_vel, lagrange_pos, grid_dx=1.0):
        lx, ly = lagrange_pos
        u_lagrange = 0.0
        for i in range(len(grid_vel)):
            for j in range(len(grid_vel[0])):
                weight = self.dirac_delta_2d((lx - i*grid_dx) / grid_dx, (ly - j*grid_dx) / grid_dx)
                u_lagrange += grid_vel[i][j] * weight
        return u_lagrange
