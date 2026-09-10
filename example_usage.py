from client import ImmersedBoundaryMethod

def main():
    print("=== Testing Immersed Boundary Method (IBM) ===")
    ibm = ImmersedBoundaryMethod()

    grid = [[2.0]*6 for _ in range(6)]
    pos = (2.5, 2.5)

    u_solid = ibm.interpolate_velocity(grid, pos)
    print(f"Interpolated velocity at submerged point {pos}: {round(u_solid, 4)}")
    assert abs(u_solid - 2.0) < 0.2
    print("=== All tests passed successfully! ===")

if __name__ == "__main__":
    main()
