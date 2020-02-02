import heterocl as hcl
import numpy as np
import time

import math


class grid:
  def __init__(self, min, max, dims ,pts_each_dim, pDim):
        self.max = max
        self.min = min
        self.dims = dims
        self.pts_each_dim = pts_each_dim
        self.pDim = pDim

        # Make some modifications to the initialized
        self.max[pDim] = self.min[pDim] + (self.max[pDim] - self.min[pDim])  * (1 - 1/self.pts_each_dim[pDim])
        self.dx = (self.max - self.min) / (self.pts_each_dim - 1.0)
        self.vs = []
        for i in range(0,3):
            tmp = np.arange(self.min[i], self.max[i] + self.dx[i], self.dx[i])
            broadcast_map = np.ones(self.dims, dtype=int)
            broadcast_map[i] = self.pts_each_dim[i]
            tmp = np.reshape(tmp, tuple(broadcast_map))
            self.vs.append(tmp)
      # Reshape

# This functino creates a cyclinderical shape
def CyclinderShape(grid, ignore_dim, center, radius):
    data = np.zeros(grid.pts_each_dim)

    for i in range (0, 3):
        if i != ignore_dim-1:
            data = data + np.power(grid.vs[i] - center[i],  2)
    data = np.sqrt(data) - radius
    return data





g = grid(np.array([-5.0, -5.0, -math.pi]), np.array([5.0, 5.0, math.pi]), 3 ,np.array([100,100,100]), 2)
shape = CyclinderShape(g, 3, np.zeros(3), 1)
print(shape)
# Declare constants
dx = np.array([0.34, 0.22, 0.21])
t_diff = 0.05




def HJ_PDE_solver(V_new, V_init, thetas,t):
    # Used for calculating time bound
    #step_bound = hcl.scalar(0, "step_bound")
    #delta_T = hcl.scalar(0, "delta_T")
    # These variables are used to dissipation calculation
    max_alpha1 = hcl.scalar(-1e9, "max_alpha1")
    max_alpha2 = hcl.scalar(-1e9, "max_alpha2")
    max_alpha3 = hcl.scalar(-1e9, "max_alpha3")

    # Custom function
    def my_abs(my_x):
        abs_value = hcl.scalar(0, "abs_value", dtype=hcl.Float())
        with hcl.if_(my_x > 0):
            abs_value.v = my_x
        with hcl.else_():
            abs_value.v = -my_x
        return abs_value[0]

    def my_sign(x):
        sign = hcl.scalar(0, "sign", dtype=hcl.Float())
        with hcl.if_(x == 0):
            sign[0] = 0
        with hcl.if_(x > 0):
            sign[0] = 1
        with hcl.if_(x < 0):
            sign[0] = -1
        return sign[0]

    # Calculate spatial derivative
    def spa_derivX(i,j,k):
        left_deriv = hcl.scalar(0, "left_deriv")
        right_deriv = hcl.scalar(0, "right_deriv")
        with hcl.if_(i == 0):
            left_boundary = hcl.scalar(0, "left_boundary")
            left_boundary[0] =  V_init[i,j,k] + my_abs(V_init[i+1,j,k] - V_init[i,j,k]) * my_sign(V_init[i,j,k])
            left_deriv[0] = (V_init[i,j,k] - left_boundary[0])/dx[0]
            right_deriv[0] = (V_init[i+1,j,k] - V_init[i,j,k])/dx[0]
        with hcl.elif_(i == V_init.shape[0] - 1):
            right_boundary = hcl.scalar(0, "right_boundary")
            right_boundary[0] = V_init[i,j,k] +my_abs(V_init[i,j,k] - V_init[i-1,j,k]) * my_sign(V_init[i,j,k])
            left_deriv[0] = (V_init[i,j,k] - V_init[i-1,j,k])/dx[0]
            right_deriv[0] = (right_boundary[0] -V_init[i,j,k])/dx[0]
        with hcl.elif_(i != 0 and i != V_init.shape[0] -  1):
            left_deriv[0] = (V_init[i, j, k] - V_init[i - 1, j, k]) / dx[0]
            right_deriv[0] = (V_init[i + 1, j, k] - V_init[i, j, k]) / dx[0]
        return left_deriv[0], right_deriv[0]

    def spa_derivY(i,j,k):
        left_deriv = hcl.scalar(0, "left_deriv")
        right_deriv = hcl.scalar(0, "right_deriv")
        with hcl.if_(j == 0):
            left_boundary = hcl.scalar(0, "left_boundary")
            left_boundary[0] = V_init[i, j, k] + my_abs(V_init[i, j + 1, k] - V_init[i, j, k]) * my_sign(
                V_init[i, j, k])
            left_deriv[0] = (V_init[i, j, k] - left_boundary[0]) / dx[1]
            right_deriv[0] = (V_init[i, j + 1, k] - V_init[i, j, k]) / dx[1]
        with hcl.elif_(j == V_init.shape[1] - 1):
            right_boundary = hcl.scalar(0, "right_boundary")
            right_boundary[0] = V_init[i, j, k] + my_abs(V_init[i, j, k] - V_init[i, j - 1, k]) * my_sign(
                V_init[i, j, k])
            left_deriv[0] = (V_init[i, j, k] - V_init[i, j - 1, k]) / dx[1]
            right_deriv[0] = (right_boundary[0] - V_init[i, j, k]) / dx[1]
        with hcl.elif_(j != 0 and j != V_init.shape[1] - 1):
            left_deriv[0] = (V_init[i, j, k] - V_init[i, j - 1, k]) / dx[1]
            right_deriv[0] = (V_init[i, j + 1, k] - V_init[i, j, k]) / dx[1]
        return left_deriv[0], right_deriv[0]


    def spa_derivT(i,j,k):
        left_deriv = hcl.scalar(0, "left_deriv")
        right_deriv = hcl.scalar(0, "right_deriv")
        with hcl.if_(k == 0):
            left_boundary = hcl.scalar(0, "left_boundary")
            left_boundary[0] = V_init[i,j,50]
            left_deriv[0] = (V_init[i, j, k] - left_boundary[0]) / dx[2]
            right_deriv[0] = (V_init[i, j, k + 1] - V_init[i, j, k]) / dx[2]
        with hcl.elif_(k == V_init.shape[2] - 1):
            right_boundary = hcl.scalar(0, "right_boundary")
            right_boundary[0] = V_init[i,j,0]
            left_deriv[0] = (V_init[i, j, k] - V_init[i, j, k - 1]) / dx[2]
            right_deriv[0] = (right_boundary[0] - V_init[i, j, k]) / dx[2]
        with hcl.elif_(k != 0 and k != V_init.shape[2] - 1):
            left_deriv[0] = (V_init[i, j, k] - V_init[i, j, k - 1]) / dx[2]
            right_deriv[0] = (V_init[i, j ,k + 1] - V_init[i, j, k]) / dx[2]
        return left_deriv[0], right_deriv[0]

    def step_bound(time): # Function to calculate time step
        stepBoundInv = hcl.scalar(0, "stepBoundInv")
        stepBound    = hcl.scalar(0, "stepBound")
        stepBoundInv[0] = max_alpha1[0]/dx[0] + max_alpha2[0]/dx[1] + max_alpha3[0]/dx[2]

        stepBound[0] = 0.8/stepBoundInv[0]
        with hcl.if_(stepBound > t_diff):
            stepBound[0] = t_diff
        time = time + stepBound[0]
        return time

    # Calculate Hamiltonian for every grid point in V_init
    with hcl.Stage("Hamiltonian"):
        with hcl.for_(0, V_init.shape[0], name="k") as k: # Plus 1 as for loop count stops at V_init.shape[0]
            with hcl.for_(0, V_init.shape[1], name="j") as j:
                with hcl.for_(0, V_init.shape[2], name="i") as i:
                    # Variables to calculate dV_dx
                    dV_dx_L = hcl.scalar(0, "dV_dx_L")
                    dV_dx_R = hcl.scalar(0, "dV_dx_R")
                    dV_dx = hcl.scalar(0, "dV_dx")
                    # Variables to calculate dV_dy
                    dV_dy_L = hcl.scalar(0, "dV_dy_L")
                    dV_dy_R = hcl.scalar(0, "dV_dy_R")
                    dV_dy = hcl.scalar(0, "dV_dy")
                    # Variables to calculate dV_dtheta
                    dV_dT_L = hcl.scalar(0, "dV_dT_L")
                    dV_dT_R = hcl.scalar(0, "dV_dT_R")
                    dV_dT = hcl.scalar(0, "dV_dT")
                    # Variables to keep track of dynamics
                    dx_dt = hcl.scalar(0, "dx_dt")
                    dy_dt = hcl.scalar(0, "dy_dt")
                    dtheta_dt = hcl.scalar(0, "dtheta_dt")

                    # No tensor slice operation
                    dV_dx_L[0], dV_dx_R[0] = spa_derivX(i, j, k)
                    dV_dy_L[0], dV_dy_R[0] = spa_derivY(i, j, k)
                    dV_dT_L[0], dV_dT_R[0] = spa_derivT(i, j, k)

                    # Calculate average gradient
                    dV_dx[0] = (dV_dx_L + dV_dx_R) / 2
                    dV_dy[0] = (dV_dy_L + dV_dy_R) / 2
                    dV_dT[0] = (dV_dT_L + dV_dT_R) / 2

                    # Declare optimal control
                    uOpt = hcl.scalar(0, "uOpt")

                    # Declare Velocity
                    vel = hcl.scalar(1,"vel")

                    # Assume that mode is min
                    with hcl.if_(dV_dT > 0):
                        uOpt.v = -uOpt.v


                    # Calculate dynamical rates of changes
                    dx_dt[0] = vel.v * hcl.cos(thetas[k])
                    dy_dt[0] = vel.v * hcl.sin(thetas[k])
                    dtheta_dt[0] = uOpt.v

                    # Calculate Hamiltonian terms:
                    V_new[i, j, k] =  -(dx_dt[0] * dV_dx[0] + dy_dt * dV_dy[0] + dtheta_dt[0] * dV_dT[0])

                    # Calculate dissipation step
                    dx_dt[0] = my_abs(dx_dt[0])
                    dy_dt[0] = my_abs(dy_dt[0])
                    dtheta_dt[0] = my_abs(dtheta_dt[0])
                    diss = hcl.scalar(0, "diss")
                    diss[0] = 0.5*((dV_dx_R[0] - dV_dx_L[0])*dx_dt[0] + (dV_dy_R[0] - dV_dy_L[0])*dy_dt[0] + (dV_dT_R[0] - dV_dT_L[0])* dtheta_dt[0])
                    V_new[i, j, k] = -(V_new[i, j, k] - diss[0])

                    # Calculate alphas
                    with hcl.if_(dx_dt[0] > max_alpha1):
                        max_alpha1[0] = dx_dt[0]
                    with hcl.if_(dy_dt[0] > max_alpha2):
                        max_alpha2[0] = dy_dt[0]
                    with hcl.if_(dtheta_dt[0] > max_alpha3):
                        max_alpha3[0] = dtheta_dt[0]

    hcl.update(t, lambda x: step_bound(t))
    result = hcl.update(V_new, lambda i,j,k: V_init[i,j,k] + V_new[i,j,k] * t[0])

    return result

def main():
    hcl.init()
    hcl.config.init_dtype = hcl.Float()
    V_f = hcl.placeholder((100, 100, 100), name="V_f", dtype = hcl.Float())
    V_init = hcl.placeholder((100, 100, 100), name="V_init", dtype=hcl.Float())
    thetas = hcl.placeholder((100,), name="thetas", dtype=hcl.Float())
    t = hcl.placeholder((1,), name="t", dtype=hcl.Float())

    # Create schedule
    s = hcl.create_schedule([V_f, V_init, thetas,t], HJ_PDE_solver)

    # Here comes the optimization

    # Accessing the hamiltonian stage
    s_H = HJ_PDE_solver.Hamiltonian
    # Split the loops
    k_out, k_in = s[s_H].split(s_H.k, 4) # These numbers are experimental, changable
    j_out, j_in = s[s_H].split(s_H.j, 4)
    i_out, i_in = s[s_H].split(s_H.i, 4)

    # Reorder the loops
    s[s_H].reorder(j_out, k_in)
    s[s_H].reorder(i_out, k_in)
    s[s_H].reorder(k_in, j_in)

    # FPGA Back end - parallel specs
    s[s_H].pipeline(k_in)
    s[s_H].unroll(i_out, 5)

    # If CPU option
    s[s_H].parallel(k_out)

    # Inspect IR
    #print(hcl.lower(s))

    # Build the code - CPU Back end
    solve_pde = hcl.build(s)
    #print(f)

    # Prepare numpy array for graph computation
    V_0 = hcl.asarray(np.random.rand(100,100,100))
    V_1=  hcl.asarray(np.random.rand(100,100,100))
    t_minh = hcl.asarray(np.zeros(1))
    list_theta = np.zeros(100)

    for i in range(1,100):
        list_theta[i] = list_theta[i-1] + dx[2]

    list_theta = hcl.asarray(list_theta)


    # Variables used for timing
    total_time = 0
    # Test the executable from heteroCL:
    for i in range(0,20):
        start = time.time()
        solve_pde(V_1, V_0, list_theta, t_minh)
        if i != 0: # Exclude first time of the computation
            total_time += time.time() - start
        print("Kernel time (s): {:.5f}".format(time.time() - start))
        V_0 = V_1

    print("Average kernel time (s): {:.5f}".format(total_time / 19))
    print("Finished solving\n")

if __name__ == '__main__':
  main()
