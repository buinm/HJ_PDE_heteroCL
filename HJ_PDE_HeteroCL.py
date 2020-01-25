import heterocl as hcl
import numpy as np

import math

hcl.config.init_dtype = hcl.Float()

class grid:
  def __init__(self, max, min, pts_each_dim, pDim):
      self.max = max
      self.min = min
      self.pts_each_dim = pts_each_dim
      self.pDim = pDim
      self.dx = (max - min)/(pts_each_dim - 1.0)


# NOTE: No information about structure of grid, dynamics of system passed in for now
# Hardcoded these information

# Global variables
def cleaner_HJ_PDE_solver(V_new, V_init, thetas, dx):

    # Calculate dV_dx
    dV_dx_L = hcl.compute((50, 50, 50), lambda i, j, k: spa_derivXL(i,j,k), name = "dV_dx_L")
    dV_dx_R = hcl.compute((50, 50, 50), lambda i, j, k: spa_derivXR(i, j, k), name = "dV_dx_R")
    dV_dx   = hcl.compute((50, 50, 50), lambda i, j, k: (dV_dx_L[i,j,k] + dV_dx_R[i,j,k])/2, name = "dV_dx")

    # Calculate dV_dy
    dV_dy_L = hcl.compute((50, 50, 50), lambda i, j, k: spa_derivYL(i, j, k), name = "dV_dy_L")
    dV_dy_R = hcl.compute((50, 50, 50), lambda i, j, k: spa_derivYR(i, j, k), name = "dV_dy_R")
    dV_dy   = hcl.compute((50, 50, 50), lambda i, j, k: (dV_dy_L[i,j,k] + dV_dy_R[i,j,k])/2, name = "dV_dy")

    # Calculate dV_dT
    dV_dT_L = hcl.compute((50, 50, 50), lambda i, j, k: spa_derivTL(i, j, k), name="dV_dT_L")
    dV_dT_R = hcl.compute((50, 50, 50), lambda i, j, k: spa_derivTR(i, j, k), name="dV_dT_R")
    dV_dT   = hcl.compute((50, 50, 50), lambda i, j, k: (dV_dy_L[i, j, k] + dV_dy_R[i, j, k]) / 2, name="dV_dT")



def HJ_PDE_solver(V_new, V_init, thetas, dx):

    # HeteroC
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
            left_deriv[0] = (V_init[i,j,k] - left_boundary[0])/dx[1]
            right_deriv[0] = (V_init[i+1,j,k] - V_init[i,j,k])/dx[1]
        with hcl.elif_(i == V_init.shape[0]):
            right_boundary = hcl.scalar(0, "right_boundary")
            right_boundary[0] = V_init[i,j,k] +my_abs(V_init[i,j,k] - V_init[i-1,j,k]) * my_sign(V_init[i,j,k])
            left_deriv[0] = (V_init[i,j,k] - V_init[i-1,j,k])/dx[1]
            right_deriv[0] = (right_boundary[0] -V_init[i,j,k])/dx[1]
        with hcl.elif_(i != 0 and i != V_init.shape[0]):
            left_deriv[0] = (V_init[i, j, k] - V_init[i - 1, j, k]) / dx[1]
            right_deriv[0] = (V_init[i + 1, j, k] - V_init[i, j, k]) / dx[1]
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
        with hcl.elif_(j == V_init.shape[1]):
            right_boundary = hcl.scalar(0, "right_boundary")
            right_boundary[0] = V_init[i, j, k] + my_abs(V_init[i, j, k] - V_init[i, j - 1, k]) * my_sign(
                V_init[i, j, k])
            left_deriv[0] = (V_init[i, j, k] - V_init[i, j - 1, k]) / dx[1]
            right_deriv[0] = (right_boundary[0] - V_init[i, j, k]) / dx[1]
        with hcl.elif_(j != 0 and j != V_init.shape[1]):
            left_deriv[0] = (V_init[i, j, k] - V_init[i, j - 1, k]) / dx[1]
            right_deriv[0] = (V_init[i, j + 1, k] - V_init[i, j, k]) / dx[1]
        return left_deriv[0], right_deriv[0]


    def spa_derivT(i,j,k):
        left_deriv = hcl.scalar(0, "left_deriv")
        right_deriv = hcl.scalar(0, "right_deriv")
        with hcl.if_(k == 0):
            left_boundary = hcl.scalar(0, "left_boundary")
            left_boundary[0] = V_init[i,j,50]
            left_deriv[0] = (V_init[i, j, k] - left_boundary[0]) / dx[1]
            right_deriv[0] = (V_init[i, j, k + 1] - V_init[i, j, k]) / dx[1]
        with hcl.elif_(k == V_init.shape[2]):
            right_boundary = hcl.scalar(0, "right_boundary")
            right_boundary[0] = V_init[i,j,0]
            left_deriv[0] = (V_init[i, j, k] - V_init[i, j, k - 1]) / dx[1]
            right_deriv[0] = (right_boundary[0] - V_init[i, j, k]) / dx[1]
        with hcl.elif_(j != 0 and j != V_init.shape[2]):
            left_deriv[0] = (V_init[i, j, k] - V_init[i, j, k - 1]) / dx[1]
            right_deriv[0] = (V_init[i, j ,k + 1] - V_init[i, j, k]) / dx[1]
        return left_deriv[0], right_deriv[0]

    # Calculate Hamiltonian for every grid point in V_init
    with hcl.Stage("Hamiltonian"):
        with hcl.for_(0, V_init.shape[0] + 1, name="i") as i: # Plus 1 as for loop count stops at V_init.shape[0]
            with hcl.for_(0, V_init.shape[1] + 1, name="j") as j:
                with hcl.for_(0, V_init.shape[2] + 1, name="k") as k:
                    # Calculate dV_dx
                    dV_dx_L = hcl.scalar(0, "dV_dx_L")
                    dV_dx_R = hcl.scalar(0, "dV_dx_R")
                    dV_dx   = hcl.scalar(0, "dV_dx")

                    # Calculate dV_dy
                    dV_dy_L = hcl.scalar(0, "dV_dy_L")
                    dV_dy_R = hcl.scalar(0, "dV_dy_R")
                    dV_dy   = hcl.scalar(0, "dV_dy")

                    # Calculate dV_dtheta
                    dV_dT_L = hcl.scalar(0, "dV_dT_L")
                    dV_dT_R = hcl.scalar(0, "dV_dT_R")
                    dV_dT   = hcl.scalar(0, "dV_dT")

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

                    # Calculate dynamical changes
                    V_new[i, j, k] =  vel.v * hcl.cos(thetas[k]) * dV_dx[0] + vel.v * hcl.sin(thetas[k]) * dV_dy[0] + uOpt.v * dV_dT[0]

def main():
    hcl.init()
    V_f = hcl.placeholder((50, 50, 50), name="V_f", dtype = hcl.Float())
    V_init = hcl.placeholder((50, 50, 50), name="V_init", dtype=hcl.Float())
    thetas = hcl.placeholder((50,), name="thetas", dtype=hcl.Float())
    dx = hcl.placeholder((3,), name="dx", dtype=hcl.Float())

    # Create schedule
    s = hcl.create_schedule([V_f, V_init, thetas, dx], HJ_PDE_solver)

    # Inspect IR
    print(hcl.lower(s))

if __name__ == '__main__':
  main()