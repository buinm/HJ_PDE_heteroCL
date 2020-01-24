import heterocl as hcl
import numpy as np

import math


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

def HJ_PDE_solver(V_new, V_init, thetas, dx):
    hcl.config.init_dtype = hcl.Float()

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

    # Calculate spatial derivative based on index and dimension number
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

            #return left_deriv[0],right_deriv[0]

            #return left_deriv[0], right_deriv[0]
        return left_deriv[0], right_deriv[0]

    def spa_derivY(i,j,k):
        return i + j +

    def spa_derivTheta(i,j,k):
        return i - j + k, k*j

    # Calculate Hamiltonian for every grid point in V_init
    with hcl.Stage("Hamiltonian"):
        with hcl.for_(0, V_init.shape[0] + 1, name="i") as i:
            with hcl.for_(0, V_init.shape[1] + 1, name="j") as j:
                with hcl.for_(0, V_init.shape[2] + 1, name="k") as k:
                    # Calculate dV_dx
                    dV_dx_L = hcl.scalar(0, "dV_dx_L")
                    dV_dx_R = hcl.scalar(0, "dV_dx_R")
                    dV_dx_C = hcl.scalar(0, "dV_dx_C")
                    dV_dx_L[0], dV_dx_R[0] = spa_derivX(i,j,k)
                    #dV_dy_L, dV_dy_R = spa_derivY(i,j,k)
                    #dV_dtheta_L, dV_dtheta_R = spa_derivTheta(i, j, k)

                    # Calculate average gradient
                    dV_dx_C[0] = (dV_dx_L[0] + dV_dx_R[0])/2
                    #dV_dy_C = (dV_dy_L + dV_dy_R) / 2
                    #dV_dtheta_C = (dV_dtheta_L + dV_dtheta_R) / 2

                    # Declare optimal control
                    uOpt = hcl.scalar(0, "uOpt")

                    # Declare Velocity
                    vel = hcl.scalar(1,"vel")

                    # Assume that mode is min
                    #with hcl.if_(dV_dtheta_C > 0):
                    #    uOpt.v = -uOpt.v

                    # Calculate dynamics function
                    #V_new[i,j,k] = 1 * cos(thetas[k]) * dV_dx_C +1 * sin(thetas[k]) * dV_dy_C +uOpt * dV_theta_C
                    #angle = hcl.scalar(thetas[k], "angle")
                    #V_new[i,j,k] = v * hcl.cos(thetas[k]) * dV_dx_C + v * hcl.sin(thetas[k]) * dV_dy_C +  dV_dtheta_C * uOpt
                    V_new[i, j, k] =  vel.v * hcl.cos(thetas[k]) * dV_dx_C[0]

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