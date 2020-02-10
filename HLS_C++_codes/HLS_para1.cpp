#include <ap_int.h>
#include <ap_fixed.h>
#include <math.h>

void default_function(float V_f[100][100][100], float V_init[100][100][100], float thetas[100], float t[1]) {
  float max_alpha2;
  max_alpha2 = -1.000000e+09f;
  float max_alpha3;
  max_alpha3 = -1.000000e+09f;
  float max_alpha1;
  max_alpha1 = -1.000000e+09f;
  float Hamiltonian;
  for (ap_int<32> k_outer = 0; k_outer < 10; ++k_outer) {
    for (ap_int<32> j_outer = 0; j_outer < 10; ++j_outer) {
      for (ap_int<32> i_outer = 0; i_outer < 10; ++i_outer) {
      #pragma HLS unroll factor=5
        for (ap_int<32> k_inner = 0; k_inner < 10; ++k_inner) {
        #pragma HLS pipeline
          for (ap_int<32> j_inner = 0; j_inner < 10; ++j_inner) {
            for (ap_int<32> i_inner = 0; i_inner < 10; ++i_inner) {
              float dV_dx_L;
              dV_dx_L = 0.000000e+00f;
              float dV_dx_R;
              dV_dx_R = 0.000000e+00f;
              float dV_dx;
              dV_dx = 0.000000e+00f;
              float dV_dy_L;
              dV_dy_L = 0.000000e+00f;
              float dV_dy_R;
              dV_dy_R = 0.000000e+00f;
              float dV_dy;
              dV_dy = 0.000000e+00f;
              float dV_dT_L;
              dV_dT_L = 0.000000e+00f;
              float dV_dT_R;
              dV_dT_R = 0.000000e+00f;
              float dV_dT;
              dV_dT = 0.000000e+00f;
              float dx_dt;
              dx_dt = 0.000000e+00f;
              float dy_dt;
              dy_dt = 0.000000e+00f;
              float dtheta_dt;
              dtheta_dt = 0.000000e+00f;
              float left_deriv;
              left_deriv = 0.000000e+00f;
              float right_deriv;
              right_deriv = 0.000000e+00f;
              if ((i_inner + (i_outer * 10)) == 0) {
                float left_boundary;
                left_boundary = 0.000000e+00f;
                float abs_value;
                abs_value = 0.000000e+00f;
                if (0.000000e+00f < (V_init[1][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))] - V_init[0][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))])) {
                  abs_value = (V_init[1][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))] - V_init[0][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))]);
                } else {
                  abs_value = (V_init[0][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))] - V_init[1][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))]);
                }
                float sign;
                sign = 0.000000e+00f;
                if (V_init[0][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))] == 0.000000e+00f) {
                  sign = 0.000000e+00f;
                }
                if (0.000000e+00f < V_init[0][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))]) {
                  sign = 1.000000e+00f;
                }
                if (V_init[0][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))] < 0.000000e+00f) {
                  sign = -1.000000e+00f;
                }
                left_boundary = (V_init[0][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))] + (abs_value * sign));
                left_deriv = ((V_init[0][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))] - left_boundary) * 9.900001e+00f);
                right_deriv = ((V_init[1][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))] - V_init[0][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))]) * 9.900001e+00f);
              } else {
                if ((i_inner + (i_outer * 10)) == 99) {
                  float right_boundary;
                  right_boundary = 0.000000e+00f;
                  float abs_value1;
                  abs_value1 = 0.000000e+00f;
                  if (0.000000e+00f < (V_init[99][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))] - V_init[98][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))])) {
                    abs_value1 = (V_init[99][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))] - V_init[98][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))]);
                  } else {
                    abs_value1 = (V_init[98][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))] - V_init[99][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))]);
                  }
                  float sign1;
                  sign1 = 0.000000e+00f;
                  if (V_init[99][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))] == 0.000000e+00f) {
                    sign1 = 0.000000e+00f;
                  }
                  if (0.000000e+00f < V_init[99][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))]) {
                    sign1 = 1.000000e+00f;
                  }
                  if (V_init[99][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))] < 0.000000e+00f) {
                    sign1 = -1.000000e+00f;
                  }
                  right_boundary = (V_init[99][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))] + (abs_value1 * sign1));
                  left_deriv = ((V_init[99][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))] - V_init[98][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))]) * 9.900001e+00f);
                  right_deriv = ((right_boundary - V_init[99][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))]) * 9.900001e+00f);
                } else {
                  if ((i_inner + (i_outer * 10)) != 99) {
                    left_deriv = ((V_init[(i_inner + (i_outer * 10))][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))] - V_init[((i_inner + (i_outer * 10)) + -1)][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))]) * 9.900001e+00f);
                    right_deriv = ((V_init[((i_inner + (i_outer * 10)) + 1)][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))] - V_init[(i_inner + (i_outer * 10))][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))]) * 9.900001e+00f);
                  }
                }
              }
              dV_dx_L = left_deriv;
              dV_dx_R = right_deriv;
              float left_deriv1;
              left_deriv1 = 0.000000e+00f;
              float right_deriv1;
              right_deriv1 = 0.000000e+00f;
              if ((j_inner + (j_outer * 10)) == 0) {
                float left_boundary1;
                left_boundary1 = 0.000000e+00f;
                float abs_value2;
                abs_value2 = 0.000000e+00f;
                if (0.000000e+00f < (V_init[(i_inner + (i_outer * 10))][1][(k_inner + (k_outer * 10))] - V_init[(i_inner + (i_outer * 10))][0][(k_inner + (k_outer * 10))])) {
                  abs_value2 = (V_init[(i_inner + (i_outer * 10))][1][(k_inner + (k_outer * 10))] - V_init[(i_inner + (i_outer * 10))][0][(k_inner + (k_outer * 10))]);
                } else {
                  abs_value2 = (V_init[(i_inner + (i_outer * 10))][0][(k_inner + (k_outer * 10))] - V_init[(i_inner + (i_outer * 10))][1][(k_inner + (k_outer * 10))]);
                }
                float sign2;
                sign2 = 0.000000e+00f;
                if (V_init[(i_inner + (i_outer * 10))][0][(k_inner + (k_outer * 10))] == 0.000000e+00f) {
                  sign2 = 0.000000e+00f;
                }
                if (0.000000e+00f < V_init[(i_inner + (i_outer * 10))][0][(k_inner + (k_outer * 10))]) {
                  sign2 = 1.000000e+00f;
                }
                if (V_init[(i_inner + (i_outer * 10))][0][(k_inner + (k_outer * 10))] < 0.000000e+00f) {
                  sign2 = -1.000000e+00f;
                }
                left_boundary1 = (V_init[(i_inner + (i_outer * 10))][0][(k_inner + (k_outer * 10))] + (abs_value2 * sign2));
                left_deriv1 = ((V_init[(i_inner + (i_outer * 10))][0][(k_inner + (k_outer * 10))] - left_boundary1) * 9.900001e+00f);
                right_deriv1 = ((V_init[(i_inner + (i_outer * 10))][1][(k_inner + (k_outer * 10))] - V_init[(i_inner + (i_outer * 10))][0][(k_inner + (k_outer * 10))]) * 9.900001e+00f);
              } else {
                if ((j_inner + (j_outer * 10)) == 99) {
                  float right_boundary1;
                  right_boundary1 = 0.000000e+00f;
                  float abs_value3;
                  abs_value3 = 0.000000e+00f;
                  if (0.000000e+00f < (V_init[(i_inner + (i_outer * 10))][99][(k_inner + (k_outer * 10))] - V_init[(i_inner + (i_outer * 10))][98][(k_inner + (k_outer * 10))])) {
                    abs_value3 = (V_init[(i_inner + (i_outer * 10))][99][(k_inner + (k_outer * 10))] - V_init[(i_inner + (i_outer * 10))][98][(k_inner + (k_outer * 10))]);
                  } else {
                    abs_value3 = (V_init[(i_inner + (i_outer * 10))][98][(k_inner + (k_outer * 10))] - V_init[(i_inner + (i_outer * 10))][99][(k_inner + (k_outer * 10))]);
                  }
                  float sign3;
                  sign3 = 0.000000e+00f;
                  if (V_init[(i_inner + (i_outer * 10))][99][(k_inner + (k_outer * 10))] == 0.000000e+00f) {
                    sign3 = 0.000000e+00f;
                  }
                  if (0.000000e+00f < V_init[(i_inner + (i_outer * 10))][99][(k_inner + (k_outer * 10))]) {
                    sign3 = 1.000000e+00f;
                  }
                  if (V_init[(i_inner + (i_outer * 10))][99][(k_inner + (k_outer * 10))] < 0.000000e+00f) {
                    sign3 = -1.000000e+00f;
                  }
                  right_boundary1 = (V_init[(i_inner + (i_outer * 10))][99][(k_inner + (k_outer * 10))] + (abs_value3 * sign3));
                  left_deriv1 = ((V_init[(i_inner + (i_outer * 10))][99][(k_inner + (k_outer * 10))] - V_init[(i_inner + (i_outer * 10))][98][(k_inner + (k_outer * 10))]) * 9.900001e+00f);
                  right_deriv1 = ((right_boundary1 - V_init[(i_inner + (i_outer * 10))][99][(k_inner + (k_outer * 10))]) * 9.900001e+00f);
                } else {
                  if ((j_inner + (j_outer * 10)) != 99) {
                    left_deriv1 = ((V_init[(i_inner + (i_outer * 10))][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))] - V_init[(i_inner + (i_outer * 10))][((j_inner + (j_outer * 10)) + -1)][(k_inner + (k_outer * 10))]) * 9.900001e+00f);
                    right_deriv1 = ((V_init[(((((j_inner - (((j_inner + (j_outer * 10)) + 1) % 100)) + (j_outer * 10)) + ((i_inner + (i_outer * 10)) * 100)) + 1) / 100)][(((j_inner + (j_outer * 10)) + 1) % 100)][(k_inner + (k_outer * 10))] - V_init[(i_inner + (i_outer * 10))][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))]) * 9.900001e+00f);
                  }
                }
              }
              dV_dy_L = left_deriv1;
              dV_dy_R = right_deriv1;
              float left_deriv2;
              left_deriv2 = 0.000000e+00f;
              float right_deriv2;
              right_deriv2 = 0.000000e+00f;
              if ((k_inner + (k_outer * 10)) == 0) {
                float left_boundary2;
                left_boundary2 = 0.000000e+00f;
                left_boundary2 = V_init[(i_inner + (i_outer * 10))][(j_inner + (j_outer * 10))][99];
                left_deriv2 = ((V_init[(i_inner + (i_outer * 10))][(j_inner + (j_outer * 10))][0] - left_boundary2) * 1.591549e+01f);
                right_deriv2 = ((V_init[(i_inner + (i_outer * 10))][(j_inner + (j_outer * 10))][1] - V_init[(i_inner + (i_outer * 10))][(j_inner + (j_outer * 10))][0]) * 1.591549e+01f);
              } else {
                if ((k_inner + (k_outer * 10)) == 99) {
                  float right_boundary2;
                  right_boundary2 = 0.000000e+00f;
                  right_boundary2 = V_init[(i_inner + (i_outer * 10))][(j_inner + (j_outer * 10))][0];
                  left_deriv2 = ((V_init[(i_inner + (i_outer * 10))][(j_inner + (j_outer * 10))][99] - V_init[(i_inner + (i_outer * 10))][(j_inner + (j_outer * 10))][98]) * 1.591549e+01f);
                  right_deriv2 = ((right_boundary2 - V_init[(i_inner + (i_outer * 10))][(j_inner + (j_outer * 10))][99]) * 1.591549e+01f);
                } else {
                  if ((k_inner + (k_outer * 10)) != 99) {
                    left_deriv2 = ((V_init[(i_inner + (i_outer * 10))][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))] - V_init[(i_inner + (i_outer * 10))][(j_inner + (j_outer * 10))][((k_inner + (k_outer * 10)) + -1)]) * 1.591549e+01f);
                    right_deriv2 = ((V_init[((((((k_inner - (((k_inner + (k_outer * 10)) + 1) % 100)) + (k_outer * 10)) + ((j_inner + (j_outer * 10)) * 100)) + ((i_inner + (i_outer * 10)) * 10000)) + 1) / 10000)][(((((((k_inner - (((k_inner + (k_outer * 10)) + 1) % 100)) + (k_outer * 10)) + ((j_inner + (j_outer * 10)) * 100)) + ((i_inner + (i_outer * 10)) * 10000)) + 1) / 100) % 100)][(((k_inner + (k_outer * 10)) + 1) % 100)] - V_init[(i_inner + (i_outer * 10))][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))]) * 1.591549e+01f);
                  }
                }
              }
              dV_dT_L = left_deriv2;
              dV_dT_R = right_deriv2;
              dV_dx = ((dV_dx_L + dV_dx_R) * 5.000000e-01f);
              dV_dy = ((dV_dy_L + dV_dy_R) * 5.000000e-01f);
              dV_dT = ((dV_dT_L + dV_dT_R) * 5.000000e-01f);
              float uOpt;
              uOpt = 1.000000e+00f;
              float vel;
              vel = 1.000000e+00f;
              if (0.000000e+00f < dV_dT) {
                uOpt = (uOpt * -1.000000e+00f);
              }
              dx_dt = ((float)(((double)vel) * cos(((double)thetas[(k_inner + (k_outer * 10))]))));
              dy_dt = ((float)(((double)vel) * sin(((double)thetas[(k_inner + (k_outer * 10))]))));
              dtheta_dt = uOpt;
              V_f[(i_inner + (i_outer * 10))][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))] = ((((dx_dt * dV_dx) + (dy_dt * dV_dy)) + (dtheta_dt * dV_dT)) * -1.000000e+00f);
              float abs_value4;
              abs_value4 = 0.000000e+00f;
              if (0.000000e+00f < dx_dt) {
                abs_value4 = dx_dt;
              } else {
                abs_value4 = (dx_dt * -1.000000e+00f);
              }
              dx_dt = abs_value4;
              float abs_value5;
              abs_value5 = 0.000000e+00f;
              if (0.000000e+00f < dy_dt) {
                abs_value5 = dy_dt;
              } else {
                abs_value5 = (dy_dt * -1.000000e+00f);
              }
              dy_dt = abs_value5;
              float abs_value6;
              abs_value6 = 0.000000e+00f;
              if (0.000000e+00f < dtheta_dt) {
                abs_value6 = dtheta_dt;
              } else {
                abs_value6 = (dtheta_dt * -1.000000e+00f);
              }
              dtheta_dt = abs_value6;
              float diss;
              diss = 0.000000e+00f;
              diss = (((((dV_dx_R - dV_dx_L) * dx_dt) + ((dV_dy_R - dV_dy_L) * dy_dt)) + ((dV_dT_R - dV_dT_L) * dtheta_dt)) * 5.000000e-01f);
              V_f[(i_inner + (i_outer * 10))][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))] = (diss - V_f[(i_inner + (i_outer * 10))][(j_inner + (j_outer * 10))][(k_inner + (k_outer * 10))]);
              if (max_alpha1 < dx_dt) {
                max_alpha1 = dx_dt;
              }
              if (max_alpha2 < dy_dt) {
                max_alpha2 = dy_dt;
              }
              if (max_alpha3 < dtheta_dt) {
                max_alpha3 = dtheta_dt;
              }
            }
          }
        }
      }
    }
  }
  float update0;
  float stepBoundInv;
  stepBoundInv = 0.000000e+00f;
  float stepBound;
  stepBound = 0.000000e+00f;
  stepBoundInv = (((max_alpha1 + max_alpha2) * 9.900001e+00f) + (max_alpha3 * 1.591549e+01f));
  stepBound = (8.000000e-01f / stepBoundInv);
  if (5.000000e-02f < stepBound) {
    stepBound = 5.000000e-02f;
  }
  t[0] = stepBound;
  float update1;
  for (ap_int<32> i = 0; i < 100; ++i) {
    for (ap_int<32> j = 0; j < 100; ++j) {
      for (ap_int<32> k = 0; k < 100; ++k) {
        V_f[i][j][k] = (V_init[i][j][k] + (V_f[i][j][k] * t[0]));
      }
    }
  }
  float update2;
  for (ap_int<32> i1 = 0; i1 < 100; ++i1) {
    for (ap_int<32> j1 = 0; j1 < 100; ++j1) {
      for (ap_int<32> k1 = 0; k1 < 100; ++k1) {
        V_init[i1][j1][k1] = V_f[i1][j1][k1];
      }
    }
  }
}


