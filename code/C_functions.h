#ifndef _C_functions_h_
#define _C_functions_h_

#ifdef __cplusplus
extern"C"{
#endif
  void SphDist(int* nrow, int* ncol, double* dat, double* sDmat);
  void Cov_mat(int* nrow, int* ncol, double* dat, int* nrow_Dmat, double* sDmat, double* a1, double* a2, double* a3, int* kappa, double*low_spherical, int* ncol_low_spherical, double* tau_spherical, double* R);
  void G_hat(int* nrow, int* ncol, double* dat, int* nrow_Dmat, double* sDmat, int* length_Hvec, double* Hvec, double* eps, int* t, double* G, double* N);
  //void svdcmp(double* U, double* D, double* V, double* R, int* nrow, int* ncol);
#ifdef __cplusplus
}
#endif
#endif
