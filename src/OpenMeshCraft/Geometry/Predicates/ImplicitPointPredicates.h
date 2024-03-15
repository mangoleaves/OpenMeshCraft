#pragma once

#include "OpenMeshCraft/NumberTypes/NumberTypesDecl.h"

#include "OpenMeshCraft/Utils/Macros.h"

namespace OMC {

#if defined(INDIRECT_PREDICATES)

inline bool lambda2d_SSI_filtered(double ea1x, double ea1y, double ea2x,
                                  double ea2y, double eb1x, double eb1y,
                                  double eb2x, double eb2y, double &lambda_x,
                                  double &lambda_y, double &lambda_det,
                                  double &max_var);

template <typename IT>
bool lambda2d_SSI_interval(IT ea1x, IT ea1y, IT ea2x, IT ea2y, IT eb1x, IT eb1y,
                           IT eb2x, IT eb2y, IT &lambda_x, IT &lambda_y,
                           IT &lambda_det);

template <typename ET>
void lambda2d_SSI_exact(const ET &ea1x, const ET &ea1y, const ET &ea2x,
                        const ET &ea2y, const ET &eb1x, const ET &eb1y,
                        const ET &eb2x, const ET &eb2y, ET &lambda_x,
                        ET &lambda_y, ET &lambda_det);

inline void lambda2d_SSI_expansion(double ea1x, double ea1y, double ea2x,
                                   double ea2y, double eb1x, double eb1y,
                                   double eb2x, double eb2y, double **lambda_x,
                                   int &lambda_x_len, double **lambda_y,
                                   int &lambda_y_len, double **lambda_det,
                                   int &lambda_det_len);

inline bool lambda3d_SSI_filtered(double xa, double ya, double za, double xb,
                                  double yb, double zb, double xp, double yp,
                                  double xq, double yq, double &lambda_x,
                                  double &lambda_y, double &lambda_z,
                                  double &lambda_d, double &max_var);

template <typename IT>
bool lambda3d_SSI_interval(IT xa, IT ya, IT za, IT xb, IT yb, IT zb, IT xp,
                           IT yp, IT xq, IT yq, IT &lambda_x, IT &lambda_y,
                           IT &lambda_z, IT &lambda_d);

template <typename ET>
void lambda3d_SSI_exact(ET xa, ET ya, ET za, ET xb, ET yb, ET zb, ET xp, ET yp,
                        ET xq, ET yq, ET &lambda_x, ET &lambda_y, ET &lambda_z,
                        ET &lambda_d);

inline void lambda3d_SSI_expansion(double xa, double ya, double za, double xb,
                                   double yb, double zb, double xp, double yp,
                                   double xq, double yq, double **lambda_x,
                                   int &lambda_x_len, double **lambda_y,
                                   int &lambda_y_len, double **lambda_z,
                                   int &lambda_z_len, double **lambda_d,
                                   int &lambda_d_len);

inline bool lambda3d_LPI_filtered(double px, double py, double pz, double qx,
                                  double qy, double qz, double rx, double ry,
                                  double rz, double sx, double sy, double sz,
                                  double tx, double ty, double tz,
                                  double &lambda_d, double &lambda_x,
                                  double &lambda_y, double &lambda_z,
                                  double &max_var);

template <typename IT>
bool lambda3d_LPI_interval(IT px, IT py, IT pz, IT qx, IT qy, IT qz, IT rx,
                           IT ry, IT rz, IT sx, IT sy, IT sz, IT tx, IT ty,
                           IT tz, IT &lambda_d, IT &lambda_x, IT &lambda_y,
                           IT &lambda_z);

template <typename ET>
void lambda3d_LPI_exact(const ET &px, const ET &py, const ET &pz, const ET &qx,
                        const ET &qy, const ET &qz, const ET &rx, const ET &ry,
                        const ET &rz, const ET &sx, const ET &sy, const ET &sz,
                        const ET &tx, const ET &ty, const ET &tz, ET &lambda_d,
                        ET &lambda_x, ET &lambda_y, ET &lambda_z);

inline void lambda3d_LPI_expansion(double px, double py, double pz, double qx,
                                   double qy, double qz, double rx, double ry,
                                   double rz, double sx, double sy, double sz,
                                   double tx, double ty, double tz,
                                   double **lambda_d, int &lambda_d_len,
                                   double **lambda_x, int &lambda_x_len,
                                   double **lambda_y, int &lambda_y_len,
                                   double **lambda_z, int &lambda_z_len);

inline bool lambda3d_TPI_filtered(
  double ov1x, double ov1y, double ov1z, double ov2x, double ov2y, double ov2z,
  double ov3x, double ov3y, double ov3z, double ow1x, double ow1y, double ow1z,
  double ow2x, double ow2y, double ow2z, double ow3x, double ow3y, double ow3z,
  double ou1x, double ou1y, double ou1z, double ou2x, double ou2y, double ou2z,
  double ou3x, double ou3y, double ou3z, double &lambda_x, double &lambda_y,
  double &lambda_z, double &lambda_d, double &max_var);

template <typename IT>
bool lambda3d_TPI_interval(IT ov1x, IT ov1y, IT ov1z, IT ov2x, IT ov2y, IT ov2z,
                           IT ov3x, IT ov3y, IT ov3z, IT ow1x, IT ow1y, IT ow1z,
                           IT ow2x, IT ow2y, IT ow2z, IT ow3x, IT ow3y, IT ow3z,
                           IT ou1x, IT ou1y, IT ou1z, IT ou2x, IT ou2y, IT ou2z,
                           IT ou3x, IT ou3y, IT ou3z, IT &lambda_x,
                           IT &lambda_y, IT &lambda_z, IT &lambda_d);

template <typename ET>
void lambda3d_TPI_exact(const ET &ov1x, const ET &ov1y, const ET &ov1z,
                        const ET &ov2x, const ET &ov2y, const ET &ov2z,
                        const ET &ov3x, const ET &ov3y, const ET &ov3z,
                        const ET &ow1x, const ET &ow1y, const ET &ow1z,
                        const ET &ow2x, const ET &ow2y, const ET &ow2z,
                        const ET &ow3x, const ET &ow3y, const ET &ow3z,
                        const ET &ou1x, const ET &ou1y, const ET &ou1z,
                        const ET &ou2x, const ET &ou2y, const ET &ou2z,
                        const ET &ou3x, const ET &ou3y, const ET &ou3z,
                        ET &lambda_x, ET &lambda_y, ET &lambda_z, ET &lambda_d);

inline void lambda3d_TPI_expansion(
  double ov1x, double ov1y, double ov1z, double ov2x, double ov2y, double ov2z,
  double ov3x, double ov3y, double ov3z, double ow1x, double ow1y, double ow1z,
  double ow2x, double ow2y, double ow2z, double ow3x, double ow3y, double ow3z,
  double ou1x, double ou1y, double ou1z, double ou2x, double ou2y, double ou2z,
  double ou3x, double ou3y, double ou3z, double **lambda_x, int &lambda_x_len,
  double **lambda_y, int &lambda_y_len, double **lambda_z, int &lambda_z_len,
  double **lambda_d, int &lambda_d_len);

template <typename IT>
bool lambda3d_LNC_interval(IT px, IT py, IT pz, IT qx, IT qy, IT qz, IT t,
                           IT &lambda_x, IT &lambda_y, IT &lambda_z,
                           IT &lambda_d);

template <typename ET>
void lambda3d_LNC_exact(ET px, ET py, ET pz, ET qx, ET qy, ET qz, ET t,
                        ET &lambda_x, ET &lambda_y, ET &lambda_z, ET &lambda_d);

inline void lambda3d_LNC_expansion(double px, double py, double pz, double qx,
                                   double qy, double qz, double t,
                                   double **lambda_x, int &lambda_x_len,
                                   double **lambda_y, int &lambda_y_len,
                                   double **lambda_z, int &lambda_z_len,
                                   double **lambda_d, int &lambda_d_len);

#elif defined(OFFSET_PREDICATES)

inline bool lambda3d_LPI_filtered(double xa, double ya, double za, double xb,
                                  double yb, double zb, double xo, double yo,
                                  double zo, double xp, double yp, double zp,
                                  double xq, double yq, double zq,
                                  double &lambda_d, double &lambda_x,
                                  double &lambda_y, double &lambda_z,
                                  double &beta_x, double &beta_y,
                                  double &beta_z, double &max_var);

template <typename IT>
bool lambda3d_LPI_interval(IT xa, IT ya, IT za, IT xb, IT yb, IT zb, IT xo,
                           IT yo, IT zo, IT xp, IT yp, IT zp, IT xq, IT yq,
                           IT zq, IT &lambda_d, IT &lambda_x, IT &lambda_y,
                           IT &lambda_z, IT &beta_x, IT &beta_y, IT &beta_z);

template <typename ET>
void lambda3d_LPI_exact(ET xa, ET ya, ET za, ET xb, ET yb, ET zb, ET xo, ET yo,
                        ET zo, ET xp, ET yp, ET zp, ET xq, ET yq, ET zq,
                        ET &lambda_d, ET &lambda_x, ET &lambda_y, ET &lambda_z,
                        ET &beta_x, ET &beta_y, ET &beta_z);

inline void lambda3d_LPI_expansion(
  double xa, double ya, double za, double xb, double yb, double zb, double xo,
  double yo, double zo, double xp, double yp, double zp, double xq, double yq,
  double zq, double **lambda_d, int &lambda_d_len, double **lambda_x,
  int &lambda_x_len, double **lambda_y, int &lambda_y_len, double **lambda_z,
  int &lambda_z_len, double &beta_x, double &beta_y, double &beta_z);

inline bool lambda3d_SSI_filtered(double xa, double ya, double za, double xb,
                                  double yb, double zb, double xp, double yp,
                                  double xq, double yq, double &lambda_d,
                                  double &lambda_x, double &lambda_y,
                                  double &lambda_z, double &beta_x,
                                  double &beta_y, double &beta_z,
                                  double &max_var);

template <typename IT>
bool lambda3d_SSI_interval(IT xa, IT ya, IT za, IT xb, IT yb, IT zb, IT xp,
                           IT yp, IT xq, IT yq, IT &lambda_d, IT &lambda_x,
                           IT &lambda_y, IT &lambda_z, IT &beta_x, IT &beta_y,
                           IT &beta_z);

template <typename ET>
void lambda3d_SSI_exact(ET xa, ET ya, ET za, ET xb, ET yb, ET zb, ET xp, ET yp,
                        ET xq, ET yq, ET &lambda_d, ET &lambda_x, ET &lambda_y,
                        ET &lambda_z, ET &beta_x, ET &beta_y, ET &beta_z);

inline void lambda3d_SSI_expansion(double xa, double ya, double za, double xb,
                                   double yb, double zb, double xp, double yp,
                                   double xq, double yq, double **lambda_d,
                                   int &lambda_d_len, double **lambda_x,
                                   int &lambda_x_len, double **lambda_y,
                                   int &lambda_y_len, double **lambda_z,
                                   int &lambda_z_len, double &beta_x,
                                   double &beta_y, double &beta_z);

inline bool lambda3d_TPI_filtered(
  double xa, double ya, double za, double xb, double yb, double zb, double xc,
  double yc, double zc, double xo, double yo, double zo, double xp, double yp,
  double zp, double xq, double yq, double zq, double xr, double yr, double zr,
  double xs, double ys, double zs, double xt, double yt, double zt,
  double &lambda_d, double &lambda_x, double &lambda_y, double &lambda_z,
  double &beta_x, double &beta_y, double &beta_z, double &max_var);

template <typename IT>
bool lambda3d_TPI_interval(IT xa, IT ya, IT za, IT xb, IT yb, IT zb, IT xc,
                           IT yc, IT zc, IT xo, IT yo, IT zo, IT xp, IT yp,
                           IT zp, IT xq, IT yq, IT zq, IT xr, IT yr, IT zr,
                           IT xs, IT ys, IT zs, IT xt, IT yt, IT zt,
                           IT &lambda_d, IT &lambda_x, IT &lambda_y,
                           IT &lambda_z, IT &beta_x, IT &beta_y, IT &beta_z);

template <typename ET>
void lambda3d_TPI_exact(ET xa, ET ya, ET za, ET xb, ET yb, ET zb, ET xc, ET yc,
                        ET zc, ET xo, ET yo, ET zo, ET xp, ET yp, ET zp, ET xq,
                        ET yq, ET zq, ET xr, ET yr, ET zr, ET xs, ET ys, ET zs,
                        ET xt, ET yt, ET zt, ET &lambda_d, ET &lambda_x,
                        ET &lambda_y, ET &lambda_z, ET &beta_x, ET &beta_y,
                        ET &beta_z);

inline void lambda3d_TPI_expansion(
  double xa, double ya, double za, double xb, double yb, double zb, double xc,
  double yc, double zc, double xo, double yo, double zo, double xp, double yp,
  double zp, double xq, double yq, double zq, double xr, double yr, double zr,
  double xs, double ys, double zs, double xt, double yt, double zt,
  double **lambda_d, int &lambda_d_len, double **lambda_x, int &lambda_x_len,
  double **lambda_y, int &lambda_y_len, double **lambda_z, int &lambda_z_len,
  double &beta_x, double &beta_y, double &beta_z);

#endif

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "ImplicitPointPredicates.inl"
#endif