#pragma once

#include "OpenMeshCraft/NumberTypes/NumberTypesDecl.h"

#include "OpenMeshCraft/Utils/Macros.h"

namespace OMC {

inline bool lambda2d_SSI_filtered(double ea1x, double ea1y, double ea2x,
                                  double ea2y, double eb1x, double eb1y,
                                  double eb2x, double eb2y, double &lambda_x,
                                  double &lambda_y, double &lambda_det,
                                  double &max_var);

template <typename _IT>
bool lambda2d_SSI_interval(_IT ea1x, _IT ea1y, _IT ea2x, _IT ea2y, _IT eb1x,
                           _IT eb1y, _IT eb2x, _IT eb2y, _IT &lambda_x,
                           _IT &lambda_y, _IT &lambda_det);

template <typename _ET>
void lambda2d_SSI_exact(const _ET &ea1x, const _ET &ea1y, const _ET &ea2x,
                        const _ET &ea2y, const _ET &eb1x, const _ET &eb1y,
                        const _ET &eb2x, const _ET &eb2y, _ET &lambda_x,
                        _ET &lambda_y, _ET &lambda_det);

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

template <typename _IT>
bool lambda3d_LPI_interval(_IT px, _IT py, _IT pz, _IT qx, _IT qy, _IT qz,
                           _IT rx, _IT ry, _IT rz, _IT sx, _IT sy, _IT sz,
                           _IT tx, _IT ty, _IT tz, _IT &lambda_d, _IT &lambda_x,
                           _IT &lambda_y, _IT &lambda_z);

template <typename _ET>
void lambda3d_LPI_exact(const _ET &px, const _ET &py, const _ET &pz,
                        const _ET &qx, const _ET &qy, const _ET &qz,
                        const _ET &rx, const _ET &ry, const _ET &rz,
                        const _ET &sx, const _ET &sy, const _ET &sz,
                        const _ET &tx, const _ET &ty, const _ET &tz,
                        _ET &lambda_d, _ET &lambda_x, _ET &lambda_y,
                        _ET &lambda_z);

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

template <typename _IT>
bool lambda3d_TPI_interval(_IT ov1x, _IT ov1y, _IT ov1z, _IT ov2x, _IT ov2y,
                           _IT ov2z, _IT ov3x, _IT ov3y, _IT ov3z, _IT ow1x,
                           _IT ow1y, _IT ow1z, _IT ow2x, _IT ow2y, _IT ow2z,
                           _IT ow3x, _IT ow3y, _IT ow3z, _IT ou1x, _IT ou1y,
                           _IT ou1z, _IT ou2x, _IT ou2y, _IT ou2z, _IT ou3x,
                           _IT ou3y, _IT ou3z, _IT &lambda_x, _IT &lambda_y,
                           _IT &lambda_z, _IT &lambda_d);

template <typename _ET>
void lambda3d_TPI_exact(const _ET &ov1x, const _ET &ov1y, const _ET &ov1z,
                        const _ET &ov2x, const _ET &ov2y, const _ET &ov2z,
                        const _ET &ov3x, const _ET &ov3y, const _ET &ov3z,
                        const _ET &ow1x, const _ET &ow1y, const _ET &ow1z,
                        const _ET &ow2x, const _ET &ow2y, const _ET &ow2z,
                        const _ET &ow3x, const _ET &ow3y, const _ET &ow3z,
                        const _ET &ou1x, const _ET &ou1y, const _ET &ou1z,
                        const _ET &ou2x, const _ET &ou2y, const _ET &ou2z,
                        const _ET &ou3x, const _ET &ou3y, const _ET &ou3z,
                        _ET &lambda_x, _ET &lambda_y, _ET &lambda_z,
                        _ET &lambda_d);

inline void lambda3d_TPI_expansion(
  double ov1x, double ov1y, double ov1z, double ov2x, double ov2y, double ov2z,
  double ov3x, double ov3y, double ov3z, double ow1x, double ow1y, double ow1z,
  double ow2x, double ow2y, double ow2z, double ow3x, double ow3y, double ow3z,
  double ou1x, double ou1y, double ou1z, double ou2x, double ou2y, double ou2z,
  double ou3x, double ou3y, double ou3z, double **lambda_x, int &lambda_x_len,
  double **lambda_y, int &lambda_y_len, double **lambda_z, int &lambda_z_len,
  double **lambda_d, int &lambda_d_len);

template <typename _IT>
bool lambda3d_LNC_interval(_IT px, _IT py, _IT pz, _IT qx, _IT qy, _IT qz,
                           _IT t, _IT &lambda_x, _IT &lambda_y, _IT &lambda_z,
                           _IT &lambda_d);

template <typename _ET>
void lambda3d_LNC_exact(_ET px, _ET py, _ET pz, _ET qx, _ET qy, _ET qz, _ET t,
                        _ET &lambda_x, _ET &lambda_y, _ET &lambda_z,
                        _ET &lambda_d);

inline void lambda3d_LNC_expansion(double px, double py, double pz, double qx,
                                   double qy, double qz, double t,
                                   double **lambda_x, int &lambda_x_len,
                                   double **lambda_y, int &lambda_y_len,
                                   double **lambda_z, int &lambda_z_len,
                                   double **lambda_d, int &lambda_d_len);

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "ImplicitPointPredicates.inl"
#endif