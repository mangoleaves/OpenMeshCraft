#pragma once

#include "OpenMeshCraft/Geometry/Primitives/GenericPoint2T.h"
#include "OpenMeshCraft/Geometry/Primitives/GenericPoint3T.h"

#include "IndirectPredicateDetails.h"

#include "OpenMeshCraft/NumberTypes/NumberUtils.h"

namespace OMC {

/*************** Orient2D *******************/

inline Sign orient2d_filtered(double p1x, double p1y, double p2x, double p2y,
                              double p3x, double p3y);

inline Sign orient2d_expansion(double p1x, double p1y, double p2x, double p2y,
                               double p3x, double p3y);

inline Sign orient2d(double p1x, double p1y, double p2x, double p2y, double p3x,
                     double p3y);

inline Sign orient2d(const double *p1, const double *p2, const double *p3);

template <typename IT, typename ET>
Sign orient2d(const GenericPoint2T<IT, ET> &p1,
              const GenericPoint2T<IT, ET> &p2,
              const GenericPoint2T<IT, ET> &p3);

/*************** Orient3D *******************/

inline Sign orient3d_filtered(double px, double py, double pz, double qx,
                              double qy, double qz, double rx, double ry,
                              double rz, double sx, double sy, double sz);

inline Sign orient3d_expansion(double px, double py, double pz, double qx,
                               double qy, double qz, double rx, double ry,
                               double rz, double sx, double sy, double sz);

inline Sign orient3d(double px, double py, double pz, double qx, double qy,
                     double qz, double rx, double ry, double rz, double sx,
                     double sy, double sz);

inline Sign orient3d(const double *p, const double *q, const double *r,
                     const double *s);

template <typename IT, typename ET>
Sign orient3d(const GenericPoint3T<IT, ET> &p, const GenericPoint3T<IT, ET> &q,
              const GenericPoint3T<IT, ET> &r, const GenericPoint3T<IT, ET> &s);

// generally, pb, pc and pd come from three vertices of a triangle.
// this convention keeps consistent with shewchuk predicates.
inline void orient3d_get_minors(const double *pa, const double *pb,
                                const double *pc, double *minor, double *perm);

// generally, pb, pc and pd come from three vertices of a triangle,
// pa is the fourth query point.
// this convention keeps consistent with shewchuk predicates.
inline Sign orient3d_with_cached_minors(const double *pa, const double *pb,
                                        const double *pc, const double *pd,
                                        const double *minor,
                                        const double *perm);

/*************** OrientOn2D *******************/

inline Sign orient2dxy(double p1x, double p1y, double p2x, double p2y,
                       double p3x, double p3y);

inline Sign orient2dxy(const double *p1, const double *p2, const double *p3);

template <typename IT, typename ET>
Sign orient2dxy(const GenericPoint3T<IT, ET> &p1,
                const GenericPoint3T<IT, ET> &p2,
                const GenericPoint3T<IT, ET> &p3);

inline Sign orient2dyz(double p1y, double p1z, double p2y, double p2z,
                       double p3y, double p3z);

inline Sign orient2dyz(const double *p1, const double *p2, const double *p3);

template <typename IT, typename ET>
Sign orient2dyz(const GenericPoint3T<IT, ET> &p1,
                const GenericPoint3T<IT, ET> &p2,
                const GenericPoint3T<IT, ET> &p3);

inline Sign orient2dzx(double p1x, double p1z, double p2x, double p2z,
                       double p3x, double p3z);

inline Sign orient2dzx(const double *p1, const double *p2, const double *p3);

template <typename IT, typename ET>
Sign orient2dzx(const GenericPoint3T<IT, ET> &p1,
                const GenericPoint3T<IT, ET> &p2,
                const GenericPoint3T<IT, ET> &p3);

/*************** Max Component In Triangle Normal *******************/

inline int maxComponentInTriangleNormal_filtered(double ov1x, double ov1y,
                                                 double ov1z, double ov2x,
                                                 double ov2y, double ov2z,
                                                 double ov3x, double ov3y,
                                                 double ov3z);

inline int maxComponentInTriangleNormal_expansion(double ov1x, double ov1y,
                                                  double ov1z, double ov2x,
                                                  double ov2y, double ov2z,
                                                  double ov3x, double ov3y,
                                                  double ov3z);

inline int maxComponentInTriangleNormal(double ov1x, double ov1y, double ov1z,
                                        double ov2x, double ov2y, double ov2z,
                                        double ov3x, double ov3y, double ov3z);

/*************** LessThanOnAll (wrap to LessThanOnX/Y/Z) *******************/

inline std::array<Sign, 3> lessThanOnAll_EE(double x1, double y1, double z1,
                                            double x2, double y2, double z2);

template <typename IT, typename ET>
std::array<Sign, 3> lessThanOnAll_EE(const GenericPoint3T<IT, ET> &a,
                                     const GenericPoint3T<IT, ET> &b);

template <typename IT, typename ET, bool WithSSFilter>
std::array<Sign, 3> lessThanOnAll_IE(const GenericPoint3T<IT, ET> &p1,
                                     double x2, double y2, double z2,
                                     PntArr3 arr);

template <typename IT, typename ET, bool WithSSFilter>
std::array<Sign, 3> lessThanOnAll_IE(const GenericPoint3T<IT, ET> &a,
                                     const GenericPoint3T<IT, ET> &b,
                                     PntArr3                       arr);

template <typename IT, typename ET, bool WithSSFilter>
std::array<Sign, 3> lessThanOnAll_II(const GenericPoint3T<IT, ET> &p1,
                                     const GenericPoint3T<IT, ET> &p2,
                                     PntArr3                       arr);

/*************** LessThan (wrap to LessThanOnX/Y/Z) *******************/

inline Sign lessThan_EE(double x1, double y1, double z1, double x2, double y2,
                        double z2);

template <typename IT, typename ET>
Sign lessThan_EE(const GenericPoint3T<IT, ET> &a,
                 const GenericPoint3T<IT, ET> &b);

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThan_IE(const GenericPoint3T<IT, ET> &p1, double x2, double y2,
                 double z2, PntArr3 arr);

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThan_IE(const GenericPoint3T<IT, ET> &a,
                 const GenericPoint3T<IT, ET> &b, PntArr3 arr);

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThan_II(const GenericPoint3T<IT, ET> &p1,
                 const GenericPoint3T<IT, ET> &p2, PntArr3 arr);

/*************** OrienOn2D xy/yz/zx  IIE/III *******************/

template <typename IT, typename ET>
Sign orientOn2Dxy_IIE_expansion(const GenericPoint3T<IT, ET> &p1,
                                const GenericPoint3T<IT, ET> &p2, double op3x,
                                double op3y);

template <typename IT, typename ET>
Sign orientOn2Dxy_III_expansion(const GenericPoint3T<IT, ET> &p1,
                                const GenericPoint3T<IT, ET> &p2,
                                const GenericPoint3T<IT, ET> &p3);

template <typename IT, typename ET>
Sign orientOn2Dyz_IIE_expansion(const GenericPoint3T<IT, ET> &p1,
                                const GenericPoint3T<IT, ET> &p2, double op3y,
                                double op3z);

template <typename IT, typename ET>
Sign orientOn2Dyz_III_expansion(const GenericPoint3T<IT, ET> &p1,
                                const GenericPoint3T<IT, ET> &p2,
                                const GenericPoint3T<IT, ET> &p3);

template <typename IT, typename ET>
Sign orientOn2Dzx_IIE_expansion(const GenericPoint3T<IT, ET> &p1,
                                const GenericPoint3T<IT, ET> &p2, double op3x,
                                double op3z);

template <typename IT, typename ET>
Sign orientOn2Dzx_III_expansion(const GenericPoint3T<IT, ET> &p1,
                                const GenericPoint3T<IT, ET> &p2,
                                const GenericPoint3T<IT, ET> &p3);



} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "IndirectPredicateDetailsHand.inl"
#endif