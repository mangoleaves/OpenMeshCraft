#pragma once

#include "OpenMeshCraft/Geometry/Primitives/GenericPoint2T.h"
#include "OpenMeshCraft/Geometry/Primitives/GenericPoint3T.h"

#include "IndirectDefinitions.h"

#include "OpenMeshCraft/NumberTypes/NumberUtils.h"

namespace OMC {

inline Sign dotProductSign2D_filtered(double px, double py, double rx,
                                      double ry, double qx, double qy);

template <typename IT>
Sign dotProductSign2D_interval(IT px, IT py, IT rx, IT ry, IT qx, IT qy);

template <typename ET>
Sign dotProductSign2D_exact(ET px, ET py, ET rx, ET ry, ET qx, ET qy);

inline Sign dotProductSign2D_expansion(double px, double py, double rx,
                                       double ry, double qx, double qy);

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSign2D(double px, double py, double rx, double ry, double qx,
                      double qy);

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSign2D(const GenericPoint2T<IT, ET> &p,
                      const GenericPoint2T<IT, ET> &r,
                      const GenericPoint2T<IT, ET> &q);

inline Sign dotProductSign2D4P_filtered(double px, double py, double rx,
                                        double ry, double qx, double qy,
                                        double sx, double sy);

template <typename IT>
Sign dotProductSign2D4P_interval(IT px, IT py, IT rx, IT ry, IT qx, IT qy,
                                 IT sx, IT sy);

template <typename ET>
Sign dotProductSign2D4P_exact(ET px, ET py, ET rx, ET ry, ET qx, ET qy, ET sx,
                              ET sy);

inline Sign dotProductSign2D4P_expansion(double px, double py, double rx,
                                         double ry, double qx, double qy,
                                         double sx, double sy);

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSign2D4P(double px, double py, double rx, double ry, double qx,
                        double qy, double sx, double sy);

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSign2D4P(const GenericPoint2T<IT, ET> &p,
                        const GenericPoint2T<IT, ET> &r,
                        const GenericPoint2T<IT, ET> &q,
                        const GenericPoint2T<IT, ET> &s);

inline Sign dotProductSign3D_filtered(double px, double py, double pz,
                                      double rx, double ry, double rz,
                                      double qx, double qy, double qz);

template <typename IT>
Sign dotProductSign3D_interval(IT px, IT py, IT pz, IT rx, IT ry, IT rz, IT qx,
                               IT qy, IT qz);

template <typename ET>
Sign dotProductSign3D_exact(ET px, ET py, ET pz, ET rx, ET ry, ET rz, ET qx,
                            ET qy, ET qz);

inline Sign dotProductSign3D_expansion(double px, double py, double pz,
                                       double rx, double ry, double rz,
                                       double qx, double qy, double qz);

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSign3D(double px, double py, double pz, double rx, double ry,
                      double rz, double qx, double qy, double qz);

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSign3D(const GenericPoint3T<IT, ET> &p,
                      const GenericPoint3T<IT, ET> &r,
                      const GenericPoint3T<IT, ET> &q);

inline Sign dotProductSign3D4P_filtered(double px, double py, double pz,
                                        double rx, double ry, double rz,
                                        double qx, double qy, double qz,
                                        double sx, double sy, double sz);

template <typename IT>
Sign dotProductSign3D4P_interval(IT px, IT py, IT pz, IT rx, IT ry, IT rz,
                                 IT qx, IT qy, IT qz, IT sx, IT sy, IT sz);

template <typename ET>
Sign dotProductSign3D4P_exact(ET px, ET py, ET pz, ET rx, ET ry, ET rz, ET qx,
                              ET qy, ET qz, ET sx, ET sy, ET sz);

inline Sign dotProductSign3D4P_expansion(double px, double py, double pz,
                                         double rx, double ry, double rz,
                                         double qx, double qy, double qz,
                                         double sx, double sy, double sz);

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSign3D4P(double px, double py, double pz, double rx, double ry,
                        double rz, double qx, double qy, double qz, double sx,
                        double sy, double sz);

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSign3D4P(const GenericPoint3T<IT, ET> &p,
                        const GenericPoint3T<IT, ET> &r,
                        const GenericPoint3T<IT, ET> &q,
                        const GenericPoint3T<IT, ET> &s);

inline Sign dotProductSignOn2Dxy4P_filtered(double px, double py, double rx,
                                            double ry, double qx, double qy,
                                            double sx, double sy);

template <typename IT>
Sign dotProductSignOn2Dxy4P_interval(IT px, IT py, IT rx, IT ry, IT qx, IT qy,
                                     IT sx, IT sy);

template <typename ET>
Sign dotProductSignOn2Dxy4P_exact(ET px, ET py, ET rx, ET ry, ET qx, ET qy,
                                  ET sx, ET sy);

inline Sign dotProductSignOn2Dxy4P_expansion(double px, double py, double rx,
                                             double ry, double qx, double qy,
                                             double sx, double sy);

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSignOn2Dxy4P(double px, double py, double rx, double ry,
                            double qx, double qy, double sx, double sy);

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSignOn2Dxy4P(const GenericPoint3T<IT, ET> &p,
                            const GenericPoint3T<IT, ET> &r,
                            const GenericPoint3T<IT, ET> &q,
                            const GenericPoint3T<IT, ET> &s);

inline Sign dotProductSignOn2Dyz4P_filtered(double py, double pz, double ry,
                                            double rz, double qy, double qz,
                                            double sy, double sz);

template <typename IT>
Sign dotProductSignOn2Dyz4P_interval(IT py, IT pz, IT ry, IT rz, IT qy, IT qz,
                                     IT sy, IT sz);

template <typename ET>
Sign dotProductSignOn2Dyz4P_exact(ET py, ET pz, ET ry, ET rz, ET qy, ET qz,
                                  ET sy, ET sz);

inline Sign dotProductSignOn2Dyz4P_expansion(double py, double pz, double ry,
                                             double rz, double qy, double qz,
                                             double sy, double sz);

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSignOn2Dyz4P(double py, double pz, double ry, double rz,
                            double qy, double qz, double sy, double sz);

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSignOn2Dyz4P(const GenericPoint3T<IT, ET> &p,
                            const GenericPoint3T<IT, ET> &r,
                            const GenericPoint3T<IT, ET> &q,
                            const GenericPoint3T<IT, ET> &s);

inline Sign dotProductSignOn2Dzx4P_filtered(double px, double pz, double rx,
                                            double rz, double qx, double qz,
                                            double sx, double sz);

template <typename IT>
Sign dotProductSignOn2Dzx4P_interval(IT px, IT pz, IT rx, IT rz, IT qx, IT qz,
                                     IT sx, IT sz);

template <typename ET>
Sign dotProductSignOn2Dzx4P_exact(ET px, ET pz, ET rx, ET rz, ET qx, ET qz,
                                  ET sx, ET sz);

inline Sign dotProductSignOn2Dzx4P_expansion(double px, double pz, double rx,
                                             double rz, double qx, double qz,
                                             double sx, double sz);

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSignOn2Dzx4P(double px, double pz, double rx, double rz,
                            double qx, double qz, double sx, double sz);

template <typename IT, typename ET, bool WithSSFilter>
Sign dotProductSignOn2Dzx4P(const GenericPoint3T<IT, ET> &p,
                            const GenericPoint3T<IT, ET> &r,
                            const GenericPoint3T<IT, ET> &q,
                            const GenericPoint3T<IT, ET> &s);

inline Sign inCircle_filtered(double pax, double pay, double pbx, double pby,
                              double pcx, double pcy, double pdx, double pdy);

template <typename IT>
Sign inCircle_interval(IT pax, IT pay, IT pbx, IT pby, IT pcx, IT pcy, IT pdx,
                       IT pdy);

template <typename ET>
Sign inCircle_exact(ET pax, ET pay, ET pbx, ET pby, ET pcx, ET pcy, ET pdx,
                    ET pdy);

inline Sign inCircle_expansion(double pax, double pay, double pbx, double pby,
                               double pcx, double pcy, double pdx, double pdy);

template <typename IT, typename ET, bool WithSSFilter>
Sign inCircle(double pax, double pay, double pbx, double pby, double pcx,
              double pcy, double pdx, double pdy);

template <typename IT, typename ET, bool WithSSFilter>
Sign inCircle(const GenericPoint2T<IT, ET> &pa,
              const GenericPoint2T<IT, ET> &pb,
              const GenericPoint2T<IT, ET> &pc,
              const GenericPoint2T<IT, ET> &pd);

inline Sign inSphere_filtered(double pax, double pay, double paz, double pbx,
                              double pby, double pbz, double pcx, double pcy,
                              double pcz, double pdx, double pdy, double pdz,
                              double pex, double pey, double pez);

template <typename IT>
Sign inSphere_interval(IT pax, IT pay, IT paz, IT pbx, IT pby, IT pbz, IT pcx,
                       IT pcy, IT pcz, IT pdx, IT pdy, IT pdz, IT pex, IT pey,
                       IT pez);

template <typename ET>
Sign inSphere_exact(ET pax, ET pay, ET paz, ET pbx, ET pby, ET pbz, ET pcx,
                    ET pcy, ET pcz, ET pdx, ET pdy, ET pdz, ET pex, ET pey,
                    ET pez);

inline Sign inSphere_expansion(double pax, double pay, double paz, double pbx,
                               double pby, double pbz, double pcx, double pcy,
                               double pcz, double pdx, double pdy, double pdz,
                               double pex, double pey, double pez);

template <typename IT, typename ET, bool WithSSFilter>
Sign inSphere(double pax, double pay, double paz, double pbx, double pby,
              double pbz, double pcx, double pcy, double pcz, double pdx,
              double pdy, double pdz, double pex, double pey, double pez);

template <typename IT, typename ET, bool WithSSFilter>
Sign inSphere(const GenericPoint3T<IT, ET> &pa,
              const GenericPoint3T<IT, ET> &pb,
              const GenericPoint3T<IT, ET> &pc,
              const GenericPoint3T<IT, ET> &pd,
              const GenericPoint3T<IT, ET> &pe);

inline Sign squareDistance2D_filtered(double px, double py, double qx,
                                      double qy, double dis);

template <typename IT>
Sign squareDistance2D_interval(IT px, IT py, IT qx, IT qy, IT dis);

template <typename ET>
Sign squareDistance2D_exact(ET px, ET py, ET qx, ET qy, ET dis);

inline Sign squareDistance2D_expansion(double px, double py, double qx,
                                       double qy, double dis);

template <typename IT, typename ET, bool WithSSFilter>
Sign squareDistance2D(double px, double py, double qx, double qy, double dis);

template <typename IT, typename ET, bool WithSSFilter>
Sign squareDistance2D(const GenericPoint2T<IT, ET> &p,
                      const GenericPoint2T<IT, ET> &q, double dis);

inline Sign squareDistance3D_filtered(double px, double py, double pz,
                                      double qx, double qy, double qz,
                                      double dis);

template <typename IT>
Sign squareDistance3D_interval(IT px, IT py, IT pz, IT qx, IT qy, IT qz,
                               IT dis);

template <typename ET>
Sign squareDistance3D_exact(ET px, ET py, ET pz, ET qx, ET qy, ET qz, ET dis);

inline Sign squareDistance3D_expansion(double px, double py, double pz,
                                       double qx, double qy, double qz,
                                       double dis);

template <typename IT, typename ET, bool WithSSFilter>
Sign squareDistance3D(double px, double py, double pz, double qx, double qy,
                      double qz, double dis);

template <typename IT, typename ET, bool WithSSFilter>
Sign squareDistance3D(const GenericPoint3T<IT, ET> &p,
                      const GenericPoint3T<IT, ET> &q, double dis);

template <typename IT, typename ET>
Sign dotProductSign2D_EEI_interval(const GenericPoint2T<IT, ET> &q, IT px,
                                   IT py, IT rx, IT ry);

template <typename IT, typename ET>
Sign dotProductSign2D_EEI_exact(const GenericPoint2T<IT, ET> &q, ET px, ET py,
                                ET rx, ET ry);

template <typename IT, typename ET>
Sign dotProductSign2D_EEI_expansion(const GenericPoint2T<IT, ET> &q, double px,
                                    double py, double rx, double ry);

template <typename IT, typename ET>
Sign dotProductSign2D_EEI(const GenericPoint2T<IT, ET> &q, double px, double py,
                          double rx, double ry);

template <typename IT, typename ET>
Sign dotProductSign2D_EEI(const GenericPoint2T<IT, ET> &q,
                          const GenericPoint2T<IT, ET> &p,
                          const GenericPoint2T<IT, ET> &r);

template <typename IT, typename ET>
Sign dotProductSign2D_IEE_interval(const GenericPoint2T<IT, ET> &p, IT rx,
                                   IT ry, IT qx, IT qy);

template <typename IT, typename ET>
Sign dotProductSign2D_IEE_exact(const GenericPoint2T<IT, ET> &p, ET rx, ET ry,
                                ET qx, ET qy);

template <typename IT, typename ET>
Sign dotProductSign2D_IEE_expansion(const GenericPoint2T<IT, ET> &p, double rx,
                                    double ry, double qx, double qy);

template <typename IT, typename ET>
Sign dotProductSign2D_IEE(const GenericPoint2T<IT, ET> &p, double rx, double ry,
                          double qx, double qy);

template <typename IT, typename ET>
Sign dotProductSign2D_IEE(const GenericPoint2T<IT, ET> &p,
                          const GenericPoint2T<IT, ET> &r,
                          const GenericPoint2T<IT, ET> &q);

template <typename IT, typename ET>
Sign dotProductSign2D_IEI_interval(const GenericPoint2T<IT, ET> &p,
                                   const GenericPoint2T<IT, ET> &q, IT rx,
                                   IT ry);

template <typename IT, typename ET>
Sign dotProductSign2D_IEI_exact(const GenericPoint2T<IT, ET> &p,
                                const GenericPoint2T<IT, ET> &q, ET rx, ET ry);

template <typename IT, typename ET>
Sign dotProductSign2D_IEI_expansion(const GenericPoint2T<IT, ET> &p,
                                    const GenericPoint2T<IT, ET> &q, double rx,
                                    double ry);

template <typename IT, typename ET>
Sign dotProductSign2D_IEI(const GenericPoint2T<IT, ET> &p,
                          const GenericPoint2T<IT, ET> &q, double rx,
                          double ry);

template <typename IT, typename ET>
Sign dotProductSign2D_IEI(const GenericPoint2T<IT, ET> &p,
                          const GenericPoint2T<IT, ET> &q,
                          const GenericPoint2T<IT, ET> &r);

template <typename IT, typename ET>
Sign dotProductSign2D_IIE_interval(const GenericPoint2T<IT, ET> &p,
                                   const GenericPoint2T<IT, ET> &r, IT qx,
                                   IT qy);

template <typename IT, typename ET>
Sign dotProductSign2D_IIE_exact(const GenericPoint2T<IT, ET> &p,
                                const GenericPoint2T<IT, ET> &r, ET qx, ET qy);

template <typename IT, typename ET>
Sign dotProductSign2D_IIE_expansion(const GenericPoint2T<IT, ET> &p,
                                    const GenericPoint2T<IT, ET> &r, double qx,
                                    double qy);

template <typename IT, typename ET>
Sign dotProductSign2D_IIE(const GenericPoint2T<IT, ET> &p,
                          const GenericPoint2T<IT, ET> &r, double qx,
                          double qy);

template <typename IT, typename ET>
Sign dotProductSign2D_IIE(const GenericPoint2T<IT, ET> &p,
                          const GenericPoint2T<IT, ET> &r,
                          const GenericPoint2T<IT, ET> &q);

template <typename IT, typename ET>
Sign dotProductSign2D_III_interval(const GenericPoint2T<IT, ET> &p,
                                   const GenericPoint2T<IT, ET> &r,
                                   const GenericPoint2T<IT, ET> &q);

template <typename IT, typename ET>
Sign dotProductSign2D_III_exact(const GenericPoint2T<IT, ET> &p,
                                const GenericPoint2T<IT, ET> &r,
                                const GenericPoint2T<IT, ET> &q);

template <typename IT, typename ET>
Sign dotProductSign2D_III_expansion(const GenericPoint2T<IT, ET> &p,
                                    const GenericPoint2T<IT, ET> &r,
                                    const GenericPoint2T<IT, ET> &q);

template <typename IT, typename ET>
Sign dotProductSign2D_III(const GenericPoint2T<IT, ET> &p,
                          const GenericPoint2T<IT, ET> &r,
                          const GenericPoint2T<IT, ET> &q);

template <typename IT, typename ET>
Sign dotProductSign3D_EEI_interval(const GenericPoint3T<IT, ET> &q, IT px,
                                   IT py, IT pz, IT rx, IT ry, IT rz);

template <typename IT, typename ET>
Sign dotProductSign3D_EEI_exact(const GenericPoint3T<IT, ET> &q, ET px, ET py,
                                ET pz, ET rx, ET ry, ET rz);

template <typename IT, typename ET>
Sign dotProductSign3D_EEI_expansion(const GenericPoint3T<IT, ET> &q, double px,
                                    double py, double pz, double rx, double ry,
                                    double rz);

template <typename IT, typename ET>
Sign dotProductSign3D_EEI(const GenericPoint3T<IT, ET> &q, double px, double py,
                          double pz, double rx, double ry, double rz);

template <typename IT, typename ET>
Sign dotProductSign3D_EEI(const GenericPoint3T<IT, ET> &q,
                          const GenericPoint3T<IT, ET> &p,
                          const GenericPoint3T<IT, ET> &r);

template <typename IT, typename ET>
Sign dotProductSign3D_IEE_interval(const GenericPoint3T<IT, ET> &p, IT rx,
                                   IT ry, IT rz, IT qx, IT qy, IT qz);

template <typename IT, typename ET>
Sign dotProductSign3D_IEE_exact(const GenericPoint3T<IT, ET> &p, ET rx, ET ry,
                                ET rz, ET qx, ET qy, ET qz);

template <typename IT, typename ET>
Sign dotProductSign3D_IEE_expansion(const GenericPoint3T<IT, ET> &p, double rx,
                                    double ry, double rz, double qx, double qy,
                                    double qz);

template <typename IT, typename ET>
Sign dotProductSign3D_IEE(const GenericPoint3T<IT, ET> &p, double rx, double ry,
                          double rz, double qx, double qy, double qz);

template <typename IT, typename ET>
Sign dotProductSign3D_IEE(const GenericPoint3T<IT, ET> &p,
                          const GenericPoint3T<IT, ET> &r,
                          const GenericPoint3T<IT, ET> &q);

template <typename IT, typename ET>
Sign dotProductSign3D_IEI_interval(const GenericPoint3T<IT, ET> &p,
                                   const GenericPoint3T<IT, ET> &q, IT rx,
                                   IT ry, IT rz);

template <typename IT, typename ET>
Sign dotProductSign3D_IEI_exact(const GenericPoint3T<IT, ET> &p,
                                const GenericPoint3T<IT, ET> &q, ET rx, ET ry,
                                ET rz);

template <typename IT, typename ET>
Sign dotProductSign3D_IEI_expansion(const GenericPoint3T<IT, ET> &p,
                                    const GenericPoint3T<IT, ET> &q, double rx,
                                    double ry, double rz);

template <typename IT, typename ET>
Sign dotProductSign3D_IEI(const GenericPoint3T<IT, ET> &p,
                          const GenericPoint3T<IT, ET> &q, double rx, double ry,
                          double rz);

template <typename IT, typename ET>
Sign dotProductSign3D_IEI(const GenericPoint3T<IT, ET> &p,
                          const GenericPoint3T<IT, ET> &q,
                          const GenericPoint3T<IT, ET> &r);

template <typename IT, typename ET>
Sign dotProductSign3D_IIE_interval(const GenericPoint3T<IT, ET> &p,
                                   const GenericPoint3T<IT, ET> &r, IT qx,
                                   IT qy, IT qz);

template <typename IT, typename ET>
Sign dotProductSign3D_IIE_exact(const GenericPoint3T<IT, ET> &p,
                                const GenericPoint3T<IT, ET> &r, ET qx, ET qy,
                                ET qz);

template <typename IT, typename ET>
Sign dotProductSign3D_IIE_expansion(const GenericPoint3T<IT, ET> &p,
                                    const GenericPoint3T<IT, ET> &r, double qx,
                                    double qy, double qz);

template <typename IT, typename ET>
Sign dotProductSign3D_IIE(const GenericPoint3T<IT, ET> &p,
                          const GenericPoint3T<IT, ET> &r, double qx, double qy,
                          double qz);

template <typename IT, typename ET>
Sign dotProductSign3D_IIE(const GenericPoint3T<IT, ET> &p,
                          const GenericPoint3T<IT, ET> &r,
                          const GenericPoint3T<IT, ET> &q);

template <typename IT, typename ET>
Sign dotProductSign3D_III_interval(const GenericPoint3T<IT, ET> &p,
                                   const GenericPoint3T<IT, ET> &r,
                                   const GenericPoint3T<IT, ET> &q);

template <typename IT, typename ET>
Sign dotProductSign3D_III_exact(const GenericPoint3T<IT, ET> &p,
                                const GenericPoint3T<IT, ET> &r,
                                const GenericPoint3T<IT, ET> &q);

template <typename IT, typename ET>
Sign dotProductSign3D_III_expansion(const GenericPoint3T<IT, ET> &p,
                                    const GenericPoint3T<IT, ET> &r,
                                    const GenericPoint3T<IT, ET> &q);

template <typename IT, typename ET>
Sign dotProductSign3D_III(const GenericPoint3T<IT, ET> &p,
                          const GenericPoint3T<IT, ET> &r,
                          const GenericPoint3T<IT, ET> &q);

template <typename IT, typename ET>
Sign inCirclexy_IEEE_interval(const GenericPoint3T<IT, ET> &p1, IT pbx, IT pby,
                              IT pcx, IT pcy, IT pdx, IT pdy);

template <typename IT, typename ET>
Sign inCirclexy_IEEE_exact(const GenericPoint3T<IT, ET> &p1, ET pbx, ET pby,
                           ET pcx, ET pcy, ET pdx, ET pdy);

template <typename IT, typename ET>
Sign inCirclexy_IEEE_expansion(const GenericPoint3T<IT, ET> &p1, double pbx,
                               double pby, double pcx, double pcy, double pdx,
                               double pdy);

template <typename IT, typename ET>
Sign inCirclexy_IEEE(const GenericPoint3T<IT, ET> &p1, double pbx, double pby,
                     double pcx, double pcy, double pdx, double pdy);

template <typename IT, typename ET>
Sign inCirclexy_IEEE(const GenericPoint3T<IT, ET> &p1,
                     const GenericPoint3T<IT, ET> &pb,
                     const GenericPoint3T<IT, ET> &pc,
                     const GenericPoint3T<IT, ET> &pd);

template <typename IT, typename ET>
Sign inCirclexy_IIEE_interval(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2, IT pcx, IT pcy,
                              IT pdx, IT pdy);

template <typename IT, typename ET>
Sign inCirclexy_IIEE_exact(const GenericPoint3T<IT, ET> &p1,
                           const GenericPoint3T<IT, ET> &p2, ET pcx, ET pcy,
                           ET pdx, ET pdy);

template <typename IT, typename ET>
Sign inCirclexy_IIEE_expansion(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2, double pcx,
                               double pcy, double pdx, double pdy);

template <typename IT, typename ET>
Sign inCirclexy_IIEE(const GenericPoint3T<IT, ET> &p1,
                     const GenericPoint3T<IT, ET> &p2, double pcx, double pcy,
                     double pdx, double pdy);

template <typename IT, typename ET>
Sign inCirclexy_IIEE(const GenericPoint3T<IT, ET> &p1,
                     const GenericPoint3T<IT, ET> &p2,
                     const GenericPoint3T<IT, ET> &pc,
                     const GenericPoint3T<IT, ET> &pd);

template <typename IT, typename ET>
Sign inCirclexy_IIIE_interval(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2,
                              const GenericPoint3T<IT, ET> &p3, IT pdx, IT pdy);

template <typename IT, typename ET>
Sign inCirclexy_IIIE_exact(const GenericPoint3T<IT, ET> &p1,
                           const GenericPoint3T<IT, ET> &p2,
                           const GenericPoint3T<IT, ET> &p3, ET pdx, ET pdy);

template <typename IT, typename ET>
Sign inCirclexy_IIIE_expansion(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3, double pdx,
                               double pdy);

template <typename IT, typename ET>
Sign inCirclexy_IIIE(const GenericPoint3T<IT, ET> &p1,
                     const GenericPoint3T<IT, ET> &p2,
                     const GenericPoint3T<IT, ET> &p3, double pdx, double pdy);

template <typename IT, typename ET>
Sign inCirclexy_IIIE(const GenericPoint3T<IT, ET> &p1,
                     const GenericPoint3T<IT, ET> &p2,
                     const GenericPoint3T<IT, ET> &p3,
                     const GenericPoint3T<IT, ET> &pd);

template <typename IT, typename ET>
Sign inCirclexy_IIII_interval(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2,
                              const GenericPoint3T<IT, ET> &p3,
                              const GenericPoint3T<IT, ET> &p4);

template <typename IT, typename ET>
Sign inCirclexy_IIII_exact(const GenericPoint3T<IT, ET> &p1,
                           const GenericPoint3T<IT, ET> &p2,
                           const GenericPoint3T<IT, ET> &p3,
                           const GenericPoint3T<IT, ET> &p4);

template <typename IT, typename ET>
Sign inCirclexy_IIII_expansion(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3,
                               const GenericPoint3T<IT, ET> &p4);

template <typename IT, typename ET>
Sign inCirclexy_IIII(const GenericPoint3T<IT, ET> &p1,
                     const GenericPoint3T<IT, ET> &p2,
                     const GenericPoint3T<IT, ET> &p3,
                     const GenericPoint3T<IT, ET> &p4);

template <typename IT, typename ET>
Sign inCircle_IEEE_interval(const GenericPoint2T<IT, ET> &p1, IT pbx, IT pby,
                            IT pcx, IT pcy, IT pdx, IT pdy);

template <typename IT, typename ET>
Sign inCircle_IEEE_exact(const GenericPoint2T<IT, ET> &p1, ET pbx, ET pby,
                         ET pcx, ET pcy, ET pdx, ET pdy);

template <typename IT, typename ET>
Sign inCircle_IEEE_expansion(const GenericPoint2T<IT, ET> &p1, double pbx,
                             double pby, double pcx, double pcy, double pdx,
                             double pdy);

template <typename IT, typename ET>
Sign inCircle_IEEE(const GenericPoint2T<IT, ET> &p1, double pbx, double pby,
                   double pcx, double pcy, double pdx, double pdy);

template <typename IT, typename ET>
Sign inCircle_IEEE(const GenericPoint2T<IT, ET> &p1,
                   const GenericPoint2T<IT, ET> &pb,
                   const GenericPoint2T<IT, ET> &pc,
                   const GenericPoint2T<IT, ET> &pd);

template <typename IT, typename ET>
Sign inCircle_IIEE_interval(const GenericPoint2T<IT, ET> &p1,
                            const GenericPoint2T<IT, ET> &p2, IT pcx, IT pcy,
                            IT pdx, IT pdy);

template <typename IT, typename ET>
Sign inCircle_IIEE_exact(const GenericPoint2T<IT, ET> &p1,
                         const GenericPoint2T<IT, ET> &p2, ET pcx, ET pcy,
                         ET pdx, ET pdy);

template <typename IT, typename ET>
Sign inCircle_IIEE_expansion(const GenericPoint2T<IT, ET> &p1,
                             const GenericPoint2T<IT, ET> &p2, double pcx,
                             double pcy, double pdx, double pdy);

template <typename IT, typename ET>
Sign inCircle_IIEE(const GenericPoint2T<IT, ET> &p1,
                   const GenericPoint2T<IT, ET> &p2, double pcx, double pcy,
                   double pdx, double pdy);

template <typename IT, typename ET>
Sign inCircle_IIEE(const GenericPoint2T<IT, ET> &p1,
                   const GenericPoint2T<IT, ET> &p2,
                   const GenericPoint2T<IT, ET> &pc,
                   const GenericPoint2T<IT, ET> &pd);

template <typename IT, typename ET>
Sign inCircle_IIIE_interval(const GenericPoint2T<IT, ET> &p1,
                            const GenericPoint2T<IT, ET> &p2,
                            const GenericPoint2T<IT, ET> &p3, IT pdx, IT pdy);

template <typename IT, typename ET>
Sign inCircle_IIIE_exact(const GenericPoint2T<IT, ET> &p1,
                         const GenericPoint2T<IT, ET> &p2,
                         const GenericPoint2T<IT, ET> &p3, ET pdx, ET pdy);

template <typename IT, typename ET>
Sign inCircle_IIIE_expansion(const GenericPoint2T<IT, ET> &p1,
                             const GenericPoint2T<IT, ET> &p2,
                             const GenericPoint2T<IT, ET> &p3, double pdx,
                             double pdy);

template <typename IT, typename ET>
Sign inCircle_IIIE(const GenericPoint2T<IT, ET> &p1,
                   const GenericPoint2T<IT, ET> &p2,
                   const GenericPoint2T<IT, ET> &p3, double pdx, double pdy);

template <typename IT, typename ET>
Sign inCircle_IIIE(const GenericPoint2T<IT, ET> &p1,
                   const GenericPoint2T<IT, ET> &p2,
                   const GenericPoint2T<IT, ET> &p3,
                   const GenericPoint2T<IT, ET> &pd);

template <typename IT, typename ET>
Sign inCircle_IIII_interval(const GenericPoint2T<IT, ET> &p1,
                            const GenericPoint2T<IT, ET> &p2,
                            const GenericPoint2T<IT, ET> &p3,
                            const GenericPoint2T<IT, ET> &p4);

template <typename IT, typename ET>
Sign inCircle_IIII_exact(const GenericPoint2T<IT, ET> &p1,
                         const GenericPoint2T<IT, ET> &p2,
                         const GenericPoint2T<IT, ET> &p3,
                         const GenericPoint2T<IT, ET> &p4);

template <typename IT, typename ET>
Sign inCircle_IIII_expansion(const GenericPoint2T<IT, ET> &p1,
                             const GenericPoint2T<IT, ET> &p2,
                             const GenericPoint2T<IT, ET> &p3,
                             const GenericPoint2T<IT, ET> &p4);

template <typename IT, typename ET>
Sign inCircle_IIII(const GenericPoint2T<IT, ET> &p1,
                   const GenericPoint2T<IT, ET> &p2,
                   const GenericPoint2T<IT, ET> &p3,
                   const GenericPoint2T<IT, ET> &p4);

template <typename IT, typename ET>
Sign inSphere_IEEEE_interval(const GenericPoint3T<IT, ET> &p1, IT pbx, IT pby,
                             IT pbz, IT pcx, IT pcy, IT pcz, IT pdx, IT pdy,
                             IT pdz, IT pex, IT pey, IT pez);

template <typename IT, typename ET>
Sign inSphere_IEEEE_exact(const GenericPoint3T<IT, ET> &p1, ET pbx, ET pby,
                          ET pbz, ET pcx, ET pcy, ET pcz, ET pdx, ET pdy,
                          ET pdz, ET pex, ET pey, ET pez);

template <typename IT, typename ET>
Sign inSphere_IEEEE_expansion(const GenericPoint3T<IT, ET> &p1, double pbx,
                              double pby, double pbz, double pcx, double pcy,
                              double pcz, double pdx, double pdy, double pdz,
                              double pex, double pey, double pez);

template <typename IT, typename ET>
Sign inSphere_IEEEE(const GenericPoint3T<IT, ET> &p1, double pbx, double pby,
                    double pbz, double pcx, double pcy, double pcz, double pdx,
                    double pdy, double pdz, double pex, double pey, double pez);

template <typename IT, typename ET>
Sign inSphere_IEEEE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &pb,
                    const GenericPoint3T<IT, ET> &pc,
                    const GenericPoint3T<IT, ET> &pd,
                    const GenericPoint3T<IT, ET> &pe);

template <typename IT, typename ET>
Sign inSphere_IIEEE_interval(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2, IT pcx, IT pcy,
                             IT pcz, IT pdx, IT pdy, IT pdz, IT pex, IT pey,
                             IT pez);

template <typename IT, typename ET>
Sign inSphere_IIEEE_exact(const GenericPoint3T<IT, ET> &p1,
                          const GenericPoint3T<IT, ET> &p2, ET pcx, ET pcy,
                          ET pcz, ET pdx, ET pdy, ET pdz, ET pex, ET pey,
                          ET pez);

template <typename IT, typename ET>
Sign inSphere_IIEEE_expansion(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2, double pcx,
                              double pcy, double pcz, double pdx, double pdy,
                              double pdz, double pex, double pey, double pez);

template <typename IT, typename ET>
Sign inSphere_IIEEE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2, double pcx, double pcy,
                    double pcz, double pdx, double pdy, double pdz, double pex,
                    double pey, double pez);

template <typename IT, typename ET>
Sign inSphere_IIEEE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2,
                    const GenericPoint3T<IT, ET> &pc,
                    const GenericPoint3T<IT, ET> &pd,
                    const GenericPoint3T<IT, ET> &pe);

template <typename IT, typename ET>
Sign inSphere_IIIEE_interval(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2,
                             const GenericPoint3T<IT, ET> &p3, IT pdx, IT pdy,
                             IT pdz, IT pex, IT pey, IT pez);

template <typename IT, typename ET>
Sign inSphere_IIIEE_exact(const GenericPoint3T<IT, ET> &p1,
                          const GenericPoint3T<IT, ET> &p2,
                          const GenericPoint3T<IT, ET> &p3, ET pdx, ET pdy,
                          ET pdz, ET pex, ET pey, ET pez);

template <typename IT, typename ET>
Sign inSphere_IIIEE_expansion(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2,
                              const GenericPoint3T<IT, ET> &p3, double pdx,
                              double pdy, double pdz, double pex, double pey,
                              double pez);

template <typename IT, typename ET>
Sign inSphere_IIIEE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2,
                    const GenericPoint3T<IT, ET> &p3, double pdx, double pdy,
                    double pdz, double pex, double pey, double pez);

template <typename IT, typename ET>
Sign inSphere_IIIEE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2,
                    const GenericPoint3T<IT, ET> &p3,
                    const GenericPoint3T<IT, ET> &pd,
                    const GenericPoint3T<IT, ET> &pe);

template <typename IT, typename ET>
Sign inSphere_IIIIE_interval(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2,
                             const GenericPoint3T<IT, ET> &p3,
                             const GenericPoint3T<IT, ET> &p4, IT pex, IT pey,
                             IT pez);

template <typename IT, typename ET>
Sign inSphere_IIIIE_exact(const GenericPoint3T<IT, ET> &p1,
                          const GenericPoint3T<IT, ET> &p2,
                          const GenericPoint3T<IT, ET> &p3,
                          const GenericPoint3T<IT, ET> &p4, ET pex, ET pey,
                          ET pez);

template <typename IT, typename ET>
Sign inSphere_IIIIE_expansion(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2,
                              const GenericPoint3T<IT, ET> &p3,
                              const GenericPoint3T<IT, ET> &p4, double pex,
                              double pey, double pez);

template <typename IT, typename ET>
Sign inSphere_IIIIE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2,
                    const GenericPoint3T<IT, ET> &p3,
                    const GenericPoint3T<IT, ET> &p4, double pex, double pey,
                    double pez);

template <typename IT, typename ET>
Sign inSphere_IIIIE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2,
                    const GenericPoint3T<IT, ET> &p3,
                    const GenericPoint3T<IT, ET> &p4,
                    const GenericPoint3T<IT, ET> &pe);

template <typename IT, typename ET>
Sign inSphere_IIIII_interval(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2,
                             const GenericPoint3T<IT, ET> &p3,
                             const GenericPoint3T<IT, ET> &p4,
                             const GenericPoint3T<IT, ET> &p5);

template <typename IT, typename ET>
Sign inSphere_IIIII_exact(const GenericPoint3T<IT, ET> &p1,
                          const GenericPoint3T<IT, ET> &p2,
                          const GenericPoint3T<IT, ET> &p3,
                          const GenericPoint3T<IT, ET> &p4,
                          const GenericPoint3T<IT, ET> &p5);

template <typename IT, typename ET>
Sign inSphere_IIIII_expansion(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2,
                              const GenericPoint3T<IT, ET> &p3,
                              const GenericPoint3T<IT, ET> &p4,
                              const GenericPoint3T<IT, ET> &p5);

template <typename IT, typename ET>
Sign inSphere_IIIII(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2,
                    const GenericPoint3T<IT, ET> &p3,
                    const GenericPoint3T<IT, ET> &p4,
                    const GenericPoint3T<IT, ET> &p5);

template <typename IT, typename ET>
Sign lessThanOnX_IE_filtered(const GenericPoint3T<IT, ET> &p1, double bx,
                             PntArr3 arr);

template <typename IT, typename ET>
Sign lessThanOnX_IE_interval(const GenericPoint3T<IT, ET> &p1, IT bx);

template <typename IT, typename ET>
Sign lessThanOnX_IE_exact(const GenericPoint3T<IT, ET> &p1, ET bx);

template <typename IT, typename ET>
Sign lessThanOnX_IE_expansion(const GenericPoint3T<IT, ET> &p1, double bx);

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnX_IE(const GenericPoint3T<IT, ET> &p1, double bx, PntArr3 arr);

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnX_IE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &b, PntArr3 arr);

template <typename IT, typename ET>
Sign lessThanOnX_II_filtered(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2, PntArr3 arr);

template <typename IT, typename ET>
Sign lessThanOnX_II_interval(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2);

template <typename IT, typename ET>
Sign lessThanOnX_II_exact(const GenericPoint3T<IT, ET> &p1,
                          const GenericPoint3T<IT, ET> &p2);

template <typename IT, typename ET>
Sign lessThanOnX_II_expansion(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2);

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnX_II(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2, PntArr3 arr);

template <typename IT, typename ET>
Sign lessThanOnY_IE_filtered(const GenericPoint3T<IT, ET> &p1, double by,
                             PntArr3 arr);

template <typename IT, typename ET>
Sign lessThanOnY_IE_interval(const GenericPoint3T<IT, ET> &p1, IT by);

template <typename IT, typename ET>
Sign lessThanOnY_IE_exact(const GenericPoint3T<IT, ET> &p1, ET by);

template <typename IT, typename ET>
Sign lessThanOnY_IE_expansion(const GenericPoint3T<IT, ET> &p1, double by);

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnY_IE(const GenericPoint3T<IT, ET> &p1, double by, PntArr3 arr);

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnY_IE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &b, PntArr3 arr);

template <typename IT, typename ET>
Sign lessThanOnY_II_filtered(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2, PntArr3 arr);

template <typename IT, typename ET>
Sign lessThanOnY_II_interval(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2);

template <typename IT, typename ET>
Sign lessThanOnY_II_exact(const GenericPoint3T<IT, ET> &p1,
                          const GenericPoint3T<IT, ET> &p2);

template <typename IT, typename ET>
Sign lessThanOnY_II_expansion(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2);

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnY_II(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2, PntArr3 arr);

template <typename IT, typename ET>
Sign lessThanOnZ_IE_filtered(const GenericPoint3T<IT, ET> &p1, double bz,
                             PntArr3 arr);

template <typename IT, typename ET>
Sign lessThanOnZ_IE_interval(const GenericPoint3T<IT, ET> &p1, IT bz);

template <typename IT, typename ET>
Sign lessThanOnZ_IE_exact(const GenericPoint3T<IT, ET> &p1, ET bz);

template <typename IT, typename ET>
Sign lessThanOnZ_IE_expansion(const GenericPoint3T<IT, ET> &p1, double bz);

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnZ_IE(const GenericPoint3T<IT, ET> &p1, double bz, PntArr3 arr);

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnZ_IE(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &b, PntArr3 arr);

template <typename IT, typename ET>
Sign lessThanOnZ_II_filtered(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2, PntArr3 arr);

template <typename IT, typename ET>
Sign lessThanOnZ_II_interval(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2);

template <typename IT, typename ET>
Sign lessThanOnZ_II_exact(const GenericPoint3T<IT, ET> &p1,
                          const GenericPoint3T<IT, ET> &p2);

template <typename IT, typename ET>
Sign lessThanOnZ_II_expansion(const GenericPoint3T<IT, ET> &p1,
                              const GenericPoint3T<IT, ET> &p2);

template <typename IT, typename ET, bool WithSSFilter>
Sign lessThanOnZ_II(const GenericPoint3T<IT, ET> &p1,
                    const GenericPoint3T<IT, ET> &p2, PntArr3 arr);

template <typename IT, typename ET>
Sign orient2D_IEE_filtered(const GenericPoint2T<IT, ET> &p1, double p2x,
                           double p2y, double p3x, double p3y, PntArr2 arr);

template <typename IT, typename ET>
Sign orient2D_IEE_interval(const GenericPoint2T<IT, ET> &p1, IT p2x, IT p2y,
                           IT p3x, IT p3y);

template <typename IT, typename ET>
Sign orient2D_IEE_exact(const GenericPoint2T<IT, ET> &p1, ET p2x, ET p2y,
                        ET p3x, ET p3y);

template <typename IT, typename ET>
Sign orient2D_IEE_expansion(const GenericPoint2T<IT, ET> &p1, double p2x,
                            double p2y, double p3x, double p3y);

template <typename IT, typename ET, bool WithSSFilter>
Sign orient2D_IEE(const GenericPoint2T<IT, ET> &p1, double p2x, double p2y,
                  double p3x, double p3y, PntArr2 arr);

template <typename IT, typename ET, bool WithSSFilter>
Sign orient2D_IEE(const GenericPoint2T<IT, ET> &p1,
                  const GenericPoint2T<IT, ET> &p2,
                  const GenericPoint2T<IT, ET> &p3, PntArr2 arr);

template <typename IT, typename ET>
Sign orient2D_IIE_filtered(const GenericPoint2T<IT, ET> &p1,
                           const GenericPoint2T<IT, ET> &p2, double p3x,
                           double p3y, PntArr2 arr);

template <typename IT, typename ET>
Sign orient2D_IIE_interval(const GenericPoint2T<IT, ET> &p1,
                           const GenericPoint2T<IT, ET> &p2, IT p3x, IT p3y);

template <typename IT, typename ET>
Sign orient2D_IIE_exact(const GenericPoint2T<IT, ET> &p1,
                        const GenericPoint2T<IT, ET> &p2, ET p3x, ET p3y);

template <typename IT, typename ET>
Sign orient2D_IIE_expansion(const GenericPoint2T<IT, ET> &p1,
                            const GenericPoint2T<IT, ET> &p2, double p3x,
                            double p3y);

template <typename IT, typename ET, bool WithSSFilter>
Sign orient2D_IIE(const GenericPoint2T<IT, ET> &p1,
                  const GenericPoint2T<IT, ET> &p2, double p3x, double p3y,
                  PntArr2 arr);

template <typename IT, typename ET, bool WithSSFilter>
Sign orient2D_IIE(const GenericPoint2T<IT, ET> &p1,
                  const GenericPoint2T<IT, ET> &p2,
                  const GenericPoint2T<IT, ET> &p3, PntArr2 arr);

template <typename IT, typename ET>
Sign orient2D_III_filtered(const GenericPoint2T<IT, ET> &p1,
                           const GenericPoint2T<IT, ET> &p2,
                           const GenericPoint2T<IT, ET> &p3, PntArr2 arr);

template <typename IT, typename ET>
Sign orient2D_III_interval(const GenericPoint2T<IT, ET> &p1,
                           const GenericPoint2T<IT, ET> &p2,
                           const GenericPoint2T<IT, ET> &p3);

template <typename IT, typename ET>
Sign orient2D_III_exact(const GenericPoint2T<IT, ET> &p1,
                        const GenericPoint2T<IT, ET> &p2,
                        const GenericPoint2T<IT, ET> &p3);

template <typename IT, typename ET>
Sign orient2D_III_expansion(const GenericPoint2T<IT, ET> &p1,
                            const GenericPoint2T<IT, ET> &p2,
                            const GenericPoint2T<IT, ET> &p3);

template <typename IT, typename ET, bool WithSSFilter>
Sign orient2D_III(const GenericPoint2T<IT, ET> &p1,
                  const GenericPoint2T<IT, ET> &p2,
                  const GenericPoint2T<IT, ET> &p3, PntArr2 arr);

template <typename IT, typename ET>
Sign orient3D_IEEE_filtered(const GenericPoint3T<IT, ET> &p1, double ax,
                            double ay, double az, double bx, double by,
                            double bz, double cx, double cy, double cz,
                            PntArr3 arr);

template <typename IT, typename ET>
Sign orient3D_IEEE_interval(const GenericPoint3T<IT, ET> &p1, IT ax, IT ay,
                            IT az, IT bx, IT by, IT bz, IT cx, IT cy, IT cz);

template <typename IT, typename ET>
Sign orient3D_IEEE_exact(const GenericPoint3T<IT, ET> &p1, ET ax, ET ay, ET az,
                         ET bx, ET by, ET bz, ET cx, ET cy, ET cz);

template <typename IT, typename ET>
Sign orient3D_IEEE_expansion(const GenericPoint3T<IT, ET> &p1, double ax,
                             double ay, double az, double bx, double by,
                             double bz, double cx, double cy, double cz);

template <typename IT, typename ET, bool WithSSFilter>
Sign orient3D_IEEE(const GenericPoint3T<IT, ET> &p1, double ax, double ay,
                   double az, double bx, double by, double bz, double cx,
                   double cy, double cz, PntArr3 arr);

template <typename IT, typename ET, bool WithSSFilter>
Sign orient3D_IEEE(const GenericPoint3T<IT, ET> &p1,
                   const GenericPoint3T<IT, ET> &a,
                   const GenericPoint3T<IT, ET> &b,
                   const GenericPoint3T<IT, ET> &c, PntArr3 arr);

template <typename IT, typename ET>
Sign orient3D_IIEE_filtered(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2, double p3x,
                            double p3y, double p3z, double p4x, double p4y,
                            double p4z, PntArr3 arr);

template <typename IT, typename ET>
Sign orient3D_IIEE_interval(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2, IT p3x, IT p3y,
                            IT p3z, IT p4x, IT p4y, IT p4z);

template <typename IT, typename ET>
Sign orient3D_IIEE_exact(const GenericPoint3T<IT, ET> &p1,
                         const GenericPoint3T<IT, ET> &p2, ET p3x, ET p3y,
                         ET p3z, ET p4x, ET p4y, ET p4z);

template <typename IT, typename ET>
Sign orient3D_IIEE_expansion(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2, double p3x,
                             double p3y, double p3z, double p4x, double p4y,
                             double p4z);

template <typename IT, typename ET, bool WithSSFilter>
Sign orient3D_IIEE(const GenericPoint3T<IT, ET> &p1,
                   const GenericPoint3T<IT, ET> &p2, double p3x, double p3y,
                   double p3z, double p4x, double p4y, double p4z, PntArr3 arr);

template <typename IT, typename ET, bool WithSSFilter>
Sign orient3D_IIEE(const GenericPoint3T<IT, ET> &p1,
                   const GenericPoint3T<IT, ET> &p2,
                   const GenericPoint3T<IT, ET> &p3,
                   const GenericPoint3T<IT, ET> &p4, PntArr3 arr);

template <typename IT, typename ET>
Sign orient3D_IIIE_filtered(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2,
                            const GenericPoint3T<IT, ET> &p3, double p4x,
                            double p4y, double p4z, PntArr3 arr);

template <typename IT, typename ET>
Sign orient3D_IIIE_interval(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2,
                            const GenericPoint3T<IT, ET> &p3, IT p4x, IT p4y,
                            IT p4z);

template <typename IT, typename ET>
Sign orient3D_IIIE_exact(const GenericPoint3T<IT, ET> &p1,
                         const GenericPoint3T<IT, ET> &p2,
                         const GenericPoint3T<IT, ET> &p3, ET p4x, ET p4y,
                         ET p4z);

template <typename IT, typename ET>
Sign orient3D_IIIE_expansion(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2,
                             const GenericPoint3T<IT, ET> &p3, double p4x,
                             double p4y, double p4z);

template <typename IT, typename ET, bool WithSSFilter>
Sign orient3D_IIIE(const GenericPoint3T<IT, ET> &p1,
                   const GenericPoint3T<IT, ET> &p2,
                   const GenericPoint3T<IT, ET> &p3, double p4x, double p4y,
                   double p4z, PntArr3 arr);

template <typename IT, typename ET, bool WithSSFilter>
Sign orient3D_IIIE(const GenericPoint3T<IT, ET> &p1,
                   const GenericPoint3T<IT, ET> &p2,
                   const GenericPoint3T<IT, ET> &p3,
                   const GenericPoint3T<IT, ET> &p4, PntArr3 arr);

template <typename IT, typename ET>
Sign orient3D_IIII_filtered(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2,
                            const GenericPoint3T<IT, ET> &p3,
                            const GenericPoint3T<IT, ET> &p4, PntArr3 arr);

template <typename IT, typename ET>
Sign orient3D_IIII_interval(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2,
                            const GenericPoint3T<IT, ET> &p3,
                            const GenericPoint3T<IT, ET> &p4);

template <typename IT, typename ET>
Sign orient3D_IIII_exact(const GenericPoint3T<IT, ET> &p1,
                         const GenericPoint3T<IT, ET> &p2,
                         const GenericPoint3T<IT, ET> &p3,
                         const GenericPoint3T<IT, ET> &p4);

template <typename IT, typename ET>
Sign orient3D_IIII_expansion(const GenericPoint3T<IT, ET> &p1,
                             const GenericPoint3T<IT, ET> &p2,
                             const GenericPoint3T<IT, ET> &p3,
                             const GenericPoint3T<IT, ET> &p4);

template <typename IT, typename ET, bool WithSSFilter>
Sign orient3D_IIII(const GenericPoint3T<IT, ET> &p1,
                   const GenericPoint3T<IT, ET> &p2,
                   const GenericPoint3T<IT, ET> &p3,
                   const GenericPoint3T<IT, ET> &p4, PntArr3 arr);

template <typename IT, typename ET>
Sign orientOn2Dxy_IEE_filtered(const GenericPoint3T<IT, ET> &p1, double p2x,
                               double p2y, double p3x, double p3y, PntArr3 arr);

template <typename IT, typename ET>
Sign orientOn2Dxy_IEE_interval(const GenericPoint3T<IT, ET> &p1, IT p2x, IT p2y,
                               IT p3x, IT p3y);

template <typename IT, typename ET>
Sign orientOn2Dxy_IEE_exact(const GenericPoint3T<IT, ET> &p1, ET p2x, ET p2y,
                            ET p3x, ET p3y);

template <typename IT, typename ET>
Sign orientOn2Dxy_IEE_expansion(const GenericPoint3T<IT, ET> &p1, double p2x,
                                double p2y, double p3x, double p3y);

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dxy_IEE(const GenericPoint3T<IT, ET> &p1, double p2x, double p2y,
                      double p3x, double p3y, PntArr3 arr);

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dxy_IEE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &p3, PntArr3 arr);

template <typename IT, typename ET>
Sign orientOn2Dxy_IIE_filtered(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2, double op3x,
                               double op3y, PntArr3 arr);

template <typename IT, typename ET>
Sign orientOn2Dxy_IIE_interval(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2, IT op3x,
                               IT op3y);

template <typename IT, typename ET>
Sign orientOn2Dxy_IIE_exact(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2, ET op3x, ET op3y);

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dxy_IIE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2, double op3x,
                      double op3y, PntArr3 arr);

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dxy_IIE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &op3, PntArr3 arr);

template <typename IT, typename ET>
Sign orientOn2Dxy_III_filtered(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3, PntArr3 arr);

template <typename IT, typename ET>
Sign orientOn2Dxy_III_interval(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3);

template <typename IT, typename ET>
Sign orientOn2Dxy_III_exact(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2,
                            const GenericPoint3T<IT, ET> &p3);

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dxy_III(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &p3, PntArr3 arr);

template <typename IT, typename ET>
Sign orientOn2Dyz_IEE_filtered(const GenericPoint3T<IT, ET> &p1, double p2y,
                               double p2z, double p3y, double p3z, PntArr3 arr);

template <typename IT, typename ET>
Sign orientOn2Dyz_IEE_interval(const GenericPoint3T<IT, ET> &p1, IT p2y, IT p2z,
                               IT p3y, IT p3z);

template <typename IT, typename ET>
Sign orientOn2Dyz_IEE_exact(const GenericPoint3T<IT, ET> &p1, ET p2y, ET p2z,
                            ET p3y, ET p3z);

template <typename IT, typename ET>
Sign orientOn2Dyz_IEE_expansion(const GenericPoint3T<IT, ET> &p1, double p2y,
                                double p2z, double p3y, double p3z);

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dyz_IEE(const GenericPoint3T<IT, ET> &p1, double p2y, double p2z,
                      double p3y, double p3z, PntArr3 arr);

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dyz_IEE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &p3, PntArr3 arr);

template <typename IT, typename ET>
Sign orientOn2Dyz_IIE_filtered(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2, double op3y,
                               double op3z, PntArr3 arr);

template <typename IT, typename ET>
Sign orientOn2Dyz_IIE_interval(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2, IT op3y,
                               IT op3z);

template <typename IT, typename ET>
Sign orientOn2Dyz_IIE_exact(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2, ET op3y, ET op3z);

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dyz_IIE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2, double op3y,
                      double op3z, PntArr3 arr);

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dyz_IIE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &op3, PntArr3 arr);

template <typename IT, typename ET>
Sign orientOn2Dyz_III_filtered(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3, PntArr3 arr);

template <typename IT, typename ET>
Sign orientOn2Dyz_III_interval(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3);

template <typename IT, typename ET>
Sign orientOn2Dyz_III_exact(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2,
                            const GenericPoint3T<IT, ET> &p3);

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dyz_III(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &p3, PntArr3 arr);

template <typename IT, typename ET>
Sign orientOn2Dzx_IEE_filtered(const GenericPoint3T<IT, ET> &p1, double p2x,
                               double p2z, double p3x, double p3z, PntArr3 arr);

template <typename IT, typename ET>
Sign orientOn2Dzx_IEE_interval(const GenericPoint3T<IT, ET> &p1, IT p2x, IT p2z,
                               IT p3x, IT p3z);

template <typename IT, typename ET>
Sign orientOn2Dzx_IEE_exact(const GenericPoint3T<IT, ET> &p1, ET p2x, ET p2z,
                            ET p3x, ET p3z);

template <typename IT, typename ET>
Sign orientOn2Dzx_IEE_expansion(const GenericPoint3T<IT, ET> &p1, double p2x,
                                double p2z, double p3x, double p3z);

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dzx_IEE(const GenericPoint3T<IT, ET> &p1, double p2x, double p2z,
                      double p3x, double p3z, PntArr3 arr);

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dzx_IEE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &p3, PntArr3 arr);

template <typename IT, typename ET>
Sign orientOn2Dzx_IIE_filtered(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2, double op3x,
                               double op3z, PntArr3 arr);

template <typename IT, typename ET>
Sign orientOn2Dzx_IIE_interval(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2, IT op3x,
                               IT op3z);

template <typename IT, typename ET>
Sign orientOn2Dzx_IIE_exact(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2, ET op3x, ET op3z);

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dzx_IIE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2, double op3x,
                      double op3z, PntArr3 arr);

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dzx_IIE(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &op3, PntArr3 arr);

template <typename IT, typename ET>
Sign orientOn2Dzx_III_filtered(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3, PntArr3 arr);

template <typename IT, typename ET>
Sign orientOn2Dzx_III_interval(const GenericPoint3T<IT, ET> &p1,
                               const GenericPoint3T<IT, ET> &p2,
                               const GenericPoint3T<IT, ET> &p3);

template <typename IT, typename ET>
Sign orientOn2Dzx_III_exact(const GenericPoint3T<IT, ET> &p1,
                            const GenericPoint3T<IT, ET> &p2,
                            const GenericPoint3T<IT, ET> &p3);

template <typename IT, typename ET, bool WithSSFilter>
Sign orientOn2Dzx_III(const GenericPoint3T<IT, ET> &p1,
                      const GenericPoint3T<IT, ET> &p2,
                      const GenericPoint3T<IT, ET> &p3, PntArr3 arr);

template <typename IT, typename ET>
Sign squareDistance2D_IE_interval(const GenericPoint2T<IT, ET> &p, IT qx, IT qy,
                                  IT dis);

template <typename IT, typename ET>
Sign squareDistance2D_IE_exact(const GenericPoint2T<IT, ET> &p, ET qx, ET qy,
                               ET dis);

template <typename IT, typename ET>
Sign squareDistance2D_IE_expansion(const GenericPoint2T<IT, ET> &p, double qx,
                                   double qy, double dis);

template <typename IT, typename ET>
Sign squareDistance2D_IE(const GenericPoint2T<IT, ET> &p, double qx, double qy,
                         double dis);

template <typename IT, typename ET>
Sign squareDistance2D_IE(const GenericPoint2T<IT, ET> &p,
                         const GenericPoint2T<IT, ET> &q, double dis);

template <typename IT, typename ET>
Sign squareDistance2D_II_interval(const GenericPoint2T<IT, ET> &p,
                                  const GenericPoint2T<IT, ET> &q, IT dis);

template <typename IT, typename ET>
Sign squareDistance2D_II_exact(const GenericPoint2T<IT, ET> &p,
                               const GenericPoint2T<IT, ET> &q, ET dis);

template <typename IT, typename ET>
Sign squareDistance2D_II_expansion(const GenericPoint2T<IT, ET> &p,
                                   const GenericPoint2T<IT, ET> &q, double dis);

template <typename IT, typename ET>
Sign squareDistance2D_II(const GenericPoint2T<IT, ET> &p,
                         const GenericPoint2T<IT, ET> &q, double dis);

template <typename IT, typename ET>
Sign squareDistance3D_IE_interval(const GenericPoint3T<IT, ET> &p, IT qx, IT qy,
                                  IT qz, IT dis);

template <typename IT, typename ET>
Sign squareDistance3D_IE_exact(const GenericPoint3T<IT, ET> &p, ET qx, ET qy,
                               ET qz, ET dis);

template <typename IT, typename ET>
Sign squareDistance3D_IE_expansion(const GenericPoint3T<IT, ET> &p, double qx,
                                   double qy, double qz, double dis);

template <typename IT, typename ET>
Sign squareDistance3D_IE(const GenericPoint3T<IT, ET> &p, double qx, double qy,
                         double qz, double dis);

template <typename IT, typename ET>
Sign squareDistance3D_IE(const GenericPoint3T<IT, ET> &p,
                         const GenericPoint3T<IT, ET> &q, double dis);

template <typename IT, typename ET>
Sign squareDistance3D_II_interval(const GenericPoint3T<IT, ET> &p,
                                  const GenericPoint3T<IT, ET> &q, IT dis);

template <typename IT, typename ET>
Sign squareDistance3D_II_exact(const GenericPoint3T<IT, ET> &p,
                               const GenericPoint3T<IT, ET> &q, ET dis);

template <typename IT, typename ET>
Sign squareDistance3D_II_expansion(const GenericPoint3T<IT, ET> &p,
                                   const GenericPoint3T<IT, ET> &q, double dis);

template <typename IT, typename ET>
Sign squareDistance3D_II(const GenericPoint3T<IT, ET> &p,
                         const GenericPoint3T<IT, ET> &q, double dis);

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "IndirectPredicateDetails.inl"
#endif