
#ifndef SMEAR_H
#define SMEAR_H

#ifndef SU3_H
#include "su3.h"
#endif

#define _su3_lcomb(z,a,u,b,v) \
  (z).c11.re=((a)*(u).c11.re+(b)*(v).c11.re);\
  (z).c11.im=((a)*(u).c11.im+(b)*(v).c11.im);\
  (z).c12.re=((a)*(u).c12.re+(b)*(v).c12.re);\
  (z).c12.im=((a)*(u).c12.im+(b)*(v).c12.im);\
  (z).c13.re=((a)*(u).c13.re+(b)*(v).c13.re);\
  (z).c13.im=((a)*(u).c13.im+(b)*(v).c13.im);\
  (z).c21.re=((a)*(u).c21.re+(b)*(v).c21.re);\
  (z).c21.im=((a)*(u).c21.im+(b)*(v).c21.im);\
  (z).c22.re=((a)*(u).c22.re+(b)*(v).c22.re);\
  (z).c22.im=((a)*(u).c22.im+(b)*(v).c22.im);\
  (z).c23.re=((a)*(u).c23.re+(b)*(v).c23.re);\
  (z).c23.im=((a)*(u).c23.im+(b)*(v).c23.im);\
  (z).c31.re=((a)*(u).c31.re+(b)*(v).c31.re);\
  (z).c31.im=((a)*(u).c31.im+(b)*(v).c31.im);\
  (z).c32.re=((a)*(u).c32.re+(b)*(v).c32.re);\
  (z).c32.im=((a)*(u).c32.im+(b)*(v).c32.im);\
  (z).c33.re=((a)*(u).c33.re+(b)*(v).c33.re);\
  (z).c33.im=((a)*(u).c33.im+(b)*(v).c33.im)

#define _su3_add(u,v,w) \
   (u).c11.re=(v).c11.re+(w).c11.re;\
   (u).c11.im=(v).c11.im+(w).c11.im;\
   (u).c12.re=(v).c12.re+(w).c12.re;\
   (u).c12.im=(v).c12.im+(w).c12.im;\
   (u).c13.re=(v).c13.re+(w).c13.re;\
   (u).c13.im=(v).c13.im+(w).c13.im;\
   (u).c21.re=(v).c21.re+(w).c21.re;\
   (u).c21.im=(v).c21.im+(w).c21.im;\
   (u).c22.re=(v).c22.re+(w).c22.re;\
   (u).c22.im=(v).c22.im+(w).c22.im;\
   (u).c23.re=(v).c23.re+(w).c23.re;\
   (u).c23.im=(v).c23.im+(w).c23.im;\
   (u).c31.re=(v).c31.re+(w).c31.re;\
   (u).c31.im=(v).c31.im+(w).c31.im;\
   (u).c32.re=(v).c32.re+(w).c32.re;\
   (u).c32.im=(v).c32.im+(w).c32.im;\
   (u).c33.re=(v).c33.re+(w).c33.re;\
   (u).c33.im=(v).c33.im+(w).c33.im

#ifndef HYP_LINKS_C
extern void free_hyp(void);
extern void hyp_time_links(double alpha1, double alpha2, double alpha3, int proj_su3, int piter, int keep);
extern void hyp_space_links(double alpha1, double alpha2, double alpha3, int proj_su3, int piter, int keep);
extern void hyp_links(double alpha1, double alpha2, double alpha3, int proj_su3, int piter, int keep);
#endif

#ifndef HYP_SPATIAL_LINKS_C
extern void free_hyp_spatial(void);
extern void hyp_spatial_links(su3_dble *ub, double alpha1, double alpha2, int proj_su3, int piter, int keep);
#endif

#ifndef APE_LINKS_C
extern void ape_links(su3_dble* ub, double a, int proj_su3, int piter, int keep);
extern void ape_spatial_links(su3_dble* ub, double a, int proj_su3, int piter, int keep);
#endif

#ifndef STAPLE_SUM_C
extern void staple_sum(su3_dble *u1,su3_dble *u2,int mu,int nu,su3_dble *vb,su3_dble *cbuf);
#endif

#ifndef PROJECTION_C
extern void approx_project_to_su3_dble(su3_dble *u, int iter);
#endif

#endif
