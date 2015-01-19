
#ifndef LNK_H
#define LNK_H

#ifndef GLOBAL_H
#include "global.h"
#endif

#ifndef SU3_H
#include "su3.h"
#endif

#define LOCAL ((4*VOLUME)+(BNDRY/2))
#define BNDRY_K (3*(BNDRY/2))

#ifndef LNKCOM_C
extern void copy_bnd(su3_dble *ub);
extern void copy_u0(su3_dble *ub);
extern void copyback_u0(su3_dble *ub);
extern void copy_u0n(su3_dble *ub);
extern void copyback_u0n(su3_dble *ub);
extern void copy_uk(su3_dble *ub);
extern void copyback_uk(su3_dble *ub);
extern void send_uk(su3_dble *ub, su3_dble *sbuf);
extern void free_lnk(su3_dble *ub, int cb);
extern void free_buf_uk(su3_dble *ub);
extern su3_dble* alloc_lnk(void);
extern su3_dble* alloc_buf_uk(void);
extern su3_dble* lnk_buf(su3_dble* ub, int ix, int mu, int face);
extern su3_dble* lnkf(su3_dble* ub, int ix, int mu, int face);
extern su3_dble* lnk(su3_dble* ub, int ix, int mu);
extern void lnk_plaq(su3_dble* ub,int ix,int mu,int nu,su3_dble **u);
extern void assign_lnk2lnk(su3_dble* ub1,su3_dble* ub2);
extern void assign_ud2lnk(su3_dble* ub);
extern void assign_lnk2ud(su3_dble* ub);
extern int* get_lnk_idx_u0(void);
extern int* get_lnk_idx_uk(void);
#endif

#ifndef LNKSTAPLE_C
extern void add_up_staple(su3_dble *u1o,su3_dble *u2o,int ix,int mu,int nu,su3_dble *v);
extern void add_dn_staple(su3_dble *u1o,su3_dble *u2o,int ix,int mu,int nu,su3_dble *v);
extern void add_up_staple_lnks(su3_dble *u1,su3_dble *u2,su3_dble *u3,su3_dble *v);
extern void add_dn_staple_lnks(su3_dble *u1,su3_dble *u2,su3_dble *u3,su3_dble *v);
extern void add_up_staple_save(su3_dble *u1o,su3_dble *u2o,int ix,int mu,int nu,su3_dble *v);
extern void add_dn_staple_save(su3_dble *u1o,su3_dble *u2o,int ix,int mu,int nu,su3_dble *v);
extern void add_up_staple_lnks_save(su3_dble *u1,su3_dble *u2,su3_dble *u3,su3_dble *v);
extern void add_dn_staple_lnks_save(su3_dble *u1,su3_dble *u2,su3_dble *u3,su3_dble *v);
#endif

#ifndef LNKPLAQ_C
extern double lnk_plaq_sum(su3_dble *ub);
#endif

#endif
