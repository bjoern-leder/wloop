
#ifndef WLOOP_H
#define WLOOP_H

#ifndef GLOBAL_H
#include "global.h"
#endif

#ifndef SU3_H
#include "su3.h"
#endif

#define MAX_SM_LEVELS (5)
#define SLICE0        (L1*L2*L3)

#define wloop_RELEASE "wloop-1.4"

typedef enum
{
   APE,HYP,
   OP_SMEARS
} op_smear_t;

typedef struct
{
   int msl,nss[MAX_SM_LEVELS];
   int mwlt,mwlr,wls;
   op_smear_t op_smear;
   double alpha_action[3];
   double alpha_ape;
   double alpha_hyp[2];
   int proj;
   int proj_iter;
   double *wl_mat;
} wloop_parms_t;


#define pindx(sl,ix,k) \
   ((k)+3*((ix)+VOLUME*(sl)))

#define pindxs(sl,ixs,x0,k) \
   ((k)+3*((ipt[(ixs)+SLICE0*(x0)])+VOLUME*(sl)))
   
   
#ifndef WLOOP_C
extern void free_wloop(void);
extern void wloop_sum(int keep);
#endif

#ifndef WLOOP_PARMS_C
extern wloop_parms_t set_wloop_parms(int msl,int nss[],int mwlt,int mwlr,op_smear_t op_smear,
                              double alpha_action[3],double alpha_ape, double alpha_hyp[2],
                              int proj_iter);
extern wloop_parms_t wloop_parms(void);
extern void read_wloop_parms(void);
extern void print_wloop_parms(void);
extern void write_wloop_parms(FILE *fdat);
extern void check_wloop_parms(FILE *fdat);
extern int wl_idx(int sl1, int sl2, int t, int r, int x0);
#endif

#ifndef SMEAR_C
extern void smear(void);
extern su3_dble** get_pus(void);
extern void free_smear(void);
#endif

#endif
