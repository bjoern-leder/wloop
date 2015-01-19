/*******************************************************************************
*
* File lnkplaq.c
*
* Copyright (C) 2009 Bjorn Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Calculation of the plaquette sum
*
* The externally accessible functions are
*
*   double lnk_plaq_sum(su3_dble *ub)
*     Returns the sum of Re(tr{U(p)}) over all 6*VOLUME*NPROC unoriented
*     plaquettes p, where U(p) is the product of the double-precision link
*     variables around p. The links are taken from the link-field ub.
*
* Notes:
*
* The return values are exactly the same on all processes and the required
* communications (including buffer allocation) are performed automatically
* using the programs in the module lnkcom.c
*
*******************************************************************************/

#define LNKPLAQ_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "lattice.h"
#include "flags.h"
#include "global.h"
#include "lnk.h"

#define MAX_LEVELS 8
#define BLK_LENGTH 8
#define N0 (NPROC0*L0)

static int cnt[MAX_LEVELS];
static double smx[MAX_LEVELS];
static su3_dble wd1,wd2,*ud[4];


static double plaq_dble(su3_dble *ub, int ix,int mu,int nu)
{
   double sm;
   su3_vector_dble *v1,*v2;
   
   lnk_plaq(ub,ix,mu,nu,ud);
   
   _su3_times_su3(wd1,(*ud[0]),(*ud[1]));
   _su3_times_su3(wd2,(*ud[2]),(*ud[3]));
   
   v1=(su3_vector_dble*)(&wd1.c11);
   v2=(su3_vector_dble*)(&wd2.c11);
   sm=_vector_prod_re((*(v1  )),(*(v2  )))+
   _vector_prod_re((*(v1+1)),(*(v2+1)))+
   _vector_prod_re((*(v1+2)),(*(v2+2)));
   
   return sm;
}


static double local_plaq_sum_dble(su3_dble *ub)
{
   int ix,mu,nu,n,bc,t;
   double pa;
   
   bc=bc_type();

   for (n=0;n<MAX_LEVELS;n++)
   {
      cnt[n]=0;
      smx[n]=0.0;
   }
   
   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);
      pa=0.0;
      
      for (mu=1;mu<4;mu++)
      {
         if ((t<(N0-1))||(bc!=0))
            pa+=plaq_dble(ub,ix,mu,0);

         for (nu=1;nu<mu;nu++)
            pa+=plaq_dble(ub,ix,mu,nu);
      }
      
      cnt[0]+=1;
      smx[0]+=pa;
      
      for (n=1;(cnt[n-1]>=BLK_LENGTH)&&(n<MAX_LEVELS);n++)
      {
         cnt[n]+=1;
         smx[n]+=smx[n-1];
         
         cnt[n-1]=0;
         smx[n-1]=0.0;
      }               
   }
   
   pa=0.0;
   
   for (n=0;n<MAX_LEVELS;n++)
      pa+=smx[n]; 
   
   return pa;
}


double lnk_plaq_sum(su3_dble *ub)
{
   double p,pa;
   
   copy_bnd(ub);
   
   p=local_plaq_sum_dble(ub);
   MPI_Reduce(&p,&pa,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   MPI_Bcast(&pa,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   
   return pa;
}
