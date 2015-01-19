/*******************************************************************************
*
* File hyp_spatial_links.c
*
* Copyright (C) 2009 Bjorn Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Calculation of HYP-smeared spatial links
*
* The externally accessible functions are
*
*   void free_hyp_spatial()
*     Frees all the allocated memory.
*
*   void hyp_spatial_links(su3_dble* ub, double alpha1, double alpha2,
*                          int proj_su3, int piter, int keep)
*     Replaces all 3-d double-precision spatial link variables by the HYP-smeared
*     ones. alpha1-2 are the two parameters of the HYP-smearing and for proj_su3
*     unequal zero a projection back to SU(3) is performed after each blocking
*     step. piter is a iteration number controlling the projection methode used.
*     For keep unequal zero the allocated memory is NOT freed at the end.
*     ub is either the base address of a link-field allocated by alloc_lnk or the
*     address returned by udfld().
*
* Notes:
*
* The definition of the HYP-smeared link can be found in hep-lat/0103029. The
* approximate projection to SU(3) used here is defined in hep-lat/0506008.
*
* The program is fully parallelised.
*
* The required communications (including buffer allocation) are performed
* automatically using the programs in the module lnk.
*
*******************************************************************************/

#define HYP_SPATIAL_LINKS_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "lattice.h"
#include "global.h"
#include "smear.h"
#include "uflds.h"
#include "flags.h"
#include "lnk.h"


typedef struct
{
   int iup;
   int idn;
   int iuk;
   int fc;
} bnd_idx_t;

static double alpha[2],beta[2];
static int proj,iter;
static const su3_dble v0={{0.0}};
static su3_dble *dlnks[3]={NULL},*buf,*pub;
static int dlnks_idx[4][4];
static bnd_idx_t *bnd_idx=NULL;


static void alloc_bnd_idx(void)
{
   int nu,ix,iy,iz;
   bnd_idx_t *bx;
   
   if (bnd_idx==NULL)
   {
      bnd_idx=amalloc((BNDRY/2)*sizeof(bnd_idx_t),3);
      error(bnd_idx==NULL,1,"alloc_bufs [hyp_spatial_links.c]","Unable to allocate index array");
      bx=bnd_idx;
      for (ix=0;ix<VOLUME;ix++)
      {
         for (nu=0;nu<4;nu++)
         {
            iy=iup[ix][nu];
            if (iy>=VOLUME)
            {
               iz=map[iy-VOLUME];
               (*bx).idn=ix;
               (*bx).iup=iz;
               (*bx).iuk=iy;
               (*bx).fc=nu;
               bx+=1;
            }
         }
      }
   }
}


static void alloc_bufs(void)
{
   int mu,nu;
   
   if (dlnks[0]==NULL)
   {
      for (mu=0;mu<3;mu++)
         dlnks[mu]=alloc_lnk();
      
      if (NPROC>1)
      {
         buf=alloc_buf_uk();
         
         alloc_bnd_idx();
      }
   }
   
   for (nu=0;nu<4;nu++)
      for (mu=0;mu<4;mu++)
         dlnks_idx[nu][mu]=0;
   dlnks_idx[1][2]=0;
   dlnks_idx[1][3]=1;
   dlnks_idx[2][1]=0;
   dlnks_idx[2][3]=1;
   dlnks_idx[3][1]=0;
   dlnks_idx[3][2]=1;
}


static void set_zero(int vol, su3_dble *u)
{
   su3_dble *v,*vm;
   
   v=u;
   vm=v+vol;
   for (;v<vm;v++)
      *v=v0;
}


static void decorate_all(int ia, su3_dble *vb)
{
   int ix,mu;
   su3_dble *v,*u;
   
   for (ix=0;ix<VOLUME;ix++)
   {
      for (mu=1;mu<4;mu++)
      {
         u=lnk(pub,ix,mu);
         v=lnk(vb,ix,mu);
         
         _su3_lcomb(*v,alpha[ia-1],*u,beta[ia-1],*v);
         if (proj)
            approx_project_to_su3_dble(v,iter);
      }
   }
}


static void decorated_links2(void)
{
   int mu,nu,i,n;
   su3_dble *vo,*u1,*u2;
   bnd_idx_t *bx,*bxm;
   
   vo=dlnks[2];
   set_zero(LOCAL,vo);
   if (NPROC>1)
      set_zero(BNDRY_K,buf);
   
   for (mu=1;mu<4;mu++)
   {
      for (i=0;i<2;i++)
      {
         nu=((mu+i)%3)+1;
         n=dlnks_idx[nu][mu];
         u1=dlnks[n];
         n=dlnks_idx[mu][nu];
         u2=dlnks[n];
         
         staple_sum(u1,u2,mu,nu,vo,buf);
      }
   }
   
   if (NPROC>1)
   { 
      send_uk(dlnks[0],buf);
   
      bx=bnd_idx;
      bxm=bnd_idx+BNDRY/2;
      for (;bx<bxm;bx++)
      {
         for (mu=1;mu<4;mu++)
         {
            if (mu!=(*bx).fc)
            {
               u1=lnk(dlnks[2],(*bx).iup,mu);
               u2=lnkf(dlnks[0],(*bx).iuk,mu,(*bx).fc);
               _su3_add(*u1,*u1,*u2);
            }
         }
      }
   }
   
   decorate_all(1,vo);
   copyback_u0(vo);
}


static void decorated_links1(void)
{
   int mu,eta,i,n;
   su3_dble *vo,*u;
   su3_dble *v,*w;
   bnd_idx_t *bx,*bxm;
   
   u=dlnks[2];
   
   for (i=0;i<2;i++)
       set_zero(LOCAL,dlnks[i]);
   
   if (NPROC>1)
      set_zero(BNDRY_K,buf);
   
   for (mu=1;mu<4;mu++)
   {
      for (i=0;i<2;i++)
      {
         eta=((mu+i)%3)+1;
         n=dlnks_idx[mu][eta];
         vo=dlnks[n];
         
         staple_sum(u,u,mu,eta,vo,buf);
      }
   }
   
   if (NPROC>1)
   { 
      send_uk(dlnks[2],buf);
      
      bx=bnd_idx;
      bxm=bnd_idx+BNDRY/2;
      for (;bx<bxm;bx++)
      {
         for (mu=1;mu<4;mu++)
         {
            if (mu!=(*bx).fc)
            {
               n=dlnks_idx[mu][(*bx).fc];
                     
               v=lnk(dlnks[n],(*bx).iup,mu);
               w=lnkf(dlnks[2],(*bx).iuk,mu,(*bx).fc);
               _su3_add(*v,*v,*w);
            }
         }
      }
   }
      
   for (i=0;i<2;i++)
   {
      vo=dlnks[i];
      decorate_all(2,vo);
      copy_uk(vo);
   }
}


void free_hyp_spatial(void)
{
   int i;
   
   for (i=0;i<3;i++)
   {
      free_lnk(dlnks[i],0);
      dlnks[i]=NULL;
   }
   free_lnk(NULL,1);
   
   if (NPROC>1)
   {
      free_buf_uk(buf);
      
      afree(bnd_idx);
      bnd_idx=NULL;
   }
}


void hyp_spatial_links(su3_dble* ub, double alpha1, double alpha2, int proj_su3,
                       int piter, int keep)
{
   int i,mu;
   double a[2];
   su3_dble *v,*u,*um,*udb;
   
   udb=udfld();
   pub=ub;
   if (pub==udb)
      copy_bnd_ud();
   else
      copy_bnd(pub);
   
   alloc_bufs();
   assign_lnk2lnk(pub,dlnks[2]);
   copy_bnd(dlnks[2]);
   
   a[0]=alpha1;
   a[1]=alpha2;

   for (i=0; i<2; i++)
   {
      alpha[i]=1.0-a[i];
      beta[i]=a[i]/(double)(2*(2-i));
   }

   proj=proj_su3;
   iter=piter;
   
   decorated_links1();
   decorated_links2();
   
   u=pub;
   um=u+4*VOLUME;
   v=dlnks[2];
   for (;u<um;)
   {
      for (mu=1;mu<4;mu++)
      {
         *(u+2*mu)=*(v+2*mu);
         *(u+2*mu+1)=*(v+2*mu+1);
      }
      u+=8;
      v+=8;
   }
   
   if (pub==udb)
   {
      set_flags(UPDATED_UD);
      set_bc();
   }
   
   if (keep==0)
      free_hyp_spatial();
}
