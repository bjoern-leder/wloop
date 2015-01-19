/*******************************************************************************
*
* File ape_links.c
*
* Copyright (C) 2009 Bjorn Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Calculation of APE-smeared links
*
* The externally accessible functions are
*
*   void free_ape_links()
*     Frees all the allocated memory.
*
*   void ape_links(su3_dble* ub, double a, int proj_su3, int piter, int keep)
*     Replaces all 3-d double-precision link variables by the APE-smeared
*     ones. alpha is the parameters of the APE-smearing and for proj_su3
*     unequal zero a projection back to SU(3) is performed. piter is a iteration
*     number controlling the projection methode used.
*     For keep unequal zero the allocated memory is NOT freed at the end.
*     ub is either the base address of a link-field allocated by alloc_lnk or
*     the address returned by udfld().
*
*   void ape_spatial_links(su3_dble* ub, double a, int proj_su3, int piter, int keep)
*     Same as above but only spatial links are used and smeared.
*
* Notes:
*
* The program is fully parallelised.
*
* The required communications (including buffer allocation) are performed
* automatically using the programs in the module lnk.
*
*******************************************************************************/

#define APE_LINKS_C

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


#define _su3_mulr(z,a,u) \
  (z).c11.re=(a)*(u).c11.re;\
  (z).c11.im=(a)*(u).c11.im;\
  (z).c12.re=(a)*(u).c12.re;\
  (z).c12.im=(a)*(u).c12.im;\
  (z).c13.re=(a)*(u).c13.re;\
  (z).c13.im=(a)*(u).c13.im;\
  (z).c21.re=(a)*(u).c21.re;\
  (z).c21.im=(a)*(u).c21.im;\
  (z).c22.re=(a)*(u).c22.re;\
  (z).c22.im=(a)*(u).c22.im;\
  (z).c23.re=(a)*(u).c23.re;\
  (z).c23.im=(a)*(u).c23.im;\
  (z).c31.re=(a)*(u).c31.re;\
  (z).c31.im=(a)*(u).c31.im;\
  (z).c32.re=(a)*(u).c32.re;\
  (z).c32.im=(a)*(u).c32.im;\
  (z).c33.re=(a)*(u).c33.re;\
  (z).c33.im=(a)*(u).c33.im


typedef struct
{
   int iup;
   int idn;
   int iuk;
   int fc;
} bnd_idx_t;

static const su3_dble v0={{0.0}};
static su3_dble *dlnks[2]={NULL},*buf,*pub;
static bnd_idx_t *bnd_idx=NULL;


static void alloc_bnd_idx(void)
{
   int nu,ix,iy,iz;
   bnd_idx_t *bx;
   
   if (bnd_idx==NULL)
   {
      bnd_idx=amalloc((BNDRY/2)*sizeof(bnd_idx_t),3);
      error(bnd_idx==NULL,1,"alloc_bufs [ape_links.c]","Unable to allocate index array");
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
   int mu;
   
   if (dlnks[0]==NULL)
   {
      for (mu=0;mu<2;mu++)
         dlnks[mu]=alloc_lnk();
      
      if (NPROC>1)
      {
         buf=alloc_buf_uk();
         
         alloc_bnd_idx();
      }
   }
}


static void set_zero(int vol, su3_dble *u)
{
   su3_dble *v,*vm;
   
   v=u;
   vm=v+vol;
   for (;v<vm;v++)
      *v=v0;
}


static void decorated_links(int time)
{
   int mu,eta;
   su3_dble *vo,*u;
   su3_dble *v,*w;
   bnd_idx_t *bx,*bxm;
   
   u=dlnks[1];
   vo=dlnks[0];
   
   set_zero(LOCAL,vo);
   
   if (NPROC>1)
      set_zero(BNDRY_K,buf);
   
   for (mu=(1-(time>0));mu<4;mu++)
   {
      for (eta=(1-(time>0));eta<4;eta++)
      {
         if (eta!=mu)
            staple_sum(u,u,mu,eta,vo,buf);
      }
   }
   
   if (NPROC>1)
   { 
      send_uk(u,buf);
      
      bx=bnd_idx;
      bxm=bnd_idx+BNDRY/2;
      for (;bx<bxm;bx++)
      {
         for (mu=(1-(time>0));mu<4;mu++)
         {
            if (mu!=(*bx).fc)
            {
               v=lnk(vo,(*bx).iup,mu);
               w=lnkf(u,(*bx).iuk,mu,(*bx).fc);
               _su3_add(*v,*v,*w);
            }
         }
      }
   }
      
   copy_uk(vo);
   copyback_u0(vo);
}


void free_ape_links(void)
{
   int i;
   
   for (i=0;i<2;i++)
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


static void ape_links_all(su3_dble* ub, double alpha, int proj_su3, int piter, int time, int keep)
{
   int mu,mu0;
   su3_dble *v,*u,*um,*udb;
   
   udb=udfld();
   pub=ub;
   if (pub==udb)
      copy_bnd_ud();
   else
      copy_bnd(pub);
   
   alloc_bufs();
   assign_lnk2lnk(pub,dlnks[1]);
   copy_bnd(dlnks[1]);
      
   decorated_links(time);
   
   u=pub;
   um=u+4*VOLUME;
   v=dlnks[0];
   mu0=(1-(time>0));
   for (;u<um;)
   {
      u+=2*mu0;
      v+=2*mu0;
      
      for (mu=mu0;mu<4;mu++)
      {
         _su3_mulr(*v,alpha,*v);
         _su3_add(*u,*u,*v);
         if (proj_su3)
            approx_project_to_su3_dble(u,piter);
         
         u+=1;
         v+=1;
         
         _su3_mulr(*v,alpha,*v);
         _su3_add(*u,*u,*v);
         if (proj_su3)
            approx_project_to_su3_dble(u,piter);
         
         u+=1;
         v+=1;
      }
   }
   
   if (pub==udb)
   {
      set_flags(UPDATED_UD);
      set_bc();
   }
      
   if (keep==0)
      free_ape_links();
}


void ape_links(su3_dble* ub, double a, int proj_su3, int piter, int keep)
{
   ape_links_all(ub,a,proj_su3,piter,1,keep);
}


void ape_spatial_links(su3_dble* ub, double a, int proj_su3, int piter, int keep)
{
   ape_links_all(ub,a,proj_su3,piter,0,keep);
}
