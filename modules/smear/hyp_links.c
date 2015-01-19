/*******************************************************************************
*
* File hyp_links.c
*
* Copyright (C) 2009 Bjorn Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Calculation of HYP-smeared links
*
* The externally accessible functions are
*
*   void free_hyp()
*     Frees all the allocated memory.
*
*   void hyp_links(double alpha1, double alpha2, double alpha3, int proj_su3,
*                  int piter, int keep)
*     Replaces all double-precision link variables by the HYP-smeared
*     ones. alpha1-3 are the three parameters of the HYP-smearing and for proj_su3
*     unequal zero a projection back to SU(3) is performed after each blocking
*     step. piter is a iteration number controlling the projection methode used.
*     For keep unequal zero the allocated memory is NOT freed at the end.
*
*   void hyp_time_links(double alpha1, double alpha2, double alpha3, int proj_su3,
*                       int piter, int keep)
*     Same as hyp_links, but only the temporal links are smeared.
*
*   void hyp_space_links(double alpha1, double alpha2, double alpha3, int proj_su3,
*                        int piter, int keep)
*     Same as hyp_links, but only the spatial links are smeared.
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

#define HYP_LINKS_C

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

static double alpha[3],beta[3];
static int proj,iter;
static const su3_dble v0={{0.0}};
static su3_dble *dlnks[6]={NULL},*buf[3],*udb;
static int dlnks_idx[4][4];
static bnd_idx_t *bnd_idx=NULL;


static void alloc_bnd_idx(void)
{
   int nu,ix,iy,iz;
   bnd_idx_t *bx;
   
   if (bnd_idx==NULL)
   {
      bnd_idx=amalloc((BNDRY/2)*sizeof(bnd_idx_t),3);
      error(bnd_idx==NULL,1,"alloc_bufs [hyp_links.c]","Unable to allocate index array");
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
      for (mu=0;mu<6;mu++)
         dlnks[mu]=alloc_lnk();
      
      if (NPROC>1)
      {
         for (mu=0;mu<3;mu++)
            buf[mu]=alloc_buf_uk();
         
         alloc_bnd_idx();
      }
   }
   
   for (nu=0;nu<4;nu++)
      for (mu=0;mu<4;mu++)
         dlnks_idx[nu][mu]=0;
      
   dlnks_idx[0][1]=0;
   dlnks_idx[0][2]=1;
   dlnks_idx[0][3]=2;
   dlnks_idx[1][0]=0;
   dlnks_idx[1][2]=2;
   dlnks_idx[1][3]=1;
   dlnks_idx[2][0]=1;
   dlnks_idx[2][1]=2;
   dlnks_idx[2][3]=0;
   dlnks_idx[3][0]=2;
   dlnks_idx[3][1]=1;
   dlnks_idx[3][2]=0;
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
      for (mu=0;mu<4;mu++)
      {
         u=lnk(udb,ix,mu);
         v=lnk(vb,ix,mu);
         
         _su3_lcomb(*v,alpha[ia-1],*u,beta[ia-1],*v);
         if (proj)
            approx_project_to_su3_dble(v,iter);
      }
   }
}


static void decorate(int ia, su3_dble *vb, int time, int space)
{
   int ix,mu;
   su3_dble *v,*u;
   
   for (ix=0;ix<VOLUME;ix++)
   {
      for (mu=(1-(time>0));mu<(1+3*(space>0));mu++)
      {
         u=lnk(udb,ix,mu);
         v=lnk(vb,ix,mu);
         
         _su3_lcomb(*v,alpha[ia-1],*u,beta[ia-1],*v);
         if (proj)
            approx_project_to_su3_dble(v,iter);
      }
   }
}


static void decorated_links3(int time, int space)
{
   int mu,nu,i,n;
   su3_dble *vo,*u1,*u2;
   bnd_idx_t *bx,*bxm;
   
   vo=dlnks[0];
   set_zero(LOCAL,vo);
   if (NPROC>1)
      set_zero(BNDRY_K,buf[0]);
   
   for (mu=(1-(time>0));mu<(1+3*(space>0));mu++)
   {
      for (i=1;i<4;i++)
      {
         nu=(mu+i)&0x3;
         n=dlnks_idx[nu][mu];
         u1=dlnks[(3+n)];
         n=dlnks_idx[mu][nu];
         u2=dlnks[(3+n)];
         
         staple_sum(u1,u2,mu,nu,vo,buf[0]);
      }
   }
   
   if (NPROC>1)
   { 
      send_uk(dlnks[1],buf[0]);
   
      bx=bnd_idx;
      bxm=bnd_idx+BNDRY/2;
      for (;bx<bxm;bx++)
      {
         for (mu=(1-(time>0));mu<(1+3*(space>0));mu++)
         {
            if (mu!=(*bx).fc)
            {
               u1=lnk(dlnks[0],(*bx).iup,mu);
               u2=lnkf(dlnks[1],(*bx).iuk,mu,(*bx).fc);
               _su3_add(*u1,*u1,*u2);
            }
         }
      }
   }
   
   decorate(1,vo,time,space);
   copyback_u0(vo);
}


static void decorated_links2(void)
{
   int mu,nu,rho,i,n,d;
   su3_dble *vo,*u1,*u2;
   su3_dble *v,*w;
   bnd_idx_t *bx,*bxm;
   
   for (i=0;i<3;i++)
   {
      set_zero(LOCAL,dlnks[(3+i)]);
      if (NPROC>1)
         set_zero(BNDRY_K,buf[i]);
   }
   
   for (mu=0;mu<4;mu++)
   {
      for (i=1;i<4;i++)
      {
         nu=(mu+i)&0x3;
         n=dlnks_idx[mu][nu];
         d=nu+mu;
         if (d>3)
            d=6-d;
         u1=dlnks[(d-1)];
         vo=dlnks[(3+n)];
         for (rho=0;rho<4;rho++)
         {
            if ((rho!=mu) && (rho!=nu))
            {
               d=nu+rho;
               if (d>3)
                  d=6-d;
               u2=dlnks[(d-1)];
         
               staple_sum(u1,u2,mu,rho,vo,buf[n]);
            }
         }
      }
   }
   
   if (NPROC>1)
   { 
      for (i=0;i<3;i++)
         send_uk(dlnks[i],buf[i]);
      
      bx=bnd_idx;
      bxm=bnd_idx+BNDRY/2;
      for (;bx<bxm;bx++)
      {
         for (i=1;i<4;i++)
         {
            mu=((*bx).fc+i)&0x3;
            for (nu=0;nu<4;nu++)
            {
               if (nu!=mu && nu!=(*bx).fc)
               {
                  n=dlnks_idx[mu][nu];
                        
                  v=lnk(dlnks[(3+n)],(*bx).iup,mu);
                  w=lnkf(dlnks[n],(*bx).iuk,mu,(*bx).fc);
                  _su3_add(*v,*v,*w);
               }
            }
         }
      }
   }
   
   for (i=0;i<3;i++)
   {
      vo=dlnks[(3+i)];
      decorate_all(2,vo);
      copy_uk(vo);
   }
}


static void decorated_links1(void)
{
   int mu,eta,i,n;
   su3_dble *vo,*u;
   su3_dble *v,*w;
   bnd_idx_t *bx,*bxm;
   
   u=dlnks[3];
   
   for (i=0;i<3;i++)
      set_zero(LOCAL,dlnks[i]);
   
   if (NPROC>1)
      set_zero(BNDRY_K,buf[0]);
   
   for (mu=0;mu<4;mu++)
   {
      for (i=1;i<4;i++)
      {
         eta=(mu+i)&0x3;
         n=dlnks_idx[mu][eta];
         vo=dlnks[n];
         
         staple_sum(u,u,mu,eta,vo,buf[0]);
      }
   }
   
   if (NPROC>1)
   { 
      send_uk(dlnks[3],buf[0]);
      
      bx=bnd_idx;
      bxm=bnd_idx+BNDRY/2;
      for (;bx<bxm;bx++)
      {
         for (i=1;i<4;i++)
         {
            mu=((*bx).fc+i)&0x3;
            n=dlnks_idx[mu][(*bx).fc];
                  
            v=lnk(dlnks[n],(*bx).iup,mu);
            w=lnkf(dlnks[3],(*bx).iuk,mu,(*bx).fc);
            _su3_add(*v,*v,*w);
         }
      }
   }
      
   for (i=0;i<3;i++)
   {
      vo=dlnks[i];
      decorate_all(3,vo);
      copy_uk(vo);
   }
}


void free_hyp(void)
{
   int i;
   
   for (i=0;i<6;i++)
   {
      free_lnk(dlnks[i],0);
      dlnks[i]=NULL;
   }
   free_lnk(NULL,1);
   
   if (NPROC>1)
   {
      for (i=0; i<3; i++)
         free_buf_uk(buf[i]);
      
      afree(bnd_idx);
      bnd_idx=NULL;
   }
}


static void decorated_links(double alpha1, double alpha2, double alpha3, \
      int proj_su3, int piter, int time, int space, int keep)
{
   int i,mu;
   double a[3];
   su3_dble *v,*u,*um;
   
   if (alpha1==0.0)
      return;
   
   udb=udfld();
   copy_bnd_ud();
   
   alloc_bufs();
   assign_ud2lnk(dlnks[3]);
   copy_uk(dlnks[3]);
   
   a[0]=alpha1;
   a[1]=alpha2;
   a[2]=alpha3;

   for (i=0; i<3; i++)
   {
      alpha[i]=1.0-a[i];
      beta[i]=a[i]/(double)(2*(3-i));
   }

   proj=proj_su3;
   iter=piter;
   
   decorated_links1();
   decorated_links2();
   decorated_links3(time,space);
   
   u=udb;
   um=u+4*VOLUME;
   v=dlnks[0];
   for (;u<um;)
   {
      for (mu=(1-(time>0));mu<(1+3*(space>0));mu++)
      {
         *(u+2*mu)=*(v+2*mu);
         *(u+2*mu+1)=*(v+2*mu+1);
      }
      u+=8;
      v+=8;
   }
   
   set_flags(UPDATED_UD);
   set_bc();
   
   if (keep==0)
      free_hyp();
}


void hyp_time_links(double alpha1, double alpha2, double alpha3, int proj_su3,
                    int piter, int keep)
{
   decorated_links(alpha1,alpha2,alpha3,proj_su3,piter,1,0,keep);
}


void hyp_space_links(double alpha1, double alpha2, double alpha3, int proj_su3,
                      int piter, int keep)
{
   decorated_links(alpha1,alpha2,alpha3,proj_su3,piter,0,1,keep);
}


void hyp_links(double alpha1, double alpha2, double alpha3, int proj_su3,
                int piter, int keep)
{
   decorated_links(alpha1,alpha2,alpha3,proj_su3,piter,1,1,keep);
}
