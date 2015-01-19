/*******************************************************************************
*
* File lnkstaple.c
*
* Copyright (C) 2009 Bjorn Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of staples
*
* The externally accessible functions are
*
*   void add_up_staple(su3_dble *u1o,su3_dble *u2o,int ix,int mu,int nu,su3_dble *v)
*     Computes the up-staple where u1 and u3 are taken from the link-field with
*     base address u1o and u2 is taken from u2o. The link matrices are expected
*     to be in SU(3). The result is added to v. For a definiton of u1-3 see the
*     notes below.
*
*   void add_dn_staple(su3_dble *u1o,su3_dble *u2o,int ix,int mu,int nu,su3_dble *v)
*     Same as add_up_staple, but computes the down-staple.
*
*   void add_up_staple_lnks(su3_dble *u1,su3_dble *u2,su3_dble *u3,su3_dble *v)
*     Same as add_up_staple, but uses the links in the argument.
*
*   void add_dn_staple_lnks(su3_dble *u1,su3_dble *u2,su3_dble *u3,su3_dble *v)
*     Same as add_dn_staple, but uses the links in the argument.
*
*   void add_up_staple_save(su3_dble *u1o,su3_dble *u2o,int ix,int mu,int nu,su3_dble *v)
*     Same as add_up_staple, but the links can be general 3x3 matrices.
*
*   void add_dn_staple_save(su3_dble *u1o,su3_dble *u2o,int ix,int mu,int nu,su3_dble *v)
*     Same as add_dn_staple, but the links can be general 3x3 matrices.
*
*   void add_up_staple_lnks_save(su3_dble *u1,su3_dble *u2,su3_dble *u3,su3_dble *v)
*     Same as add_up_staple_lnks, but the links can be general 3x3 matrices.
*
*   void add_dn_staple_lnks_save(su3_dble *u1,su3_dble *u2,su3_dble *u3,su3_dble *v)
*     Same as add_dn_staple_lnks, but the links can be general 3x3 matrices.
*
* Notes:
*
* The staple is a 3x3 matrix formed by the product of 3 links that form a path
* around the link U(x,mu):
*
* up-staple: v=u1*u2*u3^dag
*
*   u1 -> U(x,nu)
*   u2 -> U(x+nu,mu)
*   u3 -> U(x+mu,nu)
*
* down-staple: v=u1^dag*u2*u3
*
*   u1 -> U(x-nu,nu)
*   u2 -> U(x-nu,mu)
*   u3 -> U(x-nu+mu,nu)
*
* In the programs it is assumed that 0<=ix<VOLUME and mu!=nu.
*
* In the programs add_lnk_staple and add_lnk_staple_save, if ix-nu is on the
* boundary, nothing is added to v.
*
*******************************************************************************/

#define LNKSTAPLE_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "su3fcts.h"
#include "global.h"
#include "lnk.h"

static su3_vector_dble psi1,psi2,chi1,chi2;
static su3_dble *u1up,*u2up,*u3up;
static su3_dble *u1dn,*u2dn,*u3dn;
static su3_dble *u1,*u2,*u3;

#if (defined SSE2)
static su3_dble tmp __attribute__ ((aligned (16)));
#else
static su3_dble tmp;
#endif

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


static void add_to_v(su3_dble *v)
{
   (*v).c11.re+=psi1.c1.re;
   (*v).c11.im+=psi1.c1.im;
   (*v).c21.re+=psi1.c2.re;
   (*v).c21.im+=psi1.c2.im;
   (*v).c31.re+=psi1.c3.re;
   (*v).c31.im+=psi1.c3.im;

   (*v).c12.re+=psi2.c1.re;
   (*v).c12.im+=psi2.c1.im;
   (*v).c22.re+=psi2.c2.re;
   (*v).c22.im+=psi2.c2.im;
   (*v).c32.re+=psi2.c3.re;
   (*v).c32.im+=psi2.c3.im;

   (*v).c13.re+=
         (psi1.c2.re*psi2.c3.re-psi1.c2.im*psi2.c3.im)-
         (psi2.c2.re*psi1.c3.re-psi2.c2.im*psi1.c3.im);
   (*v).c13.im+=
         (psi2.c2.re*psi1.c3.im+psi2.c2.im*psi1.c3.re)-
         (psi1.c2.re*psi2.c3.im+psi1.c2.im*psi2.c3.re);
   (*v).c23.re+=
         (psi1.c3.re*psi2.c1.re-psi1.c3.im*psi2.c1.im)-
         (psi2.c3.re*psi1.c1.re-psi2.c3.im*psi1.c1.im);
   (*v).c23.im+=
         (psi2.c3.re*psi1.c1.im+psi2.c3.im*psi1.c1.re)-
         (psi1.c3.re*psi2.c1.im+psi1.c3.im*psi2.c1.re);
   (*v).c33.re+=
         (psi1.c1.re*psi2.c2.re-psi1.c1.im*psi2.c2.im)-
         (psi2.c1.re*psi1.c2.re-psi2.c1.im*psi1.c2.im);
   (*v).c33.im+=
         (psi2.c1.re*psi1.c2.im+psi2.c1.im*psi1.c2.re)-
         (psi1.c1.re*psi2.c2.im+psi1.c1.im*psi2.c2.re);
}


static void up_staple(void)
{
   chi1.c1.re=
         (*u2up).c11.re*(*u3up).c11.re+(*u2up).c11.im*(*u3up).c11.im+
         (*u2up).c12.re*(*u3up).c12.re+(*u2up).c12.im*(*u3up).c12.im+ 
         (*u2up).c13.re*(*u3up).c13.re+(*u2up).c13.im*(*u3up).c13.im;
   chi1.c1.im=
         (*u2up).c11.im*(*u3up).c11.re-(*u2up).c11.re*(*u3up).c11.im+
         (*u2up).c12.im*(*u3up).c12.re-(*u2up).c12.re*(*u3up).c12.im+ 
         (*u2up).c13.im*(*u3up).c13.re-(*u2up).c13.re*(*u3up).c13.im;
   chi1.c2.re=
         (*u2up).c21.re*(*u3up).c11.re+(*u2up).c21.im*(*u3up).c11.im+ 
         (*u2up).c22.re*(*u3up).c12.re+(*u2up).c22.im*(*u3up).c12.im+
         (*u2up).c23.re*(*u3up).c13.re+(*u2up).c23.im*(*u3up).c13.im;
   chi1.c2.im=
         (*u2up).c21.im*(*u3up).c11.re-(*u2up).c21.re*(*u3up).c11.im+
         (*u2up).c22.im*(*u3up).c12.re-(*u2up).c22.re*(*u3up).c12.im+
         (*u2up).c23.im*(*u3up).c13.re-(*u2up).c23.re*(*u3up).c13.im;
   chi1.c3.re=
         (*u2up).c31.re*(*u3up).c11.re+(*u2up).c31.im*(*u3up).c11.im+
         (*u2up).c32.re*(*u3up).c12.re+(*u2up).c32.im*(*u3up).c12.im+ 
         (*u2up).c33.re*(*u3up).c13.re+(*u2up).c33.im*(*u3up).c13.im;
   chi1.c3.im=
         (*u2up).c31.im*(*u3up).c11.re-(*u2up).c31.re*(*u3up).c11.im+ 
         (*u2up).c32.im*(*u3up).c12.re-(*u2up).c32.re*(*u3up).c12.im+
         (*u2up).c33.im*(*u3up).c13.re-(*u2up).c33.re*(*u3up).c13.im;

   chi2.c1.re=
         (*u2up).c11.re*(*u3up).c21.re+(*u2up).c11.im*(*u3up).c21.im+
         (*u2up).c12.re*(*u3up).c22.re+(*u2up).c12.im*(*u3up).c22.im+ 
         (*u2up).c13.re*(*u3up).c23.re+(*u2up).c13.im*(*u3up).c23.im;
   chi2.c1.im=
         (*u2up).c11.im*(*u3up).c21.re-(*u2up).c11.re*(*u3up).c21.im+
         (*u2up).c12.im*(*u3up).c22.re-(*u2up).c12.re*(*u3up).c22.im+ 
         (*u2up).c13.im*(*u3up).c23.re-(*u2up).c13.re*(*u3up).c23.im;
   chi2.c2.re=
         (*u2up).c21.re*(*u3up).c21.re+(*u2up).c21.im*(*u3up).c21.im+ 
         (*u2up).c22.re*(*u3up).c22.re+(*u2up).c22.im*(*u3up).c22.im+
         (*u2up).c23.re*(*u3up).c23.re+(*u2up).c23.im*(*u3up).c23.im;
   chi2.c2.im=
         (*u2up).c21.im*(*u3up).c21.re-(*u2up).c21.re*(*u3up).c21.im+
         (*u2up).c22.im*(*u3up).c22.re-(*u2up).c22.re*(*u3up).c22.im+
         (*u2up).c23.im*(*u3up).c23.re-(*u2up).c23.re*(*u3up).c23.im;
   chi2.c3.re=
         (*u2up).c31.re*(*u3up).c21.re+(*u2up).c31.im*(*u3up).c21.im+
         (*u2up).c32.re*(*u3up).c22.re+(*u2up).c32.im*(*u3up).c22.im+ 
         (*u2up).c33.re*(*u3up).c23.re+(*u2up).c33.im*(*u3up).c23.im;
   chi2.c3.im=
         (*u2up).c31.im*(*u3up).c21.re-(*u2up).c31.re*(*u3up).c21.im+ 
         (*u2up).c32.im*(*u3up).c22.re-(*u2up).c32.re*(*u3up).c22.im+
         (*u2up).c33.im*(*u3up).c23.re-(*u2up).c33.re*(*u3up).c23.im;

   _su3_multiply(psi1,(*u1up),chi1);
   _su3_multiply(psi2,(*u1up),chi2);   
}


static void dn_staple(void)
{
   chi1.c1.re=
         (*u2dn).c11.re*(*u3dn).c11.re-(*u2dn).c11.im*(*u3dn).c11.im+
         (*u2dn).c12.re*(*u3dn).c21.re-(*u2dn).c12.im*(*u3dn).c21.im+ 
         (*u2dn).c13.re*(*u3dn).c31.re-(*u2dn).c13.im*(*u3dn).c31.im;
   chi1.c1.im=
         (*u2dn).c11.im*(*u3dn).c11.re+(*u2dn).c11.re*(*u3dn).c11.im+
         (*u2dn).c12.im*(*u3dn).c21.re+(*u2dn).c12.re*(*u3dn).c21.im+ 
         (*u2dn).c13.im*(*u3dn).c31.re+(*u2dn).c13.re*(*u3dn).c31.im;
   chi1.c2.re=
         (*u2dn).c21.re*(*u3dn).c11.re-(*u2dn).c21.im*(*u3dn).c11.im+ 
         (*u2dn).c22.re*(*u3dn).c21.re-(*u2dn).c22.im*(*u3dn).c21.im+
         (*u2dn).c23.re*(*u3dn).c31.re-(*u2dn).c23.im*(*u3dn).c31.im;
   chi1.c2.im=
         (*u2dn).c21.im*(*u3dn).c11.re+(*u2dn).c21.re*(*u3dn).c11.im+
         (*u2dn).c22.im*(*u3dn).c21.re+(*u2dn).c22.re*(*u3dn).c21.im+
         (*u2dn).c23.im*(*u3dn).c31.re+(*u2dn).c23.re*(*u3dn).c31.im;
   chi1.c3.re=
         (*u2dn).c31.re*(*u3dn).c11.re-(*u2dn).c31.im*(*u3dn).c11.im+
         (*u2dn).c32.re*(*u3dn).c21.re-(*u2dn).c32.im*(*u3dn).c21.im+ 
         (*u2dn).c33.re*(*u3dn).c31.re-(*u2dn).c33.im*(*u3dn).c31.im;
   chi1.c3.im=
         (*u2dn).c31.im*(*u3dn).c11.re+(*u2dn).c31.re*(*u3dn).c11.im+ 
         (*u2dn).c32.im*(*u3dn).c21.re+(*u2dn).c32.re*(*u3dn).c21.im+
         (*u2dn).c33.im*(*u3dn).c31.re+(*u2dn).c33.re*(*u3dn).c31.im;

   chi2.c1.re=
         (*u2dn).c11.re*(*u3dn).c12.re-(*u2dn).c11.im*(*u3dn).c12.im+
         (*u2dn).c12.re*(*u3dn).c22.re-(*u2dn).c12.im*(*u3dn).c22.im+ 
         (*u2dn).c13.re*(*u3dn).c32.re-(*u2dn).c13.im*(*u3dn).c32.im;
   chi2.c1.im=
         (*u2dn).c11.im*(*u3dn).c12.re+(*u2dn).c11.re*(*u3dn).c12.im+
         (*u2dn).c12.im*(*u3dn).c22.re+(*u2dn).c12.re*(*u3dn).c22.im+ 
         (*u2dn).c13.im*(*u3dn).c32.re+(*u2dn).c13.re*(*u3dn).c32.im;
   chi2.c2.re=
         (*u2dn).c21.re*(*u3dn).c12.re-(*u2dn).c21.im*(*u3dn).c12.im+ 
         (*u2dn).c22.re*(*u3dn).c22.re-(*u2dn).c22.im*(*u3dn).c22.im+
         (*u2dn).c23.re*(*u3dn).c32.re-(*u2dn).c23.im*(*u3dn).c32.im;
   chi2.c2.im=
         (*u2dn).c21.im*(*u3dn).c12.re+(*u2dn).c21.re*(*u3dn).c12.im+
         (*u2dn).c22.im*(*u3dn).c22.re+(*u2dn).c22.re*(*u3dn).c22.im+
         (*u2dn).c23.im*(*u3dn).c32.re+(*u2dn).c23.re*(*u3dn).c32.im;
   chi2.c3.re=
         (*u2dn).c31.re*(*u3dn).c12.re-(*u2dn).c31.im*(*u3dn).c12.im+
         (*u2dn).c32.re*(*u3dn).c22.re-(*u2dn).c32.im*(*u3dn).c22.im+ 
         (*u2dn).c33.re*(*u3dn).c32.re-(*u2dn).c33.im*(*u3dn).c32.im;
   chi2.c3.im=
         (*u2dn).c31.im*(*u3dn).c12.re+(*u2dn).c31.re*(*u3dn).c12.im+ 
         (*u2dn).c32.im*(*u3dn).c22.re+(*u2dn).c32.re*(*u3dn).c22.im+
         (*u2dn).c33.im*(*u3dn).c32.re+(*u2dn).c33.re*(*u3dn).c32.im;

   _su3_inverse_multiply(psi1,(*u1dn),chi1);
   _su3_inverse_multiply(psi2,(*u1dn),chi2);      
}


void add_up_staple(su3_dble *u1o,su3_dble *u2o,int ix,int mu,int nu,su3_dble *v)
{
   int iy;
   
   u1up=lnk(u1o,ix,nu);
   
   iy=iup[ix][nu];
   u2up=lnkf(u2o,iy,mu,nu);
   
   iy=iup[ix][mu];
   u3up=lnkf(u1o,iy,nu,mu);

   up_staple();
   add_to_v(v);
}


void add_dn_staple(su3_dble *u1o,su3_dble *u2o,int ix,int mu,int nu,su3_dble *v)
{
   int iy;

   ix=idn[ix][nu];
   if (ix<VOLUME)
   {
      u1dn=lnk(u1o,ix,nu);
      
      u2dn=lnk(u2o,ix,mu);
      
      iy=iup[ix][mu];
      u3dn=lnkf(u1o,iy,nu,mu);
   
      dn_staple();
      add_to_v(v);
   }
}


void add_up_staple_lnks(su3_dble *v1,su3_dble *v2,su3_dble *v3,su3_dble *v)
{
   u1up=v1;
   u2up=v2;
   u3up=v3;

   up_staple();
   add_to_v(v);
}


void add_dn_staple_lnks(su3_dble *v1,su3_dble *v2,su3_dble *v3,su3_dble *v)
{
   u1dn=v1;
   u2dn=v2;
   u3dn=v3;
   
   dn_staple();
   add_to_v(v);
}


static void add_up_to_v_save(su3_dble *v)
{
#if (defined SSE2)
   _prefetch_su3_dble(u1);
#endif      
   su3xsu3dag(u2,u3,&tmp);
   su3xsu3(u1,&tmp,&tmp);
   _su3_add(*v,*v,tmp);
}


static void add_dn_to_v_save(su3_dble *v)
{
#if (defined SSE2)
   _prefetch_su3_dble(u1);
#endif      
   su3xsu3(u2,u3,&tmp);
   su3dagxsu3(u1,&tmp,&tmp);
   _su3_add(*v,*v,tmp);
}


void add_up_staple_save(su3_dble *u1o,su3_dble *u2o,int ix,int mu,int nu,su3_dble *v)
{
   int iy;
   
   iy=iup[ix][nu];
   u2=lnkf(u2o,iy,mu,nu);
   
   iy=iup[ix][mu];
   u3=lnkf(u1o,iy,nu,mu);

#if (defined SSE2)
   _prefetch_su3_dble(u2);
   _prefetch_su3_dble(u3);
#endif      

   u1=lnk(u1o,ix,nu);
   
   add_up_to_v_save(v);
}


void add_dn_staple_save(su3_dble *u1o,su3_dble *u2o,int ix,int mu,int nu,su3_dble *v)
{
   int iy;

   ix=idn[ix][nu];
   if (ix<VOLUME)
   {
      u2=lnk(u2o,ix,mu);
      
      iy=iup[ix][mu];
      u3=lnkf(u1o,iy,nu,mu);

#if (defined SSE2)
      _prefetch_su3_dble(u2);
      _prefetch_su3_dble(u3);
#endif      
   
      u1=lnk(u1o,ix,nu);
      
      add_dn_to_v_save(v);
   }
}


void add_up_staple_lnks_save(su3_dble *v1,su3_dble *v2,su3_dble *v3,su3_dble *v)
{
   u1=v1;
   u2=v2;
   u3=v3;

#if (defined SSE2)
   _prefetch_su3_dble(u2);
   _prefetch_su3_dble(u3);
#endif      
   add_up_to_v_save(v);
}


void add_dn_staple_lnks_save(su3_dble *v1,su3_dble *v2,su3_dble *v3,su3_dble *v)
{
   u1=v1;
   u2=v2;
   u3=v3;
   
#if (defined SSE2)
   _prefetch_su3_dble(u2);
   _prefetch_su3_dble(u3);
#endif      
   add_dn_to_v_save(v);
}
