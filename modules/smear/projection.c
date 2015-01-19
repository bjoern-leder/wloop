/*******************************************************************************
*
* File projection.c
*
* Copyright (C) 2009 Bjorn Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Projection to SU(3)
*
* The externally accessible functions are
*
*   void approx_project_to_su3_dble(su3_dble *u, int iter)
*     The approximate projection to SU(3) implemented here is defined in
*     hep-lat/0506008.
*
*******************************************************************************/

#define PROJECTION_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "su3.h"
#include "lattice.h"
#include "su3fcts.h"
#include "global.h"
#include "uflds.h"

#define SRT18 4.242640687119285


#if (defined SSE2)
static su3_dble x __attribute__ ((aligned (16)));
#else
static su3_dble x;
#endif


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

#define _su3_mulc(z,c,u) \
  (z).c11.re=(c).re*(u).c11.re-(c).im*(u).c11.im;\
  (z).c11.im=(c).re*(u).c11.im+(c).im*(u).c11.re;\
  (z).c12.re=(c).re*(u).c12.re-(c).im*(u).c12.im;\
  (z).c12.im=(c).re*(u).c12.im+(c).im*(u).c12.re;\
  (z).c13.re=(c).re*(u).c13.re-(c).im*(u).c13.im;\
  (z).c13.im=(c).re*(u).c13.im+(c).im*(u).c13.re;\
  (z).c21.re=(c).re*(u).c21.re-(c).im*(u).c21.im;\
  (z).c21.im=(c).re*(u).c21.im+(c).im*(u).c21.re;\
  (z).c22.re=(c).re*(u).c22.re-(c).im*(u).c22.im;\
  (z).c22.im=(c).re*(u).c22.im+(c).im*(u).c22.re;\
  (z).c23.re=(c).re*(u).c23.re-(c).im*(u).c23.im;\
  (z).c23.im=(c).re*(u).c23.im+(c).im*(u).c23.re;\
  (z).c31.re=(c).re*(u).c31.re-(c).im*(u).c31.im;\
  (z).c31.im=(c).re*(u).c31.im+(c).im*(u).c31.re;\
  (z).c32.re=(c).re*(u).c32.re-(c).im*(u).c32.im;\
  (z).c32.im=(c).re*(u).c32.im+(c).im*(u).c32.re;\
  (z).c33.re=(c).re*(u).c33.re-(c).im*(u).c33.im;\
  (z).c33.im=(c).re*(u).c33.im+(c).im*(u).c33.re

#define _su3_addr(z,a) \
  (z).c11.re+=a;\
  (z).c22.re+=a;\
  (z).c33.re+=a


static double det_im_dble(su3_dble *u)
{
   double detuim;
   complex_dble det1,det2,det3;

   det1.re=
         ((*u).c22.re*(*u).c33.re-(*u).c22.im*(*u).c33.im)-
         ((*u).c23.re*(*u).c32.re-(*u).c23.im*(*u).c32.im);
   det1.im=
         ((*u).c22.re*(*u).c33.im+(*u).c22.im*(*u).c33.re)-
         ((*u).c23.re*(*u).c32.im+(*u).c23.im*(*u).c32.re);
   det2.re=
         ((*u).c21.re*(*u).c33.re-(*u).c21.im*(*u).c33.im)-
         ((*u).c23.re*(*u).c31.re-(*u).c23.im*(*u).c31.im);
   det2.im=
         ((*u).c21.re*(*u).c33.im+(*u).c21.im*(*u).c33.re)-
         ((*u).c23.re*(*u).c31.im+(*u).c23.im*(*u).c31.re);
   det3.re=
         ((*u).c21.re*(*u).c32.re-(*u).c21.im*(*u).c32.im)-
         ((*u).c22.re*(*u).c31.re-(*u).c22.im*(*u).c31.im);
   det3.im=
         ((*u).c21.re*(*u).c32.im+(*u).c21.im*(*u).c32.re)-
         ((*u).c22.re*(*u).c31.im+(*u).c22.im*(*u).c31.re);

   detuim=
         ((*u).c11.re*det1.im+(*u).c11.im*det1.re)-
         ((*u).c12.re*det2.im+(*u).c12.im*det2.re)+
         ((*u).c13.re*det3.im+(*u).c13.im*det3.re);

   return detuim;
}


void approx_project_to_su3_dble(su3_dble *u, int iter)
{
   int i;
   double d;
   complex_dble z;
   
   d= (*u).c11.re*(*u).c11.re+(*u).c11.im*(*u).c11.im  \
         +(*u).c12.re*(*u).c12.re+(*u).c12.im*(*u).c12.im \
         +(*u).c13.re*(*u).c13.re+(*u).c13.im*(*u).c13.im \
         +(*u).c21.re*(*u).c21.re+(*u).c21.im*(*u).c21.im \
         +(*u).c22.re*(*u).c22.re+(*u).c22.im*(*u).c22.im \
         +(*u).c23.re*(*u).c23.re+(*u).c23.im*(*u).c23.im \
         +(*u).c31.re*(*u).c31.re+(*u).c31.im*(*u).c31.im \
         +(*u).c32.re*(*u).c32.re+(*u).c32.im*(*u).c32.im \
         +(*u).c33.re*(*u).c33.re+(*u).c33.im*(*u).c33.im;
   
   if(d<SRT18*DBL_EPSILON)
      return;

   d=1.0/sqrt(d/3.0);
   
   _su3_mulr(*u,d,*u);
   
   z.re=1.0;
   for (i=0;i<iter;i++){
      su3dagxsu3(u,u,&x);
      _su3_mulr(x,-0.5,x);
      _su3_addr(x,1.5);
      su3xsu3(u,&x,&x);

      z.im=det_im_dble(&x);
      z.im/=-3.0;
      _su3_mulc(*u,z,x);
   }
}


