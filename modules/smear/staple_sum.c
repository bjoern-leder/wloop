/*******************************************************************************
*
* File staple_sum.c
*
* Copyright (C) 2009 Bjorn Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computes the sum of up- and down-staples
*
* The externally accessible functions are
*
* void staple_sum(su3_dble *u1,su3_dble *u2,int mu,int nu,su3_dble *vb,
*                 su3_dble *cbuf)
*   Computes the sum of up- and down-staples. u1 and u2 are base addresses of
*   the global gauge field (udfld) or a lnk field  (see alloc_lnk()).
*   The links in nu direction are taken from u1, the links in mu direction
*   from u2. The result is sored in the lnk field vb and the boundary lnk
*   field vb (see alloc_buf_uk()) 
*   
* Notes:
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
* If the at compile time FSTPL is defined the staples are computed under
* the assumption that all matrices are in SU(3).
*
*******************************************************************************/

#define STAPLE_SUM_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "global.h"
#include "su3.h"
#include "lnk.h"



void staple_sum(su3_dble *u1,su3_dble *u2,int mu,int nu,su3_dble *vb,su3_dble *cbuf)
{
   int ix,iy,iz;
   su3_dble *v,*u1dn,*u2dn,*u3dn;
   
   for (ix=0;ix<VOLUME;ix++)
   {
      v=lnk(vb,ix,mu);
            
#if (defined FSTPL)
      add_up_staple(u1,u2,ix,mu,nu,v);
      add_dn_staple(u1,u2,ix,mu,nu,v);
#else
      add_up_staple_save(u1,u2,ix,mu,nu,v);
      add_dn_staple_save(u1,u2,ix,mu,nu,v);
#endif
      
      iy=iup[ix][nu];
      if (iy>=VOLUME)
      {
         v=lnk_buf(cbuf,iy,mu,nu);
         
         u1dn=lnk(u1,ix,nu);
         u2dn=lnk(u2,ix,mu);
         iz=iup[ix][mu];
         u3dn=lnkf(u1,iz,nu,mu);
         
#if (defined FSTPL)
         add_dn_staple_lnks(u1dn,u2dn,u3dn,v);
#else
         add_dn_staple_lnks_save(u1dn,u2dn,u3dn,v);
#endif
      }
   }
}
