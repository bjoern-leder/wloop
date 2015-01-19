/*******************************************************************************
*
* File smear.c
*
* Copyright (C) 2009 Bjorn Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Smearing of spatial links.
* 
* The externally accessible functions are
*
* void smear(void) 
*   Applies smearing operators to the spatial links. The number of applications
*  (smearing level) is obtainted form the wloop_parms data base.
* 
* su3_dble** get_pus(void)
*   Returns array of base addresses for smeared spatial links.
* 
* void free_smear(void)
*   Free auxillary memory.
* 
*******************************************************************************/
#define SMEAR_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "lattice.h"
#include "global.h"
#include "wloop.h"
#include "uflds.h"
#include "flags.h"
#include "smear.h"
#include "lnk.h"

static su3_dble *usb=NULL,**pus=NULL;
static su3_dble *plnk;
static wloop_parms_t wl;


static void alloc_us(void)
{
   int sl,msl,x0,ixs,k;
   su3_dble *p;

   if(usb==NULL)
   {
      plnk=alloc_lnk();
      
      usb=amalloc(3*(2*wl.msl)*VOLUME*sizeof(su3_dble),ALIGN);
      pus=amalloc(3*(2*wl.msl)*VOLUME*sizeof(su3_dble*),ALIGN);
      error_root((usb==NULL)||(pus==NULL),1,"alloc_us [smear.c]",
               "Could not allocate memory space for the smear field");

      p=usb;
      msl=2*wl.msl;
      for(sl=0;sl<msl;++sl)
      {
         for(x0=0;x0<L0;++x0)
         {
            for(ixs=0;ixs<SLICE0;++ixs)
            {
               for(k=0;k<3;++k)
               {
                  pus[pindxs(sl,ixs,x0,k)]=p;
                  ++p;
               }
            }
         }
      }
   }
}


static void copy_lnk2us(int sl)
{
   int ix,k;
   su3_dble *v,*u;

   copy_u0(plnk);
   
   for(ix=0;ix<VOLUME;++ix)
   {
      for(k=1;k<4;++k)
      {
         v=lnk(plnk,ix,k);
         u=pus[pindx(sl,ix,k-1)];
         
         *u=*v;
      }
   }
}


void smear(void)
{
   int i,sl,nstep;

   wl=wloop_parms();
      
   alloc_us();
   
   assign_ud2lnk(plnk);
   
   for(sl=0;sl<wl.msl;++sl)
   {
      nstep=wl.nss[sl];
      for(i=0;i<nstep;++i)
      {
         if (wl.op_smear==HYP)
            hyp_spatial_links(plnk,wl.alpha_hyp[0],wl.alpha_hyp[1],wl.proj,wl.proj_iter,1);
         else if (wl.op_smear==APE)
            ape_spatial_links(plnk,wl.alpha_ape,wl.proj,wl.proj_iter,1);
      }
      copy_lnk2us(sl);
   }
}


su3_dble** get_pus(void)
{
   return pus;
}


void free_smear(void)
{
   if (usb!=NULL)
   {
      afree(usb);
      usb=NULL;
      afree(pus);
      pus=NULL;
      free_lnk(plnk,1);
      free_hyp_spatial();
   }
}
