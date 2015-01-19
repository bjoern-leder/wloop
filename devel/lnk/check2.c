/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2009 Bjorn Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Staple routines provided by lnkstaple
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "flags.h"
#include "global.h"
#include "lattice.h"
#include "uflds.h"
#include "lnk.h"

static const su3_dble v0={{0.0}};
static su3_vector_dble psi1,psi2,chi1,chi2;
static su3_dble *u1up,*u2up,*u3up;
static su3_dble *u1dn,*u2dn,*u3dn;

typedef struct
{
   int ix,mu;
} link_t;

static int nlks;
static link_t *lks=NULL;



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
   

void staples(int ix,int mu,su3_dble *v)
{
   int i,nu,ixpmu,ixpnu,ixmnu,ixpmumnu;
   su3_dble *udb;

   udb=udfld();
   *v=v0;
   ixpmu=iup[ix][mu];

   for (i=1;i<4;i++)
   {
      nu=(mu+i)&0x3;
      ixpnu=iup[ix][nu];
      ixmnu=idn[ix][nu];
      ixpmumnu=idn[ixpmu][nu];
      
      u1up=lnk(udb,ix,nu);
      u2up=lnk(udb,ixpnu,mu);
      u3up=lnk(udb,ixpmu,nu);
   
      u1dn=lnk(udb,ixmnu,nu);
      u2dn=lnk(udb,ixmnu,mu);
      u3dn=lnk(udb,ixpmumnu,nu);
   
      up_staple();
      add_to_v(v);      
      dn_staple();
      add_to_v(v);
   }
}


static void my_staples(su3_dble *u,int ix,int mu,su3_dble *v)
{
   int i,nu;

   *v=v0;

   for (i=1;i<4;i++)
   {
      nu=(mu+i)&0x3;
      
      add_up_staple(u,u,ix,mu,nu,v);
      add_dn_staple(u,u,ix,mu,nu,v);
   }
}


static void my_staples_save(su3_dble *u,int ix,int mu,su3_dble *v)
{
   int i,nu;

   *v=v0;

   for (i=1;i<4;i++)
   {
      nu=(mu+i)&0x3;
      
      add_up_staple_save(u,u,ix,mu,nu,v);
      add_dn_staple_save(u,u,ix,mu,nu,v);
   }
}

static void active_links(void)
{
   int bs[4],ix,iy,mu,nu,il;
   int itx,ity;

   bs[0]=L0-2*(1-(NPROC0%2));
   bs[1]=L1-2*(1-(NPROC1%2));
   bs[2]=L2-2*(1-(NPROC2%2));
   bs[3]=L3-2*(1-(NPROC3%2));

   nlks=4*bs[0]*bs[1]*bs[2]*bs[3];

   if (NPROC0>1)
      nlks+=bs[1]*bs[2]*bs[3];
   if (NPROC1>1)
      nlks+=bs[2]*bs[3]*bs[0];
   if (NPROC2>1)
      nlks+=bs[3]*bs[0]*bs[1];
   if (NPROC3>1)
      nlks+=bs[0]*bs[1]*bs[2];

   lks=amalloc(nlks*sizeof(link_t),3);
   error(lks==NULL,1,"active_links [check1.c]",
         "Unable to allocate index array");

   il=0;

   for (ix=0;ix<VOLUME;ix++)
   {
      for (mu=0;mu<4;mu++)
      {
         iy=iup[ix][mu];

         if (iy<VOLUME)
         {
            itx=0;
            ity=0;

            for (nu=0;nu<4;nu++)
            {
               if ((iup[ix][nu]>=VOLUME)||(idn[ix][nu]>=VOLUME))
                  itx=1;
               if ((iup[iy][nu]>=VOLUME)||(idn[iy][nu]>=VOLUME))
                  ity=1;
            }

            if ((itx==0)||(ity==0))
            {
               lks[il].ix=ix;
               lks[il].mu=mu;
               il+=1;
            }
         }
      }
   }
}


int main(int argc,char *argv[])
{
   int my_rank,ix,mu,n;
   double d1,d1max,d2,d2max,dmax_all;
   double *pd,*pds,*pdm,*pf;
   su3_dble wd,vd,vds,*plnk;
   FILE *flog=NULL;   

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout);

      printf("\n");
      printf("Staple routines provided by lnkstaple\n");
      printf("----------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
      fflush(flog);
   }

   start_ranlux(0,12345);
   geometry();

   plnk=alloc_lnk();
   message("Fields allocated\n\n");
 
   random_ud();
   assign_ud2lnk(plnk);
   assign_ud2u();
   
   d1max=0.0;
   d2max=0.0;

   active_links();
   for (n=0;n<nlks;n++)
   {
      ix=lks[n].ix;
      mu=lks[n].mu;
      
      staples(ix,mu,&wd);
      
      my_staples(plnk,ix,mu,&vd);
      
      pf=(double*)&wd;
      pd=(double*)&vd;
      pdm=pd+18;
      d1=0.0;
      for (;pd<pdm;pd++)
      {
         d1+=fabs((*pf)-*pd);
         pf+=1;
      }
            
      if (d1>d1max)
         d1max=d1;
         
      my_staples_save(plnk,ix,mu,&vds);
      
      pds=(double*)&vds;
      pd=(double*)&vd;
      d2=0.0;
      for (;pd<pdm;pd++)
      {
         d2+=fabs(*pds-*pd);
         pds+=1;
      }
            
      if (d2>d2max)
         d2max=d2;
   }

   MPI_Reduce(&d1max,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
 
   if (my_rank==0)
   {   
      printf("Fast version (assuming all matrices are SU(3))\n");
      printf("Maximal absolute deviation between staples and lnkstaple: %.1e\n",dmax_all);
      printf("(Should be around 1.0e-6)\n\n");
   }   

   MPI_Reduce(&d2max,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
 
   if (my_rank==0)
   {   
      printf("Save version (no assumptions)\n");
      printf("Maximal absolute deviation between lnkstaple and lnkstaple (save): %.1e\n",dmax_all);
      printf("(Should be around 1.0e-14)\n\n");
      fclose(flog);
   }   

   MPI_Finalize();    
   exit(0);
}
