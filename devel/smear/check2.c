/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2009 Bjorn Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of the hyp_links (SU(3) normalization)
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "lattice.h"
#include "uflds.h"
#include "global.h"
#include "flags.h"
#include "smear.h"
#include "lnk.h"

static void dev(su3_dble *u,double *d1,double *d2)
{
   int i;
   double d;
   complex_dble det1,det2,det3,det;
   su3_dble v,w;

   _su3_dagger(v,(*u));
   _su3_times_su3(w,v,(*u));

   w.c11.re-=1.0;
   w.c22.re-=1.0;
   w.c33.re-=1.0;

   *d1=0.0;

   for (i=0;i<18;i++)
   {
      d=*((double*)(&w)+i);
      d=fabs(d);
      if (d>(*d1))
         *d1=d;
   }

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

   det.re=
         ((*u).c11.re*det1.re-(*u).c11.im*det1.im)-
         ((*u).c12.re*det2.re-(*u).c12.im*det2.im)+
         ((*u).c13.re*det3.re-(*u).c13.im*det3.im);
   det.im=
         ((*u).c11.re*det1.im+(*u).c11.im*det1.re)-
         ((*u).c12.re*det2.im+(*u).c12.im*det2.re)+
         ((*u).c13.re*det3.im+(*u).c13.im*det3.re);

   *d2=0.0;
   d=fabs(det.re-1.0);
   if (d>(*d2))
      *d2=d;
   d=fabs(det.im);
   if (d>(*d2))
      *d2=d;
}


static void max_dev(double *dmax1, double *dmax2)
{
   int bc,ix,mu,t;
   double d1,d2;
   su3_dble *u,*v;

   bc=bc_type();
   
   u=udfld();   
   
   *dmax1=0.0;
   *dmax2=0.0;
   
   for (ix=0;ix<VOLUME/2;ix++)
   {
      v=u+8*ix;
      t=global_time(ix);
      
      if ((t<(NPROC0*L0-1))||(bc>0))
      {
         dev(v,&d1,&d2);

         if (d1>*dmax1)
            *dmax1=d1;
         if (d2>*dmax2)
            *dmax2=d2; 
      }
      
      if ((t>0)||(bc>0))
      {
         dev(v+1,&d1,&d2);

         if (d1>*dmax1)
            *dmax1=d1;
         if (d2>*dmax2)
            *dmax2=d2; 
      }

      for (mu=2;mu<7;mu++)
      {
         dev(v+mu,&d1,&d2);

         if (d1>*dmax1)
            *dmax1=d1;
         if (d2>*dmax2)
            *dmax2=d2; 
      }
   }
}


int main(int argc,char *argv[])
{
   int my_rank,piter,bc;
   double a1,a2,a3,dmax1,dmax2,dmax1_all,dmax2_all;
   double phi[2],phi_prime[2];
   su3_dble *plnk;
   FILE *flog=NULL;   

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check2.log","w",stdout);
   
      printf("\n");
      printf("Check of the hyp_links (SU(3) normalization)\n");
      printf("-------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);      
      fflush(flog);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check2.c]",
                    "Syntax: check2 [-bc <type>]");
   }
   
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);

   error_root((bc!=0)&&(bc!=3),1,"main [check2.c]",
               "Only open and periodic boundary conditions supported");
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   set_bc_parms(bc,1.0,1.0,1.0,1.0,phi,phi_prime);
   print_bc_parms();

   start_ranlux(0,123456);   
   geometry();

   random_ud();   
   plnk=alloc_lnk();
   assign_ud2lnk(plnk);
   
   a1=0.75;
   a2=0.6;
   a3=0.3;
   piter=30;
   
   hyp_links(a1, a2, a3, 1, 0, 0);
   
   max_dev(&dmax1,&dmax2);  
   MPI_Reduce(&dmax1,&dmax1_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
   MPI_Reduce(&dmax2,&dmax2_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("After hyp_links w/o projection (parameter: %.2f, %.2f, %.2f, iter: %d): \n",
             a1,a2,a3,0); 
      printf("Maximal deviation from unitarity = %.2e\n",dmax1_all);
      printf("Maximal deviation from det{U}=1  = %.2e\n\n",dmax2_all);
   }

   assign_lnk2ud(plnk);

   hyp_links(a1, a2, a3, 1, piter, 0);
   
   max_dev(&dmax1,&dmax2);  
   MPI_Reduce(&dmax1,&dmax1_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
   MPI_Reduce(&dmax2,&dmax2_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("After hyp_links with projection (parameter: %.2f, %.2f, %.2f, iter: %d): \n",
             a1,a2,a3,piter); 
      printf("Maximal deviation from unitarity = %.2e\n",dmax1_all);
      printf("Maximal deviation from det{U}=1  = %.2e\n\n",dmax2_all);
   }

   MPI_Finalize();    
   exit(0);
}
