/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2009 Bjorn Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Check of hyp_links - ape_links equivalence
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
#include "global.h"
#include "flags.h"
#include "uflds.h"
#include "smear.h"
#include "lnk.h"
#include "archive.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

int main(int argc,char *argv[])
{
   int my_rank,bc,nplaq,i;
   double d,dmax,dmax_all,*pd1,*pd2,*pdm,p1,p2,alpha,a;
   double phi[2],phi_prime[2];
   su3_dble *v,*u,*um,*plnk;
   FILE *flog=NULL;   

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check3.log","w",stdout);
   
      printf("\n");
      printf("Check of hyp_links - ape_links equivalence\n");
      printf("-------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);      
      fflush(flog);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>]");
   }
   
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);

   error_root((bc!=0)&&(bc!=3),1,"main [check3.c]",
               "Only open and periodic boundary conditions supported");
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   set_bc_parms(bc,1.0,1.0,1.0,1.0,phi,phi_prime);
   print_bc_parms();
   
   start_ranlux(0,456);   
   geometry();
   plnk=alloc_lnk();
   
   alpha=0.5;
   a=(1.0-alpha);
   
   dmax=0.0;
   if (bc==0)
      nplaq=(double)((6*N0-3)*N1)*(double)(N2*N3);
   else if (bc==3)
      nplaq=(double)(6*N0*N1)*(double)(N2*N3);
   else
      nplaq=(double)((6*N0+3)*N1)*(double)(N2*N3);
   
   random_ud();
   assign_ud2lnk(plnk);
   
   hyp_links(alpha, 0.0, 0.0, 0, 0, 0);
   ape_links(plnk, alpha/6.0/a, 0, 0, 0);
   
   p1=plaq_sum_dble(1);
   p2=lnk_plaq_sum(plnk);
   
   u=udfld();
   um=u+4*VOLUME;
   v=plnk;
   
   for (;u<um;u++)
   {
      pd1=(double*)u;
      pd2=(double*)v;
      pdm=pd1+18;
      for (;pd1<pdm;pd1++)
      {
         d=fabs((*pd1)-a*(*pd2));
         
         if (d>dmax)
            dmax=d;
         
         pd2+=1;
      }
      
      v+=1;
   }

   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Equivalence of hyp_links and ape_links:\n");
      printf("Plaquette difference                  = %.1e\n",
             fabs(p1-a*a*a*a*p2)/nplaq);
      printf("Maximal deviation of matrix elements  = %.2e\n\n",dmax_all);
   }

   random_ud();
   assign_ud2lnk(plnk);
   
   hyp_spatial_links(udfld(),alpha, 0.0, 0, 0, 0);
   ape_spatial_links(plnk, alpha/4.0/a, 0, 0, 0);
   
   u=udfld();
   um=u+4*VOLUME;
   v=plnk;
   dmax=0.0;
   for (;u<um;)
   {
      for (i=0;i<2;i++)
      {
         pd1=(double*)u;
         pd2=(double*)v;
         pdm=pd1+18;
         for (;pd1<pdm;pd1++)
         {
            d=fabs((*pd1)-(*pd2));
            
            if (d>dmax)
               dmax=d;
            
            pd2+=1;
         }
         
         v+=1;
         u+=1;
      }
      
      for (i=0;i<6;i++)
      {
         pd1=(double*)u;
         pd2=(double*)v;
         pdm=pd1+18;
         for (;pd1<pdm;pd1++)
         {
            d=fabs((*pd1)-a*(*pd2));
            
            if (d>dmax)
               dmax=d;
            
            pd2+=1;
         }
         
         v+=1;
         u+=1;
      }
   }

   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Equivalence of hyp_spatial_links and ape_spatial_links:\n");
      printf("Maximal deviation of matrix elements  = %.2e\n\n",dmax_all);
   }

   MPI_Finalize();    
   exit(0);
}
