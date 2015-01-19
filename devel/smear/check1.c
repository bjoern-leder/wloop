/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2009 Bjorn Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Direct check of the HYP-smearing routines
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "lattice.h"
#include "random.h"
#include "uflds.h"
#include "global.h"
#include "flags.h"
#include "archive.h"
#include "smear.h"
#include "lnk.h"

#define NAME_SIZE 128

int main(int argc,char *argv[])
{
   int my_rank;
   double p0,p1,p2,nplaq;
   double d,dmax,dmax_all,*pd1,*pd2,*pdm;
   double phi[2],phi_prime[2];
   su3_dble *plnk,*v,*u,*um;
   FILE *log=NULL;
  
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  
   if (my_rank==0)
   {
      log=freopen("check1.log","w",stdout);
      
      printf("\n");
      printf("Direct check of the HYP-smearing routines\n");
      printf("---------------------------------------\n\n");
      
      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n\n",L0,L1,L2,L3);
   }

   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   set_bc_parms(3,1.0,1.0,1.0,1.0,phi,phi_prime);
   print_bc_parms();
   
   geometry();
   nplaq=(double)(6*NPROC)*(double)(VOLUME);

   import_cnfg("./8x8x8x8_cnfg"); 
   p0=plaq_sum_dble(1);

   hyp_links(0.75, 0.6, 0.3, 0, 4,0);
   
   plnk=alloc_lnk();
   assign_ud2lnk(plnk);
   p1=lnk_plaq_sum(plnk);
      
   import_cnfg("./8x8x8x8_hyp");
   
   p2=plaq_sum_dble(1);

   if (my_rank==0)
   {
      printf("Average plaquette:\n\n");
      printf("before    = %.8e\n",p0/nplaq);
      printf("after     = %.8e\n",p1/nplaq);
      printf("reference = %.8e\n\n",p2/nplaq);
   }
   
   u=udfld();
   v=plnk;
   um=u+4*VOLUME;
   for (;u<um;u++)
   {
      pd1=(double*)u;
      pd2=(double*)v;
      pdm=pd1+18;
      d=0.0;
      for (;pd1<pdm;pd1++)
      {
         d=fabs((*pd1)-(*pd2));
         
         if (d>dmax)
            dmax=d;
         
         pd2+=1;
      }
      
      v+=1;
   }

   MPI_Reduce(&dmax,&dmax_all,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

   if (my_rank==0)
   {
      printf("Plaquette difference                  = %.1e\n",
             fabs(p1-p2));
      printf("Maximal deviation of matrix elements  = %.2e\n",dmax_all);
      printf("(Should be of the order of %.2e and %.2e)\n\n",sqrt(nplaq)*FLT_EPSILON,FLT_EPSILON);
   }

   free_lnk(plnk,1);
   
   if(my_rank==0)
      fclose(log);

   MPI_Finalize(); 
   exit(0);
}
