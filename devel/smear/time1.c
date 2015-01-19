/*******************************************************************************
*
* File time1.c
*
* Copyright (C) 2009 Bjorn Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Timing of the smearing routines
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "random.h"
#include "uflds.h"
#include "global.h"
#include "lattice.h"
#include "flags.h"
#include "smear.h"
#include "lnk.h"


int main(int argc,char *argv[])
{
   int my_rank,bc;
   double wt;
   double phi[2],phi_prime[2];
   su3_dble *u;
   FILE *flog=NULL;   

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("time1.log","w",stdout);
   
      printf("\n");
      printf("Timing of the smearing routines\n");
      printf("-------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);      
      fflush(flog);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [time1.c]",
                    "Syntax: time1 [-bc <type>]");
   }
   
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);

   error_root((bc!=0)&&(bc!=3),1,"main [time1.c]",
               "Only open and periodic boundary conditions supported");
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   set_bc_parms(bc,1.0,1.0,1.0,1.0,phi,phi_prime);
   print_bc_parms();
   
   start_ranlux(0,123456);   
   geometry();
   u=udfld();
   
   message("double precision smearing:\n\n");
   
   MPI_Barrier(MPI_COMM_WORLD);
   wt=MPI_Wtime(); 
   ape_links(u,1.0,1,4,0);
   MPI_Barrier(MPI_COMM_WORLD);
   wt=(MPI_Wtime()-wt);
   
   message("ape_links:                      %10.2f sec\n\n",wt);
   
   MPI_Barrier(MPI_COMM_WORLD);
   wt=MPI_Wtime(); 
   hyp_links(1.0,1.0,1.0,1,4,0);
   MPI_Barrier(MPI_COMM_WORLD);
   wt=(MPI_Wtime()-wt);
   
   message("hyp_links:                      %10.2f sec\n\n",wt);
   
   MPI_Barrier(MPI_COMM_WORLD);
   wt=MPI_Wtime(); 
   ape_spatial_links(u,1.0,1,4,0);
   MPI_Barrier(MPI_COMM_WORLD);
   wt=(MPI_Wtime()-wt);
   
   message("ape_spatial_links:              %10.2f sec\n\n",wt);
   
   MPI_Barrier(MPI_COMM_WORLD);
   wt=MPI_Wtime(); 
   hyp_spatial_links(u,1.0,1.0,1,4,0);
   MPI_Barrier(MPI_COMM_WORLD);
   wt=(MPI_Wtime()-wt);
   
   message("hyp_spatial_links:              %10.2f sec\n\n",wt);

   MPI_Finalize();    
   exit(0);
}
