
/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2009, 2014 Bjoern Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
* 
* Simple check and timing off the Wilson loops routine.
* 
* Wloops are calculated on the fake intial gauge configuration
* in which all the links are equal to the identity.
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
#include "smear.h"
#include "wloop.h"

#define NAME_SIZE 128

int main(int argc,char *argv[])
{
   int my_rank,ii,status,bc;
   wloop_parms_t wl;
   double phi[2],phi_prime[2];
   double wt,delta,eps;
   FILE *log=NULL,*fin=NULL;   
  
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  
   
   if (my_rank==0)
   {
      log=freopen("check1.log","w",stdout);
      fin=freopen("check1.in","r",stdin);
      
      printf("\n");
      printf("Consistency checks on the wloop routine\n");
      printf("---------------------------------------\n\n");
      
      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n\n",L0,L1,L2,L3);
      fflush(log);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check1.c]",
                    "Syntax: check1 [-bc <type>]");
   }
   
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);

   error_root((bc!=0)&&(bc!=3),1,"main [check1.c]",
               "Only open and periodic boundary conditions supported");
   
   read_wloop_parms();
   if(my_rank==0)
      fclose(fin);
   wl=wloop_parms();
   
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   set_bc_parms(bc,1.0,1.0,1.0,1.0,phi,phi_prime);
   print_bc_parms();

   geometry();
   
   print_wloop_parms();

   message("HYP-smearing of temporal links\n\n");
   
   MPI_Barrier(MPI_COMM_WORLD);
   wt=MPI_Wtime(); 
   hyp_time_links(wl.alpha_action[0],wl.alpha_action[1],wl.alpha_action[2],wl.proj,wl.proj_iter,0);
   MPI_Barrier(MPI_COMM_WORLD);
   wt=(MPI_Wtime()-wt)/60.0;
  
   message("double precision hyp_time_links successfully calculated in %10.4f min.\n\n",wt);

   MPI_Barrier(MPI_COMM_WORLD);
   wt=MPI_Wtime(); 
   wloop_sum(1);
   MPI_Barrier(MPI_COMM_WORLD);
   wt=(MPI_Wtime()-wt)/60.0;

   free_wloop();
   
   status=0;
   eps=sqrt(3.0)*((double)(NPROC0*NPROC1*L0*L1/4))*DBL_EPSILON;
   for(ii=0; ii<wl.wls; ++ii)
   {   
      if (wl.wl_mat[ii]>eps)
      {
         delta=fabs(1.0-1.0/wl.wl_mat[ii]);
         if (delta>eps)
         {
            message("delta[%2d] = %.6e\n",ii,delta);
            status+=1;
         }
      }
   }
   error(status>0,1,"main [check1.c]",
         "double precision wloops have not been properly calculated!!");
   
   message("double precision wloops successfully calculated in %10.4f min.\n\n",wt);
   
   if(my_rank==0)
      fclose(log);
  
   MPI_Finalize(); 
   exit(0);
}
