
/*******************************************************************************
*
* File check3.c
*
* Copyright (C) 2009, 2014 Bjoern Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Translational invariance check on the wloop routines
* 
* Note: Translational invariance in time direction is only tested if bc==3 and
*       NPROC0==1.
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
#include "wloop.h"
#include "uflds.h"
#include "smear.h"

#define NAME_SIZE 128

int main(int argc,char *argv[])
{
   int my_rank,ii,ii2,bc,sm;
   int status,mu1,mu2,s1,s2,s[4];
   wloop_parms_t wl;
   double *wls,delta,eps,md;
   double phi[2],phi_prime[2];
   FILE *log=NULL,*fin=NULL;   
  
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  
   if (my_rank==0)
   {
      log=freopen("check3.log","w",stdout);
      fin=freopen("check3.in","r",stdin);
      
      printf("\n");
      printf("Translational invariance check on the wloops routines\n");
      printf("-----------------------------------------------------\n\n");
      
      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n\n",L0,L1,L2,L3);
      fflush(log);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check3.c]",
                    "Syntax: check3 [-bc <type>]");
   }
   
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);

   error_root((bc!=0)&&(bc!=3),1,"main [check3.c]",
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

   print_wloop_parms();

   start_ranlux(0,12345);
   geometry();

   random_ud();
  
   hyp_links(wl.alpha_action[0],wl.alpha_action[1],wl.alpha_action[2],wl.proj,wl.proj_iter,0);

   wls=amalloc(wl.wls*NPROC0*sizeof(double),ALIGN);
   error((wls==NULL),1,"main [check3.c]",
          "Unable to allocate auxiliary arrays");

   wloop_sum(1);
     
   for(ii=0; ii<wl.wls; ++ii)
      wls[ii]=wl.wl_mat[ii];
      
   status=-1;
   sm=0;
   eps=(double)(100)*sqrt((double)(1296*wl.msl*wl.msl*wl.mwlr*wl.mwlt*3*NPROC)*(double)(VOLUME))*(double)(DBL_EPSILON);
   for (mu1=0+(bc!=3);mu1<3;mu1++)
   {
      for (mu2=(mu1+1);mu2<4;mu2++)
      {
         for (s1=0;s1<2;s1++)
         {
            for (s2=0;s2<2;s2++)
            {
               s[0]=0;
               s[1]=0;
               s[2]=0;
               s[3]=0;
               
               s[mu1]=2*s1-1;
               s[mu2]=2*s2-1;
               
               sm+=s[0];

               if (my_rank==0)
                  printf("s = (%d,%d,%d,%d)\n",s[0],s[1],s[2],s[3]);
               
               shift_ud(s);
               wloop_sum(1);
  
               md=0.0;
               for(ii=0; ii<wl.wls; ++ii)
               {
                  ii2=safe_mod(ii-sm*wl.msl*wl.msl*wl.mwlr*wl.mwlt,wl.wls);
                  delta=fabs(1.0-wl.wl_mat[ii]/wls[ii2]);
                  if(delta>eps)
                     status+=2;
                  if (delta>md)
                     md=delta;
               }
               if (my_rank==0)
                  printf("Maximal deviation: %.2e (should be smaller than %.2e)\n\n",md,eps);
            }
         }
      }
   }
   
   afree(wls);
   free_wloop();
   
   error(status>0,1,"main [check3.c]",
         "wloop routine is not translationally invariant!!");
  
   message("wloops routine is translationally invariant!!\n");
  
   MPI_Finalize(); 
   exit(0);
}
