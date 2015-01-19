
/*******************************************************************************
*
* File check2.c
*
* Copyright (C) 2009, 2014 Bjoern Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Gauge invariance check on the wloop routines
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "lattice.h"
#include "su3fcts.h"
#include "random.h"
#include "uflds.h"
#include "global.h"
#include "flags.h"
#include "wloop.h"
#include "smear.h"
#include "lnk.h"

#define NAME_SIZE 128
#define N0 (NPROC0*L0)

static int bc,nfc[8],ofs[8];
static const su3_dble ud0={{0.0}};
static su3_dble *g,*gbuf;
static su3_dble wd ALIGNED16;


static void pack_gbuf(void)
{
   int ifc,ib,ix;

   nfc[0]=FACE0/2;
   nfc[1]=FACE0/2;
   nfc[2]=FACE1/2;
   nfc[3]=FACE1/2;
   nfc[4]=FACE2/2;
   nfc[5]=FACE2/2;
   nfc[6]=FACE3/2;
   nfc[7]=FACE3/2;

   ofs[0]=0;
   ofs[1]=ofs[0]+nfc[0];
   ofs[2]=ofs[1]+nfc[1];
   ofs[3]=ofs[2]+nfc[2];
   ofs[4]=ofs[3]+nfc[3];
   ofs[5]=ofs[4]+nfc[4];
   ofs[6]=ofs[5]+nfc[5];
   ofs[7]=ofs[6]+nfc[6];

   for (ifc=0;ifc<8;ifc++)
   {
      for (ib=0;ib<nfc[ifc];ib++)
      {
         ix=map[ofs[ifc]+ib];
         gbuf[ofs[ifc]+ib]=g[ix];
      }
   }
}


static void send_gbuf(void)
{
   int ifc,np,saddr,raddr;
   int nbf,tag;
   su3_dble *sbuf,*rbuf;
   MPI_Status stat;

   np=cpr[0]+cpr[1]+cpr[2]+cpr[3];

   for (ifc=0;ifc<8;ifc++)
   {
      nbf=18*nfc[ifc];

      if (nbf>0)
      {
         tag=mpi_tag();
         saddr=npr[ifc^0x1];
         raddr=npr[ifc];
         sbuf=gbuf+ofs[ifc];
         rbuf=g+VOLUME+ofs[ifc];

         if (np&0x1)
         {
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
         }
         else
         {
            MPI_Recv(rbuf,nbf,MPI_DOUBLE,raddr,tag,MPI_COMM_WORLD,&stat);
            MPI_Send(sbuf,nbf,MPI_DOUBLE,saddr,tag,MPI_COMM_WORLD);
         }
      }
   }
}


static void random_g(void)
{
   int ix,t;
   su3_dble unity,*gx;

   unity=ud0;
   unity.c11.re=1.0;
   unity.c22.re=1.0;
   unity.c33.re=1.0;
   gx=g;

   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if ((t>0)||(bc!=1))
         random_su3_dble(gx);
      else
         (*gx)=unity;

      gx+=1;
   }

   if (BNDRY>0)
   {
      pack_gbuf();
      send_gbuf();
   }
}


static void transform_ud(void)
{
   int ix,iy,t,ifc;
   su3_dble *u;

   u=udfld();

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      t=global_time(ix);

      if (t==0)
      {
         iy=iup[ix][0];
         su3xsu3dag(u,g+iy,&wd);
         su3xsu3(g+ix,&wd,u);
         u+=1;

         if (bc==3)
         {
            iy=idn[ix][0];
            su3xsu3dag(u,g+ix,&wd);
            su3xsu3(g+iy,&wd,u);
         }
         else if (bc!=0)
         {
            iy=idn[ix][0];
            su3xsu3(g+iy,u,&wd);
            (*u)=wd;
         }

         u+=1;

         for (ifc=2;ifc<8;ifc++)
         {
            if (bc!=1)
            {
               if (ifc&0x1)
               {
                  iy=idn[ix][ifc/2];
                  su3xsu3dag(u,g+ix,&wd);
                  su3xsu3(g+iy,&wd,u);
               }
               else
               {
                  iy=iup[ix][ifc/2];
                  su3xsu3dag(u,g+iy,&wd);
                  su3xsu3(g+ix,&wd,u);
               }
            }

            u+=1;
         }
      }
      else if (t==(N0-1))
      {
         if (bc==3)
         {
            iy=iup[ix][0];
            su3xsu3dag(u,g+iy,&wd);
            su3xsu3(g+ix,&wd,u);
         }
         else if (bc!=0)
         {
            su3xsu3(g+ix,u,&wd);
            (*u)=wd;
         }

         u+=1;

         for (ifc=1;ifc<8;ifc++)
         {
            if (ifc&0x1)
            {
               iy=idn[ix][ifc/2];
               su3xsu3dag(u,g+ix,&wd);
               su3xsu3(g+iy,&wd,u);
            }
            else
            {
               iy=iup[ix][ifc/2];
               su3xsu3dag(u,g+iy,&wd);
               su3xsu3(g+ix,&wd,u);
            }

            u+=1;
         }
      }
      else
      {
         for (ifc=0;ifc<8;ifc++)
         {
            if (ifc&0x1)
            {
               iy=idn[ix][ifc/2];
               su3xsu3dag(u,g+ix,&wd);
               su3xsu3(g+iy,&wd,u);
            }
            else
            {
               iy=iup[ix][ifc/2];
               su3xsu3dag(u,g+iy,&wd);
               su3xsu3(g+ix,&wd,u);
            }

            u+=1;
         }
      }
   }

   set_flags(UPDATED_UD);
}


int main(int argc,char *argv[])
{
   int my_rank,ii,status;
   wloop_parms_t wl;
   double *wls,delta,eps,dmax;
   double phi[2],phi_prime[2];
   su3_dble *plnk;
   FILE *log=NULL,*fin=NULL;   
  
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  
   if (my_rank==0)
   {
      log=freopen("check2.log","w",stdout);
      fin=freopen("check2.in","r",stdin);
      
      printf("\n");
      printf("Gauge invariance check on the wloops routines\n");
      printf("---------------------------------------------\n\n");
      
      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n\n",L0,L1,L2,L3);
      fflush(log);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check2.c]",
                    "Syntax: check2 [-bc <type>]");
   }
   
   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);

   error_root((bc!=0)&&(bc!=3),1,"main [check2.c]",
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

   g=amalloc(NSPIN*sizeof(*g),4);
   if (BNDRY!=0)
      gbuf=amalloc((BNDRY/2)*sizeof(*gbuf),4);

   error((g==NULL)||((BNDRY!=0)&&(gbuf==NULL)),1,"main [check1.c]",
         "Unable to allocate auxiliary arrays");

   start_ranlux(0,12345);
   geometry();

   random_ud();
   plnk=alloc_lnk();
   assign_ud2lnk(plnk);

   hyp_links(wl.alpha_action[0],wl.alpha_action[1],wl.alpha_action[2],wl.proj,wl.proj_iter,0);
  
   wls=amalloc(wl.wls*sizeof(double),ALIGN);
   error((wls==NULL),1,"main [check2.c]",
          "Unable to allocate auxiliary arrays");
      
   message("Measure Wloops matrix\n");
   
   wloop_sum(0);
      
   for(ii=0; ii<wl.wls; ++ii)
      wls[ii]=wl.wl_mat[ii];
      
   message("Perform gauge transformation\n");
   
   assign_lnk2ud(plnk);
   random_g();
   transform_ud();
  
   hyp_links(wl.alpha_action[0],wl.alpha_action[1],wl.alpha_action[2],wl.proj,wl.proj_iter,0);
  
   message("Measure Wloops matrix again\n");
   
   wloop_sum(0);
      
   status=-1;
   eps=(double)(100)*sqrt((double)(wl.msl*wl.msl*wl.mwlr*wl.mwlt*3*NPROC)*(double)(VOLUME))*(double)(DBL_EPSILON);
   dmax=(double)(DBL_EPSILON);
   for (ii=0; ii<wl.wls; ++ii)
   {
      delta=fabs(1.0-wl.wl_mat[ii]/wls[ii]);
      if (delta>=eps)
      {
         message("delta[%2d] = %.6e\n",ii,delta);
         status+=2;
      }
      if (delta>dmax)
         dmax=delta;
   }
   
   if (status>0)
      message("wloop routine is NOT gauge invariant!!");
   else   
      message("wloops routine is gauge invariant (max. deviation: %.1e)!!\n\n",dmax);

   afree(wls);
   afree(g);
   afree(gbuf);
   free_lnk(plnk,0);
  
   MPI_Finalize(); 
   exit(0);
}
