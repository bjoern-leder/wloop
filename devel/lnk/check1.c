/*******************************************************************************
*
* File check1.c
*
* Copyright (C) 2009 Bjorn Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Communication routines provided by lnkcom
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
#include "lattice.h"
#include "su3fcts.h"
#include "uflds.h"
#include "global.h"
#include "lnk.h"

static int nfc[8],ofs[8];
static su3_dble *g,*gbuf,*plnk;
static const su3_dble v0={{0.0}};

#define MAX_LEVELS 8
#define BLK_LENGTH 8
#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static int cnt[MAX_LEVELS];
static double smx[MAX_LEVELS];
static su3_dble wd1,wd2,*ud[4];


static double plaq_dble(su3_dble *ub, int ix,int mu,int nu)
{
   double sm;
   su3_vector_dble *v1,*v2;
   
   lnk_plaq(ub,ix,mu,nu,ud);
   
   _su3_times_su3(wd1,(*ud[0]),(*ud[1]));
   _su3_times_su3(wd2,(*ud[2]),(*ud[3]));
   
   v1=(su3_vector_dble*)(&wd1.c11);
   v2=(su3_vector_dble*)(&wd2.c11);
   sm=_vector_prod_re((*(v1  )),(*(v2  )))+
   _vector_prod_re((*(v1+1)),(*(v2+1)))+
   _vector_prod_re((*(v1+2)),(*(v2+2)));
   
   return sm;
}


static double local_plaq_sum_dble(su3_dble *ub)
{
   int ix,mu,nu,n,bc,t;
   double pa;
   
   bc=bc_type();

   for (n=0;n<MAX_LEVELS;n++)
   {
      cnt[n]=0;
      smx[n]=0.0;
   }
   
   for (ix=0;ix<VOLUME;ix++)
   {
      t=global_time(ix);
      pa=0.0;
      
      for (mu=1;mu<4;mu++)
      {
         if ((t<(N0-1))||(bc!=0))
            pa+=plaq_dble(ub,ix,mu,0);

         for (nu=1;nu<mu;nu++)
            pa+=plaq_dble(ub,ix,mu,nu);
      }
      
      cnt[0]+=1;
      smx[0]+=pa;
      
      for (n=1;(cnt[n-1]>=BLK_LENGTH)&&(n<MAX_LEVELS);n++)
      {
         cnt[n]+=1;
         smx[n]+=smx[n-1];
         
         cnt[n-1]=0;
         smx[n-1]=0.0;
      }               
   }
   
   pa=0.0;
   
   for (n=0;n<MAX_LEVELS;n++)
      pa+=smx[n]; 
   
   return pa;
}


static double my_plaq_sum_dble(void)
{
   double p,pa;
   
   p=local_plaq_sum_dble(plnk);
   MPI_Reduce(&p,&pa,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   MPI_Bcast(&pa,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   return pa;
}


static void pack_gbuf(void)
{
   int n,ix,iy,io;

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

   for (n=0;n<8;n++)
   {
      io=ofs[n];

      for (ix=0;ix<nfc[n];ix++)
      {
         iy=map[io+ix];
         gbuf[io+ix]=g[iy];
      }
   }
}


static void send_gbuf(void)
{
   int n,mu,np,saddr,raddr;
   int nbf,tag;
   double *sbuf,*rbuf;
   MPI_Status stat;

   for (n=0;n<8;n++)
   {
      nbf=18*nfc[n];

      if (nbf>0)
      {
         tag=mpi_tag();
         mu=n/2;
         np=cpr[mu];

         if (n==(2*mu))
         {
            saddr=npr[n+1];
            raddr=npr[n];
         }
         else
         {
            saddr=npr[n-1];
            raddr=npr[n];
         }

         sbuf=(double*)(gbuf+ofs[n]);
         rbuf=(double*)(g+ofs[n]+VOLUME);

         if ((np|0x1)!=np)
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
   int ix;
   float *ru,*rm;
   double *rg;
   su3 u;

   for (ix=0;ix<VOLUME;ix++)
   {
      random_su3(&u);

      ru=(float*)(&u);
      rm=ru+18;
      rg=(double*)(g+ix);

      for (;ru<rm;ru++)
      {
         *rg=(double)(*ru);
         rg+=1;
      }

      project_to_su3_dble(g+ix);
   }

   pack_gbuf();
   send_gbuf();
}


static void transform_lnk(void)
{
   int ix,iy,mu;
   su3_dble *ub,u,v,w,gx,gxi,gy,gyi;

   for (ix=(VOLUME/2);ix<VOLUME;ix++)
   {
      ub=lnk(plnk,ix,0);
      gx=g[ix];

      for (mu=0;mu<4;mu++)
      {
         iy=iup[ix][mu];
         gy=g[iy];
         u=ub[2*mu];
         _su3_dagger(gyi,gy);
         _su3_times_su3(v,u,gyi);
         _su3_times_su3(w,gx,v);
         ub[2*mu]=w;

         iy=idn[ix][mu];
         gy=g[iy];
         u=ub[2*mu+1];
         _su3_dagger(gxi,gx);
         _su3_times_su3(v,u,gxi);
         _su3_times_su3(w,gy,v);
         ub[2*mu+1]=w;
      }
   }
}


static void copy_ud2lnk(void)
{
   int ix,mu;
   su3_dble *v,*u,*udb;
   
   udb=udfld();
   for (ix=0;ix<VOLUME;ix++)
   {
      for (mu=0;mu<4;mu++)
      {
         v=lnk(plnk,ix,mu);
         u=lnk(udb,ix,mu);
         *v=*u;
      }
   }
}


static void set_zero(void)
{
   int ix,iy,mu,fc;
   su3_dble *v;
   
   for (ix=0;ix<VOLUME;ix++)
   {
      fc=0;
      for (mu=0;mu<4;mu++)
      {
         iy=idn[ix][mu];
         if (iy>=VOLUME)
         {
            if (fc)
            {
               fc=5;
               break;
            }
            else
               fc=mu+1;
         }
      }
      if (fc)
      {
         for (mu=0;mu<4;mu++)
         {
            if (fc!=(mu+1))
            {
               v=lnk(plnk,ix,mu);
               *v=v0;
            }
         }
      }
   }
}


static void set_zero_u0o(void)
{
   int ix,iy,mu;
   su3_dble *v;
   
   for (ix=0;ix<VOLUME;ix++)
   {
      for (mu=0;mu<4;mu++)
      {
         iy=iup[ix][mu];
         if ((iy>=VOLUME) && (iy<(VOLUME+BNDRY/2)))
         {
            v=lnk(plnk,ix,mu);
            *v=v0;
         }
      }
   }
}


int main(int argc,char *argv[])
{
   int my_rank,bc;
   double p1,p2,eps,nplaq1;
   double phi[2],phi_prime[2];
   FILE *flog=NULL;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      flog=freopen("check1.log","w",stdout);

      printf("\n");
      printf("Communication routines provided by lnkcom\n");
      printf("----------------------------------------------------\n\n");

      printf("%dx%dx%dx%d lattice, ",NPROC0*L0,NPROC1*L1,NPROC2*L2,NPROC3*L3);
      printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
      printf("%dx%dx%dx%d local lattice\n\n",L0,L1,L2,L3);
      fflush(flog);

      bc=find_opt(argc,argv,"-bc");

      if (bc!=0)
         error_root(sscanf(argv[bc+1],"%d",&bc)!=1,1,"main [check1.c]",
                    "Syntax: check1 [-bc <type>]");
   }

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);

   error_root((bc!=0)&&(bc!=3),1,"main [check1.c]",
               "Only open and periodic boundary conditions supported");
   phi[0]=0.0;
   phi[1]=0.0;
   phi_prime[0]=0.0;
   phi_prime[1]=0.0;
   set_bc_parms(bc,1.0,1.0,1.0,1.0,phi,phi_prime);
   print_bc_parms();
   
   start_ranlux(0,12345);
   geometry();
   plnk=alloc_lnk();
   message("Fields allocated\n");

   if (bc==0)
      nplaq1=(double)((6*N0-3)*N1)*(double)(N2*N3);
   else if (bc==3)
      nplaq1=(double)(6*N0*N1)*(double)(N2*N3);
   else
      nplaq1=(double)((6*N0+3)*N1)*(double)(N2*N3);

   eps=sqrt(nplaq1)*(double)(DBL_EPSILON);   
   
   random_ud();
   p1=plaq_sum_dble(1);

   if (my_rank==0)
   {
      printf("Average plaquette after initialization (random) = %.8e\n\n",
             p1/nplaq1);
      fflush(flog);
   }
   
   assign_ud2lnk(plnk);
   p2=lnk_plaq_sum(plnk);
      
   if (my_rank==0)
   {
      printf("lnk_plaq_sum:                  relative difference = %.1e\n",
              fabs(p1-p2));
      fflush(flog);
   }
   
   random_ud();
   p1=plaq_sum_dble(1);
   assign_ud2lnk(plnk);
   copy_bnd(plnk);
   p2=my_plaq_sum_dble();
   
   if (my_rank==0)
   {
      printf("assign_ud2lnk, copy_bnd:       relative difference = %.1e\n",
              fabs(p1-p2));
      fflush(flog);
   }
   
   random_ud();
   assign_lnk2ud(plnk);
   p2=plaq_sum_dble(1);
   
   if (my_rank==0)
   {
      printf("assign_ud2lnk, assign_lnk2ud:  relative difference = %.1e\n",
              fabs(p1-p2));
      fflush(flog);
   }
   
   random_ud();
   p1=plaq_sum_dble(1);
   copy_ud2lnk();
   copyback_u0(plnk);
   copy_u0(plnk);
   copy_uk(plnk);
   p2=my_plaq_sum_dble();
   
   if (my_rank==0)
   {
      printf("copyback_u0, copy_u0, copy_uk: relative difference = %.1e\n",
             fabs(p1-p2));
      fflush(flog);
   }
   
   set_zero();
   copyback_uk(plnk);
   p2=my_plaq_sum_dble();
   
   if (my_rank==0)
   {
      printf("copyback_uk:                   relative difference = %.1e\n",
             fabs(p1-p2));
      fflush(flog);
   }
   
   copy_u0n(plnk);
   set_zero_u0o();
   copyback_u0n(plnk);
   p2=my_plaq_sum_dble();
          
   if (my_rank==0)
   {
      printf("copy_u0n, copyback_u0n:        relative difference = %.1e\n",
             fabs(p1-p2));
      fflush(flog);
   }
   
   g=amalloc(NSPIN*sizeof(su3_dble),4);
   if (BNDRY!=0)
      gbuf=amalloc((BNDRY/2)*sizeof(su3_dble),4);

   error((g==NULL)||((BNDRY!=0)&&(gbuf==NULL)),1,"main [check1.c]",
          "Unable to allocate auxiliary arrays");

   random_g();
   transform_lnk();
   copy_bnd(plnk);
   p2=my_plaq_sum_dble();

   
   if (my_rank==0)
   {
      printf("Gauge invariance:              relative difference = %.1e\n\n",
             fabs(p1-p2));
      printf("(should all be on the order of %.1e)\n\n",eps);
      fflush(flog);
   }
   
   free_lnk(plnk,1);
      
   MPI_Finalize();
   exit(0);
}
