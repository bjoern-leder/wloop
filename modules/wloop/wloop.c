/*******************************************************************************
*
* File wloop.c
*
* Copyright (C) 2009 Bjorn Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of Wilson loops
*
* The externally accessible functions are
*
* void wloop_sum(int keep)
*   Calculates the Wilson loops correlation matrix. This matrix is stored
*   in the array wl.wl_mat[ii]; the index ii is calculated by using the
*   function wl_idx(sl1,sl2,t,r,x0) defined in the file wloop_parms.c.
*   If keep equal zero the auxillary memory will be freed at the end,
*   otherwise it is kept allocated.
*
* void free_wloop(void)
*   Free auxillary memory.
* 
* Notes:
*   The Wloops are calculated in double precision.
*   Wilson loops extend along the time direction up to wld.mwlt and along
*   the speatial directions up to wld.mwlr. 
*   The spatial long links are smearead independently up to wld.msl.
*   At fixed ii the results are obtained by an average over the 3 spatial
*   directions.
*   The "structural parameters" wl.mwlt, wl.mwlr, wl.msl are set
*   through the function set_wloop_parms defined in to the file wloop_parms.c
*   and are stored in the structure wl that is accessible through the
*   function wloop_parms();
*
*******************************************************************************/

#define WLOOP_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "lattice.h"
#include "wloop.h"
#include "global.h"
#include "su3fcts.h"
#include "uflds.h"
#include "flags.h"
#include "lnk.h"

typedef struct {
   int sl1,sl2;
   int x[4];
   su3_dble *ul[6*L0],*um[6*L0],*ub[L0],*ubR[3*L0];
} wloop_t;

#define idx_u(x0,k,s) \
         (6*x0+2*k+s)

static su3_dble *ud_long=NULL,*ud_short,**pus,*ubuf,*mat,*u0,*uk;
static int T,R,bs[4]={L0,L1,L2,L3},nproc[4]={NPROC0,NPROC1,NPROC2,NPROC3};
static wloop_parms_t wl;
static wloop_t w;
static int icpr[4],ipr;
static int debug1=0,debug2=0;


static void alloc_buf(void)
{
   int mu;
   
   if (ud_long==NULL)
   {
      ud_long=alloc_lnk();
      ud_short=alloc_lnk();
      
      ubuf=amalloc((10*L0+6)*sizeof(su3_dble),ALIGN);
      error_root(ubuf==NULL,1,"alloc_buf [wloop.c]",
                 "Unable to allocate memory for link field");
         
      uk=ubuf;
      u0=&ubuf[6*L0];
      mat=&ubuf[10*L0];
   }
   
   for(mu=0;mu<4;++mu)
      icpr[mu]=cpr[mu];
   ipr=ipr_global(cpr);
}


static void assign_uk2pus(su3_dble *uf, int sl)
{
   int ix,k;
   su3_dble *v,*u;

   copy_u0(uf);
   
   for(ix=0;ix<VOLUME;++ix)
   {
      for(k=1;k<4;++k)
      {
         v=lnk(uf,ix,k);
         u=pus[pindx(sl,ix,k-1)];
         
         *u=*v;
      }
   }
}


static void assign_pus2uk(int sl, su3_dble *ut)
{
   int ix,k;
   su3_dble *v,*u;

   for(ix=0;ix<VOLUME;++ix)
   {
      for(k=1;k<4;++k)
      {
         v=lnk(ut,ix,k);
         u=pus[pindx(sl,ix,k-1)];
         
         *v=*u;
      }
   }
   
   copyback_u0(ut);
}


static void assign_u02lnk(su3_dble *ub)
{
   int ix;
   su3_dble *v,*u,*udb;

   udb=udfld();
   
   for(ix=0;ix<VOLUME;++ix)
   {
      v=lnk(ub,ix,0);
      u=lnk(udb,ix,0);
      
      *v=*u;
   }
   
   copyback_u0(ub);
}


static void compute_u0_long(void)
{
   int ixs,ix,iy,x0;
   su3_dble *u1,*u2,*ul;
   
   copy_u0(ud_long);
   copy_u0n(ud_long);
   
   for (ixs=0;ixs<SLICE0;ixs++)
   {
      ix=ipt[ixs];
      iy=idn[ix][0];
      
      u1=lnkf(ud_short,iy,0,0);
#if (defined SSE2)
      _prefetch_su3_dble(u1);
#endif      
      ul=lnkf(ud_long,iy,0,0);
      if (NPROC0==1)
         mat[0]=*ul;
      
/*      if (ixs==0)
      {
         printf("ul, ixs=%d, cpr[0]=%d:\n",ixs,cpr[0]);
         printf("%.5f + %.5fi\t %.5f + %.5fi\t %.5f + %.5fi\n",(*ul).c11.re,(*ul).c11.im,(*ul).c12.re,(*ul).c12.im,(*ul).c13.re,(*ul).c13.im);
         printf("%.5f + %.5fi\t %.5f + %.5fi\t %.5f + %.5fi\n",(*ul).c21.re,(*ul).c21.im,(*ul).c22.re,(*ul).c22.im,(*ul).c23.re,(*ul).c23.im);
         printf("%.5f + %.5fi\t %.5f + %.5fi\t %.5f + %.5fi\n",(*ul).c31.re,(*ul).c31.im,(*ul).c32.re,(*ul).c32.im,(*ul).c33.re,(*ul).c33.im);
      }*/
      
      for (x0=0;x0<L0;x0++)
      {
         if ((NPROC0>1)||(x0<(L0-1)))
            u2=lnk(ud_long,ix,0);
         else
            u2=&mat[0];

         su3xsu3(u1,u2,ul);
         
         iy=ix;
         ix=iup[ix][0];
         
         u1=lnk(ud_short,iy,0);
#if (defined SSE2)
         _prefetch_su3_dble(u1);
#endif      
         ul=lnk(ud_long,iy,0);
      }
   }

   copy_u0(ud_long);
   copyback_u0n(ud_long);
}


static void compute_uk_long(int sl)
{
   int ix,iy,mu,sll;
   su3_dble *u1,*u2,*ul;
   
   assign_pus2uk(sl,ud_short);
   
   copy_u0(ud_short);
   copy_u0n(ud_short);
   copy_u0n(ud_long);
   
   sll=sl;
   if (R>2)
      sll=wl.msl+sl;
   
   for (ix=0;ix<VOLUME;ix++)
   {
      for (mu=1;mu<4;mu++)
      {
         iy=iup[ix][mu];
         if (iy<VOLUME)
         {
            u1=lnk(ud_short,ix,mu);
            u2=pus[pindx(sll,iy,mu-1)];
#if (defined SSE2)
            _prefetch_su3_dble(u1);
            _prefetch_su3_dble(u2);
#endif      
            ul=lnk(ud_long,ix,mu);
            
            su3xsu3(u1,u2,ul);
         }
         
         iy=idn[ix][mu];
         if (iy>=VOLUME)
         {
            u1=lnkf(ud_short,iy,mu,mu);
            u2=pus[pindx(sll,ix,mu-1)];
#if (defined SSE2)
            _prefetch_su3_dble(u1);
            _prefetch_su3_dble(u2);
#endif      
            ul=lnkf(ud_long,iy,mu,mu);
            
            su3xsu3(u1,u2,ul);
         }
      }
   }

   copyback_u0n(ud_long);
   assign_uk2pus(ud_long,wl.msl+sl);
}


static void get_wloop(void)
{
   int iu,ix,x0,xo,xp,dpr,spr,rpr,k,m,tag,io;
   double *sb,*rb;
   su3_dble *u;
   MPI_Status stat;

   for (x0=0;x0<L0;x0++)
   {
      ix=ipt[(w.x[3])+L3*(w.x[2])+L2*L3*(w.x[1])+L1*L2*L3*x0];
      w.ub[x0]=lnk(ud_long,ix,0);
      for (k=0;k<3;k++)
      {
         iu=idx_u(x0,k,0);
         w.ul[iu]=pus[pindx(w.sl1,ix,k)];
         w.ul[iu+1]=pus[pindx(w.sl2,ix,k)];
      }
   
      xp=x0+T;
      xo=xp%L0;
      dpr=xp/L0;
      
      ix=ipt[(w.x[3])+L3*(w.x[2])+L2*L3*(w.x[1])+L1*L2*L3*xo];
      for (k=0;k<3;k++)
      {
         iu=idx_u(x0,k,0);
         
         w.um[iu]=pus[pindx(w.sl2,ix,k)];
         w.um[iu+1]=pus[pindx(w.sl1,ix,k)];
      }
      
      if (NPROC0>1 && dpr>0)
      {
         icpr[0]-=1;
         spr=ipr_global(icpr);
         icpr[0]+=2;
         rpr=ipr_global(icpr);
         icpr[0]=cpr[0];
         
         iu=idx_u(x0,0,0);
         
/*         if (x0==L0-T && w.x[1]==0 && w.x[2]==0 && w.x[3]==0)
         {
            u=w.um[iu];
            printf("x0: %d, x1: %d, x2: %d, x3: %d\n",x0,w.x[1],w.x[2],w.x[3]);
            printf("%d send w.um[%d]:\n",cpr[0],iu);
            printf("%.5f + %.5fi\t %.5f + %.5fi\t %.5f + %.5fi\n",(*u).c11.re,(*u).c11.im,(*u).c12.re,(*u).c12.im,(*u).c13.re,(*u).c13.im);
            printf("%.5f + %.5fi\t %.5f + %.5fi\t %.5f + %.5fi\n",(*u).c21.re,(*u).c21.im,(*u).c22.re,(*u).c22.im,(*u).c23.re,(*u).c23.im);
            printf("%.5f + %.5fi\t %.5f + %.5fi\t %.5f + %.5fi\n",(*u).c31.re,(*u).c31.im,(*u).c32.re,(*u).c32.im,(*u).c33.re,(*u).c33.im);
         }*/
         
         io=dpr%2;
         if (io==0)
            u=&uk[iu];
         else
            u=mat;
         
         for (k=0;k<6;k++)
            *(u+k)=*(w.um[iu+k]);
         
         debug1+=1;
         
         for (k=io;k<(dpr+io);k++)
         {
            tag=mpi_tag();
            
            if (k%2==0)
            {
               sb=(double*)&uk[iu];
               rb=(double*)mat;
            }
            else
            {
               sb=(double*)mat;
               rb=(double*)&uk[iu];
            }
                  
            if (icpr[0]%2==0)
            {
               MPI_Send(sb,108,MPI_DOUBLE,spr,tag,MPI_COMM_WORLD);
               MPI_Recv(rb,108,MPI_DOUBLE,rpr,tag,MPI_COMM_WORLD,&stat);
            }
            else
            {
               MPI_Recv(rb,108,MPI_DOUBLE,rpr,tag,MPI_COMM_WORLD,&stat);
               MPI_Send(sb,108,MPI_DOUBLE,spr,tag,MPI_COMM_WORLD);
            }
         }
                  
         for (k=0;k<6;k++)            
            w.um[iu+k]=&uk[iu+k];         
         
/*         if (x0==L0-T && w.x[1]==0 && w.x[2]==0 && w.x[3]==0)
         {
            u=w.um[iu];
            printf("x0: %d, x1: %d, x2: %d, x3: %d\n",x0,w.x[1],w.x[2],w.x[3]);
            printf("%d rcvd w.um[%d]:\n",cpr[0],iu);
            printf("%.5f + %.5fi\t %.5f + %.5fi\t %.5f + %.5fi\n",(*u).c11.re,(*u).c11.im,(*u).c12.re,(*u).c12.im,(*u).c13.re,(*u).c13.im);
            printf("%.5f + %.5fi\t %.5f + %.5fi\t %.5f + %.5fi\n",(*u).c21.re,(*u).c21.im,(*u).c22.re,(*u).c22.im,(*u).c23.re,(*u).c23.im);
            printf("%.5f + %.5fi\t %.5f + %.5fi\t %.5f + %.5fi\n",(*u).c31.re,(*u).c31.im,(*u).c32.re,(*u).c32.im,(*u).c33.re,(*u).c33.im);
         }*/
      }
   }
      
   for (k=0;k<3;k++)
   {
      xp=w.x[k+1]+R;
      xo=w.x[k+1];
      w.x[k+1]=xp%bs[k+1];
      dpr=xp/bs[k+1];
      
      for (x0=0;x0<L0;x0++)
      {
         ix=ipt[(w.x[3])+L3*(w.x[2])+L2*L3*(w.x[1])+L1*L2*L3*x0];
         w.ubR[L0*k+x0]=lnk(ud_long,ix,0);
      }
      
      w.x[k+1]=xo;
      
      if (nproc[k+1]>1 && dpr>0)
      {
         icpr[k+1]-=1;
         spr=ipr_global(icpr);
         icpr[k+1]+=2;
         rpr=ipr_global(icpr);
         icpr[k+1]=cpr[k+1];
            
         io=dpr%2;
         if (io==0)
            u=&u0[L0*k];
         else
            u=&u0[3*L0];
         
         for (x0=0;x0<L0;x0++)
            *(u+x0)=*(w.ubR[L0*k+x0]);
         
         debug2+=1;
            
         for (m=io;m<(dpr+io);m++)
         {
            tag=mpi_tag();
            
            if (m%2==0)
            {
               sb=(double*)&u0[L0*k];
               rb=(double*)&u0[3*L0];
            }
            else
            {
               sb=(double*)&u0[3*L0];
               rb=(double*)&u0[L0*k];
            }
            
            if (icpr[k+1]%2==0)
            {
               MPI_Send(sb,L0*18,MPI_DOUBLE,spr,tag,MPI_COMM_WORLD);
               MPI_Recv(rb,L0*18,MPI_DOUBLE,rpr,tag,MPI_COMM_WORLD,&stat);
            }
            else
            {
               MPI_Recv(rb,L0*18,MPI_DOUBLE,rpr,tag,MPI_COMM_WORLD,&stat);
               MPI_Send(sb,L0*18,MPI_DOUBLE,spr,tag,MPI_COMM_WORLD);
            }
         }
         
         for (x0=0;x0<L0;x0++)
         {
            iu=L0*k+x0;
            w.ubR[iu]=&u0[iu];         
         }
      }
   }
}


static void wloop_fixed_t_r_sl(int sl1,int sl2,double *res)
{
   int k,x0,x1,x2,x3,iu,i;
   double trc;

	for (i=0;i<L0;i++)
	{
		res[2*i]=0.0;
		res[2*i+1]=0.0;
	}
   
   w.sl1=sl1;
   w.sl2=sl2;
   
   for (x1=0;x1<L1;x1++)
   {
      w.x[1]=x1;
      for (x2=0;x2<L2;x2++)
      {
         w.x[2]=x2;
         for (x3=0;x3<L3;x3++)
         {
            w.x[3]=x3;
            
            get_wloop();
            
            for (x0=0;x0<L0;x0++)
            {
#if (defined SSE2)
               _prefetch_su3_dble(w.ub[x0]);
#endif      
               for (k=0;k<3;k++)
               {
                  iu=idx_u(x0,k,0);
#if (defined SSE2)
                  _prefetch_su3_dble(w.um[iu]);
                  _prefetch_su3_dble(w.ul[iu]);
                  _prefetch_su3_dble(w.ubR[L0*k+x0]);
#endif      
                  
                  su3xsu3(w.ub[x0],w.um[iu],&mat[0]);
#if (defined SSE2)
                  if (sl1!=sl2)
                     _prefetch_su3_dble(w.um[iu+1]);
#endif      
                  su3xsu3(w.ul[iu],w.ubR[L0*k+x0],&mat[1]);
                  
                  trc=  _vector_prod_re(*((su3_vector_dble*)&(mat[0].c11)),*((su3_vector_dble*)&(mat[1].c11)))+
                        _vector_prod_re(*((su3_vector_dble*)&(mat[0].c21)),*((su3_vector_dble*)&(mat[1].c21)))+
                        _vector_prod_re(*((su3_vector_dble*)&(mat[0].c31)),*((su3_vector_dble*)&(mat[1].c31)));
                  
                  res[x0]+=trc;
                  
                  if (sl1!=sl2)
                  {
#if (defined SSE2)
                     _prefetch_su3_dble(w.ul[iu+1]);
#endif      
                     su3xsu3(w.ub[x0],w.um[iu+1],&mat[0]);
                     su3xsu3(w.ul[iu+1],w.ubR[L0*k+x0],&mat[1]);
                     
                     trc=  _vector_prod_re(*((su3_vector_dble*)&(mat[0].c11)),*((su3_vector_dble*)&(mat[1].c11)))+
                           _vector_prod_re(*((su3_vector_dble*)&(mat[0].c21)),*((su3_vector_dble*)&(mat[1].c21)))+
                           _vector_prod_re(*((su3_vector_dble*)&(mat[0].c31)),*((su3_vector_dble*)&(mat[1].c31)));
                  
							res[x0+L0]+=trc;
						}                                    
               }
            }
         }
      }
   }
}


void free_wloop(void)
{
   free_lnk(ud_long,0);
   ud_long=NULL;
   free_lnk(ud_short,1);
   ud_short=NULL;
   afree(ubuf);
   ubuf=NULL;
   
   /*free_smear();*/
}


static void wloops(int keep)
{
   int t,r,sl1,sl2,ii,i;
   double norm,tmp[2*L0];

   copy_bnd_ud();
   
   alloc_buf();
   assign_u02lnk(ud_short);
   copy_u0n(ud_short);
   
   norm=(double)9*L1*L2*L3;
    
   wl=wloop_parms();
   
   smear();

   pus=get_pus();

   R=0;
         
   for(r=0;r<wl.mwlr;++r)
   {
      R+=1;
      
      T=0;
      assign_u02lnk(ud_long);
      
      for(t=0;t<wl.mwlt;++t)
      {
         T+=1;
         if (T>1)
            compute_u0_long();
         
         for(sl1=0;sl1<wl.msl;++sl1)
         {
            for(sl2=sl1;sl2<wl.msl;++sl2)
            {
               if ((R>1) && (T==1) && (sl1==0))
                  compute_uk_long(sl2);

               debug1=0;
               debug2=0;
               
               if (R>1)
                  wloop_fixed_t_r_sl(wl.msl+sl1,wl.msl+sl2,tmp);
               else
                  wloop_fixed_t_r_sl(sl1,sl2,tmp);
         
               /*if (NPROCS>1)
               {
                  MPI_Reduce(&tmp,&glob,2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
                  MPI_Bcast(&glob,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
               }
               else
               {
                  glob[0]=tmp[0];
                  glob[1]=tmp[1];
               }*/

               for (i=0;i<L0;i++)
               {
                  ii=wl_idx(sl1,sl2,t,r,i);
                  wl.wl_mat[ii]=*(tmp+i)/norm;
               }
         
               if (sl1!=sl2)
               {
                  for (i=0;i<L0;i++)
                  {
                     ii=wl_idx(sl2,sl1,t,r,i);
                     wl.wl_mat[ii]=*(tmp+i+L0)/norm;
                  }
               }
         
/*               printf("wl_mat[%d] = %.5e (sl1: %d, sl2: %d, t: %d, r: %d)\n",ii,glob[0],sl1,sl2,t,r);
               printf("debug1: %d, debug2: %d\n",debug1,debug2);*/
            }
         }
      }
   }
   
   if (keep==0)
      free_wloop();
}


void wloop_sum(int keep)
{
   int i,it,ix,iy,iz,tag,n[4],my_rank,wls_local;
   int sedr,recr;
   double *pwl,*tmp;
   MPI_Status stat;

   wloops(keep);
      
   if (NPROC>1)
   {   
      MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

      n[0]=cpr[0];
      n[1]=0;
      n[2]=0;
      n[3]=0;
      recr=ipr_global(n);
      
      pwl=wl.wl_mat + wl_idx(0, 0, 0, 0, 0);
      wls_local=wl.wls/NPROC0;
      
      tmp=amalloc(wls_local*sizeof(double),ALIGN);
      error_root(tmp==NULL,1,"wloop_sum [wloop.c]",
               "Unable to allocate memory for wloop sum");

      for (ix=0;ix<NPROC1;ix++)
      {
         n[1]=ix;
         for (iy=0;iy<NPROC2;iy++)
         {
            n[2]=iy;
            for (iz=0;iz<NPROC3;iz++)
            {
               n[3]=iz;
               sedr=ipr_global(n);
               
               tag=mpi_tag();
               
               if ((my_rank==sedr)&&(recr!=sedr))
                  MPI_Send(pwl,wls_local,MPI_DOUBLE,recr,tag,MPI_COMM_WORLD);
               
               if ((my_rank==recr)&&(recr!=sedr))
               {
                  MPI_Recv(tmp,wls_local,MPI_DOUBLE,sedr,tag,MPI_COMM_WORLD,&stat);
                  for (i=0;i<wls_local;i++)
                     pwl[i]+=tmp[i];
               }
               
            }
         }  
      }
    
      n[1]=0;
      n[2]=0;
      n[3]=0;
      for (it=1;it<NPROC0;it++)
      {
         n[0]=it;
         sedr=ipr_global(n);
         
         tag=mpi_tag();
         
         if (my_rank==sedr)
            MPI_Send(pwl,wls_local,MPI_DOUBLE,0,tag,MPI_COMM_WORLD);
         
         if (my_rank==0)
            MPI_Recv(wl.wl_mat+it*wls_local,wls_local,MPI_DOUBLE,sedr,tag,MPI_COMM_WORLD,&stat);
      }
      
      MPI_Bcast(wl.wl_mat,wl.wls,MPI_DOUBLE,0,MPI_COMM_WORLD);

      for (i=0;i<wl.wls;i++)
         *(wl.wl_mat+i)/=(double)(NPROC1*NPROC2*NPROC3);
      
      afree(tmp);
   }
}

