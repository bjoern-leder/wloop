/*******************************************************************************
*
* File lnkcom.c
*
* Copyright (C) 2009 Bjorn Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Allocation of link-fields and commuincation routines.
*
*   void copy_bnd(su3_dble *ub)
*     Fetches the double-precision link variables on the boundaries of the
*     local lattice from the neighbouring processes and stores them in a
*     buffer. The links that stick out of the local lattice at even points
*     are also updated. The pointer ub is the base address of a link-field
*     allocated with alloc_lnk.
*
*   void copy_u0(su3_dble *ub)
*     The links that stick out of the local lattice at even points
*     are updated. The pointer ub is the base address of a link-field
*     allocated with alloc_lnk.
*
*   void copyback_u0(su3_dble *ub)
*     The links that stick out of the local lattice at even points
*     are copied back to the neighbouring processes. The pointer ub is the base
*     address of a link-field allocated with alloc_lnk.
*
*   void copy_uk(su3_dble *ub)
*     Fetches the double-precision link variables on the boundaries of the
*     local lattice from the neighbouring processes and stores them in a
*     buffer. The pointer ub is the base address of a link-field
*     allocated with alloc_lnk.
*
*   void copyback_uk(su3_dble *ub)
*     Reads the double-precision link variables on the boundaries of the
*     local lattice from the buffer and copies them back to the neighbouring
*     processes. The pointer ub is the base address of a link-field
*     allocated with alloc_lnk.
*
*   void send_uk(su3_dble *ub, su3_dble *sbuf)
*     Same as copyback_uk, but read the links from sbuf. The pointer sbuf is the
*     base address of a link-field allocated with alloc_lnk_uk.
*
*   void free_lnk(su3_dble *ub, int cb)
*     Frees the memory starting at the base adress ub, that was allocated by
*     alloc_lnk. For cb unequal zero the internal auxillary buffers are also
*     freed.
*
*   void free_buf_uk(su3_dble *ub)
*     Frees the memory starting at the base adress ub, that was allocated by
*     alloc_lnk_uk.
*
*   su3_dble* alloc_lnk()
*     Allocates memory for a double-precision link-field. The layout is the same
*     as for the field allocated by alloc_ud, but additional memory is allocated
*     for the boundary links. The base address is returned. The size is
*     (4*VOLUME+BNDRY/4+3*BNDRY/2)*sizeof(su3_dble).
*
*   su3_dble* alloc_buf_uk()
*     Allocates memory for double-precision boundary links. The base address is
*     returned. The size is (3*BNDRY/2)*sizeof(su3_dble).
*
*   su3_dble* lnk_buf(su3_dble* ub, int ix, int mu, int face)
*     Returns the address of the a link in a buffer for boundary links. The base
*     adress ub was allocated by alloc_lnk_uk or u+(4*VOLUME+BNDRY/4), where
*     u is a base address of memory allocated with alloc_lnk.
*     ix, mu define the link. And face specifies the boundary face. 
*
*   su3_dble* lnkf(su3_dble* ub, int ix, int mu, int face)
*     Returns the address of the a link in a link-field. The base
*     adress ub was allocated by alloc_lnk.
*     ix, mu define the link. ix>=VOLUME is allowed and then face specifies the
*     boundary face. If ix<VOLUME the value of face is ignored.
*
*   su3_dble* lnk(su3_dble* ub, int ix, int mu)
*     Returns the address of the a link in a link-field. The base
*     adress ub was allocated by alloc_lnk.
*     ix, mu define the link. ix<VOLUME is assumed here.
*
*   void lnk_plaq(su3_dble* ub, int ix,int mu,int nu,su3_dble **u)
*     Calculates the pointers u[4] to the four double-precision link
*     variables in the (mu,nu)-plaquette at the point ix on the local
*     lattice. The values stored at these memory locations are correct
*     only after copy_bnd(ub) is called. The base adress ub was allocated by
*     alloc_lnk.
*
*   void assign_lnk2lnk(su3_dble* ub,su3_dble* ub)
*     Copy the links of the link-field ub1 to the link-field ub2.
*     Also the links that stick out of the local lattice at even points are
*     copied. The base addresses ub1,ub2 were allocated by alloc_lnk or 
*     are the base address of the local gauge field pud[VOLUME/2][0].
*
*   void assign_ud2lnk(su3_dble* ub)
*     Copy the links of the local lattice from pud[VOLUME/2][0] to ub.
*     Also the links that stick out of the local lattice at even points are
*     copied. The base address ub was allocated by alloc_lnk.
*
*   void assign_lnk2ud(su3_dble* ub)
*     Copy the links of the local lattice from ud to pud[VOLUME/2][0] .
*     Also the links that stick out of the local lattice at even points are
*     copied. The base address ub was allocated by alloc_lnk.
*
*   int* get_lnk_idx_u0()
*     Returns the offsets of the links that stick out of the
*     local lattice in the negative directions. These are the links that are sent
*     to the neighbouring processes in copy_u0.
*
*   int* get_lnk_idx_uk()
*     Returns the offsets of the links that that are sent
*     to the neighbouring processes in copy_uk.
*
*
* Notes:
*
* It is possible to allocate as many link-fields as desired. The programs keep
* track of the base addresses. They should be always freed using free_lnk.
*
* The boundary links that are communicated by the routines copy* are only those
* on the boundary in the positive directions.
*
* The link variables in the (mu,nu)-plaquette at the point x are ordered
* according to
*
*   u[0] -> U(x,mu)
*   u[1] -> U(x+mu,nu)
*   u[2] -> U(x,nu)
*   u[3] -> U(x+nu,mu)
*
* In the program lnk_plaq_u() it is taken for granted that the
* arguments satisfy mu!=nu and 0<=ix<VOLUME
*
*******************************************************************************/

#define LNKCOM_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "su3.h"
#include "global.h"
#include "flags.h"
#include <uflds.h>
#include "lnk.h"

typedef struct
{
   int nu0,nuk;
   int saddr,raddr;
   int iu0,iuk;
} comdat_t;

static int nfc[4],ofs[4],snu[4],ofs_uke[4],ofs_uko[4],ofs_u0n[4];
static int *idx_u0=NULL,*idx_uk,*idx_u0n;
static su3_dble *sdbuf_u0=NULL,*sdbuf_uk,*rdbuf_u0,*rdbuf_uk;
static comdat_t comdat[4];
static su3_dble *lnks=NULL;
static su3_dble *tbuf_u0,*tbuf_uk;
static int nu=0,nb=0;
static su3_dble **us,**bufks;


#define UP   1
#define DOWN 2


static void alloc_idx(void)
{
   int n,mu,iu0,iuk;
   comdat_t *c;

   nfc[0]=FACE0;
   nfc[1]=FACE1;
   nfc[2]=FACE2;
   nfc[3]=FACE3;

   ofs[0]=FACE0/2;
   ofs[1]=ofs[0]+(FACE0+FACE1)/2;
   ofs[2]=ofs[1]+(FACE1+FACE2)/2;
   ofs[3]=ofs[2]+(FACE2+FACE3)/2;

   snu[0]=0;
   snu[1]=snu[0]+(FACE0/2);
   snu[2]=snu[1]+(FACE1/2);
   snu[3]=snu[2]+(FACE2/2);

   n=BNDRY/4;
   idx_u0=amalloc(7*n*sizeof(int),3);
   error(idx_u0==NULL,1,"alloc_idx [lnkcom.c]","Unable to allocate index array");

   idx_u0n=amalloc(n*sizeof(int),3);
   error(idx_u0n==NULL,1,"alloc_idx [lnkcom.c]","Unable to allocate index array");

   idx_uk=idx_u0+n;
   iu0=0;
   iuk=0;

   for (mu=0;mu<4;mu++)
   {
      c=comdat+mu;

      (*c).nu0=nfc[mu]/2;
      (*c).nuk=3*nfc[mu];
      (*c).saddr=npr[2*mu];
      (*c).raddr=npr[2*mu+1];
      (*c).iu0=iu0;
      (*c).iuk=iuk;

      ofs_uke[mu]=iuk-3*(VOLUME+ofs[mu]);
      ofs_uko[mu]=iuk+3*(nfc[mu]/2)-3*(VOLUME+(BNDRY/2)+ofs[mu]);
      
      ofs_u0n[mu]=iu0-VOLUME-BNDRY/2-ofs[mu]+nfc[mu]/2;

      iu0+=(*c).nu0;
      iuk+=(*c).nuk;
   }
}


static void set_idx(void)
{
   int mu,nu;
   int ioe,ioo,ix,iye,iyo,izo;
   int nu0,*u0,*uke,*uko,*u0n;

   u0=idx_u0;
   uke=idx_uk;
   
   u0n=idx_u0n;

   for (mu=0;mu<4;mu++)
   {
      nu0=nfc[mu]/2;
      ioe=ofs[mu];
      ioo=ioe+(BNDRY/2);
      uko=uke+3*nu0;

      for (ix=0;ix<nu0;ix++)
      {
         iyo=map[ioo-nu0+ix];
         
         *u0n=8*(iyo-(VOLUME/2))+2*mu;
         u0n+=1;
         
         iye=map[ioe+ix];
         iyo=map[ioo+ix];

         *u0=8*(iyo-(VOLUME/2))+2*mu+1;
         u0+=1;

         for (nu=0;nu<4;nu++)
         {
            if (nu!=mu)
            {
               izo=iup[iye][nu];

               if (izo<VOLUME)
                  *uke=8*(izo-(VOLUME/2))+2*nu+1;
               else
                  *uke=4*VOLUME+comdat[nu].iu0+(izo-VOLUME-(BNDRY/2)-ofs[nu]);

               *uko=8*(iyo-(VOLUME/2))+2*nu;

               uke+=1;
               uko+=1;
            }
         }
      }

      uke=uko;
   }
}


static void alloc_bufs(void)
{
   int n;

   if (idx_u0==NULL)
   {
      alloc_idx();
      set_idx();
   }
   
   n=BNDRY/4;
   sdbuf_u0=amalloc(14*n*sizeof(su3_dble),ALIGN);
   error(sdbuf_u0==NULL,1,"alloc_bufs [lnkcom.c]","Unable to allocate buffers");
   
   sdbuf_uk=sdbuf_u0+n;
   tbuf_u0=sdbuf_uk+6*n;
   tbuf_uk=tbuf_u0+n;
}


static int is_in_list(su3_dble **l, int count, su3_dble *ub)
{
   int i,ret;
   
   ret=-1;
   for (i=0;i<count;i++)
   {
      if (l[i]==ub)
      {
         ret=i;
         break;
      }
   }
   
   return ret;
}


static void add_to_list(su3_dble ***pl, int *count, su3_dble *ub)
{
   int nc,i;
   su3_dble **lnew;
   
   nc=(*count)+1;
   lnew=amalloc(nc*sizeof(su3_dble*),3);
   error(lnew==NULL,1,"add_to_list [lnkcom.c]","Unable to allocate buffer");
   
   for (i=0;i<nc-1;i++)
      lnew[i]=(*pl)[i];
   
   afree(*pl);
   
   lnew[nc-1]=ub;
   *pl=lnew;
   *count=nc;
}


static void del_from_list(su3_dble ***pl, int *count, su3_dble *ub)
{
   int nc,n,i;
   su3_dble **lnew=NULL;
   
   if (*count)
   {
      nc=(*count)-1;
      if (nc>0)
      {
         lnew=amalloc(nc*sizeof(su3_dble*),3);
         error(lnew==NULL,1,"del_from_list [lnkcom.c]","Unable to allocate buffer");
      
         n=0;
         for (i=0;i<nc+1;i++)
         {
            if ((*pl)[i]!=ub)
            {
               lnew[n]=(*pl)[i];
               n+=1;
            }
         }
      }
      
      afree(*pl);
      
      if (nc>0)
         *pl=lnew;
      *count=nc;
   }
}


static void pack_ud0(void)
{
   int *iu,*ium;
   su3_dble *u,*ub;

   u=sdbuf_u0;
   ub=lnks;
   iu=idx_u0;
   ium=idx_u0+(BNDRY/4);

   for (;iu<ium;iu++)
   {
      *u=*(ub+(*iu));
      u+=1;
   }
}


static void pack_ud0n(void)
{
   int *iu,*ium;
   su3_dble *u,*ub;

   u=sdbuf_u0;
   ub=lnks;
   iu=idx_u0n;
   ium=idx_u0n+(BNDRY/4);

   for (;iu<ium;iu++)
   {
      *u=*(ub+(*iu));
      u+=1;
   }
}


static void unpack_ud0(void)
{
   int *iu,*ium;
   su3_dble *u,*ub;

   u=rdbuf_u0;
   ub=lnks;
   iu=idx_u0;
   ium=idx_u0+(BNDRY/4);

   for (;iu<ium;iu++)
   {
      *(ub+(*iu))=*u;
      u+=1;
   }
}


static void unpack_ud0n(void)
{
   int *iu,*ium;
   su3_dble *u,*ub;

   u=rdbuf_u0;
   ub=lnks;
   iu=idx_u0n;
   ium=idx_u0n+(BNDRY/4);

   for (;iu<ium;iu++)
   {
      *(ub+(*iu))=*u;
      u+=1;
   }
}


static void pack_udk(void)
{
   int *iu,*ium;
   su3_dble *u,*ub;

   u=sdbuf_uk;
   ub=lnks;
   iu=idx_uk;
   ium=idx_uk+BNDRY_K;

   for (;iu<ium;iu++)
   {
      *u=*(ub+(*iu));
      u+=1;
   }
}


static void unpack_udk(void)
{
   int *iu,*ium;
   su3_dble *u,*ub;

   u=rdbuf_uk;
   ub=lnks;
   iu=idx_uk;
   ium=idx_uk+BNDRY_K;

   for (;iu<ium;iu++)
   {
      *(ub+(*iu))=*u;
      u+=1;
   }
}


static void set_ud0(su3_dble *buf)
{
   su3_dble *u,*um,*v;
   
   u=sdbuf_u0;
   um=u+(BNDRY/4);
   v=buf;
   for (;u<um;u++)
   {
      *u=*v;
      v+=1;
   }
}


static void set_udk(su3_dble *buf)
{
   su3_dble *u,*um,*v;
   
   u=sdbuf_uk;
   um=u+BNDRY_K;
   v=buf;
   for (;u<um;u++)
   {
      *u=*v;
      v+=1;
   }
}


static void send_ud0(int dir)
{
   int tag,nbf,saddr,raddr,np;
   double *sbuf,*rbuf;
   comdat_t *c,*cm;
   MPI_Status stat;

   np=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;
   c=comdat;
   cm=c+4;

   for (;c<cm;c++)
   {
      nbf=18*(*c).nu0;

      if (nbf>0)
      {
         tag=mpi_tag();
         if (dir==DOWN)
         {
            saddr=(*c).saddr;
            raddr=(*c).raddr;
         }
         else
         {
            saddr=(*c).raddr;
            raddr=(*c).saddr;
         }
         sbuf=(double*)(sdbuf_u0+(*c).iu0);
         rbuf=(double*)(rdbuf_u0+(*c).iu0);

         if (np==0)
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


static void send_udk(int dir)
{
   int tag,nbf,saddr,raddr,np;
   double *sbuf,*rbuf;
   comdat_t *c,*cm;
   MPI_Status stat;

   np=(cpr[0]+cpr[1]+cpr[2]+cpr[3])&0x1;
   c=comdat;
   cm=c+4;

   for (;c<cm;c++)
   {
      nbf=18*(*c).nuk;

      if (nbf>0)
      {
         tag=mpi_tag();
         if (dir==DOWN)
         {
            saddr=(*c).saddr;
            raddr=(*c).raddr;
         }
         else
         {
            saddr=(*c).raddr;
            raddr=(*c).saddr;
         }
         sbuf=(double*)(sdbuf_uk+(*c).iuk);
         rbuf=(double*)(rdbuf_uk+(*c).iuk);

         if (np==0)
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


void copy_u0(su3_dble *ub)
{
   if (NPROC>1)
   {
      if (sdbuf_u0==NULL)
         alloc_bufs();
      
      lnks=ub;
      rdbuf_u0=ub+4*VOLUME;

      pack_ud0();
      send_ud0(DOWN);
   }
}


void copy_u0n(su3_dble *ub)
{
   if (NPROC>1)
   {
      if (sdbuf_u0==NULL)
         alloc_bufs();
      
      lnks=ub;
      rdbuf_u0=ub+(4*VOLUME+BNDRY/4);

      pack_ud0n();
      send_ud0(UP);
   }
}


void copyback_u0(su3_dble *ub)
{
   if (NPROC>1)
   {
      if (sdbuf_u0==NULL)
         alloc_bufs();
      
      lnks=ub;
      rdbuf_u0=tbuf_u0;

      set_ud0(ub+(4*VOLUME));
      send_ud0(UP);
      unpack_ud0();
   }
}


void copyback_u0n(su3_dble *ub)
{
   if (NPROC>1)
   {
      if (sdbuf_u0==NULL)
         alloc_bufs();
      
      lnks=ub;
      rdbuf_u0=tbuf_u0;

      set_ud0(ub+(4*VOLUME+BNDRY/4));
      send_ud0(DOWN);
      unpack_ud0n();
   }
}


void copy_uk(su3_dble *ub)
{
   if (NPROC>1)
   {
      if (sdbuf_u0==NULL)
         alloc_bufs();
      
      lnks=ub;
      rdbuf_uk=ub+LOCAL;

      pack_udk();
      send_udk(DOWN);
   }
}


void copyback_uk(su3_dble *ub)
{
   if (NPROC>1)
   {
      if (sdbuf_u0==NULL)
         alloc_bufs();
      
      lnks=ub;
      rdbuf_uk=tbuf_uk;

      set_udk(ub+LOCAL);
      send_udk(UP);
      unpack_udk();
   }
}


void send_uk(su3_dble *ub, su3_dble *sbuf)
{
   if (NPROC>1)
   {
      if (sdbuf_u0==NULL)
         alloc_bufs();
      
      lnks=ub;
      rdbuf_uk=ub+LOCAL;

      set_udk(sbuf);
      send_udk(UP);
   }
}


void copy_bnd(su3_dble *ub)
{
   copy_u0(ub);
   copy_uk(ub);
}


void free_lnk(su3_dble *ub, int cb)
{
   if ((cb) && (NPROC>1))
   {
      afree(sdbuf_u0);
      sdbuf_u0=NULL;
   
      afree(idx_u0);
      idx_u0=NULL;
   }
   
   if (nu>0)
   {
      if (is_in_list(us,nu,ub)>=0)
      {
         del_from_list(&us,&nu,ub);
         afree(ub);
      }
   }
   
   if (nu==0)
   {
      afree(us);
      us=NULL;
   }
}


void free_buf_uk(su3_dble *ub)
{
   if (nu>0)
   {
      if (is_in_list(bufks,nb,ub)>=0)
      {
         del_from_list(&bufks,&nb,ub);
         afree(ub);
      }
   }
   
   if (nb==0)
   {
      afree(bufks);
      bufks=NULL;
   }
}


su3_dble* alloc_lnk(void)
{
   su3_dble *ub;
   
   if (NPROC>1)
   {
      if (sdbuf_u0==NULL)
         alloc_bufs();
   }
   
   ub=amalloc((LOCAL+BNDRY_K)*sizeof(su3_dble),ALIGN);
   error(ub==NULL,1,"alloc_lnk [lnkcom.c]","Unable to allocate memory for link field");
   
   add_to_list(&us,&nu,ub);
      
   return ub;
}


su3_dble* alloc_buf_uk(void)
{
   su3_dble *ub;
   
   ub=NULL;
   if (NPROC>1)
   {
      if (sdbuf_u0==NULL)
         alloc_bufs();
      
      ub=amalloc(BNDRY_K*sizeof(su3_dble),ALIGN);
      error(ub==NULL,1,"alloc_lnk [lnkcom.c]","Unable to allocate buffers");
      
      add_to_list(&bufks,&nb,ub);
   }
      
   return ub;
}


su3_dble* lnk_buf(su3_dble* ub, int ix, int mu, int face)
{
   su3_dble *u;
   
   if (ix<(VOLUME+(BNDRY/2)))
      u=ub+ofs_uke[face];
   else
      u=ub+ofs_uko[face];

   u+=(3*ix+mu-(mu>face));

   return u;
}


static su3_dble* lnk_u0n(su3_dble* ub, int ix, int mu)
{
   su3_dble *u;
   int iy;
   
   if (ix<(VOLUME+(BNDRY/2)))
   {
      iy=map[ix-VOLUME];
      iy=iup[iy][mu];
      iy=map[iy-VOLUME];
      
      u=ub+8*(iy-VOLUME/2)+2*mu+1;
   }
   else
      u=ub+(4*VOLUME)+(BNDRY/4)+ix+ofs_u0n[mu];
   
   return u;
}


static int ofs_lnk(int ix,int mu)
{
   int iy,ib;

   if (ix<(VOLUME/2))
   {
      iy=iup[ix][mu];

      if (iy<VOLUME)
         return 8*(iy-(VOLUME/2))+2*mu+1;
      else
      {
         ib=iy-VOLUME-ofs[mu]-(BNDRY/2);

         return 4*VOLUME+snu[mu]+ib;
      }
   }
   else
      return 8*(ix-(VOLUME/2))+2*mu;
}


su3_dble* lnkf(su3_dble* ub, int ix, int mu, int face)
{
   su3_dble *u;
   
   if (ix<VOLUME)
      u=ub+ofs_lnk(ix,mu);
   else
   {
      if (mu!=face)
         u=lnk_buf(ub+LOCAL,ix,mu,face);
      else
         u=lnk_u0n(ub,ix,mu);
   }
   
   return u;
}


su3_dble* lnk(su3_dble* ub, int ix, int mu)
{
   return (ub+ofs_lnk(ix,mu));
}


void lnk_plaq(su3_dble* ub,int ix,int mu,int nu,su3_dble **u)
{
   int iy;

   u[0]=lnk(ub,ix,mu);
   iy=iup[ix][mu];
   u[1]=lnkf(ub,iy,nu,mu);
   u[2]=lnk(ub,ix,nu);
   iy=iup[ix][nu];
   u[3]=lnkf(ub,iy,mu,nu);
}


void assign_lnk2lnk(su3_dble* ub1, su3_dble* ub2)
{
   su3_dble *u,*v,*vm;
   
   u=ub2;
   v=ub1;
   vm=v+(4*VOLUME)+(BNDRY/4);
   for (;v<vm;v++)
   {
      *u=*v;
      u+=1;
   }
}


void assign_ud2lnk(su3_dble* ub)
{
   su3_dble *u,*v,*vm;
   
   u=ub;
   v=udfld();
   vm=v+(4*VOLUME)+(BNDRY/4);
   for (;v<vm;v++)
   {
      *u=*v;
      u+=1;
   }
}


void assign_lnk2ud(su3_dble* ub)
{
   su3_dble *u,*v,*vm;
   
   u=ub;
   v=udfld();
   vm=v+(4*VOLUME)+(BNDRY/4);
   for (;v<vm;v++)
   {
      *v=*u;
      u+=1;
   }

   set_flags(UPDATED_UD);
}


int* get_lnk_idx_u0(void)
{
   if (idx_u0==NULL)
      alloc_idx();
   
   return idx_u0;
}


int* get_lnk_idx_uk(void)
{
   if (idx_u0==NULL)
      alloc_idx();
   
   return idx_uk;
}
