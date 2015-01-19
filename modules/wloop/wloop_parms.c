
/*******************************************************************************
*
* File wloop_parms.c
*
* Copyright (C) 2009, 2014 Bjoern Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Wilson loops parameter data base
*
* The externally accessible functions are
*
*   wloop_parms_t set_wloop_parms(int msl,int nss[],int mwlt,int mwlr,op_smear_t op_smear,
*                                 double alpha_action[3],double alpha_ape, double alpha_hyp[2],
*                                 int proj_iter)
*     Sets the wloop parameters and returns a structure containing them. The
*     parameters are
* 
*       msl          number of smearing levels
*       nss[msl]     number of smearing steps per levels
*       mwlt         maximal wloop temporal extension
*       mwlr         maximal wloop spacial extension
*       op_smear     maximal wloop spacial extension
*       alpha_action HYP smearing parameters for static action
*       alpha_ape    APE smearing parameter for spatial links
*       alpha_hyp    HYP smearing parameters for spatial links
*       proj_iter    number of iterations of the SU(3)-projection         
*
*   wloop_parms_t wloop_parms()
*     Returns a structure containing the wloop parameter set.
*
*   void read_wloop_parms()
*     On process 0, this program scans stdin for a line starting with the
*     string "[Wloop]" (after any number of blanks). An error occurs if no such
*     line or more than one is found. The lines 
*
*       msl          <int>
*       nss          <int> <int> ...
*       mwlt         <int>
*       mwlr         <int>
*       alpha_action <double> <double> <double>
*       proj_iter    <int>
*       op_smear     <op_smear_t>
*       ...
*     
*     if APE operator smearing is selected:
*       alpha        <double>
* 
*     if HYP operator smearing is selected:
*       alpha        <double> <double>
*
*     are then read one by one using read_line() [utils/mutils.c].
*     The data are then added to the  data base by calling set_wloop_parms(...).
*
*   void print_wloop_parms()
*     Prints the wloop parameters to stdout on MPI process 0.
*     On MPI processes other than 0, the program does nothing.
*
*   void write_wloop_parms(FILE *fdat)
*     Writes the wloop parameters to the file fdat on MPI process 0.
*
*   void check_wloop_parms(FILE *fdat)
*     Compares the wloop parameters with those stored on the file fdat on
*     MPI process 0, assuming the latter were written to the file by the
*     program write_wloop_parms().
*
* Notes:
*
* For each fixed temporal and spatial extension the msl x msl Wilson loop
* correlation matrix will be measured.
*
* Except for wloop_parms(), the programs in this module perform global
* operations and must be called simultaneously on all MPI processes.
*
*******************************************************************************/

#define WLOOP_PARMS_C

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "utils.h"
#include "flags.h"
#include "global.h"
#include "wloop.h"

static op_smear_t op_smear[]={APE,HYP,OP_SMEARS};
static wloop_parms_t wl={1,{0},0,0,0,OP_SMEARS,{0.0},0.0,{0.0},0,0,NULL};
static double *wlm=NULL;


wloop_parms_t set_wloop_parms(int msl,int nss[],int mwlt,int mwlr,op_smear_t op_smear,
                              double alpha_action[3],double alpha_ape, double alpha_hyp[2],
                              int proj_iter)
{
   int k,*iprms,iprm1,ifail;
   int i;
   double dprms[6];

   error_root((op_smear!=APE)&&(op_smear!=HYP),1,
              "set_wloop_parms [wloop_parms.c]","Unknown type of operator smearing");

   if (NPROC>1)
   {
      iprm1=msl;
      
      MPI_Bcast(&iprm1,1,MPI_INT,0,MPI_COMM_WORLD);
      
      error_root((iprm1!=msl),1,
             "set_wloop_parms [wloop_parms.c]","Parameters are not global");
      
      iprms=amalloc((msl+4)*sizeof(int),3);
      
      for(k=0;k<msl;k++)
         iprms[k]=nss[k];
      iprms[msl+0]=mwlt;
      iprms[msl+1]=mwlr;
      iprms[msl+2]=(int)(op_smear);
      iprms[msl+3]=proj_iter;
      
      dprms[0]=alpha_ape;
      dprms[1]=alpha_hyp[0];
      dprms[2]=alpha_hyp[1];
      dprms[3]=alpha_action[0];
      dprms[4]=alpha_action[1];
      dprms[5]=alpha_action[2];

      MPI_Bcast(iprms,msl+4,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(dprms,6,MPI_DOUBLE,0,MPI_COMM_WORLD);
      
      ifail=0;
      
      ifail|=(dprms[0]!=alpha_ape);      
      ifail|=(dprms[1]!=alpha_hyp[0]);      
      ifail|=(dprms[2]!=alpha_hyp[1]);      
      ifail|=(dprms[3]!=alpha_action[0]);      
      ifail|=(dprms[4]!=alpha_action[1]);      
      ifail|=(dprms[5]!=alpha_action[2]);      
      
      for(k=0;k<msl;k++)
         ifail|=(iprms[k]!=nss[k]);

      ifail|=(iprms[msl]!=mwlt);
      ifail|=(iprms[msl+1]!=mwlr);
      ifail|=(iprms[msl+2]!=(int)(op_smear));
      ifail|=(iprms[msl+3]!=proj_iter);

      error_root(ifail!=0,1,
             "set_wloop_parms [wloop_parms.c]","Parameters are not global");
   }
   
   error_root((msl<1)||(msl>MAX_SM_LEVELS),1,"set_wloop_parms [wloop_parms.c]",
          "max smearing level out of range");

   wl.msl=msl;

   for(k=0;k<wl.msl;++k) 
      wl.nss[k]=nss[k];

   wl.alpha_ape=alpha_ape;
  
   wl.alpha_hyp[0]=alpha_hyp[0];
   wl.alpha_hyp[1]=alpha_hyp[1];
   wl.alpha_action[0]=alpha_action[0];
   wl.alpha_action[1]=alpha_action[1];
   wl.alpha_action[2]=alpha_action[2];
   wl.proj=(proj_iter>0);
   wl.proj_iter=proj_iter;

   error_root((mwlt<1)||(mwlt>NPROC0*L0/2+1),1,"set_wloop_parms [wloop_parms.c]",
          "max wloop time extension out of range");
   error_root((mwlr<1)||(mwlr>NPROC1*L1/2+1)||(mwlr>NPROC2*L2/2+1)||(mwlr>NPROC3*L3/2+1),1,
          "set_wloop_parms [wloop_parms.c]","max wloop spatial extension out of range");

   wl.mwlt=mwlt;
   wl.mwlr=mwlr;
   wl.op_smear=op_smear;
   wl.wls=wl.msl*wl.msl*wl.mwlt*wl.mwlr*L0*NPROC0;
   
   if (wlm!=NULL)
      afree(wlm);
   wlm=amalloc(wl.wls*sizeof(double),ALIGN);
   error_root(wlm==NULL,1,"set_wloop_parms [wloop_parms.c]","Could not allocate wloop matrix");
   for (i=0;i<wl.wls;i++)
      *(wlm+i)=1.0;
   
   wl.wl_mat=wlm;
   
   return wl;
}


wloop_parms_t wloop_parms(void)
{
   return wl;
}


void read_wloop_parms(void)
{
   int my_rank,msl,*nss=NULL,mwlt,mwlr,proj_iter,idr,i;
   double alpha_ape, alpha_hyp[2], alpha_action[3];
   char line[NAME_SIZE];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    
   alpha_ape=0.0;
   alpha_hyp[0]=0.0;
   alpha_hyp[1]=0.0;
   idr=0;
   
   if (my_rank==0)
   {
      find_section("Wloop");

      read_line("op_smear","%s",line);

      if (strcmp(line,"APE")==0)
         idr=0;
      else if (strcmp(line,"HYP")==0)
         idr=1;
      else
      {
         idr=2;
         error_root(1,1,"read_wloop_parms [wloop_parms.c]",
                    "Unknown operator smearing %s",line);
      }

      read_line("mwlt","%d",&mwlt);
      read_line("mwlr","%d",&mwlr);
      read_line("proj_iter","%d",&proj_iter);
      read_dprms("alpha_action",3,alpha_action);
      read_line("msl","%d",&msl);
   }
   
   if (NPROC>1)
   {
      MPI_Bcast(&idr,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&mwlt,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&mwlr,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&proj_iter,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&msl,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(alpha_action,3,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }

   if (idr==0 || idr==1)
   {
      nss=malloc(msl*sizeof(int));
      error((nss==NULL),1,"read_wloop_parms [wloop_parms.c]",
            "Unable to allocated data arrays");
      for (i=0;i<msl;i++)
         nss[i]=0;
   }   
   
   if (my_rank==0)
   {
      
      if (idr==0 || idr==1)
         read_iprms("nss",msl,nss);

      if (idr==0)
         read_line("alpha","%lf",&alpha_ape);

      if (idr==1)
         read_dprms("alpha",2,alpha_hyp);    
   }

   if (NPROC>1)
   {     
      if (idr==0 || idr==1)
         MPI_Bcast(nss,msl,MPI_INT,0,MPI_COMM_WORLD);

      if (idr==0)
         MPI_Bcast(&alpha_ape,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      if (idr==1)
         MPI_Bcast(alpha_hyp,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
   }
   
   set_wloop_parms(msl,nss,mwlt,mwlr,op_smear[idr],alpha_action,alpha_ape,alpha_hyp,proj_iter);
}


void print_wloop_parms(void)
{
   int my_rank,i,n;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      printf("Wloop parameters:\n");

      printf("msl = %d\n",wl.msl);
      
      printf("nss =");
      for (i=0;i<wl.msl;i++)
         printf(" %d",wl.nss[i]);
      printf("\n");

      printf("mwlt = %d\n",wl.mwlt);
      printf("mwlr = %d\n",wl.mwlr);
      
      printf("alpha_action =");
      for (i=0;i<3;i++)
      {
         n=fdigits(wl.alpha_action[i]);
         printf(" %.*f",IMAX(n,1),wl.alpha_action[i]);
      }
      printf("\n");
            
      printf("proj_iter = %d\n",wl.proj_iter);

      printf("operator smearing = ");
      if (wl.op_smear==APE)
      {
         printf("APE\n");
         n=fdigits(wl.alpha_ape);
         printf("alpha = %.*f\n",IMAX(n,1),wl.alpha_ape);
      }
      else if (wl.op_smear==HYP)
      {
         printf("HYP\n");
         
         printf("alpha =");
         for (i=0;i<2;i++)
         {
            n=fdigits(wl.alpha_hyp[i]);
            printf(" %.*f",IMAX(n,1),wl.alpha_hyp[i]);
         }
         printf("\n");
      }
      
      printf("\n");
   }
}


void write_wloop_parms(FILE *fdat)
{
   int my_rank,endian;
   int iw,i;
   stdint_t *istd;
   double dstd[6];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();
   
   istd=amalloc((wl.msl+4)*sizeof(int),3);
        
   if (my_rank==0)
   {
      for(i=0;i<wl.msl;i++)
         istd[i]=(stdint_t)(wl.nss[i]);
      istd[wl.msl+0]=(stdint_t)(wl.mwlt);
      istd[wl.msl+1]=(stdint_t)(wl.mwlr);
      istd[wl.msl+2]=(stdint_t)(wl.op_smear);
      istd[wl.msl+3]=(stdint_t)(wl.proj_iter);
      
      dstd[0]=wl.alpha_ape;
      dstd[1]=wl.alpha_hyp[0];
      dstd[2]=wl.alpha_hyp[1];
      dstd[3]=wl.alpha_action[0];
      dstd[4]=wl.alpha_action[1];
      dstd[5]=wl.alpha_action[2];

      if (endian==BIG_ENDIAN)
      {
         bswap_int(wl.msl+4,istd);
         bswap_double(6,dstd);
      }
      
      iw=fwrite(istd,sizeof(stdint_t),wl.msl+4,fdat);
      iw+=fwrite(dstd,sizeof(double),6,fdat);
      error_root(iw!=wl.msl+10,1,"write_wloop_parms [wloop_parms.c]",
                  "Incorrect write count");
   }
}


void check_wloop_parms(FILE *fdat)
{
   int my_rank,endian;
   int ir,ie,i;
   stdint_t *istd;
   double dstd[6];

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
   endian=endianness();
   
   istd=amalloc((wl.msl+4)*sizeof(int),3);

   if (my_rank==0)
   {
      ie=0;
      
      ir=fread(istd,sizeof(stdint_t),wl.msl+4,fdat);
      ir+=fread(dstd,sizeof(double),6,fdat);
      error_root(ir!=wl.msl+10,1,"check_wloop_parms [wloop_parms.c]",
                  "Incorrect read count");

      if (endian==BIG_ENDIAN)
      {
         bswap_int(wl.msl+4,istd);
         bswap_double(6,dstd);
      }
            
      for(i=0;i<wl.msl;i++)
         ie|=(istd[i]!=(stdint_t)(wl.nss[i]));
      ie|=(istd[wl.msl+0]!=(stdint_t)(wl.mwlt));
      ie|=(istd[wl.msl+1]!=(stdint_t)(wl.mwlr));
      ie|=(istd[wl.msl+2]!=(stdint_t)(wl.op_smear));
      ie|=(istd[wl.msl+3]!=(stdint_t)(wl.proj_iter));
      
      ie|=(dstd[0]!=wl.alpha_ape);
      ie|=(dstd[1]!=wl.alpha_hyp[0]);
      ie|=(dstd[2]!=wl.alpha_hyp[1]);
      ie|=(dstd[3]!=wl.alpha_action[0]);
      ie|=(dstd[4]!=wl.alpha_action[1]);
      ie|=(dstd[5]!=wl.alpha_action[2]);
         
      error_root(ie!=0,1,"check_wloop_parms [wloop_parms.c]",
                 "Parameters do not match");         
   }
}

int wl_idx(int sl1, int sl2, int t, int r, int x0)
{
   return (sl1+wl.msl*(sl2+wl.msl*(t+wl.mwlt*(r+wl.mwlr*(cpr[0]*L0+x0)))));
}
