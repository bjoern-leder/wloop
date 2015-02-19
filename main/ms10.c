
/*******************************************************************************
*
* File ms10.c
*
* Copyright (C) 2012, 2013 Martin Luescher, 2014 Bjoern Leder
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Computation of Wilson loops.
*
* Syntax: ms10 -i <input file> [-noexp] [-a]
*
* For usage instructions see the file README.ms10.
*
*******************************************************************************/

#define MAIN_PROGRAM

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
#include "flags.h"
#include "random.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "archive.h"
#include "smear.h"
#include "wloop.h"
#include "version.h"
#include "global.h"

#define N0 (NPROC0*L0)
#define N1 (NPROC1*L1)
#define N2 (NPROC2*L2)
#define N3 (NPROC3*L3)

static struct
{
   int wls,msl,wlmt,wlmr,tmax;
} file_head;

static struct
{
   int nc;
   double *wl_mat;
} data;

static int my_rank,noexp,append,endian;
static int first,last,step;
static int ipgrd[2],flint;

static char line[NAME_SIZE];
static char log_dir[NAME_SIZE],dat_dir[NAME_SIZE];
static char loc_dir[NAME_SIZE],cnfg_dir[NAME_SIZE];
static char log_file[NAME_SIZE],log_save[NAME_SIZE],end_file[NAME_SIZE];
static char par_file[NAME_SIZE],par_save[NAME_SIZE];
static char dat_file[NAME_SIZE],dat_save[NAME_SIZE];
static char cnfg_file[NAME_SIZE],nbase[NAME_SIZE];
static FILE *fin=NULL,*flog=NULL,*fdat=NULL,*fend=NULL;

static bc_parms_t bcp;
static wloop_parms_t wl;


static void alloc_data(void)
{
   int wls;
   double *p;

   wls=file_head.wls;

   p=amalloc(wls*sizeof(*p),4);

   error((p==NULL),1,"alloc_data [ms10.c]",
         "Unable to allocate data arrays");

   data.wl_mat=p;
}


static void write_file_head(void)
{
   int iw;
   stdint_t istd[5];

   istd[0]=(stdint_t)(file_head.wls);
   istd[1]=(stdint_t)(file_head.msl);
   istd[2]=(stdint_t)(file_head.wlmt);
   istd[3]=(stdint_t)(file_head.wlmr);
   istd[4]=(stdint_t)(file_head.tmax);

   if (endian==BIG_ENDIAN)
      bswap_int(5,istd);

   iw=fwrite(istd,sizeof(stdint_t),5,fdat);

   error_root(iw!=5,1,"write_file_head [ms10.c]",
              "Incorrect write count");
}


static void check_file_head(void)
{
   int ir;
   stdint_t istd[5];

   ir=fread(istd,sizeof(stdint_t),5,fdat);

   error_root(ir!=5,1,"check_file_head [ms10.c]",
              "Incorrect read count");

   if (endian==BIG_ENDIAN)
      bswap_int(5,istd);

   error_root(((int)(istd[0])!=file_head.wls)||
              ((int)(istd[1])!=file_head.msl)||
              ((int)(istd[2])!=file_head.wlmt)||
              ((int)(istd[3])!=file_head.wlmr)||
              ((int)(istd[4])!=file_head.tmax),1,"check_file_head [ms10.c]",
              "Unexpected value of dn,nn,tmax or eps");
}


static void write_data(void)
{
   int iw,wls,i;
   stdint_t istd[1];
   double dstd[1];

   istd[0]=(stdint_t)(data.nc);

   if (endian==BIG_ENDIAN)
      bswap_int(1,istd);

   iw=fwrite(istd,sizeof(stdint_t),1,fdat);

   wls=file_head.wls;

   for (i=0;i<wls;i++)
   {
      dstd[0]=data.wl_mat[i];

      if (endian==BIG_ENDIAN)
         bswap_double(1,dstd);

      iw+=fwrite(dstd,sizeof(double),1,fdat);
   }

   error_root(iw!=(1+wls),1,"write_data [ms10.c]",
              "Incorrect write count");
}


static int read_data(void)
{
   int ir,wls,i;
   stdint_t istd[1];
   double dstd[1];

   ir=fread(istd,sizeof(stdint_t),1,fdat);

   if (ir!=1)
      return 0;

   if (endian==BIG_ENDIAN)
      bswap_int(1,istd);

   data.nc=(int)(istd[0]);

   wls=file_head.wls;

   for (i=0;i<wls;i++)
   {
      ir+=fread(dstd,sizeof(double),1,fdat);

      if (endian==BIG_ENDIAN)
         bswap_double(1,dstd);

      data.wl_mat[i]=dstd[0];
   }

   error_root(ir!=(1+wls),1,"read_data [ms10.c]",
              "Read error or incomplete data record");

   return 1;
}


static void read_dirs(void)
{
   if (my_rank==0)
   {
      find_section("Run name");
      read_line("name","%s",nbase);

      find_section("Directories");
      read_line("log_dir","%s",log_dir);
      read_line("dat_dir","%s",dat_dir);

      if (noexp)
      {
         read_line("loc_dir","%s",loc_dir);
         cnfg_dir[0]='\0';
      }
      else
      {
         read_line("cnfg_dir","%s",cnfg_dir);
         loc_dir[0]='\0';
      }

      find_section("Configurations");
      read_line("first","%d",&first);
      read_line("last","%d",&last);
      read_line("step","%d",&step);

      error_root((last<first)||(step<1)||(((last-first)%step)!=0),1,
                 "read_dirs [ms10.c]","Improper configuration range");
   }

   MPI_Bcast(nbase,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

   MPI_Bcast(log_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(dat_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(loc_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);
   MPI_Bcast(cnfg_dir,NAME_SIZE,MPI_CHAR,0,MPI_COMM_WORLD);

   MPI_Bcast(&first,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&last,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&step,1,MPI_INT,0,MPI_COMM_WORLD);
}


static void setup_files(void)
{
   if (noexp)
      error_root(name_size("%s/%sn%d_%d",loc_dir,nbase,last,NPROC-1)>=NAME_SIZE,
                 1,"setup_files [ms10.c]","loc_dir name is too long");
   else
      error_root(name_size("%s/%sn%d",cnfg_dir,nbase,last)>=NAME_SIZE,
                 1,"setup_files [ms10.c]","cnfg_dir name is too long");

   check_dir_root(log_dir);
   check_dir_root(dat_dir);
   error_root(name_size("%s/%s.ms10.log~",log_dir,nbase)>=NAME_SIZE,
              1,"setup_files [ms10.c]","log_dir name is too long");
   error_root(name_size("%s/%s.ms10.dat~",dat_dir,nbase)>=NAME_SIZE,
              1,"setup_files [ms10.c]","dat_dir name is too long");

   sprintf(log_file,"%s/%s.ms10.log",log_dir,nbase);
   sprintf(par_file,"%s/%s.ms10.par",dat_dir,nbase);
   sprintf(dat_file,"%s/%s.ms10.dat",dat_dir,nbase);
   sprintf(end_file,"%s/%s.ms10.end",log_dir,nbase);
   sprintf(log_save,"%s~",log_file);
   sprintf(par_save,"%s~",par_file);
   sprintf(dat_save,"%s~",dat_file);
}


static void read_bc_parms(void)
{
   int bc;
   double phi[2],phi_prime[2];

   if (my_rank==0)
   {
      find_section("Boundary conditions");
      read_line("type","%d",&bc);

      phi[0]=0.0;
      phi[1]=0.0;
      phi_prime[0]=0.0;
      phi_prime[1]=0.0;

      if (bc==1)
         read_dprms("phi",2,phi);

      if ((bc==1)||(bc==2))
         read_dprms("phi'",2,phi_prime);
   }

   MPI_Bcast(&bc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(phi,2,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Bcast(phi_prime,2,MPI_DOUBLE,0,MPI_COMM_WORLD);

   bcp=set_bc_parms(bc,1.0,1.0,1.0,1.0,phi,phi_prime);

   if (append)
      check_bc_parms(fdat);
   else
      write_bc_parms(fdat);
}


static void read_wl_parms(void)
{
   read_wloop_parms();
   wl=wloop_parms();

   file_head.msl=wl.msl;
   file_head.wls=wl.wls;
   file_head.tmax=N0;
   file_head.wlmr=wl.mwlr;
   file_head.wlmt=wl.mwlt;

   if (append)
      check_wloop_parms(fdat);
   else
      write_wloop_parms(fdat);
}


static void read_infile(int argc,char *argv[])
{
   int ifile;

   if (my_rank==0)
   {
      flog=freopen("STARTUP_ERROR","w",stdout);

      ifile=find_opt(argc,argv,"-i");
      endian=endianness();

      error_root((ifile==0)||(ifile==(argc-1)),1,"read_infile [ms10.c]",
                 "Syntax: ms10 -i <input file> [-noexp] [-a]");

      error_root(endian==UNKNOWN_ENDIAN,1,"read_infile [ms10.c]",
                 "Machine has unknown endianness");

      noexp=find_opt(argc,argv,"-noexp");
      append=find_opt(argc,argv,"-a");

      fin=freopen(argv[ifile+1],"r",stdin);
      error_root(fin==NULL,1,"read_infile [ms10.c]",
                 "Unable to open input file");
   }

   MPI_Bcast(&endian,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&noexp,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&append,1,MPI_INT,0,MPI_COMM_WORLD);

   read_dirs();
   setup_files();

   if (my_rank==0)
   {
      if (append)
         fdat=fopen(par_file,"rb");
      else
         fdat=fopen(par_file,"wb");

      error_root(fdat==NULL,1,"read_infile [ms10.c]",
                 "Unable to open parameter file");
   }

   read_bc_parms();
   read_wl_parms();

   if (my_rank==0)
   {
      fclose(fin);
      fclose(fdat);

      if (append==0)
         copy_file(par_file,par_save);
   }
}


static void check_old_log(int *fst,int *lst,int *stp)
{
   int ie,ic,isv;
   int fc,lc,dc,pc;
   int np[4],bp[4];

   fend=fopen(log_file,"r");
   error_root(fend==NULL,1,"check_old_log [ms10.c]",
              "Unable to open log file");

   fc=0;
   lc=0;
   dc=0;
   pc=0;

   ie=0x0;
   ic=0;
   isv=0;

   while (fgets(line,NAME_SIZE,fend)!=NULL)
   {
      if (strstr(line,"process grid")!=NULL)
      {
         if (sscanf(line,"%dx%dx%dx%d process grid, %dx%dx%dx%d",
                    np,np+1,np+2,np+3,bp,bp+1,bp+2,bp+3)==8)
         {
            ipgrd[0]=((np[0]!=NPROC0)||(np[1]!=NPROC1)||
                      (np[2]!=NPROC2)||(np[3]!=NPROC3));
            ipgrd[1]=((bp[0]!=NPROC0_BLK)||(bp[1]!=NPROC1_BLK)||
                      (bp[2]!=NPROC2_BLK)||(bp[3]!=NPROC3_BLK));
         }
         else
            ie|=0x1;
      }
      else if (strstr(line,"fully processed")!=NULL)
      {
         pc=lc;

         if (sscanf(line,"Configuration no %d",&lc)==1)
         {
            ic+=1;
            isv=1;
         }
         else
            ie|=0x1;

         if (ic==1)
            fc=lc;
         else if (ic==2)
            dc=lc-fc;
         else if ((ic>2)&&(lc!=(pc+dc)))
            ie|=0x2;
      }
      else if (strstr(line,"Configuration no")!=NULL)
         isv=0;
   }

   fclose(fend);

   error_root((ie&0x1)!=0x0,1,"check_old_log [ms10.c]",
              "Incorrect read count");
   error_root((ie&0x2)!=0x0,1,"check_old_log [ms10.c]",
              "Configuration numbers are not equally spaced");
   error_root(isv==0,1,"check_old_log [ms10.c]",
              "Log file extends beyond the last configuration save");

   (*fst)=fc;
   (*lst)=lc;
   (*stp)=dc;
}


static void check_old_dat(int fst,int lst,int stp)
{
   int ie,ic;
   int fc,lc,dc,pc;

   fdat=fopen(dat_file,"rb");
   error_root(fdat==NULL,1,"check_old_dat [ms10.c]",
              "Unable to open data file");

   check_file_head();

   fc=0;
   lc=0;
   dc=0;
   pc=0;

   ie=0x0;
   ic=0;

   while (read_data()==1)
   {
      pc=lc;
      lc=data.nc;
      ic+=1;

      if (ic==1)
         fc=lc;
      else if (ic==2)
         dc=lc-fc;
      else if ((ic>2)&&(lc!=(pc+dc)))
         ie|=0x1;
   }

   fclose(fdat);

   error_root(ic==0,1,"check_old_dat [ms10.c]",
              "No data records found");
   error_root((ie&0x1)!=0x0,1,"check_old_dat [ms10.c]",
              "Configuration numbers are not equally spaced");
   error_root((fst!=fc)||(lst!=lc)||(stp!=dc),1,"check_old_dat [ms10.c]",
              "Configuration range is not as reported in the log file");
}


static void check_files(void)
{
   int fst,lst,stp;

   ipgrd[0]=0;
   ipgrd[1]=0;

   if (my_rank==0)
   {
      if (append)
      {
         check_old_log(&fst,&lst,&stp);
         check_old_dat(fst,lst,stp);

         error_root((fst!=lst)&&(stp!=step),1,"check_files [ms10.c]",
                    "Continuation run:\n"
                    "Previous run had a different configuration separation");
         error_root(first!=lst+step,1,"check_files [ms10.c]",
                    "Continuation run:\n"
                    "Configuration range does not continue the previous one");
      }
      else
      {
         fin=fopen(log_file,"r");
         fdat=fopen(dat_file,"rb");

         error_root((fin!=NULL)||(fdat!=NULL),1,"check_files [ms10.c]",
                    "Attempt to overwrite old *.log or *.dat file");

         fdat=fopen(dat_file,"wb");
         error_root(fdat==NULL,1,"check_files [ms10.c]",
                    "Unable to open data file");
         write_file_head();
         fclose(fdat);
      }
   }
}


static void print_info(void)
{
   int n[3];
   long ip;

   if (my_rank==0)
   {
      ip=ftell(flog);
      fclose(flog);

      if (ip==0L)
         remove("STARTUP_ERROR");

      if (append)
         flog=freopen(log_file,"a",stdout);
      else
         flog=freopen(log_file,"w",stdout);

      error_root(flog==NULL,1,"print_info [ms10.c]","Unable to open log file");
      printf("\n");

      if (append)
         printf("Continuation run\n\n");
      else
      {
         printf("Computation of Wilson loops\n");
         printf("--------------------------------------\n\n");
      }

      printf("Wloop version %s\n",wloop_RELEASE);
      printf("openQCD version %s\n",openQCD_RELEASE);

      if (endian==LITTLE_ENDIAN)
         printf("The machine is little endian\n");
      else
         printf("The machine is big endian\n");
      if (noexp)
         printf("Configurations are read in imported file format\n\n");
      else
         printf("Configurations are read in exported file format\n\n");

      if ((ipgrd[0]!=0)&&(ipgrd[1]!=0))
         printf("Process grid and process block size changed:\n");
      else if (ipgrd[0]!=0)
         printf("Process grid changed:\n");
      else if (ipgrd[1]!=0)
         printf("Process block size changed:\n");

      if ((append==0)||(ipgrd[0]!=0)||(ipgrd[1]!=0))
      {
         printf("%dx%dx%dx%d lattice, ",N0,N1,N2,N3);
         printf("%dx%dx%dx%d local lattice\n",L0,L1,L2,L3);
         printf("%dx%dx%dx%d process grid, ",NPROC0,NPROC1,NPROC2,NPROC3);
         printf("%dx%dx%dx%d process block size\n\n",
                NPROC0_BLK,NPROC1_BLK,NPROC2_BLK,NPROC3_BLK);
      }

      if (append==0)
      {
         if (bcp.type==0)
            printf("Open boundary conditions\n\n");
         else if (bcp.type==1)
         {
            printf("SF boundary conditions\n");

            n[0]=fdigits(bcp.phi[0][0]);
            n[1]=fdigits(bcp.phi[0][1]);
            n[2]=fdigits(bcp.phi[0][2]);
            printf("phi = %.*f,%.*f,%.*f\n",IMAX(n[0],1),bcp.phi[0][0],
                   IMAX(n[1],1),bcp.phi[0][1],IMAX(n[2],1),bcp.phi[0][2]);

            n[0]=fdigits(bcp.phi[1][0]);
            n[1]=fdigits(bcp.phi[1][1]);
            n[2]=fdigits(bcp.phi[1][2]);
            printf("phi' = %.*f,%.*f,%.*f\n\n",IMAX(n[0],1),bcp.phi[1][0],
                   IMAX(n[1],1),bcp.phi[1][1],IMAX(n[2],1),bcp.phi[1][2]);
         }
         else if (bcp.type==2)
         {
            printf("Open-SF boundary conditions\n");

            n[0]=fdigits(bcp.phi[1][0]);
            n[1]=fdigits(bcp.phi[1][1]);
            n[2]=fdigits(bcp.phi[1][2]);
            printf("phi' = %.*f,%.*f,%.*f\n\n",IMAX(n[0],1),bcp.phi[1][0],
                   IMAX(n[1],1),bcp.phi[1][1],IMAX(n[2],1),bcp.phi[1][2]);
         }
         else
            printf("Periodic boundary conditions\n\n");
      }
   }
   
   if (append==0)
      print_wloop_parms();

   if (my_rank==0)
   {
      printf("Configurations no %d -> %d in steps of %d\n\n",
             first,last,step);
      fflush(flog);
   }
}


static void set_data(int nc)
{
   int i,wls;

   data.nc=nc;
   wls=file_head.wls;

   hyp_links(wl.alpha_action[0],wl.alpha_action[1],wl.alpha_action[2],
                  wl.proj,wl.proj_iter,0);
   wloop_sum(1);
   
   for (i=0;i<wls;i++)
      data.wl_mat[i]=wl.wl_mat[i];
}


static void save_data(void)
{
   if (my_rank==0)
   {
      fdat=fopen(dat_file,"ab");
      error_root(fdat==NULL,1,"save_data [ms10.c]",
                 "Unable to open data file");
      write_data();
      fclose(fdat);
   }
}

/*
static void print_log(void)
{
   int in,dn,nn,din;
   double eps;

   if (my_rank==0)
   {
      dn=file_head.dn;
      nn=file_head.nn;
      eps=file_head.eps;

      din=nn/10;
      if (din<1)
         din=1;

      for (in=0;in<=nn;in+=din)
         printf("n = %3d, t = %.2e, Wact = %.6e, Yact = %.6e, Q = % .2e\n",
                in*dn,eps*(double)(in*dn),Wact[in],Yact[in],Qtop[in]);
   }
}
*/

static void check_endflag(int *iend)
{
   if (my_rank==0)
   {
      fend=fopen(end_file,"r");

      if (fend!=NULL)
      {
         fclose(fend);
         remove(end_file);
         (*iend)=1;
         printf("End flag set, run stopped\n\n");
      }
      else
         (*iend)=0;
   }

   MPI_Bcast(iend,1,MPI_INT,0,MPI_COMM_WORLD);
}


int main(int argc,char *argv[])
{
   int nc,iend;
   double wt1,wt2,wtavg;

   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   read_infile(argc,argv);
   alloc_data();
   check_files();
   print_info();

   geometry();
   if (flint)
      alloc_wfd(1);

   iend=0;
   wtavg=0.0;

   for (nc=first;(iend==0)&&(nc<=last);nc+=step)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      wt1=MPI_Wtime();

      if (my_rank==0)
         printf("Configuration no %d\n",nc);

      if (noexp)
      {
         sprintf(cnfg_file,"%s/%sn%d_%d",loc_dir,nbase,nc,my_rank);
         read_cnfg(cnfg_file);
      }
      else
      {
         sprintf(cnfg_file,"%s/%sn%d",cnfg_dir,nbase,nc);
         import_cnfg(cnfg_file);
      }

      set_data(nc);
      save_data();
      /*print_log();*/

      MPI_Barrier(MPI_COMM_WORLD);
      wt2=MPI_Wtime();
      wtavg+=(wt2-wt1);
      error_chk();

      if (my_rank==0)
      {
         printf("Configuration no %d fully processed in %.2e sec ",
                nc,wt2-wt1);
         printf("(average = %.2e sec)\n\n",
                wtavg/(double)((nc-first)/step+1));
         fflush(flog);

         copy_file(log_file,log_save);
         copy_file(dat_file,dat_save);
      }

      check_endflag(&iend);
   }

   error_chk();

   if (my_rank==0)
   {
      fflush(flog);
      copy_file(log_file,log_save);
      fclose(flog);
   }

   MPI_Finalize();
   exit(0);
}
