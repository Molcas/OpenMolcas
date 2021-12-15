/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2012, Victor P. Vysotskiy                              *
***********************************************************************/
/******************************************************************************/
/*                   Molcas Memory Allocator: Sys V Shared memory             */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*    written by:                                                             */
/*    V. P. Vysotskiy                                                         */
/*    University of Lund, Sweden, 2012                                        */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*    history: Initial revision                                               */
/*                                                                            */
/******************************************************************************/
#ifdef  _CYGWIN_
#define _DARWIN_
#endif

#include <sys/mman.h>
#include <stddef.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <sys/ipc.h>
#include <sys/shm.h>

#include "molcastype.h"

#include "mma.h"

#ifdef _CAPITALS_
#define c_getshmem C_GETSHMEM
#else
#ifndef ADD_
#define c_getshmem c_getshmem_
#endif
#endif

#define MXLINE         9
#define MAXREC     32768   /* 8192 */
#define UNDEF_MEM      0

char* getenvc(const char*);

typedef struct mentry mentry;
struct  mentry{
    char  elbl[MXLINE];       /* label                                      */
    char  etyp[MXLINE];       /* datatype                                   */
    INT   offset;             /* offset relative to Fortran  xWrkSpc arrays */
    INT   len;                /* len = numel*dsize(datatype)                */
    INT   atime;              /* an unique ID based on the creation time    */
    void *addr;               /* memory adress returned by malloc           */
};

typedef struct  mstat  mstat;
struct  mstat{
    INT   nmentry;             /* Total number of entries                    */
    INT   naccess;             /* Total number of allocs                     */
    INT   mxmem;               /* MOLCASMEM                                  */
    INT   avmem;               /* maximal available memory for MOLCAS        */
};

extern mstat   MlM;
extern double *dptr;
extern float  *sptr;
extern INT    *iptr;
extern char   *cptr;


INT dsize(char datatype[]);
INT memop(char *op);
void string2UC(char *src, char *dest);
void dump_mentry(char *tag, mentry *curr);
void set_mentry(mentry *tmp, char *lbl, char *typ, INT *offset, INT *len);
void print_params(char *func, char *name, char* Op, char *dtyp, INT *offset, INT *len);
INT c_getmem_kern(INT *op, mentry *tmp, INT *offset, INT *len);
char *woff2cptr(char typ[], INT offset);

/*-----------------------------------------------------------------------------*/
INT add_shmentry(mstat *MM, mentry *tmp, char *path, INT *Id) {
    int   shmid;
    key_t SHM_key;
    char *defpath=NULL,*wpath;
#ifdef _DEBUGPRINT_MEM_
    int    rc;
    struct shmid_ds shm_stat;
#endif

    ptrdiff_t dist=UNDEF_MEM;
    char   *wrkspc=NULL;

    if(defpath==NULL) defpath=getenvc("MOLCAS");
#ifdef _DEBUGPRINT_MEM_
    printf("defpath=%s\n",defpath);
#endif
    wpath=(strlen(path))?path:defpath;
    SHM_key=ftok(wpath,(int) (*Id));
/*
    if(strlen(path)==0) {
        SHM_key=ftok(defpath,(int) (*Id));
    } else {
        SHM_key=ftok(path,(int) (*Id));
    }
*/
/* Linux part is here */
#ifdef _DEBUGPRINT_MEM_
    printf("++++++++ Adding new shared memory entry %s of type = %s with length=%ld and the key=%ld\n",tmp->elbl,tmp->etyp,LIFMT(tmp->len),LIFMT(SHM_key));
#endif
#ifdef _HUGE_PAGES_
#ifdef _AIX_
    shmid = shmget(SHM_key, (size_t) tmp->len, 0644|IPC_CREAT|IPC_EXCL|SHM_LGPAGE|SHM_PIN);
#else
    shmid = shmget(SHM_key, (size_t) tmp->len, 0644|SHM_HUGETLB|IPC_CREAT|IPC_EXCL|SHM_LOCKED|SHM_NORESERVE);
#endif
#else
#ifdef _DARWIN_
      shmid = shmget(SHM_key, (size_t) tmp->len, 0644|IPC_CREAT|IPC_EXCL);
#elif _AIX_
      shmid = shmget(SHM_key, (size_t) tmp->len, 0644|IPC_CREAT|IPC_EXCL|SHM_PIN);
#else
      shmid = shmget(SHM_key, (size_t) tmp->len, 0644|IPC_CREAT|IPC_EXCL|SHM_LOCKED);
/*    shmid = shmget(SHM_key, (size_t) tmp->len, 0644|IPC_CREAT|IPC_EXCL|SHM_LOCKED|SHM_NORESERVE);*/
#endif
#endif

/* cat /proc/sysvipc/shm ?*/

    if(shmid==-1) {
      shmid = shmget(SHM_key, 0,0644|IPC_CREAT);
#ifdef _DEBUGPRINT_MEM_
      rc=shmctl(shmid,IPC_STAT,&shm_stat);
      if(rc==0) {
         printf("Shared memory segment with the key '%x' is already exist and its size = %ld bytes\n",SHM_key,LIFMT(shm_stat.shm_segsz));
      } else {
         printf("Something is terrible wrong, please contact author!\n");
         perror("shmget failed to get id");
         return(-2);
      }
#endif
      tmp->len=0;
    } else {
/* Memory pinning
      rc=shmctl(shmid,SHM_LOCK,&shm_stat);
      if(rc==-1) perror("Failed to pinn a shared memory");
*/
      MM->avmem-=tmp->len;
    }

   *Id=(INT) shmid;

    tmp->atime=(INT) (-shmid);

    wrkspc = shmat(shmid, NULL, 0);
    if (wrkspc == (char *)-1) {
        perror("shmop: shmat failed");
        return(-2);
    }
#ifdef _DEBUGPRINT_MEM_
    bzero((void *) wrkspc, (size_t) tmp->len);
#endif
    tmp->addr=(void *) wrkspc;

    switch (tmp->etyp[0]) {
        case 'R':
            dist=(double *) wrkspc-dptr;
            break;
        case 'S':
            dist=(float *)  wrkspc-sptr;
            break;
        case 'I':
            dist=(INT *)    wrkspc-iptr;
            break;
        case 'C':
            dist=           wrkspc-cptr;
            break;
        default: printf("MMA: not supported datatype %s\n",tmp->etyp);
    }
#ifdef _DEBUGPRINT_MEM_
    printf("dist=%ld\n",LIFMT(dist));
    printf("entry %s has been allocated at %p. The corresponding <key> and [id] are <%x> and [%ld], respectively\n",tmp->elbl,tmp->addr,SHM_key,LIFMT(*Id));
#endif
    tmp->offset=(INT) dist;
#ifdef _DEBUGPRINT_MEM_
    dump_mentry("NEW SHARED MEMORY ENTRY", tmp);
#endif

    return(tmp->offset);
}

/*-----------------------------------------------------------------------------*/
INT del_shmentry(mentry *tmp, INT shmid) {
    INT     rc;
    char   *wrkspc=NULL;

#ifdef _DEBUGPRINT_MEM_
    printf("-------- Deleting shared memory entry '%9s' of type = '%9s' with length=%12ld (offset=%12ld)\n",tmp->elbl,tmp->etyp,LIFMT(tmp->len),LIFMT(tmp->offset));
#endif

#ifdef _DEBUGPRINT_MEM_
    dump_mentry("DELETING SHARED MEMORY ENTRY",tmp);
    printf("requested shmid=%ld\n",shmid);
#endif

    wrkspc=woff2cptr(tmp->etyp, tmp->offset);

#ifdef _DEBUGPRINT_MEM_
    printf("Deallocating memory %s at adress %p\n",tmp->elbl, wrkspc);
#endif
    rc=shmdt(wrkspc);
    if(rc==-1) {
       perror("Cannot detache the  shared memory segment");
       return(-2);
    }
    shmctl(shmid, IPC_RMID, 0);

    return(0);
}


INT c_getshmem(char *name, char* Op, char *dtyp, INT *offset, INT *len, char *path, INT *SHM_Id) {

    enum memops {ALLO,FREE,LENG,CHAN,CHEC,MAX,LIST,TERM,FLUS,PINN,RGST,EXCL};

    mentry  tmp;
    INT     op,rc=1,blen;

    char e_name[MXLINE];
    char e_dtyp[MXLINE];
    char action[MXLINE];

    string2UC(name,e_name);strcat(e_name, "_s");
    string2UC(dtyp,e_dtyp);
    string2UC(Op,action);

    op=memop(action);
#ifdef _DEBUGPRINT_MEM_
    print_params("C_Get_SHMEM",name, Op, dtyp, offset,len);
    printf("Op=%s (%ld)\n",Op,LIFMT(memop(Op)));
#endif
    blen=(e_dtyp[0]=='C')?1:0;
    blen+=dsize(e_dtyp)*(*len);

    bzero(&tmp,sizeof(mentry));
    set_mentry(&tmp, e_name, e_dtyp, offset, &blen);

    switch (op) {
        case ALLO:
/* Checking for memory leaks     */
            if(MlM.nmentry==MAXREC) {
               printf("MEMORY ERROR: Possible memory leak detected: The number of memory blocks exceeds the limit of %d entries\n",MAXREC);
               exit(-3);
            }
/* Checking for free memory available */
            if(MlM.avmem<blen) {
              if(MlM.mxmem<blen) {
                printf("MEMORY ERROR: Memory is exhausted!\n");
                exit(-3);
              } else {
                printf("MEMORY WARNING: MOLCAS_MEM has been increased by MOLCAS_MAXMEM (%ld) !\n",LIFMT(MlM.mxmem));
                MlM.avmem+=MlM.mxmem;
                MlM.mxmem=0;
              }
            }
            op=RGST;
           *offset=add_shmentry(&MlM, &tmp, path, SHM_Id);
            rc=c_getmem_kern(&op, &tmp, offset, len);
            break;
        case FREE:
            if(MlM.nmentry==0) {
              printf("WARNING: Attempt to operate on zero allocated memory blocks\n");
              exit(-3);
            }
            rc=del_shmentry(&tmp, *SHM_Id);
            op=EXCL;
            rc=c_getmem_kern(&op, &tmp, offset, len);
            if(rc<0) print_params("C_Get_SHMEM",name, Op, dtyp, offset,len);
            break;

        default: rc=-1;printf("Unsupported shared memory operation %s\n",Op);
    }
    return(rc);
}
