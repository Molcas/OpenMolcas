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
* Copyright (C) 2012,2013, Victor P. Vysotskiy                         *
***********************************************************************/
/******************************************************************************/
/*                         Molcas Memory Allocator                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*    written by:                                                             */
/*    V. P. Vysotskiy                                                         */
/*    University of Lund, Sweden, 2012-2013                                   */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*                                                                            */
/*    history: Initial revision                                               */
/*                                                                            */
/******************************************************************************/
#ifdef  _DEBUG_MEM_
#define _MEMORY_TRACE_
#endif
#ifdef  _CYGWIN_
#define _DARWIN_
#undef  _MEMORY_TRACE_
#endif
#define _POSIX_C_SOURCE 200112L
#include <sys/mman.h>
#include <stddef.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "molcastype.h"
#ifdef _MEMORY_TRACE_
#include <mcheck.h>
#endif
#ifdef _GARBLE_
#include <limits.h>
#include <float.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#include "mma.h"

#define SFCTR(X)       (X-X/10)
#define MXLINE         9
#define PGSIZE      4096
#define MB       1000000
#define EMPTYE     -1000
#define PINNED_MEM   999
#define UNDEF_MEM      0
#define SYS_ATIME      0
#define ALLOC_FLD     -2
#define MINMEMPTR -577777000306848070
#ifdef _DEBUG_MEM_
#define MAXREC    32768
#else
#define MAXREC    32768 /* 8192 */
#endif

#if defined (_MEM_PROF_) || defined (_GARBLE_)

#define FCNST 133679
#define MEM_THRSHLD 1

#ifdef _CAPITALS_
#define dcopy DCOPY
#define scopy SCOPY
#define icopy ICOPY
#else
#ifndef ADD_
#define dcopy dcopy_
#define scopy scopy_
#define icopy icopy_
#endif
#endif

void dcopy(INT *N, double *X, INT *INCX, double *Y, INT *INCY);
void scopy(INT *N, float  *X, INT *INCX, float  *Y, INT *INCY);
void icopy(INT *N, INT    *X, INT *INCX, INT    *Y, INT *INCY);
#endif

char* getenvc(const char*);
INT dsize(char datatype[]);

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
    INT   totmem;              /* initial MOLCAS_MEM                         */
};

mstat   MlM={0,SYS_ATIME+1,0,0,0};

double *dptr;
float  *sptr;
INT    *iptr;
char   *cptr;

#ifdef _OPENMP
omp_lock_t mma_lock;
#endif
#ifdef _TRACK_
INT   trc_mentry(mstat *MM, mentry mentries[], mentry *tmp) {
      INT       i,j,len,nentry;
      char     *p,*base;
      ptrdiff_t match=-1,dist=0;

      p=tmp->addr;

      for(nentry=MM->nmentry,j=0;j<nentry;j++) {
          for(base=mentries[j].addr,len=mentries[j].len,i=0;i<len;i++)if((match=(&base[i]-p))==0)  break;
          if(match==0) break;
      }

      if(match==0) {
        printf("Tracked %p address belongs to the '%s' memory entry %p of %ld bytes\n",tmp->addr,mentries[j].elbl,mentries[j].addr,LIFMT(mentries[j].len));
        dist=(mentries[j].addr+mentries[j].len)-tmp->addr;
      } else {
        printf("Could not track the %p address\n",tmp->addr);
      }

      return(dist);
}
#endif
#ifdef _MEM_PROF_
void anlz_mem(mentry *tgt) {
    INT    n,i,nusd=0;
    INT    flIcnst=(INT)    FCNST, *iw;
    double flRcnst=(double) FCNST, *rw;
    float  flScnst=(float)  FCNST, *sw;
    float  mfrc=100;

    if(tgt->len<MEM_THRSHLD) return;
    n=tgt->len/dsize(tgt->etyp);

    switch (tgt->etyp[0]) {
        case 'R':
            rw=(double *) tgt->addr;
            for(i=0;i<n;i++) if(rw[i]==flRcnst) nusd++;
            break;
        case 'S':
            sw=(float  *) tgt->addr;
            for(i=0;i<n;i++) if(sw[i]==flScnst) nusd++;
            break;
        case 'I':
            iw=(INT    *) tgt->addr;
            for(i=0;i<n;i++) if(iw[i]==flIcnst) nusd++;
            break;
        default: return;
    }
    mfrc*=(float) nusd; mfrc/=n;
    printf("Fraction of memory usage for entry '%s' (offset=%ld) = %4.3f \%\n",tgt->elbl, LIFMT(tgt->offset),mfrc);
    }
#endif


INT testmem(INT *MOLCASMEM) {
    INT     rc;
    char   *twrkspc=NULL;

    twrkspc=malloc(sizeof(char)*(*MOLCASMEM));

    if(twrkspc) {
        rc=1;
        free(twrkspc);
    } else {
        rc=-1;
#ifdef _DEBUG_MEM_
        printf("MOLCAS cannot get %ld bytes of memory !\n",LIFMT(*MOLCASMEM));
#endif
    }

    return(rc);
}


INT allocmem(double ref[],char cref[],INT *intof,INT *dblof,INT *sglof, INT *chrof,INT *size) {
    INT     rc,MOLCASMEM,MAXMEM;
    char c;
    char *ptr;
    int factor=1;
#ifdef _DEMO_
    char    memsize[]="128";
#else
    char   *memsize=NULL;
#endif
/* I. THE 'MOLCASMEM' PART */
#ifndef _DEMO_
    memsize=getenvc("MOLCAS_MEM");
#endif
    if(memsize==NULL) {
         printf("MOLCAS_MEM is not defined!\n");
         return(-1);
    } else {
         ptr=strstr(memsize,"b");
         if(ptr==NULL) ptr=strstr(memsize,"B");
         if(ptr!=NULL)
           {
             c=*(ptr-1);
             switch (c)
              {
               case 'G':
               case 'g': factor=1024; *(ptr-1)=0; break;
               case 'M':
               case 'm': *(ptr-1)=0; break;
               case 'T':
               case 't': factor=1024*1024;*(ptr-1)=0; break;
               default : printf("Unknown units for MOLCAS_MEM\n"); break;
              }
           }
         MOLCASMEM=strtol(memsize,NULL,10);
         MOLCASMEM*=MB*factor;
         factor=1;
#ifdef _MEMORY_TRACE_
         mcheck(NULL);
         mcheck_pedantic(NULL);
         mtrace();
#endif
#ifdef _DEBUG_MEM_
         printf("MOLCAS_MEM=%ld byte\n",LIFMT(MOLCASMEM));
#endif
    }

    /* when malloc is intercepted, we should not test if the memory
    * can be allocated (testmem) as that might result in very expensive
    * memset calls from various memory inspection tools */
#ifndef _MALLOC_INTERCEPT_
     rc=(testmem(&MOLCASMEM))?0:-1;
#else
     rc=0;
#endif

    *size=MOLCASMEM/sizeof(double);

    *dblof=*sglof=*intof=*chrof=1;

     dptr=          ref;
     sptr=(float *) ref;
     iptr=(INT *)   ref;
     cptr=(char *) cref;

     MlM.avmem=MOLCASMEM;
     MlM.totmem=MOLCASMEM;
#ifndef _DEMO_
     free(memsize);
#endif
/* II. THE 'MOLCAS_MAXMEM' PART */
#ifndef _DEMO_
     memsize=getenvc("MOLCAS_MAXMEM");
#endif
     if(memsize==NULL) {
#ifdef _DEBUG_MEM_
         printf("WARNING: MOLCAS_MAXMEM is not defined!\n");
#endif
     } else {
         ptr=strstr(memsize,"b");
         if(ptr==NULL) ptr=strstr(memsize,"B");
         if(ptr!=NULL)
           {
             c=*(ptr-1);
             switch (c)
              {
               case 'G':
               case 'g': factor=1024; *(ptr-1)=0; break;
               case 'M':
               case 'm': factor=1; *(ptr-1)=0; break;
               case 'T':
               case 't': factor=1024*1024;*(ptr-1)=0; break;
               default : printf("Unknown units for MOLCAS_MEM\n"); break;
              }
           }
         MAXMEM=strtol(memsize,NULL,10);
         MAXMEM*=MB*factor;
#ifdef _DEBUG_MEM_
         printf("MOLCAS_MAXMEM=%ld byte\n",LIFMT(MAXMEM));
#endif
         MlM.mxmem=MAXMEM-MOLCASMEM;
         if(MlM.mxmem<0) {
             printf("WARNING: MOLCAS_MAXMEM (%ld) < MOLCAS_MEM (%ld)\n",LIFMT(MAXMEM),LIFMT(MOLCASMEM));
             MlM.mxmem=0;
         }
#ifndef _DEMO_
         free(memsize);
#endif
     }
#ifdef _DEBUG_MEM_
     printf("ref=%p\n",ref);
     printf("cref=%p\n",cref);
     setvbuf(stdout, NULL, _IOLBF,0);
#endif
#ifdef _OPENMP
     omp_init_lock(&mma_lock);
#endif
     return(rc);
}


/*-----------------------------------------------------------------------------*/
INT dsize(char datatype[]) {
    INT bsize=-1;
    switch (datatype[0]) {
        case 'R':
            bsize=sizeof(double);
            break;
        case 'S':
            bsize=sizeof(float);
            break;
        case 'I':
            bsize=sizeof(INT);
            break;
        case 'C':
            bsize=sizeof(char);
            break;
        default: printf("MMA: not supported datatype '%s'\n",datatype);
#ifdef _DEBUG_MEM_
                 abort();
#endif
    }
    return(bsize);
}

/*-----------------------------------------------------------------------------*/
void l2u(char *str) {

    INT i,n=0;

    if(str!=NULL) n=strlen(str);

    for(i=0; i<n; i++) {
        str[i]=toupper(str[i]);
        if(str[i]==32) break;
    }
    str[i]='\0';

    return;
}
/*-----------------------------------------------------------------------------*/
void string2UC(char *src, char *dest) {
    strcpy(dest,src);l2u(dest);
    return;
}

INT memop(char *op) {
#ifdef _TRACK_
    enum memops {ALLO,FREE,LENG,CHEC,MAX,LIST,TERM,FLUS,PINN,RGST,EXCL,TRCK};
#else
    enum memops {ALLO,FREE,LENG,CHEC,MAX,LIST,TERM,FLUS,PINN,RGST,EXCL};
#endif
    if(strstr(op,"ALLO")) return(ALLO);
    if(strstr(op,"FREE")) return(FREE);
    if(strstr(op,"LENG")) return(LENG);
    if(strstr(op,"CHEC")) return(CHEC);
    if(strstr(op,"MAX"))  return(MAX);
    if(strstr(op,"LIST")) return(LIST);
    if(strstr(op,"TERM")) return(TERM);
    if(strstr(op,"FLUS")) return(FLUS);
    if(strstr(op,"PINN")) return(PINN);
    if(strstr(op,"RGST")) return(RGST);
    if(strstr(op,"EXCL")) return(EXCL);
#ifdef _TRACK_
    if(strstr(op,"TRCK")) return(TRCK);
#endif
    return(-1);
}

/*-----------------------------------------------------------------------------*/
void dump_mentry(char *tag, mentry *curr) {
    if(curr) {
        printf("MA_DUMP_INFO < %s > name=%s, datatype=%s, offset=%ld (adress=%p), len=%ld\n",tag,curr->elbl,curr->etyp,LIFMT(curr->offset),curr->addr,LIFMT(curr->len));
    }else{
        printf("MA_DUMP_INFO < %s >  EMPTY RECORD!\n",tag);
    }
    return;
}
/*-----------------------------------------------------------------------------*/
INT find_mentry(mentry mentries[], mentry *tmp) {
    INT i;
#ifdef _DEBUG_MEM_
    mentry *tgt=NULL;
#endif
    for (i=0;i<MAXREC;i++) if(mentries[i].offset==tmp->offset) break;
#ifdef _DEBUG_MEM_
    printf("#i=%ld\n",LIFMT(i));
    if(mentries[i].offset==tmp->offset)  tgt=&mentries[i];
    dump_mentry("TARGET",tgt);
#endif

#ifdef _DEBUG_MEM_
    if(tgt) {
        if(strcmp(tgt->elbl,tmp->elbl)!=0) printf("WARNING: Data  labels    are not matching: %s (stored) and %s (requested)\n", tgt->elbl, tmp->elbl);
        if(strcmp(tgt->etyp,tmp->etyp)!=0) printf("WARNING: Data  types     are not matching: %s (stored) and %s (requested)\n", tgt->etyp, tmp->etyp);
        if(       tgt->len!=tmp->len     ) printf("WARNING: Entry lengths   are not matching: %ld (stored) and %ld (requested)\n", LIFMT(tgt->len),  LIFMT(tmp->len));
    }
#endif
    return(i);
}
/*-----------------------------------------------------------------------------*/
INT ismax_mentry(INT i) {
    if(i==MAXREC) {
#ifdef _DEBUG_MEM_
        abort();
#else
        return(1);
#endif
    } else {
        return 0;
    }
}
/*-----------------------------------------------------------------------------*/
char *woff2cptr(char etyp[], INT offset) {

    char   *wrkspc=NULL;

    switch (etyp[0]) {
        case 'R':
            wrkspc=(char *) &dptr[offset];
            break;
        case 'S':
            wrkspc=(char *) &sptr[offset];
            break;
        case 'I':
            wrkspc=(char *) &iptr[offset];
            break;
        case 'C':
            wrkspc=         &cptr[offset];
            break;
        default: printf("MMA: not supported datatype %s\n",etyp);
    }
    return(wrkspc);
}

/*-----------------------------------------------------------------------------*/
INT cptr2woff(char etyp[], void *c_ptr) {

    ptrdiff_t dist=UNDEF_MEM;

    switch (etyp[0]) {
        case 'R':
            dist=(double *) c_ptr-dptr;
            break;
        case 'S':
            dist=(float *)  c_ptr-sptr;
            break;
        case 'I':
            dist=(INT *)    c_ptr-iptr;
            break;
        case 'C':
            dist=(char *)   c_ptr-cptr;
            break;
        default: printf("MMA: not supported datatype %s\n",etyp);
    }
    return(dist);
}

/*-----------------------------------------------------------------------------*/
INT mma_avmem(void) {
        return MlM.avmem;
}

/*-----------------------------------------------------------------------------*/
INT add_mentry(mstat *MM, mentry mentries[], mentry *tmp) {
    ptrdiff_t dist=UNDEF_MEM;
    char   *wrkspc=NULL;
    mentry *newe;
#ifdef _MEM_PROF_
    INT    n,incx=0,incy=1;
    INT    flIcnst=(INT)    FCNST;
    double flRcnst=(double) FCNST;
    float  flScnst=(float)  FCNST;
#elif defined (_GARBLE_)
    INT    n,incx=0,incy=1;
#ifdef _I8_
    INT    flIcnst=(INT) LONG_MAX;
#else
    INT    flIcnst=(INT) INT_MAX;
#endif
    double flRcnst= DBL_MAX;
    float  flScnst= FLT_MAX;
#endif
#ifdef _DARWIN_
#else
    int    rc;
#endif

#ifdef _DEBUG_MEM_
    printf("++++++++ Adding new entry %s of type = %s with length=%ld\n",tmp->elbl,tmp->etyp,LIFMT(tmp->len));
#endif

    newe=&mentries[MM->nmentry++];

    memcpy(newe,tmp,sizeof(mentry));
    MM->naccess++;

    if(newe->atime) newe->atime=MM->naccess;

    if(tmp->len==0) {
#ifdef _DEBUG_MEM_
       printf("Request for a zero allocation detected!\n");
#endif
       newe->offset=MINMEMPTR+MM->naccess;
       return(newe->offset);
    }

    if(tmp->offset==UNDEF_MEM) {
       wrkspc=(char *) malloc(tmp->len);
    } else {
#ifdef _DARWIN_
       wrkspc=(char *) malloc((size_t) tmp->len);
#else
       rc=posix_memalign((void **) &wrkspc, (size_t) sysconf(_SC_PAGESIZE), (size_t) tmp->len);
       (void)rc;
#endif
#ifdef _DEBUG_MEM_
       if(mlock((void *) wrkspc, (size_t) tmp->len)) dump_mentry("Cannot lock memory:",tmp);
#else
       mlock((void *) wrkspc, (size_t) tmp->len);
#endif
    }
    if(!wrkspc&&tmp->len) return(ALLOC_FLD); /* This is an exceptional case, i.e. wrkspc==NULL */

    MM->avmem-=tmp->len;

#ifdef _DEBUG_MEM_
    bzero((void *) wrkspc, (size_t) tmp->len);
#endif
    newe->addr=(void *) wrkspc;

#if defined (_MEM_PROF_) || defined (_GARBLE_)
    n=tmp->len/dsize(tmp->etyp);
#endif
    switch (tmp->etyp[0]) {
        case 'R':
            dist=(double *) wrkspc-dptr;
#if defined (_MEM_PROF_) || defined (_GARBLE_)
            dcopy(&n,&flRcnst,&incx,(double *) wrkspc,&incy);
#endif
            break;
        case 'S':
            dist=(float *)  wrkspc-sptr;
#if defined (_MEM_PROF_) || defined (_GARBLE_)
            scopy(&n,&flScnst,&incx,(float *) wrkspc,&incy);
#endif
            break;
        case 'I':
            dist=(INT *)    wrkspc-iptr;
#if defined (_MEM_PROF_) || defined (_GARBLE_)
            icopy(&n,&flIcnst,&incx,(INT *) wrkspc,&incy);
#endif
            break;
        case 'C':
            dist=           wrkspc-cptr;
            break;
        default: printf("MMA: not supported datatype %s\n",tmp->etyp);
    }
#ifdef _DEBUG_MEM_
    printf("dist=%ld\n",LIFMT(dist));
    printf("entry %s has been allocated at %p\n",newe->elbl,newe->addr);
#endif
    newe->offset=(INT) dist;
#ifdef _DEBUG_MEM_
    dump_mentry("NEW ENTRY", newe);
#endif
    return(newe->offset);
}

/*-----------------------------------------------------------------------------*/
INT reg_mentry(mstat *MM, mentry mentries[], mentry *tmp) {
    mentry *newe;
#ifdef _DEBUG_MEM_
    printf("++++++++ Registering new entry %s of type = %s with length=%ld\n",tmp->elbl,tmp->etyp,LIFMT(tmp->len));
#endif

    newe=&mentries[MM->nmentry++];
    memcpy(newe,tmp,sizeof(mentry));
    MM->naccess++;

    if(MlM.mxmem<tmp->len) {
       MlM.avmem-=tmp->len;
    } else {
       MM->mxmem-=tmp->len;
    }
    newe->addr=woff2cptr(tmp->etyp, tmp->offset);
    newe->atime=MM->naccess;

#ifdef _DEBUG_MEM_
    dump_mentry("NEW ENTRY HAS BEEN REGISTERED", newe);
#endif
    return(newe->atime);
}

/*-----------------------------------------------------------------------------*/
INT del_mentry(mstat *MM, mentry mentries[], mentry *tmp, INT i) {

    char   *wrkspc=NULL;
    mentry *tgt,*lst;

#ifdef _DEBUG_MEM_
    printf("-------- Deleting entry '%9s' of type = '%9s' with length=%12ld (offset=%12ld)\n",tmp->elbl,tmp->etyp,LIFMT(tmp->len),LIFMT(tmp->offset));
#endif

    if(i==0) i=find_mentry(mentries,tmp);

    if(ismax_mentry(i))
#ifdef _DEBUG_MEM_
    {
       abort();
    }
#else
    {
       return(-1);
    }
#endif
    tgt=&mentries[i];
    lst=&mentries[--MM->nmentry];

#ifdef _DEBUG_MEM_
    dump_mentry("DELETING ENTRY",tgt);
#endif

#ifdef _MEM_PROF_
    anlz_mem(tgt);
#endif

    MM->avmem+=tgt->len;

    wrkspc=tgt->addr;
#ifdef _DEBUG_MEM_
    printf("Deallocating memory %s at adress %p\n",tgt->elbl, wrkspc);
#endif

//    printf("Deallocating memory %s at adress %p\n",tgt->elbl, wrkspc);
//    printf("Could you see me? - 00\n"); //yma

    if(tgt->len) free(wrkspc);

    if(tgt!=lst) memcpy(tgt,lst,sizeof(mentry));
    bzero((void *) lst, (size_t) sizeof(mentry));

    lst->len=EMPTYE;

    return(0);
}

/*-----------------------------------------------------------------------------*/
INT exc_mentry(mstat *MM, mentry mentries[], mentry *tmp) {
    INT     i;
    mentry *tgt,*lst;

#ifdef _DEBUG_MEM_
    printf("-------- Deleting entry '%9s' of type = '%9s' with length=%12ld (offset=%12ld)\n",tmp->elbl,tmp->etyp,LIFMT(tmp->len),LIFMT(tmp->offset));
#endif

    i=find_mentry(mentries,tmp);

    if(ismax_mentry(i)) return(-1);

    tgt=&mentries[i];
    lst=&mentries[--MM->nmentry];

#ifdef _DEBUG_MEM_
    dump_mentry("REMOVING ENTRY",tgt);
#endif

    MM->avmem+=tgt->len;

    memcpy(tgt,lst,sizeof(mentry));
    bzero((void *) lst, (size_t) sizeof(mentry));

    lst->len=EMPTYE;

    return(0);
}

/*-----------------------------------------------------------------------------*/
void flushMM(mstat *MM, mentry mentries[], mentry *tmp) {
    mentry *tgt;
    INT     latime,i;

    if(MM->nmentry==0) return;

    i=find_mentry(mentries,tmp);
#ifdef _DEBUG_MEM_
    if((ismax_mentry(i))||(mentries[i].len==EMPTYE)) {
       printf("It should never happen!\n");
       abort();
    }
#else
    if(ismax_mentry(i)) return;
#endif
    tgt=&mentries[i];
    latime=tgt->atime;

#ifdef _DEBUG_MEM_
    dump_mentry("FLUSH", tgt);
#endif
#ifdef _DEBUG_MEM_
    printf("Going to delete all memory entries older than atime=%ld\n",LIFMT(latime));
#endif
    for(i=MM->nmentry-1;i>0;i--) if(mentries[i].atime>latime) del_mentry( MM, mentries, &mentries[i], i);

    return;
}

/*-----------------------------------------------------------------------------*/
void list_MlM(mstat *MM, mentry mentries[]) {

    INT     i=1;
    if(MM->nmentry==0) return;
    printf("---------------------------------------------------------------------------------------------\n");
    printf("  Nr.\t Label\t\tType\t\tOffset\t\tLength\t   Atime\t  Address\n");
    printf("---------------------------------------------------------------------------------------------\n");
    for (i=0; i<MM->nmentry; i++) printf("%3ld\t%-12s\t%4s\t%14ld\t%12ld   %9ld\t[%p]\n", LIFMT(i+1),mentries[i].elbl,mentries[i].etyp,LIFMT(mentries[i].offset),LIFMT(mentries[i].len),LIFMT(mentries[i].atime),mentries[i].addr);
    printf("---------------------------------------------------------------------------------------------\n");
    printf("Maximal available memory for Molcas = %ld\n",LIFMT(MM->avmem));
    return;
}

/*-----------------------------------------------------------------------------*/
void set_mentry(mentry *tmp, char *lbl, char *typ, INT *offset, INT *len) {
     strcpy(tmp->elbl,lbl);
     strcpy(tmp->etyp,typ);
     tmp->len=(*len);
     tmp->offset=(*offset);
     tmp->atime=1;
     tmp->addr=NULL;
     return;
}

void print_params(char *func, char *name, char* Op, char *dtyp, INT *offset, INT *len) {
     printf("%s Calling parameters: ('%s','%s','%s',%ld,%ld)\n", func, name, Op, dtyp, LIFMT(*offset), LIFMT(*len));
     return;
}


INT c_getmem_kern(INT *op, mentry *tmp, INT *offset, INT *len) {
#ifdef _TRACK_
    enum memops {ALLO,FREE,LENG,CHEC,MAX,LIST,TERM,FLUS,PINN,RGST,EXCL,TRCK};
#else
    enum memops {ALLO,FREE,LENG,CHEC,MAX,LIST,TERM,FLUS,PINN,RGST,EXCL};
#endif

    static mentry MDATA[MAXREC]={{"\0","\0",0,EMPTYE,0,NULL}};

    mentry *tmpp;
    INT     rc=1,i,allo=1,maxMM;

#ifdef  _MEMORY_MTRACE_
    mcheck_check_all();
#endif
    switch (*op) {
        case PINN:
            allo=0;
            /* fall through */
        case ALLO:
            tmp->offset=(allo)?UNDEF_MEM:PINNED_MEM; /* Just to be sure that we are not going to allocate the PINNED memory for wrong reasons */
/* Checking for memory leaks     */
            if(MlM.nmentry==MAXREC) {
               list_MlM(&MlM,MDATA);
               printf("MEMORY ERROR: Possible memory leak detected: The number of memory blocks exceeds the limit of %d entries\n",MAXREC);
               return(-3);
            }
/* Checking for free memory available */
            if(MlM.avmem<tmp->len) {
              if((MlM.avmem+MlM.mxmem)<tmp->len) {
                list_MlM(&MlM,MDATA);
                printf("MEMORY ERROR: Memory is exhausted!\n");
                printf("MEMORY ERROR: Available memory = %ld ( %ld Mb ) !\n",LIFMT(MlM.avmem+MlM.mxmem),LIFMT((MlM.avmem+MlM.mxmem)/MB));
                printf("MEMORY ERROR: Requested memory = %ld ( %ld Mb ) !\n",LIFMT(tmp->len),LIFMT(tmp->len/MB));
                printf("MEMORY ERROR: The suggested MOLCAS_MEM=%ld !\n",LIFMT((tmp->len-MlM.avmem+MlM.totmem)/MB+1));
                return(-4);
              } else {
#ifdef _DEBUG_MEM_
                printf("MEMORY WARNING: MOLCAS_MEM has been increased by MOLCAS_MAXMEM (%ld) !\n",LIFMT(MlM.mxmem));
#endif
                MlM.avmem+=tmp->len;
                MlM.mxmem-=tmp->len;
              }
            }
           *offset=add_mentry(&MlM, MDATA, tmp);
            if(*offset==ALLOC_FLD) { /* Failed to allocate the requested memory */
               list_MlM(&MlM,MDATA);
               return(-5);
            }

            break;
        case FREE:
            if(MlM.nmentry==0) {
              printf("WARNING: Attempt to operate on zero allocated memory blocks\n");
              exit(-3);
            }
            rc=del_mentry(&MlM, MDATA, tmp, 0);

            if(rc<0) list_MlM(&MlM,MDATA);

            break;
        case MAX:
#ifdef _DEBUG_MEM_
            printf("MOLCAS_MEM=%ld byte \n",LIFMT(SFCTR(MlM.avmem)));
#endif
            maxMM=SFCTR(MlM.avmem);

            /* if malloc is intercepted, we should not test if the memory
            * can be allocated (testmem) as that might result in very expensive
            * memset calls from various memory inspection tools */
#ifndef _MALLOC_INTERCEPT_
            while(maxMM>0 && testmem(&maxMM)<0) maxMM=SFCTR(maxMM);
#endif
#ifdef _DEBUG_MEM_
            if(maxMM!=SFCTR(MlM.avmem)) printf("MOLCAS MAXMEM: initial = %ld, allocatable =%ld\n",LIFMT(SFCTR(MlM.avmem)),LIFMT(maxMM));
#endif
            if(maxMM<=0) {
               printf("MEMORY ERROR: the memory limit has been reached. No window for further memory allocation.\n");
               rc=-1;
            }

           *len=maxMM/dsize(tmp->etyp);

            break;
        case CHEC:
#ifdef _MEMORY_TRACE_
            mcheck_check_all();
#endif
            break;
        case LIST:
            list_MlM(&MlM,MDATA);
            break;
        case LENG:
            i=find_mentry(MDATA,tmp);
            tmpp=&MDATA[i];
           *len=(tmpp)?(tmpp->len/dsize(tmpp->etyp)):0;
            break;
        case TERM:
#ifdef  _MEMORY_MTRACE_
            muntrace();
            mcheck_check_all();
#endif
            tmp->offset=SYS_ATIME;
            if(MlM.nmentry==0) {
               break;
            } else {
#if defined(_DEBUG_MEM_) || defined(_BIGOT_)
               printf("MEMORY ERROR: some memory allocations are not released!\n");
               abort();
#else
               printf("MEMORY WARNING: some memory allocations are not released!\n");
#endif
             break;
            }
        case FLUS:
            printf("**************************************************\n");
            printf("MEMORY WARNING: use of FLUSH operation deprecated!\n");
            printf("please contact the developer of this module and\n");
            printf("ask him/her to fix this!\n");
            printf("**************************************************\n");
            flushMM(&MlM,MDATA,tmp);
            break;
        case RGST:
             rc=reg_mentry(&MlM, MDATA, tmp);
             break;
        case EXCL:
             rc=exc_mentry(&MlM, MDATA, tmp);
             break;
#ifdef _TRACK_
        case TRCK:
             rc=trc_mentry(&MlM, MDATA, tmp);
            *len=rc/(dsize(tmp->etyp));
             break;
#endif
        default: rc=-1;printf("Unsupported memory operation !\n");
    }
#ifdef _DEBUG_MEM_
    list_MlM(&MlM,MDATA);
    if(MlM.avmem<0)  {
      printf("Integer Overflow has been detected!\n");
      rc=-1;
    }
#endif
    return(rc);
}

#ifdef _TRACK_
void track_mem(char *address) {

    mentry  tmp;
    char e_name[]="TRCK",Op[]="TRCK",e_dtyp[]="CHAR";
    INT offset=0,blen=1,op;

    bzero(&tmp,sizeof(mentry));
    set_mentry(&tmp, e_name, e_dtyp, &offset, &blen);
    tmp.addr=address;

    op=memop(Op);

    c_getmem_kern(&op, &tmp, &offset, &blen);

    return;
}
#endif



INT c_getmem(char *name, char* Op, char *dtyp, INT *offset, INT *len) {
#ifdef _TRACK_
    enum memops {ALLO,FREE,LENG,CHEC,MAX,LIST,TERM,FLUS,PINN,RGST,EXCL,TRCK};
#else
    enum memops {ALLO,FREE,LENG,CHEC,MAX,LIST,TERM,FLUS,PINN,RGST,EXCL};
#endif

    mentry  tmp;

    INT rc,op,blen;
    char e_name[MXLINE];
    char e_dtyp[MXLINE];
    char action[MXLINE];

    string2UC(name,e_name);
    string2UC(dtyp,e_dtyp);
    string2UC(Op,action);

    op=memop(action);
#ifdef _DEBUG_MEM_
    print_params("C_GetMem", name, Op, dtyp, offset,len);
    printf("Op=%s (%ld)\n",Op,LIFMT(memop(Op)));
#endif

    blen=(e_dtyp[0]=='C')?1:0;
    blen+=dsize(e_dtyp)*(*len);

    bzero(&tmp,sizeof(mentry));
    set_mentry(&tmp, e_name, e_dtyp, offset, &blen);
#ifdef _OPENMP
    omp_set_lock(&mma_lock);
#endif
    rc=c_getmem_kern(&op, &tmp, offset, len);
#ifdef _OPENMP
    omp_unset_lock(&mma_lock);
    if(op==TERM) omp_destroy_lock(&mma_lock);
#endif
    if(rc<0) {
        print_params("C_GetMem",name, Op, dtyp, offset,len);
    }
    return(rc);
}


char *allomblck(char *name,  INT *len) {

      mentry  tmp;
      char    Op[]="ALLO", ctyp[]="CHAR", *mblck=NULL, e_name[MXLINE];
      INT     rc,op,blen,offset=0;

      op=memop(Op);
      string2UC(name,e_name);
#ifdef _DEBUG_MEM_
      print_params("C_GetMem", name, Op, ctyp, &offset,len);
      printf("Op=%s (%ld)\n",Op,LIFMT(memop(Op)));
#endif

      blen=(*len)+1;
      bzero(&tmp,sizeof(mentry));
      set_mentry(&tmp, e_name, ctyp, &offset, &blen);
      tmp.atime=SYS_ATIME;

      rc=c_getmem_kern(&op, &tmp, &offset, len);

      if(rc>=0) {
         mblck=woff2cptr(ctyp, offset);
      } else {
         print_params("C_GetMem",name, Op, ctyp, &offset,len);
      }

      return(mblck);
}

char *pinnmblck(char *name,  INT *len) {

      char Op[]="PINN", ctyp[]="CHAR", *mblck=NULL;
      INT rc,offset;

      rc=c_getmem(name, Op, ctyp, &offset, len);
      if(rc>=0) mblck=woff2cptr(ctyp, offset);

      return(mblck);
}


INT freemblck(char *mblck) {

    char Op[]="FREE", ctyp[]="CHAR", name[]="DELMEM";
    INT offset,len,rc;

    offset=cptr2woff(ctyp, mblck);
    rc=c_getmem(name, Op, ctyp, &offset, &len);
    mblck=NULL;
    return(rc);

}

INT lengmblck(char *mblck) {

    char Op[]="LENG", ctyp[]="CHAR", name[]="LENMEM";
    INT offset,len;

    offset=cptr2woff(ctyp, mblck);
    c_getmem(name, Op, ctyp, &offset, &len);

    return(len);

}

INT trckmblck(char *mblck) {

    char Op[]="TRCK", ctyp[]="CHAR", name[]="TRACK";
    INT offset,len,rc;

    offset=cptr2woff(ctyp, mblck);
    rc=c_getmem(name, Op, ctyp, &offset, &len);

    return(rc);

}
