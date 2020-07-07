************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
c
c     this file contains following routines:
c     daopen
c     daread
c     dawrite
c     mkabpqmap
c     mkampqmap
c     abpack
c     ampack
c     addintab
c     mkampq
c     vf
c
c     ----------------------------------------
c
       subroutine daopen (name,lun,recl,nrec)
c
c     this routine open direct acces file
c
c     name  - name of the file A8 (I)
c     lun   - logical unit number (I)
c     recl - record lenght in R8 (I)
c     nrec  - number of records (if needed) (I)
c
       integer lun,recl,nrec
       character*8 name
c
c     help variables
c
       integer recln,f_iostat
       logical is_error
c
#ifdef _DECAXP_
       recln=recl*2
#else
       recln=recl*8
#endif
c
       call molcas_open_ext2(lun,name,'direct','unformatted',
     &  f_iostat,.true.,recln,'unknown',is_error)
c       open (unit=lun,file=name,form='unformatted',access='direct',
c     & recl=recln)
c
       return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(nrec)
       end
c
c     ------------------------------
c
       subroutine daread (lun,irec0,vector,lenght,recl)
c
c     this routine read vector with required lenght from
c     opened direct access file lun starting from record number
c     irec0
c     lun   - logical unit of direct access file (I)
c     irec0 - initial recored number (I)
c     vector- vector (O)
c     lenght- number of R8 data to be readed (I)
c     recl  - lenght of one record in lun  in R8 (I)
c
       real*8 vector(1:lenght)
       integer lun,irec0,lenght,recl
c
c     help variables
c
       integer ilow,iup,need,irec,i
c
       if (lenght.eq.0) then
       return
       end if
c
c*    def need,ilow,iup,irec
c
       need=lenght
       ilow=1
       iup=0
       irec=irec0
c
 1      if (recl.ge.need) then
       iup=iup+need
       else
       iup=iup+recl
       end if
c
       read (lun,rec=irec) (vector(i),i=ilow,iup)
c
       need=need-(iup-ilow+1)
       irec=irec+1
       ilow=ilow+recl
c
       if (need.gt.0) goto 1
c
       return
       end
c
c     ------------------------
c
       subroutine dawrite (lun,irec0,vector,lenght,recl)
c
c     this routine write vector with required lenght to
c     opened direct access file lun starting from record number
c     irec0
c
c     lun   - logical unit of direct access file (I)
c     irec0 - initial recored number (I)
c     vector- vector (I)
c     lenght- number of R8 data to be readed (I)
c     recl  - lenght of one record in lun in R8 (I)
c
       real*8 vector(1:lenght)
       integer lun,irec0,lenght,recl
c
c     help variables
c
       integer ilow,iup,need,irec,i
c
       if (lenght.eq.0) then
       return
       end if
c
c*    def need,ilow,iup,irec
c
       need=lenght
       ilow=1
       irec=irec0
       iup=0
c
 1     if (recl.ge.need) then
       iup=iup+need
       else
       iup=iup+recl
       end if
c
       write (lun,rec=irec) (vector(i),i=ilow,iup)
c
       need=need-(iup-ilow+1)
       irec=irec+1
       ilow=ilow+recl
c
       if (need.gt.0) goto 1
c
       return
       end
c
c     ------------------------
c
       subroutine mkabpqmap (abmap,syma,symb,rc)
c
c     this routine prepair abmap
c
#include "reorg.fh"
#include "ccsort.fh"
c
       integer abmap(1:mbas,1:mbas,1:8)
       integer syma,symb,rc
c
c     help variables
c
       integer a,b,bup,symp,symq,symab
       integer lenghtpq,nrecc,nrest,irec
c
cT    test, if there are any ab pair
c
       if (nvb(syma)*nvb(symb).eq.0) then
       rc=1
c     RC=1 : there are no ab pair in this symmetry
       return
       else
       rc=0
       end if
c
c*    def initial address
c
       irec=1
       symab=mul(syma,symb)
c
c*    loop over all combinations
c
       do 100 symp=1,nsym
       symq=mul(symab,symp)
c
c*    define number of records, required to store this block
c     and determine shift in initial possitions
c
       lenghtpq=norb(symp)*norb(symq)
       nrecc=int(lenghtpq/recl)
       nrest=lenghtpq-nrecc*recl
       if (nrest.gt.0) then
       nrecc=nrecc+1
       end if
c
       do 101 a=1,nvb(syma)
c
       if (syma.eq.symb) then
       bup=a
       else
       bup=nvb(symb)
       end if
c
       do 102 b=1,bup
c
       abmap(a,b,symp)=irec
       irec=irec+nrecc
c
 102    continue
 101    continue
 100    continue
c
       return
       end
c
c     ------------------------
c
       subroutine mkampqmap (ammap,syma,rc)
c
c     this routine prepair ammap
c
#include "reorg.fh"
#include "ccsort.fh"
c
       integer syma,rc
       integer ammap(1:mbas,1:8,1:8)
c
c     help variables
c
       integer a,symp,symq,symm,symam
       integer lenghtmpq,nrecc,nrest,irec
c
cT    test, if there are any a in this symmtry
c
       if (nvb(syma).eq.0) then
       rc=1
c     RC=1 : there are no a in this symmetry
       return
       else
       rc=0
       end if
c
c*    def initial address
c
       irec=1
c
c*    loop over all combinations
c
       do 100 symm=1,nsym
       symam=mul(syma,symm)
       do 101 symp=1,nsym
       symq=mul(symam,symp)
c
c*    define number of records, required to store this block
c     and determine shift in initial possitions
c
       lenghtmpq=noa(symm)*norb(symp)*norb(symq)
       nrecc=int(lenghtmpq/recl)
       nrest=lenghtmpq-nrecc*recl
       if (nrest.gt.0) then
       nrecc=nrecc+1
       end if
c
       do 102 a=1,nvb(syma)
c
       ammap(a,symm,symp)=irec
       irec=irec+nrecc
c
 102    continue
 101    continue
 100    continue
c
       return
       end
c
c     --------------------------
c
       subroutine abpack (wrk,wrksize,
     & syma,symb,symp,symq,a,vint,ndimv1,ndimv2,
     & ndimv3,abmap)
c
c     this routine pack corresponding parts to ab direct acc. file
c     from given integrals <_a,bb,p,q> readed in vint
c
c     syma  - irrep of first index (I)
c     symb  - irrep of 2.nd index (I)
c     symp  - irrep of p (I)
c     symq  - irrep of q (I)
c     a- pivot virtual index (counted in nvb set) (I)
c     vint  - array of integrals <_a,bb,p,q> (I)
c     ndimv1- 1.st dimension of vint - norb(symb) (I)
c     ndimv2- 2.nd dimension of vint - norb(symp) (I)
c     ndimv3- 3.rd dimension of vint - norb(symq) (I)
c     abmap - map for storing of addresses in DA file TEMPDA1 (I)
c
#include "wrk.fh"
#include "ccsort.fh"
#include "reorg.fh"
c
       integer syma,symb,symp,symq,a,ndimv1,ndimv2,ndimv3
       real*8 vint(1:ndimv1,1:ndimv2,1:ndimv3)
       integer abmap(1:mbas,1:mbas,1:8)
c
c     help variables
c
       integer p,q,pq,irec0,lenght,b,bup,bvint
c
cT    if there are no ab pair, or no integrals in _a_bpq block return
c
       if (nvb(syma)*nvb(symb)*norb(symp)*norb(symq).eq.0) then
       return
       end if
c
c*    def lenght of _a_b(p,q) block
c
       lenght=norb(symp)*norb(symq)
c
       if (syma.eq.symb) then
       bup=a
       else
       bup=nvb(symb)
       end if
c
c*    cycle over b for given a
c
       do 1000 b=1,bup
       bvint=nob(symb)+b
c
c*    map _a_b(pq) block into #v3
c
       pq=poss30-1
       do 100 q=1,norb(symq)
       do 101 p=1,norb(symp)
       pq=pq+1
       wrk(pq)=vint(bvint,p,q)
 101    continue
 100    continue
c
c*    put this block to iappropriate possition in direct acces file
c
       irec0=abmap(a,b,symp)
       call dawrite (lunda1,irec0,wrk(poss30),lenght,recl)
c
 1000   continue
c
       return
       end
c
c     --------------------------
c
       subroutine ampack (wrk,wrksize,
     & syma,symm,symp,symq,a,vint,ndimv1,ndimv2,ndimv3,
     & ammap)
c
c     this routine pack corresponding parts to ab direct acc. file
c     from given integrals <_a,bb,p,q> readed in vint
c
c     syma  - irrep of first index (I)
c     symm  - irrep of 2.nd index - m (I)
c     symp  - irrep of p (I)
c     symq  - irrep of q (I)
c     a- pivot virtual index (counted in nvb set) (I)
c     vint  - array of integrals <_a,mm,p,q> (I)
c     ndimv1- 1.st dimension of vint - norb(symb) (I)
c     ndimv2- 2.nd dimension of vint - norb(symp) (I)
c     ndimv3- 3.rd dimension of vint - norb(symq) (I)
c     ammap - map for storing of addresses in DA file TEMPDA2 (I)
c
#include "wrk.fh"
#include "ccsort.fh"
#include "reorg.fh"
c
       integer syma,symm,symp,symq,a,ndimv1,ndimv2,ndimv3
       real*8 vint(1:ndimv1,1:ndimv2,1:ndimv3)
       integer ammap(1:mbas,1:8,1:8)
c
c     help variables
c
       integer m,p,q,pq,irec0,lenght
c
cT    if there are no a, or no integrals in _a_mpq block return
c
       if (nvb(syma)*noa(symm)*norb(symp)*norb(symq).eq.0) then
       return
       end if
c
c*    def lenght of _a(m,p,q) block
c
       lenght=noa(symm)*norb(symp)*norb(symq)
c
c*    map _a(mpq) block into #v3
c
       pq=poss30-1
c
       do 100 q=1,norb(symq)
       do 101 p=1,norb(symp)
       do 102 m=1,noa(symm)
       pq=pq+1
       wrk(pq)=vint(m,p,q)
 102    continue
 101    continue
 100    continue
c
c*    put this block to iappropriate possition in direct acces file
c
       irec0=ammap(a,symm,symp)
       call dawrite (lunda2,irec0,wrk(poss30),lenght,recl)
c
       return
       end
c
c     ------------------------
c
       subroutine addintab (wrk,wrksize,
     & syma,symb,abmap)
c
c     this routine add contribution to opened INTAB1 file,
c     comming from ab syma,symb
c
#include "wrk.fh"
#include "reorg.fh"
#include "ccsort.fh"
       integer syma,symb
       integer abmap(1:mbas,1:mbas,1:8)
c
c     help variables
c
       integer nhelp,lenght,symp,symq,symab,irec0,poss3
       integer poss,a,b,bup,ii,rc
c
c*    def symab
       symab=mul(syma,symb)
c
c*    make mapd3,mapi3 for <_a_b|pq>
c
c**   set mapi3=0 (partly)
c
       do 100 nhelp=1,nsym
       do 101 symq=1,nsym
       do 102 symp=1,nsym
       mapi3(symp,symq,nhelp)=0
 102    continue
 101    continue
 100    continue
c
c**   def 0-th row
c
       mapd3(0,1)=5
       mapd3(0,2)=5
       mapd3(0,3)=0
       mapd3(0,4)=0
       mapd3(0,5)=nsym
       mapd3(0,6)=0
c
c**   def other rows
c
       poss=poss30
       do 200 ii=1,nsym
c
       symp=ii
       symq=mul(symab,symp)
       lenght=norb(symp)*norb(symq)
       mapd3(ii,1)=poss
       mapd3(ii,2)=lenght
       mapd3(ii,3)=symp
       mapd3(ii,4)=symq
       mapd3(ii,5)=1
       mapd3(ii,6)=1
       mapi3(symp,1,1)=ii
       poss=poss+lenght
c
 200    continue
c
c*    write mapd,mapi to INTAB
       call dawrtmap (lunab,mapd3,mapi3,rc)
c
cT    if there are no _a_b,pq integrals in this symab,
c     skip sumation over ab
c
       if ((mapd3(nsym,1)+mapd3(nsym,2)).eq.poss30) then
       return
       end if
c
c*    loop over a,b
c
       do 1000 a=1,nvb(syma)
c
       if (syma.eq.symb) then
       bup=a
       else
       bup=nvb(symb)
       end if
c
       do 1001 b=1,bup
c
c**   loop over symp
c
       do 500 symp=1,nsym
c
c***  def irec0 for this a,b,symp in TEMPDA1
       irec0=abmap(a,b,symp)
c
c***  def corresponding possition and lenght in #3
       ii=mapi3(symp,1,1)
       poss3=mapd3(ii,1)
       lenght=mapd3(ii,2)
c
c***  read this block to #3
       if (lenght.gt.0) then
       call daread (lunda1,irec0,wrk(poss3),lenght,recl)
       end if
c
 500    continue
c
c**   since there must be some integrals, write them to TEMPAB
c
       call deflenght (mapd3,lenght)
       call dawri (lunab,lenght,wrk(poss30))
c
 1001   continue
 1000   continue
c
       return
       end
c
c     ------------------------
c
       subroutine mkampq (wrk,wrksize,
     & a,syma,ammap)
c
c     this routine reconstruct #2 V2<_a,m|p,q> from corresponding TEMPDA2 file
c
#include "wrk.fh"
#include "reorg.fh"
#include "ccsort.fh"
       integer a,syma
       integer ammap(1:mbas,1:8,1:8)
c
c     help variables
c
       integer symm,symp,symq,symam
       integer iiv2,lenght,poss,irec0
c
c*    loops over symmetry combinations
       do 100 symm=1,nsym
       symam=mul(syma,symm)
       do 101 symp=1,nsym
       symq=mul(symam,symp)
c
c*    def initioal record possition in TEMPDA2
c     and corresponding possition and lenght in wrk (#2)
c
       irec0=ammap(a,symm,symp)
       iiv2=mapi2(symm,symp,1)
       poss=mapd2(iiv2,1)
       lenght=mapd2(iiv2,2)
c
       if (lenght.gt.0) then
       call daread (lunda2,irec0,wrk(poss),lenght,recl)
       end if
c
 101    continue
 100    continue
c
       return
       end
c
c     ------------------------
c
       subroutine vf (name,lun)
c
c      this routine open file vanisf file with a given name
c      name - name of the vanished file (I)
c      lun  - lun number with which file will be opened (I)
c
       character*8 name
       integer lun
c
       call molcas_open(lun,name)
c       open (unit=lun,file=name)
       write (lun,*) ' File scratched'
       close (lun)
c
       return
       end
c
c     ------------------------
c
