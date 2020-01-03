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
c     ----------------------
c
c
c     this file contains following routines
c     exppsb
c     zasun
c     inittemp
c     closetemp
c     mktempnam
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine exppsb (symp,symq,symr,syms,
     &                    valn,jn,kn,ln)
c
c     this routine realize expansion of symmetry block
c     symp,symq,symr,syms <p,q|r,s>  <-> (IJ|KL), provided such integrals exists
c     It found corresponding (IJ|KL) and expand it to opened
c     NORB(symp) TEMP files with a structure
c     indq,indr,inds,value, each TEMP for one p
c
c     N.B. This process can be accelerated, if exppbs would be
c     divided into exppsb1-8, each for given typ, since this
c     routine is common for all types.
c
c     types of (ij|kl) NI,J,K,L defined in III
c
c     1 - si=sk, si=sj, sk=sl
c     2 - si=sk, si=sj, sk>sl
c     3 - si=sk, si>sj, sk=sl
c     4 - si=sk, si>sj, sk>sl
c     5 - si>sk, si=sj, sk=sl
c     6 - si>sk, si=sj, sk>sl
c     7 - si>sk, si>sj, sk=sl
c     8 - si>sk, si>sj, sk>sl
c
c
       implicit real*8 (a-h,o-z)


#include "SysDef.fh"
#include "reorg.fh"
#include "ccsort.fh"
c
       integer symp,symq,symr,syms
       real*8 valn(1:nsize,1:mbas)
       integer jn(1:nsize,1:mbas)
       integer kn(1:nsize,1:mbas)
       integer ln(1:nsize,1:mbas)
c
c     help variables
c
       integer idis13,indtemp
       integer ni,nj,nk,nl,nsi,nsj,nsk,nsl,i1,j1,k1,l1
       integer iup,ilow,jup,jlow,kup,lup,iold,jold,kold,lold
c
       integer nhelp1,nhelp2,m3
       integer yes1,yes234,yes5,yes678
       integer typp
       integer ind(1:4)
#include "tratoc.fh"
       integer INDMAX
       parameter (INDMAX=nTraBuf)
       REAL*8    TWO(INDMAX)
c
cI    get  adress
       idis13=idis(symp,symq,symr)
c
cII   prepairing nshow vector
c
       do nhelp1=1,norb(symp)
       nshow(nhelp1)=0
       end do
c
cIII.1define order of indices
c
       ni=np(symp,symq,symr)
       nj=nq(symp,symq,symr)
       nk=nr(symp,symq,symr)
       nl=ns(symp,symq,symr)
c
cIII.2def yes1-8
c
       typp=typ(symp,symq,symr)
c
c:1   combination (ij|kl) -> (ij|kl)
c     used in types: 1,2,3,4,5,6,7,8 (all)
       yes1=1
c
c:2   combination (ij|kl) -> (ji|kl)
c:3   combination (ij|kl) -> (ij|lk)
c:4   combination (ij|kl) -> (ji|lk)
c     used in types: 1,5 since 2,3,6,7 never appear
       if ((typp.eq.1).or.(typp.eq.5)) then
       yes234=1
       else
       yes234=0
       end if
c
c:5   combination (ij|kl) -> (kl|ij)
c     used in types: 1,2,3,4
       if ((typp.ge.1).and.(typp.le.4)) then
       yes5=1
       else
       yes5=0
       end if
c
c:6   combination (ij|kl) -> (lk|ij)
c:7   combination (ij|kl) -> (kl|ji)
c:8   combination (ij|kl) -> (lk|ji)
c     used in types: 1 (since 2,3 never appeard)
       if (typp.eq.1) then
       yes678=1
       else
       yes678=0
       end if
c
c
c     define NSI,NSJ,NSK,NSL
       ind(ni)=symp
       ind(nj)=symq
       ind(nk)=symr
       ind(nl)=syms
       NSI=ind(1)
       NSJ=ind(2)
       NSK=ind(3)
       NSL=ind(4)
C
       indtemp=indmax+1
       KUP=NORB(NSK)
       DO 401 KOLD=1,KUP
         if (fullprint.ge.3) write (6,*) ' * K ind ',KOLD
C
       LUP=NORB(NSL)
       IF (NSK.EQ.NSL) LUP=KOLD
       DO 402 LOLD=1,LUP
         if (fullprint.ge.3) write (6,*) ' ** L ind ',LOLD
C
       ILOW=1
       IF (NSI.EQ.NSK) ILOW=KOLD
       IUP=NORB(NSI)
       DO 403 IOLD=ILOW,IUP
         if (fullprint.ge.3) write (6,*) ' *** I ind ',IOLD
C
       JLOW=1
       IF (NSI.EQ.NSK.AND.IOLD.EQ.KOLD) JLOW=LOLD
       JUP=NORB(NSJ)
       IF (NSI.EQ.NSJ) JUP=IOLD
       DO 404 JOLD=JLOW,JUP
         if (fullprint.ge.3) write (6,*) ' **** J ind ',JOLD
C
c
c*    read block of integrals if necessary
c
       if (indtemp.eq.(indmax+1)) then
       indtemp=1
c     read block
       CALL dDAFILE(LUINTM,2,TWO,INDMAX,IDIS13)
       end if
c
c*    write integrals to appropriate positions
c
       val1=TWO(indtemp)
c
c:1   combination (ij|kl) -> (ij|kl)
c     since yes1 is always 1, if structure is skipped
       ind(1)=iold
       ind(2)=jold
       ind(3)=kold
       ind(4)=lold
       j1=ind(nj)
       l1=ind(nl)
       if (symq.eq.syms) then
       if (l1.gt.j1) then
       goto 21
       end if
       end if
       i1=ind(ni)
       k1=ind(nk)
c
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
c
       if (m3.eq.nsize) then
       call zasun (i1,nsize,
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
c
 21     if (yes234.eq.1) then
c
c:2   combination (ij|kl) -> (ji|kl)
       ind(1)=jold
       ind(2)=iold
       ind(3)=kold
       ind(4)=lold
       j1=ind(nj)
       l1=ind(nl)
       if (symq.eq.syms) then
       if (l1.gt.j1) then
       goto 31
       end if
       end if
       i1=ind(ni)
       k1=ind(nk)
c
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
c
       if (m3.eq.nsize) then
       call zasun (i1,nsize,
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
c
c:3   combination (ij|kl) -> (ij|lk)
 31     ind(1)=iold
       ind(2)=jold
       ind(3)=lold
       ind(4)=kold
       j1=ind(nj)
       l1=ind(nl)
       if (symq.eq.syms) then
       if (l1.gt.j1) then
       goto 41
       end if
       end if
       i1=ind(ni)
       k1=ind(nk)
c
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
c
       if (m3.eq.nsize) then
       call zasun (i1,nsize,
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
c
c:4   combination (ij|kl) -> (ji|lk)
 41     ind(1)=jold
       ind(2)=iold
       ind(3)=lold
       ind(4)=kold
       j1=ind(nj)
       l1=ind(nl)
       if (symq.eq.syms) then
       if (l1.gt.j1) then
       goto 51
       end if
       end if
       i1=ind(ni)
       k1=ind(nk)
c
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
c
       if (m3.eq.nsize) then
       call zasun (i1,nsize,
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
c
       end if
c
c:5   combination (ij|kl) -> (kl|ij)
 51     if (yes5.eq.1) then
       ind(1)=kold
       ind(2)=lold
       ind(3)=iold
       ind(4)=jold
       j1=ind(nj)
       l1=ind(nl)
       if (symq.eq.syms) then
       if (l1.gt.j1) then
       goto 61
       end if
       end if
       i1=ind(ni)
       k1=ind(nk)
c
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
c
       if (m3.eq.nsize) then
       call zasun (i1,nsize,
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
       end if
c
 61     if (yes678.eq.1) then
c
c:6   combination (ij|kl) -> (lk|ij)
       ind(1)=lold
       ind(2)=kold
       ind(3)=iold
       ind(4)=jold
       j1=ind(nj)
       l1=ind(nl)
       if (symq.eq.syms) then
       if (l1.gt.j1) then
       goto 71
       end if
       end if
       i1=ind(ni)
       k1=ind(nk)
c
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
c
       if (m3.eq.nsize) then
       call zasun (i1,nsize,
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
c
c:7   combination (ij|kl) -> (kl|ji)
 71     ind(1)=kold
       ind(2)=lold
       ind(3)=jold
       ind(4)=iold
       j1=ind(nj)
       l1=ind(nl)
       if (symq.eq.syms) then
       if (l1.gt.j1) then
       goto 81
       end if
       end if
       i1=ind(ni)
       k1=ind(nk)
c
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
c
       if (m3.eq.nsize) then
       call zasun (i1,nsize,
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
c
c:8   combination (ij|kl) -> (lk|ji)
 81     ind(1)=lold
       ind(2)=kold
       ind(3)=jold
       ind(4)=iold
       j1=ind(nj)
       l1=ind(nl)
       if (symq.eq.syms) then
       if (l1.gt.j1) then
       goto 100
       end if
       end if
       i1=ind(ni)
       k1=ind(nk)
c
       m3=nshow(i1)+1
       jn(m3,i1)=j1
       kn(m3,i1)=k1
       ln(m3,i1)=l1
       valn(m3,i1)=val1
       nshow(i1)=m3
c
       if (m3.eq.nsize) then
       call zasun (i1,nsize,
     &             valn,jn,kn,ln)
       nshow(i1)=0
       end if
c
       end if
c
 100    indtemp=indtemp+1
c
 404    CONTINUE
 403    CONTINUE
 402    CONTINUE
 401    CONTINUE
C
c
cIV   write the rest integrals if needed
c
       do nhelp1=1,norb(symp)
       nhelp2=nshow(nhelp1)
       if (nhelp2.gt.0) then
       call zasun (nhelp1,nhelp2,
     &             valn,jn,kn,ln)
       end if
       end do
c
       return
       end
c
c     ----------------------------------------
c
       subroutine zasun  (i1,lenght,
     &                    valn,jn,kn,ln)
c
c      control routine over zasun process
c
c     i1 - number of pivot index (I)
c     lenght - number of valid integrals in block (I)
c
       integer lenght,i1
#include "reorg.fh"
       real*8 valn(1:nsize,1:mbas)
       integer jn(1:nsize,1:mbas)
       integer kn(1:nsize,1:mbas)
       integer ln(1:nsize,1:mbas)
c
       if (zrkey.eq.1) then
         call zasun_zr  (i1,lenght,
     &                   valn,jn,kn,ln)
       else
         call zasun_pck  (i1,lenght,
     &                    valn,jn,kn,ln)
       end if
c
       return
       end
c
c     ----------------------------------------
c
       subroutine zasun_zr  (i1,lenght,
     &                       valn,jn,kn,ln)
c
c     this routine write one block of 3-indices and appropriate
c     values of integrals into an opened TEMP-file
c
c     i1 - number of pivot index (I)
c     lenght - number of valid integrals in block (I)
c     this routine has also jn,kn,ln,valn
c     and stattemp and tmpnam as inputs, but they are
c     transpotred through commons  in reorg.fh
c
       implicit real*8 (a-h,o-z)
       integer lenght,i1
#include "reorg.fh"

#include "SysDef.fh"
       real*8 valn(1:nsize,1:mbas)
       integer jn(1:nsize,1:mbas)
       integer kn(1:nsize,1:mbas)
       integer ln(1:nsize,1:mbas)
c
c     help variable
c
       integer m2 ! ,iRec
c
       integer jkl(1:nsize)
       integer constj
       parameter (constj=1048576)
       integer constk
       parameter (constk=1024)
       integer f_iostat
       logical is_error
c
c*     pack indexes
c
       do m2=1,lenght
          jkl(m2)=ln(m2,i1)+constj*jn(m2,i1)
          jkl(m2)=jkl(m2)+constk*kn(m2,i1)
       end do
c
c*     open corresponding TEMP file in corresponding form
c
       if (iokey.eq.1) then
c
c      Fortran IO
c
       if (stattemp(i1).eq.0) then
c        file will be opened first time, it must be opened
c        whith the pointer at then first possition
         call molcas_binaryopen_vanilla(lunpublic,tmpnam(i1))
c         open (unit=lunpublic,
c     &         file=tmpnam(i1),
c     &         form='unformatted',
c     &         status='unknown')
         stattemp(i1)=1
c
       else
c        file was alredy used in expansion of this block, it must
c        be opened whith the pointer at the end of the file
c@#ifdef _DECAXP_
       call molcas_open_ext2(lunpublic,tmpnam(i1),'append',
     &   'unformatted',f_iostat,.false.,1,'unknown',is_error)
cvv         open (unit=lunpublic,
cvv     &         file=tmpnam(i1),
cvv     &         form='unformatted',
cvv     &         status='unknown',
cvv     &         access='append')
c@#else
c@      call molcas_binaryopen_vanilla(lunpublic,tmpnam(i1))
c         open (unit=lunpublic,
c     &         file=tmpnam(i1),
c     &         form='unformatted',
c     &         status='unknown')
c@       Do iRec = 1,nrectemp(i1)
c@         Read (lunpublic) m2
c@       End Do
c@#endif
c
       end if
c
       write (lunpublic) (valn(i,i1),i=1,lenght),
     &                   (jkl(i),i=1,lenght)
       close (lunpublic)
c
       else
c
c      MOLCAS IO
c
       call daname (lunpublic,tmpnam(i1))
       call ddafile (lunpublic,1,valn(1,i1),lenght,stattemp(i1))
       call idafile (lunpublic,1,jkl,       lenght,stattemp(i1))
       call daclos (lunpublic)
c
       end if
c
       nrectemp(i1)=nrectemp(i1)+1
       lrectemp(i1)=lenght
c
       return
       end
c
c     ----------------------------------------
c
       subroutine zasun_pck  (i1,lenght,
     &                        valn,jn,kn,ln)
c
c     this routine write one block of 3-indices and appropriate
c     values of integrals into an opened TEMP-file
c
c     i1 - number of pivot index (I)
c     lenght - number of valid integrals in block (I)
c     this routine has also jn,kn,ln,valn
c     and stattemp and tmpnam as inputs, but they are
c     transpotred through commons  in reorg.fh
c
       implicit real*8 (a-h,o-z)
       integer lenght,i1
#include "reorg.fh"

#include "SysDef.fh"
       real*8 valn(1:nsize,1:mbas)
       integer jn(1:nsize,1:mbas)
       integer kn(1:nsize,1:mbas)
       integer ln(1:nsize,1:mbas)
c
c     help variable
c
       integer m2,iRec
c
       integer jkl(1:nsize)
       integer constj
       parameter (constj=1048576)
       integer constk
       parameter (constk=1024)
c
       character*(RtoB+ItoB) pp(1:nsize),pphelp
       character*1 p1help(1:(RtoB+ItoB))
       real*8 rhelp
       integer ihelp
       equivalence (p1help(1),pphelp)
       equivalence (p1help(1),rhelp)
       equivalence (p1help(1+RtoB),ihelp)
c
c*     pack indexes and integral values
c
       do m2=1,lenght
       jkl(m2)=ln(m2,i1)+constj*jn(m2,i1)
       end do
       do m2=1,lenght
       jkl(m2)=jkl(m2)+constk*kn(m2,i1)
       end do
c
       do m2=1,lenght
       rhelp=valn(m2,i1)
       ihelp=jkl(m2)
       pp(m2)=pphelp
       end do
c
c*     open corresponding TEMP file in corresponding form
c
       if (iokey.eq.1) then
c
c      Fortran IO
c
       if (stattemp(i1).eq.0) then
c        file will be opened first time, it must be opened
c        whith the pointer at then first possition
         call molcas_binaryopen_vanilla(lunpublic,tmpnam(i1))
c         open (unit=lunpublic,
c     &         file=tmpnam(i1),
c     &         form='unformatted',
c     &         status='unknown')
         stattemp(i1)=1
c
       else
c        file was alredy used in expansion of this block, it must
c        be opened whith the pointer at the end of the file
#ifdef _DECAXP_
         call molcas_open_ext2(lunpublic,tmpnam(i1),'append',
     &   'unformatted',f_iostat,.false.,1,'unknown',is_error)

c         open (unit=lunpublic,
c     &         file=tmpnam(i1),
c     &         form='unformatted',
c     &         status='unknown',
c     &         access='append')
#else
         call molcas_binaryopen_vanilla(lunpublic,tmpnam(i1))
c         open (unit=lunpublic,
c     &         file=tmpnam(i1),
c     &         form='unformatted',
c     &         status='unknown')
         Do iRec = 1,nrectemp(i1)
           Read (lunpublic) m2
         End Do
#endif
c
       end if
c
       call zashlp1 (lunpublic,pp,lenght)
       close (lunpublic)
c
       else
c
c      MOLCAS IO
c
       call daname (lunpublic,tmpnam(i1))
       call cdafile (lunpublic,1,pp,(RtoB+ItoB)*lenght,stattemp(i1))
       call daclos (lunpublic)
c
       end if
c
       nrectemp(i1)=nrectemp(i1)+1
       lrectemp(i1)=lenght
c
       return
       end
c
c     ---------
c
        subroutine zashlp1 (lunpublic,pp,lenght)

#include "SysDef.fh"
        integer lunpublic,lenght
        character*(RtoB+ItoB) pp(1:lenght)
        write (lunpublic) pp
        return
        end
c
c     ----------------------
c
       subroutine inittemp (num)
c
c     this routine initialize status matrix
c     num - number of files to be used (I)
c
       implicit real*8 (a-h,o-z)
#include "reorg.fh"
       integer num
c
c     help variables
c
       integer nhelp
c
       do nhelp=1,num
       stattemp(nhelp)=0
       nrectemp(nhelp)=0
       lrectemp(nhelp)=0
       end do
c
       return
       end
c
c     ----------------------
c
#ifdef _NOTUSED_
       subroutine closetemp (num)
c
c     this routine close TEMP1 - TEMPnum Temporarry files
c     This routine scratch all files
c     num - number of files to be closed (I)
c
       implicit real*8 (a-h,o-z)
#include "reorg.fh"
       integer num
c
c     help variables
c
       integer nhelp
       real*8 scal
c
       scal=0.0d0
c
       do nhelp=1,num
c
       if (iokey.eq.1) then
c      Fortran IO
         call molcas_binaryopen_vanilla(lunpublic,tmpnam(nhelp))
c         open (unit=lunpublic,file=tmpnam(nhelp),form='unformatted')
         write (lunpublic) scal
         close (lunpublic)
c
       else
c      MOLCAS IO
         call daname (lunpublic,tmpnam(nhelp))
         call daeras (lunpublic)
       end if
c
       end do
c
       return
       end
#endif
c
c     ----------------------------------------
c
       subroutine mktempanam
c
c     this routine prepare names for TEMP and files as
c     TEMP001 - TEMPmbas and store them into
c     tmpnam and tmanam arrays (mbas-maximum number of basis functions)
c
c     variables used:
c     tmp-anam - array of TEMP file names (Transported through common /tmnames/)
c     this routine (I)
c
       implicit real*8 (a-h,o-z)
#include "reorg.fh"
       integer lun,itemp,k1
c
       lun=lunpublic
       call molcas_open(lun,'TEMP000')
c       open (unit=lun,file='TEMP000')
c
       itemp=0
       do 100 k1=1,9
       itemp=itemp+1
       if (itemp.gt.mbas) goto 500
       write (lun,99) k1
 99     format (6hTEMP00,i1)
 100    continue
c
       do 200 k1=10,99
       itemp=itemp+1
       if (itemp.gt.mbas) goto 500
       write (lun,199) k1
 199    format (5hTEMP0,i2)
 200    continue
c
       do 300 k1=100,mbas
       itemp=itemp+1
       if (itemp.gt.mbas) goto 500
       write (lun,299) k1
 299    format (4hTEMP,i3)
 300    continue
c
 500    rewind (lun)
c
       do 600 itemp=1,mbas
       read (lun,599) tmpnam(itemp)
 599    format (a7)
 600    continue
c
       rewind (lun)
       write (lun,*) ' File scratched'
       close (lun)
c
       return
       end
c
c     ----------------------------------------
c
       subroutine esb_ic_1 (symp,symq,symr,syms,
     c                      Vic,dimp,dimq,dimr,dims)
c
c     this routine realize expansion of symmetry block
c     symp,symq,symr,syms <p,q|r,s>  <-> (IJ|KL), provided such integrals exists
c     It found corresponding (IJ|KL) and expand it to
c     matrix vic (np,nq,nr,ns)
c

#include "SysDef.fh"
#include "reorg.fh"
#include "ccsort.fh"
c
        integer symp,symq,symr,syms
        integer dimp,dimq,dimr,dims
               real*8 Vic(1:dimp,1:dimq,1:dimr,1:dims)
        real*8 val1
c

       integer idis13,indtemp
       integer ni,nj,nk,nl,nsi,nsj,nsk,nsl,i1,j1,k1,l1
       integer iup,ilow,jup,jlow,kup,lup,iold,jold,kold,lold
c
c     help variables
       integer yes1,yes234,yes5,yes678
       integer typp
       integer ind(1:4)
#include "tratoc.fh"
       integer INDMAX
       parameter (INDMAX=nTraBuf)
       REAL*8    TWO(INDMAX)
c
cI    get  adress
       idis13=idis(symp,symq,symr)
c
cIII.1define order of indices
c
       ni=np(symp,symq,symr)
       nj=nq(symp,symq,symr)
       nk=nr(symp,symq,symr)
       nl=ns(symp,symq,symr)
c
cIII.2def yes1-8
c
       typp=typ(symp,symq,symr)
c
c:1   combination (ij|kl) -> (ij|kl)
c     used in types: 1,2,3,4,5,6,7,8 (all)
       yes1=1
c
c:2   combination (ij|kl) -> (ji|kl)
c:3   combination (ij|kl) -> (ij|lk)
c:4   combination (ij|kl) -> (ji|lk)
c     used in types: 1,5 since 2,3,6,7 never appear
       if ((typp.eq.1).or.(typp.eq.5)) then
       yes234=1
       else
       yes234=0
       end if
c
c:5   combination (ij|kl) -> (kl|ij)
c     used in types: 1,2,3,4
       if ((typp.ge.1).and.(typp.le.4)) then
       yes5=1
       else
       yes5=0
       end if
c
c:6   combination (ij|kl) -> (lk|ij)
c:7   combination (ij|kl) -> (kl|ji)
c:8   combination (ij|kl) -> (lk|ji)
c     used in types: 1 (since 2,3 never appeard)
       if (typp.eq.1) then
       yes678=1
       else
       yes678=0
       end if
c
c
c     define NSI,NSJ,NSK,NSL
       ind(ni)=symp
       ind(nj)=symq
       ind(nk)=symr
       ind(nl)=syms
       NSI=ind(1)
       NSJ=ind(2)
       NSK=ind(3)
       NSL=ind(4)
C
       indtemp=indmax+1
       KUP=NORB(NSK)
       DO 401 KOLD=1,KUP
C
       LUP=NORB(NSL)
       IF (NSK.EQ.NSL) LUP=KOLD
       DO 402 LOLD=1,LUP
C
       ILOW=1
       IF (NSI.EQ.NSK) ILOW=KOLD
       IUP=NORB(NSI)
       DO 403 IOLD=ILOW,IUP
C
       JLOW=1
       IF (NSI.EQ.NSK.AND.IOLD.EQ.KOLD) JLOW=LOLD
       JUP=NORB(NSJ)
       IF (NSI.EQ.NSJ) JUP=IOLD
       DO 404 JOLD=JLOW,JUP
C
c
c*    read block of integrals if neccesarry
c
       if (indtemp.eq.(indmax+1)) then
       indtemp=1
c     read block
       CALL dDAFILE(LUINTM,2,TWO,INDMAX,IDIS13)
       end if
c
c*    write integrals to appropriate possitions
c
       val1=TWO(indtemp)
c
c:1   combination (ij|kl) -> (ij|kl)
c     since yes1 is always 1, if structure is skipped
       ind(1)=iold
       ind(2)=jold
       ind(3)=kold
       ind(4)=lold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimr) then
             if (l1.le.dims) then
             Vic(i1,j1,k1,l1)=val1
             end if
           end if
         end if
       end if
c
       if (yes234.eq.1) then
c
c:2   combination (ij|kl) -> (ji|kl)
       ind(1)=jold
       ind(2)=iold
       ind(3)=kold
       ind(4)=lold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimr) then
             if (l1.le.dims) then
             Vic(i1,j1,k1,l1)=val1
             end if
           end if
         end if
       end if
c
c:3   combination (ij|kl) -> (ij|lk)
       ind(1)=iold
       ind(2)=jold
       ind(3)=lold
       ind(4)=kold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimr) then
             if (l1.le.dims) then
             Vic(i1,j1,k1,l1)=val1
             end if
           end if
         end if
       end if
c
c:4   combination (ij|kl) -> (ji|lk)
       ind(1)=jold
       ind(2)=iold
       ind(3)=lold
       ind(4)=kold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimr) then
             if (l1.le.dims) then
             Vic(i1,j1,k1,l1)=val1
             end if
           end if
         end if
       end if
c
       end if
c
c:5   combination (ij|kl) -> (kl|ij)
       if (yes5.eq.1) then
       ind(1)=kold
       ind(2)=lold
       ind(3)=iold
       ind(4)=jold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimr) then
             if (l1.le.dims) then
             Vic(i1,j1,k1,l1)=val1
             end if
           end if
         end if
       end if
c
       end if
c
       if (yes678.eq.1) then
c
c:6   combination (ij|kl) -> (lk|ij)
       ind(1)=lold
       ind(2)=kold
       ind(3)=iold
       ind(4)=jold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimr) then
             if (l1.le.dims) then
             Vic(i1,j1,k1,l1)=val1
             end if
           end if
         end if
       end if
c
c:7   combination (ij|kl) -> (kl|ji)
       ind(1)=kold
       ind(2)=lold
       ind(3)=jold
       ind(4)=iold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimr) then
             if (l1.le.dims) then
             Vic(i1,j1,k1,l1)=val1
             end if
           end if
         end if
       end if
c
c:8   combination (ij|kl) -> (lk|ji)
       ind(1)=lold
       ind(2)=kold
       ind(3)=jold
       ind(4)=iold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimr) then
             if (l1.le.dims) then
             Vic(i1,j1,k1,l1)=val1
             end if
           end if
         end if
       end if
c
       end if
c
        indtemp=indtemp+1
c
 404    CONTINUE
 403    CONTINUE
 402    CONTINUE
 401    CONTINUE
C
c
       return
       end
c
c     ----------------------------------------
c
       subroutine esb_ic_2 (symp,symq,Vic,dimp,dimq,pqind)
c
c     this routine realize expansion of symmetry block
c     symp,symq,symr,syms <p,q|r,s>  <-> (IJ|KL),
c     (for case symp=symr,symq=syms)
c     provided such integrals exists
c     It found corresponding (IJ|KL) and expand it to
c     matrix vic (pr,qs)
c

#include "SysDef.fh"
#include "reorg.fh"
#include "ccsort.fh"
c
        integer symp,symq
        integer dimp,dimq
               real*8 Vic(1:(dimp*(dimp+1)/2),1:(dimq*(dimq+1)/2))
        real*8 val1
c

       integer idis13,indtemp
       integer ni,nj,nk,nl,nsi,nsj,nsk,nsl,i1,j1,k1,l1
       integer iup,ilow,jup,jlow,kup,lup,iold,jold,kold,lold
       integer pqind(1:mbas,1:mbas)
c
c     help variables
       integer i,j,maxx
       integer yes1,yes234,yes5,yes678
       integer typp
       integer ind(1:4)
#include "tratoc.fh"
       integer INDMAX
       parameter (INDMAX=nTraBuf)
       REAL*8    TWO(INDMAX)
c
cI        calc pqind
c
        if (dimp.ge.dimq) then
          maxx=dimp
        else
          maxx=dimq
        end if
c
        do i=1,maxx
        do j=1,maxx
          if (i.ge.j) then
            pqind(i,j)=i*(i-1)/2+j
          else
            pqind(i,j)=j*(j-1)/2+i
          end if
        end do
        end do
c
cII    get  adress
       idis13=idis(symp,symq,symp)
c
cIII.1define order of indices
c
       ni=np(symp,symq,symp)
       nj=nq(symp,symq,symp)
       nk=nr(symp,symq,symp)
       nl=ns(symp,symq,symp)
c
cIII.2def yes1-8
c
       typp=typ(symp,symq,symp)
c
c:1   combination (ij|kl) -> (ij|kl)
c     used in types: 1,2,3,4,5,6,7,8 (all)
       yes1=1
c
c:2   combination (ij|kl) -> (ji|kl)
c:3   combination (ij|kl) -> (ij|lk)
c:4   combination (ij|kl) -> (ji|lk)
c     used in types: 1,5 since 2,3,6,7 never appear
       if ((typp.eq.1).or.(typp.eq.5)) then
       yes234=1
       else
       yes234=0
       end if
c
c:5   combination (ij|kl) -> (kl|ij)
c     used in types: 1,2,3,4
       if ((typp.ge.1).and.(typp.le.4)) then
       yes5=1
       else
       yes5=0
       end if
c
c:6   combination (ij|kl) -> (lk|ij)
c:7   combination (ij|kl) -> (kl|ji)
c:8   combination (ij|kl) -> (lk|ji)
c     used in types: 1 (since 2,3 never appeard)
       if (typp.eq.1) then
       yes678=1
       else
       yes678=0
       end if
c
c
c     define NSI,NSJ,NSK,NSL
       ind(ni)=symp
       ind(nj)=symq
       ind(nk)=symp
       ind(nl)=symq
       NSI=ind(1)
       NSJ=ind(2)
       NSK=ind(3)
       NSL=ind(4)
C
       indtemp=indmax+1
       KUP=NORB(NSK)
       DO 401 KOLD=1,KUP
C
       LUP=NORB(NSL)
       IF (NSK.EQ.NSL) LUP=KOLD
       DO 402 LOLD=1,LUP
C
       ILOW=1
       IF (NSI.EQ.NSK) ILOW=KOLD
       IUP=NORB(NSI)
       DO 403 IOLD=ILOW,IUP
C
       JLOW=1
       IF (NSI.EQ.NSK.AND.IOLD.EQ.KOLD) JLOW=LOLD
       JUP=NORB(NSJ)
       IF (NSI.EQ.NSJ) JUP=IOLD
       DO 404 JOLD=JLOW,JUP
C
c
c*    read block of integrals if neccesarry
c
       if (indtemp.eq.(indmax+1)) then
       indtemp=1
c     read block
       CALL dDAFILE(LUINTM,2,TWO,INDMAX,IDIS13)
       end if
c
c*    write integrals to appropriate possitions
c
       val1=TWO(indtemp)
c
c:1   combination (ij|kl) -> (ij|kl)
c     since yes1 is always 1, if structure is skipped
       ind(1)=iold
       ind(2)=jold
       ind(3)=kold
       ind(4)=lold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimp) then
             if (l1.le.dimq) then
             Vic(pqind(i1,k1),pqind(j1,l1))=val1
             end if
           end if
         end if
       end if
c
       if (yes234.eq.1) then
c
c:2   combination (ij|kl) -> (ji|kl)
       ind(1)=jold
       ind(2)=iold
       ind(3)=kold
       ind(4)=lold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimp) then
             if (l1.le.dimq) then
             Vic(pqind(i1,k1),pqind(j1,l1))=val1
             end if
           end if
         end if
       end if
c
c:3   combination (ij|kl) -> (ij|lk)
       ind(1)=iold
       ind(2)=jold
       ind(3)=lold
       ind(4)=kold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimp) then
             if (l1.le.dimq) then
             Vic(pqind(i1,k1),pqind(j1,l1))=val1
             end if
           end if
         end if
       end if
c
c:4   combination (ij|kl) -> (ji|lk)
       ind(1)=jold
       ind(2)=iold
       ind(3)=lold
       ind(4)=kold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimp) then
             if (l1.le.dimq) then
             Vic(pqind(i1,k1),pqind(j1,l1))=val1
             end if
           end if
         end if
       end if
c
       end if
c
c:5   combination (ij|kl) -> (kl|ij)
       if (yes5.eq.1) then
       ind(1)=kold
       ind(2)=lold
       ind(3)=iold
       ind(4)=jold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimp) then
             if (l1.le.dimq) then
             Vic(pqind(i1,k1),pqind(j1,l1))=val1
             end if
           end if
         end if
       end if
c
       end if
c
       if (yes678.eq.1) then
c
c:6   combination (ij|kl) -> (lk|ij)
       ind(1)=lold
       ind(2)=kold
       ind(3)=iold
       ind(4)=jold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimp) then
             if (l1.le.dimq) then
             Vic(pqind(i1,k1),pqind(j1,l1))=val1
             end if
           end if
         end if
       end if
c
c:7   combination (ij|kl) -> (kl|ji)
       ind(1)=kold
       ind(2)=lold
       ind(3)=jold
       ind(4)=iold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimp) then
             if (l1.le.dimq) then
             Vic(pqind(i1,k1),pqind(j1,l1))=val1
             end if
           end if
         end if
       end if
c
c:8   combination (ij|kl) -> (lk|ji)
       ind(1)=lold
       ind(2)=kold
       ind(3)=jold
       ind(4)=iold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
       if (i1.le.dimp) then
         if (j1.le.dimq) then
           if (k1.le.dimp) then
             if (l1.le.dimq) then
             Vic(pqind(i1,k1),pqind(j1,l1))=val1
             end if
           end if
         end if
       end if
c
       end if
c
        indtemp=indtemp+1
c
 404    CONTINUE
 403    CONTINUE
 402    CONTINUE
 401    CONTINUE
C
c
       return
       end
c
c     ----------------------------------------
c
       subroutine esb_ic_3 (symp,Vic,dimp,pqind)
c
c     this routine realize expansion of symmetry block
c     symp,symq,symr,syms <p,q|r,s>  <-> (IJ|KL),
c     (for case symp=symr=symq=syms)
c     provided such integrals exists
c     It found corresponding (IJ|KL) and expand it to
c     matrix vic (prqs)
c

#include "SysDef.fh"
#include "reorg.fh"
#include "ccsort.fh"
c
           integer symp,dimp
CLD           integer symp,dimp,dimq
               real*8 Vic(1:(dimp*(dimp+1)/2)*((dimp*(dimp+1)/2)+1)/2)
        real*8 val1
c

       integer idis13,indtemp
       integer ni,nj,nk,nl,nsi,nsj,nsk,nsl,i1,j1,k1,l1
       integer iup,ilow,jup,jlow,kup,lup,iold,jold,kold,lold
       integer pqind(1:mbas,1:mbas)
c
c     help variables
       integer i,j,maxx,ik,jl,ikjl
       integer ind(1:4)
#include "tratoc.fh"
       integer INDMAX
       parameter (INDMAX=nTraBuf)
       REAL*8    TWO(INDMAX)
c
cI        calc pqind
c
CLD        if (dimp.ge.dimp) then
          maxx=dimp
CLD        else
CLD          maxx=dimq
CLD        end if
c
        do i=1,maxx
        do j=1,maxx
          if (i.ge.j) then
            pqind(i,j)=i*(i-1)/2+j
          else
            pqind(i,j)=j*(j-1)/2+i
          end if
        end do
        end do
c
cII    get  adress
       idis13=idis(symp,symp,symp)
c
cIII.1define order of indices
c
       ni=np(symp,symp,symp)
       nj=nq(symp,symp,symp)
       nk=nr(symp,symp,symp)
       nl=ns(symp,symp,symp)
c
c
c     define NSI,NSJ,NSK,NSL
       ind(ni)=symp
       ind(nj)=symp
       ind(nk)=symp
       ind(nl)=symp
       NSI=ind(1)
       NSJ=ind(2)
       NSK=ind(3)
       NSL=ind(4)
C
       indtemp=indmax+1
       KUP=NORB(NSK)
       DO 401 KOLD=1,KUP
         if (fullprint.ge.3) write (6,*) ' * K ind ',KOLD
C
       LUP=NORB(NSL)
       IF (NSK.EQ.NSL) LUP=KOLD
       DO 402 LOLD=1,LUP
         if (fullprint.ge.3) write (6,*) ' ** L ind ',LOLD
C
       ILOW=1
       IF (NSI.EQ.NSK) ILOW=KOLD
       IUP=NORB(NSI)
       DO 403 IOLD=ILOW,IUP
         if (fullprint.ge.3) write (6,*) ' *** I ind ',IOLD
C
       JLOW=1
       IF (NSI.EQ.NSK.AND.IOLD.EQ.KOLD) JLOW=LOLD
       JUP=NORB(NSJ)
       IF (NSI.EQ.NSJ) JUP=IOLD
       DO 404 JOLD=JLOW,JUP
         if (fullprint.ge.3) write (6,*) ' **** J ind ',JOLD
C
c
c*    read block of integrals if neccesarry
c
       if (indtemp.eq.(indmax+1)) then
       indtemp=1
c     read block
       CALL dDAFILE(LUINTM,2,TWO,INDMAX,IDIS13)
       end if
c
c*    write integrals to appropriate possitions
c
       val1=TWO(indtemp)
c
c:1   combination (ij|kl) -> (ij|kl)
c     since yes1 is always 1, if structure is skipped
       ind(1)=iold
       ind(2)=jold
       ind(3)=kold
       ind(4)=lold
       j1=ind(nj)
       l1=ind(nl)
       i1=ind(ni)
       k1=ind(nk)
c
c:2    def iklj
       ik=pqind(i1,k1)
       jl=pqind(j1,l1)
       if (ik.ge.jl) then
         ikjl=ik*(ik-1)/2+jl
       else
         ikjl=jl*(jl-1)/2+ik
       end if
c
       Vic(ikjl)=val1
c
c
        indtemp=indtemp+1
c
 404    CONTINUE
 403    CONTINUE
 402    CONTINUE
 401    CONTINUE
C
c
       return
       end
c
c     ----------------------------------------
c
