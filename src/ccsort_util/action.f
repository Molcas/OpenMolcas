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
c     action
c     mkmapampq
c     mkmappqij
c     addpqij
c     mkintsta
c     expandfok
c     addinta
c     deflenght
c     fokupdate1
c     fokupdate2
c     unpackk
c     getpp
c     initintabc1
c     addintabc1
c     initwrk
c
c
c     -------------------------------------
c
       subroutine action (foka, fokb, fi,eps)
c
       implicit real*8 (a-h,o-z)
c
c      work file declaration
       integer wrksize
#include "WrkSpc.fh"
c
#include "ccsort.fh"
#include "reorg.fh"
#include "files_ccsd.fh"
       real*8 fi(*)
       real*8 eps(mbas)
c
c     help variables
c
       integer symp,symq,symr,syms,sympq,sympqr
       integer ndimv1,ndimv2,ndimv3,ndimvi
       integer p,a,posst,rc
       integer keyred,vsize,freespace,ickey,iOff_Vic,iOff
       integer t3help1,t3help2,t3help3,t3help4
       integer ipAMMAP,ipABMAP,ipJN,ipKN,ipLN,iOff_valn
       real*8 foka((mbas**2+mbas)/2)
       real*8 fokb((mbas**2+mbas)/2)

       Call qEnter('Action')

c
c*    distribute memory
c
c*.1   calc. work space requirements
       call initwrk (wrksize)
c
c*.2   test if allocation of required memory is possible
c
c*.2.1 allocate work space
c
       Call GetMem('CCSORT','Max','Real',maxspace,maxspace)
       maxspace=maxspace-4
       if (maxspace.lt.wrksize) then
         write(6,*) ' Allocation of work space failed!'
         write(6,*) ' Increase the size of the variable MOLCAS_MEM'
         Call Abend
       end if
       Call GetMem('CCSORT','Allo','Real',iOff,wrksize)
c
c*.3   set wrk = 0
       call ccsort_mv0zero (wrksize,wrksize,Work(iOff))
c
c
c*    def foka,fokb
c
       ndimv1=0
       do 100 symp=1,nsym
       ndimv1=ndimv1+(norb(symp)+1)*norb(symp)/2
 100    continue
c
       do 200 ndimv2=1,ndimv1
       foka(ndimv2)=fi(ndimv2)
       fokb(ndimv2)=fi(ndimv2)
 200    continue
c
c*    make names
c
       call mktempanam
c
c*    if T3 are requited, define T3IntPoss and calc T3Off
       if (t3key.eq.1) then
         call DefT3par (noa,nsym)
         do symp=1,4
         end do
       end if
c
c*    open files INTA1,INTA2,INTA3,INTA4 and INTAB
c
       if (iokey.eq.1) then
c      Fortran IO
       call molcas_binaryopen_vanilla(luna1,'INTA1')
       call molcas_binaryopen_vanilla(luna2,'INTA2')
       call molcas_binaryopen_vanilla(luna3,'INTA3')
       call molcas_binaryopen_vanilla(luna4,'INTA4')
       call molcas_binaryopen_vanilla(lunab,'INTAB')
c       open (unit=luna1,file='INTA1',form='unformatted')
c       open (unit=luna2,file='INTA2',form='unformatted')
c       open (unit=luna3,file='INTA3',form='unformatted')
c       open (unit=luna4,file='INTA4',form='unformatted')
c       open (unit=lunab,file='INTAB',form='unformatted')
c
       else
c      MOLCAS IO
       call daname (luna1,'INTA1')
       call daname (luna2,'INTA2')
       call daname (luna3,'INTA3')
       call daname (luna4,'INTA4')
       call daname (lunab,'INTAB')
       daddr(luna1)=0
       daddr(luna2)=0
       daddr(luna3)=0
       daddr(luna4)=0
       daddr(lunab)=0
       end if
c
c*    define #V1 (for <p,q,i,j>
       call mkmappqij
c
c*    make head of INTAB file for nonsymmetrical (C1) case
       if (nsym.eq.1) then
       call initintabc1
       end if
c
c     allocate space for ammap,abmap
       Call GetMem('AMMAP','ALLO','INTE',ipAMMAP,mbas*64)
       Call GetMem('ABMAP','ALLO','INTE',ipABMAP,mbas*mbas*8)
c
       do 1000 symp=1,nsym
c
       if (fullprint.gt.0) then
       write(6,'(6X,A,2X,I2)') 'Symmetry of the pivot index',symp
       end if
c
c*    define #V2 (for <_a,m,p,q>)
       call mkmapampq (symp)
c
c*    if T3 are required, make maps for R_i
       if ((t3key.eq.1).and.(noa(symp).gt.0)) then
c*            get mapd and mapi for  R_i(a,bc)
        call ccsort_t3grc0(3,8,4,4,4,0,symp,possri0,posst,mapdri,mapiri)
       end if
c
c*    open TEMPDA2 fils for <am|rs> integrals, if there are some virtiuals
c     in symp symmetry
       if (nvb(symp).gt.0) then
       call mkampqmap (iWork(ipAMMAP),symp,rc)
       call daopen ('TEMPDA2 ',lunda2,recl,1)
       end if
c
       do 900 symq=1,nsym
       sympq=mul(symp,symq)
c
       if ((nsym.gt.1).and.(symp.ge.symq)) then
c*    open direct acces file here, to enable exact specification of
c     the number of records (only for symmetrical cases; syma>=symb)
c     N.B. nrec is not needed now
       call daopen ('TEMPDA1 ',lunda1,recl,1)
       end if
c
c*    make abmap for syma>=symb
       if (symp.ge.symq) then
       call mkabpqmap (iWork(ipABMAP),symp,symq,rc)
       end if
c
       do 800 symr=1,nsym
       sympqr=mul(sympq,symr)
       syms=sympqr
c
c*     calc size of the integral file
       if (symp.eq.symr) then
         if (symp.eq.symq) then
          vsize=(norb(symp)*(norb(symp)+1))/2
          vsize=vsize*(vsize+1)/2
          ickey=3
         else
          vsize=(norb(symp)*(norb(symp)+1)*norb(symq)*(norb(symq)+1))/4
          ickey=2
         end if
       else
       vsize=norb(symp)*norb(symq)*norb(symr)*norb(syms)
       ickey=1
       end if
c
       if (vsize.eq.0) goto 800
c
       if (fullprint.gt.1) then
         write(6,'(6X,A,I4,4X,4I2)') 'Block',typ(symp,symq,symr),
     &                                symp,symq,symr,syms
       end if
c
c*     test for incore expansion
       freespace=0
       Call GetMem('CCSORT','Max','Real',freespace,freespace)
       freespace=freespace-4
       if (fullprint.ge.2) then
       write(6,*)
       Write(6,'(6X,A,I10)') 'Available freespace   ',freespace
       Write(6,'(6X,A,I10)') 'Available freespace/MB',freespace/131072
       write(6,'(6X,A,I10)') 'Incore expansion      ',vsize+mbas*mbas
       write(6,'(6X,A,I10)') 'Out of core expansion ',4*nsize*mbas
       endif
c
       if (freespace.ge.(vsize+mbas*mbas)) then
c      INCORE EXPANSION
         if (fullprint.ge.1) then
         write(6,*)
         write(6,'(6X,A)') 'Incore expansion      '
         end if
         Call GetMem('CCSORT','Allo','Real',iOff_Vic,vsize)
c
         if (ickey.eq.1) then
c:         case V(p,q,r,s)
           call esb_ic_1 (symp,symq,symr,syms,
     &     Work(iOff_Vic),norb(symp),norb(symq),norb(symr),norb(syms))
c
         else if (ickey.eq.2) then
c:       case V(pr,qs)
         Call GetMem('PQIND','ALLO','INTE',ipPQIND,mbas*mbas)
           call  esb_ic_2 (symp,symq,Work(iOff_Vic),
     &     norb(symp),norb(symq),iWork(ipPQIND))
         Call GetMem('PQIND','FREE','INTE',ipPQIND,mbas*mbas)
c
         else
c:       case V(prqs)
         Call GetMem('PQIND','ALLO','INTE',ipPQIND,mbas*mbas)
           call  esb_ic_3 (symp,Work(iOff_Vic),norb(symp),
     &                     iWork(ipPQIND))
         Call GetMem('PQIND','FREE','INTE',ipPQIND,mbas*mbas)
c
         end if
c
       else
c      OUT OF CORE EXPANSION
c*     init temp files and realize expansion of this block
         ickey=0
         if (fullprint.ge.1) then
         write(6,'(6X,A)') 'Out of core expansion '
         end if
       call inittemp (norb(symp))
       if (freespace.lt.(4*nsize*mbas)) then
         write(6,*) ' Allocation of work space for Out-of-core failed!'
         write(6,*) ' Increase the size of the variable MOLCAS_MEM'
         Call Abend
       endif
c
c     allocate space for valn,jn,kn,ln
       Call GetMem('VALN','ALLO','REAL',iOff_valn,nsize*mbas)
       Call GetMem('JN','ALLO','INTE',ipJN,nsize*mbas)
       Call GetMem('KN','ALLO','INTE',ipKN,nsize*mbas)
       Call GetMem('LN','ALLO','INTE',ipLN,nsize*mbas)
c
       call exppsb (symp,symq,symr,syms,Work(iOff_valn),
     &              iWork(ipJN),iWork(ipKN),iWork(ipLN))
c
c     relese space for valn,jn,kn,ln
       Call GetMem('VALN','FREE','REAL',iOff_valn,nsize*mbas)
       Call GetMem('JN','FREE','INTE',ipJN,nsize*mbas)
       Call GetMem('KN','FREE','INTE',ipKN,nsize*mbas)
       Call GetMem('LN','FREE','INTE',ipLN,nsize*mbas)
c
       end if
c
c
c*    run over all pivot indexes
c
       do 500 p=1,norb(symp)
c
c**   def dimensions of vint and read this block of integrals into vint
       ndimvi=norb(symp)
       ndimv1=norb(symq)
       ndimv2=norb(symr)
       ndimv3=norb(syms)
c
       if (ickey.eq.0) then
c      Out of core expansion
c
       if (symq.eq.syms) then
       keyred=1
       else
       keyred=0
       end if
       call unpackk (p,Work(iOff),ndimv1,ndimv2,ndimv3,keyred)
c
       else if (ickey.eq.1) then
c      else Incore expansions
c
c      case V(p,q,r,s)
         call unpackk_ic_1 (p,Work(iOff),ndimv1,ndimv2,ndimv3,
     c                      Work(iOff_Vic),ndimvi)
       else if (ickey.eq.2) then
c      case V(pr,qs)
         call unpackk_ic_2 (p,Work(iOff),ndimvi,ndimv1,Work(iOff_Vic))
       else
c      case V(prqs)
         call unpackk_ic_3 (p,Work(iOff),ndimvi,Work(iOff_Vic))
       end if
c
c
c        goto 500
c
c
c**   add integrals to T3nam if needed (v nacechranej forme)
c
       if (t3key.eq.1) then
       if (p.le.noa(symp)) then
       if (symq.gt.syms) then
c
c***  calc proper address in t3nam file
       t3help4=0
       do t3help1=1,symp-1
       t3help4=t3help4+noa(t3help1)
       end do
       t3help4=t3help4+p
       t3help1=mapiri(symr,symq,1)
       daddr(lunt3)=T3IntPoss(t3help4)+T3Off(t3help1,symp)
c
c***  def required parameters
       t3help1=nvb(symr)
       t3help2=nvb(symq)
       t3help3=nvb(syms)
c
c***  do packing
       call t3intpck2(Work(iOff),Work(iOff+possri0-1),ndimv1,ndimv2,
     & ndimv3,t3help1,t3help2,t3help3,
     & symq,symr,syms,nob,nvb)
c
       else if (symq.eq.syms) then
c
c***  calc proper address in t3nam file
       t3help4=0
       do t3help1=1,symp-1
       t3help4=t3help4+noa(t3help1)
       end do
       t3help4=t3help4+p
       t3help1=mapiri(symr,symq,1)
       daddr(lunt3)=T3IntPoss(t3help4)+T3Off(t3help1,symp)
c
c***  def required parameters
       t3help1=nvb(symr)
       t3help2=nvb(symq)*(nvb(symq)+1)/2
c
c***  do packing
       call t3intpck1(Work(iOff),Work(iOff+possri0-1),ndimv1,ndimv2,
     & ndimv3,t3help1,t3help2,
     & symq,symr,syms,nob,nvb)
c
       end if
       end if
       end if
c
c**   add integrals to #1 <pq|ij> if needed
c
       if (symr.ge.syms) then
c     contributions only for symi(r)>=symj(s)
       call addpqij (Work(iOff),wrksize,
     & symp,symq,symr,syms,p,Work(iOff),ndimv1,ndimv2,ndimv3)
       end if
c
c**   updete fok if neccesarry (only for open shell case)
c
       if (clopkey.eq.1) then
c
       if ((symp.eq.symr).and.(symq.eq.syms)
     & .and.(p.gt.nob(symp)).and.(p.le.(noa(symp)))) then
       call fokupdate1 (foka,fokb,symq,p,Work(iOff),
     &                  ndimv1,ndimv2,ndimv3)
       end if
c
       if ((symp.eq.syms).and.(symq.eq.symr)
     & .and.(p.gt.nob(symp)).and.(p.le.(noa(symp)))) then
       call fokupdate2 (foka,symq,p,Work(iOff),
     &                  ndimv1,ndimv2,ndimv3)
       end if
c
       end if
c
c**   add corresponding <am|pq> interals to TEMPDA2
c     and pack _a_brs to direct acces file TEMPDA1 if need and symm in not C1
       if (p.gt.nob(symp)) then
       a=p-nob(symp)
       call ampack (Work(iOff),wrksize,
     & symp,symq,symr,syms,a,Work(iOff),ndimv1,ndimv2,ndimv3,
     & iWork(ipAMMAP))
       if ((nsym.gt.1).and.(symp.ge.symq)) then
       call abpack (Work(iOff),wrksize,
     & symp,symq,symr,syms,a,Work(iOff),ndimv1,ndimv2,ndimv3,
     & iWork(ipABMAP))
       end if
       end if
c
c**   add INTAB file (for nonsymmetrycal (C1) state)
c
       if ((nsym.eq.1).and.(p.gt.nob(1))) then
       call addintabc1 (Work(iOff),wrksize,
     & p-nob(1),Work(iOff),ndimv1)
       end if
c
 500    continue
c
c        goto 800
c
c
c*    close temp files
c    call closetemp (norb(symp))
c
       if (ickey.ge.1) then
       Call GetMem('CCSORT','Free','Real',iOff_Vic,vsize)
       end if
c
 800    continue
c
c        goto 900
c
c
       if ((nsym.gt.1).and.(symp.ge.symq)) then
c*    add contributions to INTAB comming from symp,sumq and close TEMPDA1 file
c     only for symmetrical cases; only for syma>=symb
       call addintab (Work(iOff),wrksize,
     & symp,symq,iWork(ipABMAP))
       close (lunda1)
       call vf ('TEMPDA1 ',lunda1)
       end if
c
 900    continue
c
c        goto 1000
c
c
c*    add contributions to INTA1-4 if there are some virtuals in symp symmetry
c     and close TEMPDA2 files
c
c     if (nvb(symp).gt.0) then
       call addinta (Work(iOff),wrksize,
     & symp,iWork(ipAMMAP))
       close (lunda2)
       call vf ('TEMPDA2 ',lunda2)
c     end if
c
 1000   continue
c
c*    if T3 are required, reorganize T3nam file
       if (t3key.eq.1) then
       call t3reorg (Work(iOff),wrksize,
     &               noa,nsym)
         do symp=1,4
         end do
       end if
c
c        relese space for ammap,abmap
        Call GetMem('AMMAP','FREE','INTE',ipAMMAP,mbas*64)
        Call GetMem('ABMAP','FREE','INTE',ipABMAP,mbas*mbas*8)
c
c
c*    close files INTA1,INTA2,INTA3 and INTA4, INTAB1
c
       if (iokey.eq.1) then
c      Fortran IO
       close (luna1)
       close (luna2)
       close (luna3)
       close (luna4)
       close (lunab)
c
       else
c      MOLCAS IO
       call daclos (luna1)
       call daclos (luna2)
       call daclos (luna3)
       call daclos (luna4)
       call daclos (lunab)
       end if
c
c        return
c
c
c*    def static integrals (file INTSTA)
c
       call mkintsta (Work(iOff),wrksize,
     & foka,fokb)
c
c*    write general informations to INPDAT
c
       call molcas_binaryopen_vanilla(1,'INPDAT')
c       open (unit=1,file='INPDAT',form='unformatted')
       write (1) NACTEL,ISPIN,NSYM,LSYM,mul,
     &           noa,nob,nva,nvb,norb,eps,Escf
       close (1)
c
c      Release the memory
       Call GetMem('CCSORT','Free','Real',iOff,wrksize)
       Call qExit('Action')

       return
       end
c
c     ----------------------
c
       subroutine mkmapampq (syma)
c
c     this routine prepair mapd,mapi
c     for <am|pq> for given syma, m, p,q to mapd2,mapi2
c
#include "ccsort.fh"
#include "reorg.fh"
       integer syma
c
c     help variables
c
       integer symm,symp,symq,symmp
       integer nhelp,possition,lenght
c
c*    set mapi1 to zero
c
       do 1 symq=1,nsym
       do 2 symp=1,nsym
       do 3 symm=1,nsym
       mapi2(symm,symp,symq)=0
 3      continue
 2      continue
 1      continue
c
c     def zero-th row
c
       mapd2(0,1)=1
       mapd2(0,2)=5
       mapd2(0,3)=5
       mapd2(0,4)=0
       mapd2(0,6)=0

       nhelp=0
       possition=poss20
       do 100 symm=1,nsym
       do 101 symp=1,nsym
       symmp=mul(symm,symp)
       symq=mul(syma,symmp)
       nhelp=nhelp+1
c
c     calc. lenght
       lenght=noa(symm)*NORB(symp)*NORB(symq)
c
       mapd2(nhelp,1)=possition
       mapd2(nhelp,2)=lenght
       mapd2(nhelp,3)=symm
       mapd2(nhelp,4)=symp
       mapd2(nhelp,5)=symq
       mapd2(nhelp,6)=1
       possition=possition+lenght
c
       mapi2(symm,symp,1)=nhelp
c
 101    continue
 100    continue
c
       mapd2(0,5)=nhelp
c
       return
       end
c
c     ----------------------
c
       subroutine mkmappqij
c
c     this routine prepair mapd,mapi
c     for <pq|ij> for p,q, i>=j to mapd1,mapi1
c
#include "ccsort.fh"
#include "reorg.fh"
c
c     help variables
c
       integer symi,symj,symp,symq,sympq,sympqi
       integer nhelp,possition,lenght
c
c*    set mapi1 to zero
c
       do 1 symi=1,nsym
       do 2 symq=1,nsym
       do 3 symp=1,nsym
       mapi1(symp,symq,symi)=0
 3      continue
 2      continue
 1      continue
c
c     def zero-th row
c
       mapd1(0,1)=5
       mapd1(0,2)=5
       mapd1(0,3)=1
       mapd1(0,4)=1
       mapd1(0,6)=3

       nhelp=0
       possition=poss10
       do 100 symp=1,nsym
       do 101 symq=1,nsym
       sympq=mul(symp,symq)
       do 102 symi=1,nsym
       sympqi=mul(sympq,symi)
       symj=sympqi
       if (symj.gt.symi) goto 102
       nhelp=nhelp+1
c
c     calc. lenght
       lenght=noa(symi)*noa(symj)*NORB(symp)*NORB(symq)
c
       mapd1(nhelp,1)=possition
       mapd1(nhelp,2)=lenght
       mapd1(nhelp,3)=symp
       mapd1(nhelp,4)=symq
       mapd1(nhelp,5)=symi
       mapd1(nhelp,6)=symj
       possition=possition+lenght
c
       mapi1(symp,symq,symi)=nhelp
c
 102    continue
 101    continue
 100    continue
c
       mapd1(0,5)=nhelp
c
       return
       end
c
c     ----------------------
c
       subroutine addpqij (wrk,wrksize,
     & symp,symq,symi,symj,p,vint,ndimv1,ndimv2,
     &                     ndimv3)
c
c     this routine add corresponding part to <pq|ij> record (#1)
c     comming from readed integrals with pivot index p vint_p(q,i,j)
c
#include "wrk.fh"
#include "reorg.fh"
#include "ccsort.fh"
       integer symi,symj,symp,symq,p,ndimv1,ndimv2,ndimv3
       real*8 vint(1:ndimv1,1:ndimv2,1:ndimv3)
c
c     help variables
c
       integer ii,ij,i,j,poss0,possij0,q,pqij
c
c*    find number of this symmetry combination
c     and initial possition of this symmetry block in (1)
c
       ii=mapi1(symp,symq,symi)
       poss0=mapd1(ii,1)
c
cT0   if symi<symj return
       if (symi.lt.symj) then
       return
       end if
c
cT1   return, if lenght is 0
       if (mapd1(ii,2).eq.0) then
       return
       end if
c
       do 1000 j=1,noa(symj)
       do 1001 i=1,noa(symi)
c
c*    def ij index and initial possition for <p,q,i,j> integral
c
       ij=(j-1)*noa(symi)+i
       possij0=poss0+(norb(symp)*norb(symq))*(ij-1)
c
       do 200 q=1,norb(symq)
       pqij=possij0-1+norb(symp)*(q-1)+p
       wrk(pqij)=vint(q,i,j)
 200    continue
c
 1001   continue
 1000   continue
c
       return
       end
c
c     -----------------------------------------------
c     -----------------------------------------------
c
       subroutine mkintsta (wrk,wrksize,
     & foka,fokb)
c
c     this routine produces integral file INTSTA, which contains
c     following integrals: foka,fokb,
c     <kl||ij>aaaa,<kl||ij>bbbb,<kl||ij>abab
c     <ka||ij>aaaa,<ka||ij>bbbb,<ka||ij>abab,<ka||ij>baab
c     <ab||ij>aaaa,<ab||ij>bbbb,<ab||ij>abab
c
c     N.B. 1. work file #1 is used for <ij|pq> integrals, #2,3,4
c     must be free. possb0 must be defined
c     N.B. 2. this routine can be used only after definition of <ij|pq>
c     N.B. 3. this routine use followuing help routuines:
c     expandfok
c     wrtmediate (from SYMM)
c
c
#include "wrk.fh"
#include "reorg.fh"
#include "files_ccsd.fh"
       real*8 foka(*)
       real*8 fokb(*)
c
c     help variables
c
       integer rc
c
c*    open INTSTA file
       if (iokey.eq.1) then
c      Fortarn IO
       call molcas_binaryopen_vanilla(lunsta,'INTSTA')
c       open (unit=lunsta,file='INTSTA',form='unformatted')
c
       else
c      MOLCAS IO
       call daname (lunsta,'INTSTA')
       daddr(lunsta)=0
       end if
c
c*    expand foka into work #2 and write to INTSTA
       call expandfok (wrk,wrksize,
     & foka)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c*    expand fokb into work #2 and write to INTSTA
       call expandfok (wrk,wrksize,
     & fokb)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c
c*    get #2 <kl||ij>aaaa from #1 <p,q|i,j> and write to INTSTA
       call exppqij (wrk,wrksize,
     & 4,1,1,1,1,1,1)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c*    get #2 <kl||ij>bbbb from #1 <p,q|i,j> and write to INTSTA
       call exppqij (wrk,wrksize,
     & 4,2,2,2,2,1,1)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c*    get #2 <kl| ij>abab from #1 <p,q|i,j> and write to INTSTA
       call exppqij (wrk,wrksize,
     & 0,1,2,1,2,1,0)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c
c*    get #2 <ka||ij>aaaa from #1 <p,q|i,j> and write to INTSTA
       call exppqij (wrk,wrksize,
     & 3,1,3,1,1,1,1)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c*    get #2 <ka||ij>bbbb from #1 <p,q|i,j> and write to INTSTA
       call exppqij (wrk,wrksize,
     & 3,2,4,2,2,1,1)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c*    get #2 <ka| ij>abab from #1 <p,q|i,j> and write to INTSTA
       call exppqij (wrk,wrksize,
     & 0,1,4,1,2,1,0)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c*    get #2 <ka||ij>baab from #1 <p,q|i,j> and write to INTSTA
       call exppqij (wrk,wrksize,
     & 0,2,3,1,2,0,1)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c
c*    get #2 <ab||ij>aaaa from #1 <p,q|i,j> and write to INTSTA
       call exppqij (wrk,wrksize,
     & 4,3,3,1,1,1,1)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c*    get #2 <ab||ij>bbbb from #1 <p,q|i,j> and write to INTSTA
       call exppqij (wrk,wrksize,
     & 4,4,4,2,2,1,1)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c*    get #2 <ab| ij>abab from #1 <p,q|i,j> and write to INTSTA
       call exppqij (wrk,wrksize,
     & 0,3,4,1,2,1,0)
       call dawrtmediate (wrk,wrksize,
     & lunsta,mapd2,mapi2,rc)
c
c*    close INTSTA file
c
       if (iokey.eq.1) then
c      Fortran IO
       close (lunsta)
c
       else
c      MOLCAS IO
       call daclos (lunsta)
       end if
c
       return
       end
c
c     ---------
c
       subroutine expandfok (wrk,wrksize,
     & fok)
c
c     This routine expand fok operator to #2
c     it also defines new mapd2,mapi2
c
#include "wrk.fh"
#include "reorg.fh"
#include "ccsort.fh"
       real*8 fok(*)
c
c     help variables
c
       integer symp,symq,symr,posstemp,pqwrk,qpwrk,pqfok,p,q
c
c*    set mapi zero
c
       do 1 symr=1,nsym
       do 2 symq=1,nsym
       do 3 symp=1,nsym
       mapi2(symp,symq,symr)=0
 3      continue
 2      continue
 1      continue

c
c*    def zeroth row of mapd
c
       mapd2(0,1)=5
       mapd2(0,2)=5
       mapd2(0,3)=0
       mapd2(0,4)=0
       mapd2(0,5)=nsym
       mapd2(0,6)=0
c
       posstemp=poss20
       pqfok=0
       do 1000 symp=1,nsym
c
c*    def mapd,mapi
c
       mapd2(symp,1)=posstemp
       mapd2(symp,2)=norb(symp)*norb(symp)
       mapd2(symp,3)=symp
       mapd2(symp,4)=symp
       mapd2(symp,5)=1
       mapd2(symp,6)=1
       mapi2(symp,1,1)=symp
c
c*    expand
c
       do 100 p=1,norb(symp)
       do 101 q=1,p
c
c*    calc pq and qp possition in work and fok
c     and write integrals to this possitions
c
       pqwrk=posstemp+(norb(symp)*(p-1)+q)-1
       qpwrk=posstemp+(norb(symp)*(q-1)+p)-1
       pqfok=pqfok+1
       wrk(pqwrk)=fok(pqfok)
       wrk(qpwrk)=fok(pqfok)
c
 101    continue
 100    continue
c
       posstemp=posstemp+mapd2(symp,2)
c
 1000   continue
c
       return
       end
c
c     -----------------------------------------------
c     -----------------------------------------------
c
       subroutine addinta (wrk,wrksize,
     & syma,ammap)
c
c     this routine do for all a in syma
c     1- reconstruct #2 <_a,m,p,q> from TEMPDA2 file
c     2- prepair corresponding <_am p q> (like <amef>aaaa) to #3
c     and wrirte it to opened INTA1-4
c     N.B.  this routine use followuing foreign routuines:
c     wrtmap
c     wri
c
#include "wrk.fh"
#include "reorg.fh"
#include "ccsort.fh"
       integer syma
       integer ammap(1:mbas,1:8,1:8)
c
c     help variables
c
       integer lenefaaaa,lenefbaab,lenefbbbb,lenefabab
       integer lenejaaaa,lenejbaab,lenejbaba,lenejbbbb,lenejabab,
     & lenejabba
       integer posst,rc,a
c
c*    mapd2 and mapi2 of #2 <_a,m|p,q> are prepaired
c
c*    make required mapd3 and mapi3 and write them to INTA1-4
c     define lenghts of this mediates
c
c*1   to INTA1 <m,_a||ef>aaaa, <m,_a||ef>baab
       call ccsort_grc0(3,2,1,3,3,0,syma,poss30,posst,mapd3,mapi3)
       call deflenght (mapd3,lenefaaaa)
       call dawrtmap (luna1,mapd3,mapi3,rc)
       call ccsort_grc0(3,0,2,3,4,0,syma,poss30,posst,mapd3,mapi3)
       call deflenght (mapd3,lenefbaab)
       call dawrtmap (luna1,mapd3,mapi3,rc)
c
c*2   to INTA2 <m,_a||ef>bbbb, <m,_a||ef>abab
       call ccsort_grc0(3,2,2,4,4,0,syma,poss30,posst,mapd3,mapi3)
       call deflenght (mapd3,lenefbbbb)
       call dawrtmap (luna2,mapd3,mapi3,rc)
       call ccsort_grc0(3,0,1,3,4,0,syma,poss30,posst,mapd3,mapi3)
       call deflenght (mapd3,lenefabab)
       call dawrtmap (luna2,mapd3,mapi3,rc)
c
c*3   to INTA3 <m,_a||ej>aaaa, <m,_a||ej>baab, <m,_a||ej>baba
       call ccsort_grc0(3,0,1,3,1,0,syma,poss30,posst,mapd3,mapi3)
       call deflenght (mapd3,lenejaaaa)
       call dawrtmap (luna3,mapd3,mapi3,rc)
       call ccsort_grc0(3,0,2,3,2,0,syma,poss30,posst,mapd3,mapi3)
       call deflenght (mapd3,lenejbaab)
       call dawrtmap (luna3,mapd3,mapi3,rc)
       call ccsort_grc0(3,0,2,4,1,0,syma,poss30,posst,mapd3,mapi3)
       call deflenght (mapd3,lenejbaba)
       call dawrtmap (luna3,mapd3,mapi3,rc)
c
c*4   to INTA4 <m,_a||ej>bbbb, <m,_a||ej>abba, <m,_a||ej>abab
       call ccsort_grc0(3,0,2,4,2,0,syma,poss30,posst,mapd3,mapi3)
       call deflenght (mapd3,lenejbbbb)
       call dawrtmap (luna4,mapd3,mapi3,rc)
       call ccsort_grc0(3,0,1,4,1,0,syma,poss30,posst,mapd3,mapi3)
       call deflenght (mapd3,lenejabba)
       call dawrtmap (luna4,mapd3,mapi3,rc)
       call ccsort_grc0(3,0,1,3,2,0,syma,poss30,posst,mapd3,mapi3)
       call deflenght (mapd3,lenejabab)
       call dawrtmap (luna4,mapd3,mapi3,rc)
c
c
c*    cycle over a
c
       do 1000 a=1,nvb(syma)
c
c*    reconstruct #2 <_a,m,p,q> for given _a
       call mkampq (wrk,wrksize,
     & a,syma,ammap)
c
c*    get contributions to INTA2 <m,_a||ef>bbbb, <m,_a||ef>abab
c     and wtite it there
c
       if (lenefbbbb.gt.0) then
       call expmpq (wrk,wrksize,
     & syma,2,2,4,4,1,1)
       call dawri (luna2,lenefbbbb,wrk(mapd3(1,1)))
       end if
c
       if (lenefabab.gt.0) then
       call expmpq (wrk,wrksize,
     & syma,0,1,3,4,1,0)
       call dawri (luna2,lenefabab,wrk(mapd3(1,1)))
       end if
c
c*    get contributions to INTA4 <m,_a||ej>bbbb, <m,_a||ej>abba, <m,_a||ej>abab
c     and wtite it there
c
       if (lenejbbbb.gt.0) then
       call expmpq (wrk,wrksize,
     & syma,0,2,4,2,1,1)
       call dawri (luna4,lenejbbbb,wrk(mapd3(1,1)))
       end if
c
       if (lenejabba.gt.0) then
       call expmpq (wrk,wrksize,
     & syma,0,1,4,1,0,1)
       call dawri (luna4,lenejabba,wrk(mapd3(1,1)))
       end if
c
       if (lenejabab.gt.0) then
       call expmpq (wrk,wrksize,
     & syma,0,1,3,2,1,0)
       call dawri (luna4,lenejabab,wrk(mapd3(1,1)))
       end if
c
       if (a.gt.(nvb(syma)-nva(syma))) then
c     contributions to INTA1 and INTA3 only for a-alfa
c
c*    get contributions to INTA1 <m,_a||ef>aaaa, <m,_a||ef>baab if any
c     and wtite it there
c
       if (lenefaaaa.gt.0) then
       call expmpq (wrk,wrksize,
     & syma,2,1,3,3,1,1)
       call dawri (luna1,lenefaaaa,wrk(mapd3(1,1)))
       end if
c
       if (lenefbaab.gt.0) then
       call expmpq (wrk,wrksize,
     & syma,0,2,3,4,0,1)
       call dawri (luna1,lenefbaab,wrk(mapd3(1,1)))
       end if
c
c*    get contributions to INTA3 <m,_a||ej>aaaa, <m,_a||ej>baab, <m,_a||ej>baba
c     and wtite it there
c
       if (lenejaaaa.gt.0) then
       call expmpq (wrk,wrksize,
     & syma,0,1,3,1,1,1)
       call dawri (luna3,lenejaaaa,wrk(mapd3(1,1)))
       end if
c
       if (lenejbaab.gt.0) then
       call expmpq (wrk,wrksize,
     & syma,0,2,3,2,0,1)
       call dawri (luna3,lenejbaab,wrk(mapd3(1,1)))
       end if
c
       if (lenejbaba.gt.0) then
       call expmpq (wrk,wrksize,
     & syma,0,2,4,1,1,0)
       call dawri (luna3,lenejbaba,wrk(mapd3(1,1)))
       end if
c
       end if
c
 1000   continue
c
       return
       end
c
c     -------------------------
c
       subroutine deflenght (mapd,lenght)
c
c     this routine defines lenght of mediate, described by mapd
c
       integer mapd(0:512,1:6)
       integer lenght
c
c     help variable
c
       integer ii
c
       ii=mapd(0,5)
       lenght=mapd(ii,1)+mapd(ii,2)-mapd(1,1)
c
       return
       end
c
c     -------------------------------------------------
c
       subroutine fokupdate1 (foka,fokb,symp,i,vint,ndimv1,ndimv2,
     &                        ndimv3)
c
c     this routine realize update
c     foka(p,q) = foka(p,q) + <ip|iq>
c     fokb(p,q) = fokb(p,q) + <ip|iq>
c
c     N.B. integrals are of type <symi, symp| symi,symp>
c
c     foka    - packed Fokaa matrix (I,O)
c     fokb    - packed Fokbb matrix (I,O)
c     symp    - irrep or p (and also q) index (I)
c     i       - value of i, (I)
c     vint    - array of integrals <ip|iq> for given i (I)
c     ndimv1  - first dimension (norb(symp)) (I)
c     ndimv2  - second dimension (norb(symi)) (I)
c     ndimv3  - third dimension (norb(symp)) (I)
c
#include "ccsort.fh"
       real*8 foka(*)
       real*8 fokb(*)
       real*8 vint(1:ndimv1,1:ndimv2,1:ndimv3)
       integer symp,i,ndimv1,ndimv2,ndimv3
c
c     help variables
c
       integer nhelp1,nhelp2,p,q,pq
c
c*    calculate shift
c
       nhelp1=0
       if (symp.gt.1) then
       do 100 nhelp2=1,symp-1
       nhelp1=nhelp1+(norb(nhelp2)**2+norb(nhelp2))/2
 100    continue
       end if
c
c*    add integral
c
       pq=nhelp1
       do 200 p=1,norb(symp)
       do 201 q=1,p
       pq=pq+1
       foka(pq)=foka(pq)+vint(p,i,q)
       fokb(pq)=fokb(pq)+vint(p,i,q)
 201    continue
 200    continue
c
       return
       end
c
c     -------------------------------------------------
c
       subroutine fokupdate2 (foka,symp,i,vint,ndimv1,ndimv2,ndimv3)
c
c     this routine realize update
c     foka(p,q) = foka(p,q) - <ip|qi>
c
c     N.B. integrals are of type <symi, symp| symp, symi>
c
c     foka    - packed Fokaa matrix (I,O)
c     symp    - irrep or p (and also q) index (I)
c     i       - value of i, (I)
c     vint    - array of integrals <ip|iq> for given i (I)
c     ndimv1  - first dimension (norb(symp)) (I)
c     ndimv2  - second dimension (norb(symi)) (I)
c     ndimv3  - third dimension (norb(symp)) (I)
c
#include "ccsort.fh"
       real*8 foka(*)
       real*8 vint(1:ndimv1,1:ndimv2,1:ndimv3)
       integer symp,i,ndimv1,ndimv2,ndimv3
c
c     help variables
c
       integer nhelp1,nhelp2,p,q,pq
c
c*    calculate shift
c
       nhelp1=0
       if (symp.gt.1) then
       do 100 nhelp2=1,symp-1
       nhelp1=nhelp1+(norb(nhelp2)**2+norb(nhelp2))/2
 100    continue
       end if
c
c*    add integral
c
       pq=nhelp1
       do 200 p=1,norb(symp)
       do 201 q=1,p
       pq=pq+1
       foka(pq)=foka(pq)-vint(p,q,i)
 201    continue
 200    continue
c
       return
       end
c
c     -------------------------------------------------
c
       subroutine unpackk_ic_3 (i,vint,ndimvi,Vic)
c
c     this routine vint(j,k,l) = <i,j|k,l>
c     for given i from incore (reduced) expanded block Vic
c     ie. symp=symq=symr=syms
c
c     i      - value of pivot index (I)
c     vint   - array of integrals (O)
c     ndimvi - (norb(symi)) (I)
c     Vic    - incore expanded block of integrals (I)
c
#include "reorg.fh"

#include "SysDef.fh"
       integer i,ndimvi
       real*8 vint(1:ndimvi,1:ndimvi,1:ndimvi)
       real*8 Vic(1:(ndimvi*(ndimvi+1)/2)*(1+ndimvi*(ndimvi+1)/2)/2)
c
c     help variables
c
      integer j,k,l,ik,jl,ikjl
c
c
        do k=1,ndimvi
c
c        def ik
        if (i.ge.k) then
          ik=i*(i-1)/2+k
        else
          ik=k*(k-1)/2+i
        end if
c
          jl=0
        do j=1,ndimvi
          do l=1,j
c
c         def jl
          jl=jl+1
c
c           def ikjl
            if (ik.ge.jl) then
            ikjl=ik*(ik-1)/2+jl
            else
            ikjl=jl*(jl-1)/2+ik
            end if
c
            vint(j,k,l)=Vic(ikjl)
            vint(l,k,j)=Vic(ikjl)
c
          end do
        end do
        end do
c
       return
       end
c
c     -------------------------------------------------
c
       subroutine unpackk_ic_2 (i,vint,ndimvi,ndimvj,Vic)
c
c     this routine vint(j,k,l) = <i,j|k,l>
c     for given i from incore (reduced) expanded block Vic
c     ie. symp=symr,symq=syms
c
c     i      - value of pivot index (I)
c     vint   - array of integrals (O)
c     ndimvi - (norb(symi)) (I)
c     ndimvj - (norb(symj)) (I)
c     Vic    - incore expanded block of integrals (I)
c
#include "reorg.fh"

#include "SysDef.fh"
       integer i,ndimvi,ndimvj
       real*8 vint(1:ndimvj,1:ndimvi,1:ndimvj)
       real*8 Vic(1:(ndimvi*(ndimvi+1)/2),1:(ndimvj*(ndimvj+1)/2))
c
c     help variables
c
      integer j,k,l,ik,jl
c
c
        do k=1,ndimvi
        if (i.ge.k) then
          ik=i*(i-1)/2+k
        else
          ik=k*(k-1)/2+i
        end if
          jl=0
          do j=1,ndimvj
          do l=1,j
            jl=jl+1
            vint(j,k,l)=Vic(ik,jl)
            vint(l,k,j)=Vic(ik,jl)
          end do
          end do
        end do
c
       return
       end
c
c     -------------------------------------------------
c
       subroutine unpackk_ic_1 (i,vint,ndimv1,ndimv2,ndimv3,
     c                       Vic,ndimvi)
c
c     this routine vint(j,k,l) = <i,j|k,l>
c     for given i from incore (nonreduced) expanded block Vic
c
c     i      - value of pivot index (I)
c     vint   - array of integrals (O)
c     ndimv1 - first dimension of vint (norb(symj)) (I)
c     ndimv2 - second dimension of vint (norb(symk)) (I)
c     ndimv3 - third dimension of vint (norb(syml)) (I)
c     Vic    - incore expanded block of integrals (I)
c     ndimvi - first dimension of Vic norb(symi) (I)
c
#include "reorg.fh"

#include "SysDef.fh"
       integer i,ndimv1,ndimv2,ndimv3,ndimvi
       real*8 vint(1:ndimv1*ndimv2*ndimv3)
       real*8 Vic(1:ndimvi,1:ndimv1*ndimv2*ndimv3)
c
c     help variables
c
      integer jkl
c
c
        do jkl=1,ndimv1*ndimv2*ndimv3
        vint(jkl)=Vic(i,jkl)
        end do
c
       return
       end
c
c     -------------------------------------------------
c
       subroutine unpackk (i,vint,ndimv1,ndimv2,ndimv3,key)
c
c      unpackk process control routine
c
c     i      - value of pivot index (I)
c     vint   - array of integrals (O)
c     ndimv1 - first dimension of vint (norb(symj)) (I)
c     ndimv2 - second dimension of vint (norb(symk)) (I)
c     ndimv3 - third dimension of vint (norb(syml)) (I)
c     key    - reduced storing key (I)
c     = 0 if symj is not syml
c     = 1 if symj = syml
c
#include "reorg.fh"
       integer i,ndimv1,ndimv2,ndimv3,key
       real*8 vint(1:ndimv1,1:ndimv2,1:ndimv3)
c
       if (zrkey.eq.1) then
         call unpackk_zr (i,vint,ndimv1,ndimv2,ndimv3,key)
       else
         call unpackk_pck (i,vint,ndimv1,ndimv2,ndimv3,key)
       end if
c
       return
       end
c
c     -------------------------------------------------
c
       subroutine unpackk_zr (i,vint,ndimv1,ndimv2,ndimv3,key)
c
c     this routine expand integrals packed in i-th TEMP file
c     to vint(j,k,l) = <i,j|k,l>
c
c     i      - value of pivot index (I)
c     vint   - array of integrals (O)
c     ndimv1 - first dimension of vint (norb(symj)) (I)
c     ndimv2 - second dimension of vint (norb(symk)) (I)
c     ndimv3 - third dimension of vint (norb(syml)) (I)
c     key    - reduced storing key (I)
c     = 0 if symj is not syml
c     = 1 if symj = syml
c
#include "reorg.fh"

#include "SysDef.fh"
       integer i,ndimv1,ndimv2,ndimv3,key
       real*8 vint(1:ndimv1,1:ndimv2,1:ndimv3)
c
c     help variables
c
       integer nhelp,lenght,daddr,nrec
c
       integer constj
       parameter (constj=1048576)
       integer constk
       parameter (constk=1024)
c
       integer ihelp,ires
       Integer iBuf(1:nsize)
c
c*    set vint=0
c
       nhelp=ndimv1*ndimv2*ndimv3
       call ccsort_mv0zero (nhelp,nhelp,vint)
c
c*     open corresponding TEMP file
c
       if (iokey.eq.1) then
c      Fortran IO
       call molcas_binaryopen_vanilla(lunpublic,tmpnam(i))
c       open (unit=lunpublic,file=tmpnam(i),form='unformatted')
       else
c      MOLCAS IO
       call daname (lunpublic,tmpnam(i))
       daddr=0
       end if
c
      do nrec=1,nrectemp(i)
c
        if (nrec.ne.nrectemp(i)) then
        lenght=nsize
        else
        lenght=lrectemp(i)
        end if
c
      if (iokey.eq.1) then
c     Fortran IO
      call getpp_zr (lunpublic,valh,iBuf,lenght)
      else
c     MOLCAS IO
      call ddafile (lunpublic,2,valh,lenght,daddr)
      call idafile (lunpublic,2,iBuf,lenght,daddr)
      end if
c
c*     get indexes jh,kh,lh and value valh from packed form
c
      do nhelp=1,lenght
         ihelp=iBuf(nhelp)
         jh(nhelp)=int(ihelp/constj)
         ires=ihelp-constj*jh(nhelp)
         kh(nhelp)=int(ires/constk)
         lh(nhelp)=ires-constk*kh(nhelp)
      end do
c
      if (key.eq.0) then
      do 100 nhelp=1,lenght
      vint(jh(nhelp),kh(nhelp),lh(nhelp))=valh(nhelp)
 100  continue
      else
      do 200 nhelp=1,lenght
      vint(jh(nhelp),kh(nhelp),lh(nhelp))=valh(nhelp)
      vint(lh(nhelp),kh(nhelp),jh(nhelp))=valh(nhelp)
 200  continue
      end if
c
       end do
c
      if (iokey.eq.1) then
c     Fortran IO
       close (lunpublic)
      else
c     Molcas IO
       call daclos (lunpublic)
      end if
c
       return
       end
c
c     -------------------------------------------------
c
       subroutine unpackk_pck (i,vint,ndimv1,ndimv2,ndimv3,key)
c
c     this routine expand integrals packed in i-th TEMP file
c     to vint(j,k,l) = <i,j|k,l>
c
c     i      - value of pivot index (I)
c     vint   - array of integrals (O)
c     ndimv1 - first dimension of vint (norb(symj)) (I)
c     ndimv2 - second dimension of vint (norb(symk)) (I)
c     ndimv3 - third dimension of vint (norb(syml)) (I)
c     key    - reduced storing key (I)
c     = 0 if symj is not syml
c     = 1 if symj = syml
c
#include "reorg.fh"

#include "SysDef.fh"
       integer i,ndimv1,ndimv2,ndimv3,key
       real*8 vint(1:ndimv1,1:ndimv2,1:ndimv3)
c
c     help variables
c
       integer nhelp,lenght,daddr,nrec
c
       integer constj
       parameter (constj=1048576)
       integer constk
       parameter (constk=1024)
c
       character*(RtoB+ItoB) pp(1:nsize),pphelp
       character*1 p1help(1:(RtoB+ItoB))
       real*8 rhelp
       integer ihelp,ires
       equivalence (p1help(1),pphelp)
       equivalence (p1help(1),rhelp)
       equivalence (p1help(1+RtoB),ihelp)
c
c*    set vint=0
c
       nhelp=ndimv1*ndimv2*ndimv3
       call ccsort_mv0zero (nhelp,nhelp,vint)
c
c*     open corresponding TEMP file
c
       if (iokey.eq.1) then
c      Fortran IO
       call molcas_binaryopen_vanilla(lunpublic,tmpnam(i))
c       open (unit=lunpublic,file=tmpnam(i),form='unformatted')
       else
c      MOLCAS IO
       call daname (lunpublic,tmpnam(i))
       daddr=0
       end if
c
      do nrec=1,nrectemp(i)
c
        if (nrec.ne.nrectemp(i)) then
        lenght=nsize
        else
        lenght=lrectemp(i)
        end if
c
      if (iokey.eq.1) then
c     Fortran IO
      call getpp_pck (lunpublic,pp,lenght)
      else
c     MOLCAS IO
      call cdafile (lunpublic,2,pp,(RtoB+ItoB)*lenght,daddr)
      end if

c
c*     get indexes jh,kh,lh and value valh from packed form
c
      do nhelp=1,lenght
      pphelp=pp(nhelp)
      valh(nhelp)=rhelp
      jh(nhelp)=int(ihelp/constj)
      ires=ihelp-constj*jh(nhelp)
      kh(nhelp)=int(ires/constk)
      lh(nhelp)=ires-constk*kh(nhelp)
      end do
c
      if (key.eq.0) then
      do 100 nhelp=1,lenght
      vint(jh(nhelp),kh(nhelp),lh(nhelp))=valh(nhelp)
 100  continue
      else
      do 200 nhelp=1,lenght
      vint(jh(nhelp),kh(nhelp),lh(nhelp))=valh(nhelp)
      vint(lh(nhelp),kh(nhelp),jh(nhelp))=valh(nhelp)
 200  continue
      end if
c
       end do
c
      if (iokey.eq.1) then
c     Fortran IO
       close (lunpublic)
      else
c     Molcas IO
       call daclos (lunpublic)
      end if
c
       return
       end
c
c     ----------------------
c
       subroutine getpp_zr (lunpublic,pp,ipp,lenght)
c
#include "SysDef.fh"
       integer lunpublic,lenght
       Real*8 pp(1:lenght)
       Integer ipp(1:lenght)
c
       read (lunpublic) pp,ipp
c
       return
       end
c
c     ----------------------
c
       subroutine getpp_pck (lunpublic,pp,lenght)
c

#include "SysDef.fh"
       integer lunpublic,lenght
       character*(RtoB+ItoB) pp(1:lenght)
c
       read (lunpublic) pp
c
       return
       end
c
c     -------------------------------------------------
c
       subroutine initintabc1
c     this routine write corresponding mapd and mapi to INTAB
c     for nonsymetrical (C1) case
c
#include "reorg.fh"
#include "ccsort.fh"
c
c     help variables
c
       integer nhelp,lenght,symp,symq,symab
       integer poss,ii,syma,symb,rc
c
c*    def symab
       syma=1
       symb=1
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
c
       return
       end
c
c     ------------------------------------
c
       subroutine addintabc1 (wrk,wrksize,
     & a,vint,ndimv)
c
c     this routine add integrals <_a,_b|p,q> for given a
c     for nonsymmetrical (C1) case
c     from integrals vv _a(u,p,q)
c
#include "wrk.fh"
#include "reorg.fh"
#include "ccsort.fh"
       integer a,ndimv
       real*8 vint(1:ndimv,1:ndimv,1:ndimv)
c
c     help variables
c
       integer poss,b,bvint,p,q,lenght
c
c
cT    if there are no _a_b,pq integrals in this symab,
c     skip sumation over ab
c
       if (nvb(1).eq.0) then
       return
       end if
c
c*    loop over b
c
       do 1000 b=1,a
       bvint=b+nob(1)
c
c     map <_a,b|p,q> to wrk in #3
       poss=poss30
       do 1010 q=1,norb(1)
       do 1011 p=1,norb(1)
       wrk(poss)=vint(bvint,p,q)
       poss=poss+1
 1011   continue
 1010   continue
c
c
c**   since there must be some integrals, write them to TEMPAB
c
       lenght=poss-poss30
       call dawri (lunab,lenght,wrk(poss30))
c
 1000   continue
c
       return
       end
c
c     -------------------------
c
       subroutine initwrk (lenght)
c
c      this routine calculate required size of work space and
c      definie initial possitions of work vectors
c
#include "ccsort.fh"
#include "reorg.fh"
c
c     help variables
c
       integer n
       integer sizevint,sizev1,sizev2,sizempq,lenght,norbmax,sizeri
       integer symp,symq,symi,symj,sympq,sympqi,symm,symmp,syma
c
c1*   def maxzie of vint
c
       norbmax=norb(1)
       do 10 n=1,nsym
       if (norb(n).gt.norbmax) then
       norbmax=norb(n)
       end if
 10     continue
c
       sizevint=norbmax*norbmax*norbmax
c
c2*   def size of <pq|i>=j>, <pq|i,j>
c
       sizev1=0
       sizev2=0
       do 20 symp=1,nsym
       do 21 symq=1,nsym
       sympq=mul(symp,symq)
       do 22 symi=1,nsym
       sympqi=mul(sympq,symi)
       symj=sympqi
c      calc. lenght
       if (symj.gt.symi) then
         sizev2=sizev2+noa(symi)*noa(symj)*NORB(symp)*NORB(symq)
       else
         sizev1=sizev1+noa(symi)*noa(symj)*NORB(symp)*NORB(symq)
         sizev2=sizev2+noa(symi)*noa(symj)*NORB(symp)*NORB(symq)
       end if
22     continue
21     continue
20     continue
c
c3*   def maxsize of <_am|pq>
c

       sizempq=0
       do 50 syma=1,nsym
c
       lenght=0
       do 30 symm=1,nsym
       do 31 symp=1,nsym
       symmp=mul(symm,symp)
       symq=mul(syma,symmp)
c     calc. lenght
       lenght=lenght+noa(symm)*NORB(symp)*NORB(symq)
 31     continue
 30     continue
c
       if (sizempq.lt.lenght) then
       sizempq=lenght
       end if
c
 50     continue
c
c4*   def maxsize of R_i if needed
c
       sizeri=0
c
       if (t3key.eq.1) then
       do 60 symi=1,nsym
       call ccsort_t3grc0 (3,8,4,4,4,0,symi,1,lenght,mapdri,mapiri)
       lenght=lenght-1
       if (lenght.gt.sizeri) then
       sizeri=lenght
       end if
 60     continue
       end if

c     ******* distribution of memory ******
c
       poss10=1+sizevint
       poss20=poss10+sizev1
       poss30=poss20+sizev2
       possri0=poss30+sizempq
       lenght=possri0+sizeri-1
c
       if (fullprint.gt.1) then
       write(6,*)
       write(6,'(6X,A)') 'size of help (work) vectors:'
       write(6,'(6X,A)') '----------------------------'
       write(6,*)
       write(6,'(6X,A,I8)') 'Vints     V0 required : ',sizevint
       write(6,'(6X,A,I8)') 'PQIJ ints V1 required : ',sizev1
       write(6,'(6X,A,I8)') '          V2 required : ',sizev2
       write(6,'(6X,A,I8)') 'AMIJ ints V3 required : ',sizempq
       write(6,'(6X,A,I8)') 'R_i mtx   Ri required : ',sizeri
       end if
c
       if (fullprint.ge.0)
     &     write(6,'(6X,A,I20)') 'Required WRK size-sum : ',lenght
c
       return
       end
