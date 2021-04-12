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
       subroutine action (foka, fokb, fi,eps)
c
       implicit real*8 (a-h,o-z)
c
c      work file declaration
       integer wrksize
#include "WrkSpc.fh"
#include "stdalloc.fh"
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
       integer ipJN,ipKN,ipLN,iOff_valn
       integer, allocatable :: AMMAP(:,:,:),ABMAP(:,:,:)
       real*8 foka((mbas**2+mbas)/2)
       real*8 fokb((mbas**2+mbas)/2)


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
       Call mma_Allocate(AMMAP,mbas,8,8,Label='AMMAP')
       Call mma_Allocate(ABMAP,mbas,mbas,8,Label='ABMAP')
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
       call mkampqmap (AMMAP,symp,rc)
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
       call mkabpqmap (ABMAP,symp,symq,rc)
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
c     release space for valn,jn,kn,ln
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
     & symp,symq,symr,syms,a,Work(iOff),ndimv1,ndimv2,ndimv3,AMMAP)
       if ((nsym.gt.1).and.(symp.ge.symq)) then
       call abpack (Work(iOff),wrksize,
     & symp,symq,symr,syms,a,Work(iOff),ndimv1,ndimv2,ndimv3,ABMAP)
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
       call addintab (Work(iOff),wrksize,symp,symq,ABMAP)
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
       call addinta (Work(iOff),wrksize,symp,AMMAP)
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
c        release space for ammap,abmap
        Call mma_Deallocate(AMMAP)
        Call mma_Deallocate(ABMAP)
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

       return
       end
