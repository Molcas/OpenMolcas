!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
       subroutine action (foka, fokb, fi,eps)
!
       implicit real*8 (a-h,o-z)
!
!      work file declaration
       integer wrksize
#include "WrkSpc.fh"
#include "stdalloc.fh"
!
#include "ccsort.fh"
#include "reorg.fh"
#include "files_ccsd.fh"
       real*8 fi(*)
       real*8 eps(mbas)
!
!     help variables
!
       integer symp,symq,symr,syms,sympq,sympqr
       integer ndimv1,ndimv2,ndimv3,ndimvi
       integer p,a,posst,rc
       integer keyred,vsize,freespace,ickey,iOff_Vic,iOff
       integer t3help1,t3help2,t3help3,t3help4
       integer ipJN,ipKN,ipLN,iOff_valn
       integer, allocatable :: AMMAP(:,:,:),ABMAP(:,:,:)
       real*8 foka((mbas**2+mbas)/2)
       real*8 fokb((mbas**2+mbas)/2)


!
!*    distribute memory
!
!*.1   calc. work space requirements
       call initwrk (wrksize)
!
!*.2   test if allocation of required memory is possible
!
!*.2.1 allocate work space
!
       Call GetMem('CCSORT','Max','Real',maxspace,maxspace)
       maxspace=maxspace-4
       if (maxspace.lt.wrksize) then
         write(6,*) ' Allocation of work space failed!'
         write(6,*) ' Increase the size of the variable MOLCAS_MEM'
         Call Abend
       end if
       Call GetMem('CCSORT','Allo','Real',iOff,wrksize)
!
!*.3   set wrk = 0
       call ccsort_mv0zero (wrksize,wrksize,Work(iOff))
!
!
!*    def foka,fokb
!
       ndimv1=0
       do 100 symp=1,nsym
       ndimv1=ndimv1+(norb(symp)+1)*norb(symp)/2
 100    continue
!
       do 200 ndimv2=1,ndimv1
       foka(ndimv2)=fi(ndimv2)
       fokb(ndimv2)=fi(ndimv2)
 200    continue
!
!*    make names
!
       call mktempanam
!
!*    if T3 are requited, define T3IntPoss and calc T3Off
       if (t3key.eq.1) then
         call DefT3par (noa,nsym)
         do symp=1,4
         end do
       end if
!
!*    open files INTA1,INTA2,INTA3,INTA4 and INTAB
!
       if (iokey.eq.1) then
!      Fortran IO
       call molcas_binaryopen_vanilla(luna1,'INTA1')
       call molcas_binaryopen_vanilla(luna2,'INTA2')
       call molcas_binaryopen_vanilla(luna3,'INTA3')
       call molcas_binaryopen_vanilla(luna4,'INTA4')
       call molcas_binaryopen_vanilla(lunab,'INTAB')
!       open (unit=luna1,file='INTA1',form='unformatted')
!       open (unit=luna2,file='INTA2',form='unformatted')
!       open (unit=luna3,file='INTA3',form='unformatted')
!       open (unit=luna4,file='INTA4',form='unformatted')
!       open (unit=lunab,file='INTAB',form='unformatted')
!
       else
!      MOLCAS IO
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
!
!*    define #V1 (for <p,q,i,j>
       call mkmappqij
!
!*    make head of INTAB file for nonsymmetrical (C1) case
       if (nsym.eq.1) then
       call initintabc1
       end if
!
!     allocate space for ammap,abmap
       Call mma_Allocate(AMMAP,mbas,8,8,Label='AMMAP')
       Call mma_Allocate(ABMAP,mbas,mbas,8,Label='ABMAP')
!
       do 1000 symp=1,nsym
!
       if (fullprint.gt.0) then
       write(6,'(6X,A,2X,I2)') 'Symmetry of the pivot index',symp
       end if
!
!*    define #V2 (for <_a,m,p,q>)
       call mkmapampq (symp)
!
!*    if T3 are required, make maps for R_i
       if ((t3key.eq.1).and.(noa(symp).gt.0)) then
!*            get mapd and mapi for  R_i(a,bc)
        call ccsort_t3grc0(3,8,4,4,4,0,symp,possri0,posst,mapdri,mapiri)
       end if
!
!*    open TEMPDA2 fils for <am|rs> integrals, if there are some virtiuals
!     in symp symmetry
       if (nvb(symp).gt.0) then
       call mkampqmap (AMMAP,symp,rc)
       call daopen ('TEMPDA2 ',lunda2,recl,1)
       end if
!
       do 900 symq=1,nsym
       sympq=mul(symp,symq)
!
       if ((nsym.gt.1).and.(symp.ge.symq)) then
!*    open direct acces file here, to enable exact specification of
!     the number of records (only for symmetrical cases; syma>=symb)
!     N.B. nrec is not needed now
       call daopen ('TEMPDA1 ',lunda1,recl,1)
       end if
!
!*    make abmap for syma>=symb
       if (symp.ge.symq) then
       call mkabpqmap (ABMAP,symp,symq,rc)
       end if
!
       do 800 symr=1,nsym
       sympqr=mul(sympq,symr)
       syms=sympqr
!
!*     calc size of the integral file
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
!
       if (vsize.eq.0) goto 800
!
       if (fullprint.gt.1) then
         write(6,'(6X,A,I4,4X,4I2)') 'Block',typ(symp,symq,symr),       &
     &                                symp,symq,symr,syms
       end if
!
!*     test for incore expansion
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
!
       if (freespace.ge.(vsize+mbas*mbas)) then
!      INCORE EXPANSION
         if (fullprint.ge.1) then
         write(6,*)
         write(6,'(6X,A)') 'Incore expansion      '
         end if
         Call GetMem('CCSORT','Allo','Real',iOff_Vic,vsize)
!
         if (ickey.eq.1) then
!:         case V(p,q,r,s)
           call esb_ic_1 (symp,symq,symr,syms,                          &
     &     Work(iOff_Vic),norb(symp),norb(symq),norb(symr),norb(syms))
!
         else if (ickey.eq.2) then
!:       case V(pr,qs)
         Call GetMem('PQIND','ALLO','INTE',ipPQIND,mbas*mbas)
           call  esb_ic_2 (symp,symq,Work(iOff_Vic),                    &
     &     norb(symp),norb(symq),iWork(ipPQIND))
         Call GetMem('PQIND','FREE','INTE',ipPQIND,mbas*mbas)
!
         else
!:       case V(prqs)
         Call GetMem('PQIND','ALLO','INTE',ipPQIND,mbas*mbas)
           call  esb_ic_3 (symp,Work(iOff_Vic),norb(symp),              &
     &                     iWork(ipPQIND))
         Call GetMem('PQIND','FREE','INTE',ipPQIND,mbas*mbas)
!
         end if
!
       else
!      OUT OF CORE EXPANSION
!*     init temp files and realize expansion of this block
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
!
!     allocate space for valn,jn,kn,ln
       Call GetMem('VALN','ALLO','REAL',iOff_valn,nsize*mbas)
       Call GetMem('JN','ALLO','INTE',ipJN,nsize*mbas)
       Call GetMem('KN','ALLO','INTE',ipKN,nsize*mbas)
       Call GetMem('LN','ALLO','INTE',ipLN,nsize*mbas)
!
       call exppsb (symp,symq,symr,syms,Work(iOff_valn),                &
     &              iWork(ipJN),iWork(ipKN),iWork(ipLN))
!
!     release space for valn,jn,kn,ln
       Call GetMem('VALN','FREE','REAL',iOff_valn,nsize*mbas)
       Call GetMem('JN','FREE','INTE',ipJN,nsize*mbas)
       Call GetMem('KN','FREE','INTE',ipKN,nsize*mbas)
       Call GetMem('LN','FREE','INTE',ipLN,nsize*mbas)
!
       end if
!
!
!*    run over all pivot indexes
!
       do 500 p=1,norb(symp)
!
!**   def dimensions of vint and read this block of integrals into vint
       ndimvi=norb(symp)
       ndimv1=norb(symq)
       ndimv2=norb(symr)
       ndimv3=norb(syms)
!
       if (ickey.eq.0) then
!      Out of core expansion
!
       if (symq.eq.syms) then
       keyred=1
       else
       keyred=0
       end if
       call unpackk (p,Work(iOff),ndimv1,ndimv2,ndimv3,keyred)
!
       else if (ickey.eq.1) then
!      else Incore expansions
!
!      case V(p,q,r,s)
         call unpackk_ic_1 (p,Work(iOff),ndimv1,ndimv2,ndimv3,          &
     &                      Work(iOff_Vic),ndimvi)
       else if (ickey.eq.2) then
!      case V(pr,qs)
         call unpackk_ic_2 (p,Work(iOff),ndimvi,ndimv1,Work(iOff_Vic))
       else
!      case V(prqs)
         call unpackk_ic_3 (p,Work(iOff),ndimvi,Work(iOff_Vic))
       end if
!
!
!        goto 500
!
!
!**   add integrals to T3nam if needed (v nacechranej forme)
!
       if (t3key.eq.1) then
       if (p.le.noa(symp)) then
       if (symq.gt.syms) then
!
!***  calc proper address in t3nam file
       t3help4=0
       do t3help1=1,symp-1
       t3help4=t3help4+noa(t3help1)
       end do
       t3help4=t3help4+p
       t3help1=mapiri(symr,symq,1)
       daddr(lunt3)=T3IntPoss(t3help4)+T3Off(t3help1,symp)
!
!***  def required parameters
       t3help1=nvb(symr)
       t3help2=nvb(symq)
       t3help3=nvb(syms)
!
!***  do packing
       call t3intpck2(Work(iOff),Work(iOff+possri0-1),ndimv1,ndimv2,    &
     & ndimv3,t3help1,t3help2,t3help3,                                  &
     & symq,symr,syms,nob,nvb)
!
       else if (symq.eq.syms) then
!
!***  calc proper address in t3nam file
       t3help4=0
       do t3help1=1,symp-1
       t3help4=t3help4+noa(t3help1)
       end do
       t3help4=t3help4+p
       t3help1=mapiri(symr,symq,1)
       daddr(lunt3)=T3IntPoss(t3help4)+T3Off(t3help1,symp)
!
!***  def required parameters
       t3help1=nvb(symr)
       t3help2=nvb(symq)*(nvb(symq)+1)/2
!
!***  do packing
       call t3intpck1(Work(iOff),Work(iOff+possri0-1),ndimv1,ndimv2,    &
     & ndimv3,t3help1,t3help2,                                          &
     & symq,symr,syms,nob,nvb)
!
       end if
       end if
       end if
!
!**   add integrals to #1 <pq|ij> if needed
!
       if (symr.ge.syms) then
!     contributions only for symi(r)>=symj(s)
       call addpqij (Work(iOff),wrksize,                                &
     & symp,symq,symr,syms,p,Work(iOff),ndimv1,ndimv2,ndimv3)
       end if
!
!**   updete fok if neccesarry (only for open shell case)
!
       if (clopkey.eq.1) then
!
       if ((symp.eq.symr).and.(symq.eq.syms)                            &
     & .and.(p.gt.nob(symp)).and.(p.le.(noa(symp)))) then
       call fokupdate1 (foka,fokb,symq,p,Work(iOff),                    &
     &                  ndimv1,ndimv2,ndimv3)
       end if
!
       if ((symp.eq.syms).and.(symq.eq.symr)                            &
     & .and.(p.gt.nob(symp)).and.(p.le.(noa(symp)))) then
       call fokupdate2 (foka,symq,p,Work(iOff),                         &
     &                  ndimv1,ndimv2,ndimv3)
       end if
!
       end if
!
!**   add corresponding <am|pq> interals to TEMPDA2
!     and pack _a_brs to direct acces file TEMPDA1 if need and symm in not C1
       if (p.gt.nob(symp)) then
       a=p-nob(symp)
       call ampack (Work(iOff),wrksize,                                 &
     & symp,symq,symr,syms,a,Work(iOff),ndimv1,ndimv2,ndimv3,AMMAP)
       if ((nsym.gt.1).and.(symp.ge.symq)) then
       call abpack (Work(iOff),wrksize,                                 &
     & symp,symq,symr,syms,a,Work(iOff),ndimv1,ndimv2,ndimv3,ABMAP)
       end if
       end if
!
!**   add INTAB file (for nonsymmetrycal (C1) state)
!
       if ((nsym.eq.1).and.(p.gt.nob(1))) then
       call addintabc1 (Work(iOff),wrksize,                             &
     & p-nob(1),Work(iOff),ndimv1)
       end if
!
 500    continue
!
!        goto 800
!
!
!*    close temp files
!    call closetemp (norb(symp))
!
       if (ickey.ge.1) then
       Call GetMem('CCSORT','Free','Real',iOff_Vic,vsize)
       end if
!
 800    continue
!
!        goto 900
!
!
       if ((nsym.gt.1).and.(symp.ge.symq)) then
!*    add contributions to INTAB comming from symp,sumq and close TEMPDA1 file
!     only for symmetrical cases; only for syma>=symb
       call addintab (Work(iOff),wrksize,symp,symq,ABMAP)
       close (lunda1)
       call vf ('TEMPDA1 ',lunda1)
       end if
!
 900    continue
!
!        goto 1000
!
!
!*    add contributions to INTA1-4 if there are some virtuals in symp symmetry
!     and close TEMPDA2 files
!
!     if (nvb(symp).gt.0) then
       call addinta (Work(iOff),wrksize,symp,AMMAP)
       close (lunda2)
       call vf ('TEMPDA2 ',lunda2)
!     end if
!
 1000   continue
!
!*    if T3 are required, reorganize T3nam file
       if (t3key.eq.1) then
       call t3reorg (Work(iOff),wrksize,                                &
     &               noa,nsym)
         do symp=1,4
         end do
       end if
!
!        release space for ammap,abmap
        Call mma_Deallocate(AMMAP)
        Call mma_Deallocate(ABMAP)
!
!
!*    close files INTA1,INTA2,INTA3 and INTA4, INTAB1
!
       if (iokey.eq.1) then
!      Fortran IO
       close (luna1)
       close (luna2)
       close (luna3)
       close (luna4)
       close (lunab)
!
       else
!      MOLCAS IO
       call daclos (luna1)
       call daclos (luna2)
       call daclos (luna3)
       call daclos (luna4)
       call daclos (lunab)
       end if
!
!        return
!
!
!*    def static integrals (file INTSTA)
!
       call mkintsta (Work(iOff),wrksize,                               &
     & foka,fokb)
!
!*    write general informations to INPDAT
!
       call molcas_binaryopen_vanilla(1,'INPDAT')
!       open (unit=1,file='INPDAT',form='unformatted')
       write (1) NACTEL,ISPIN,NSYM,LSYM,mul,                            &
     &           noa,nob,nva,nvb,norb,eps,Escf
       close (1)
!
!      Release the memory
       Call GetMem('CCSORT','Free','Real',iOff,wrksize)

       return
       end
