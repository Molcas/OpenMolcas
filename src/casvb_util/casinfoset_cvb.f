************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1996-2006, T. Thorsteinsson and D. L. Cooper           *
************************************************************************
      subroutine casinfoset_cvb()
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "inpmod_cvb.fh"
#include "malloc_cvb.fh"
#include "casinfo_cvb.fh"
      logical hadinput
      logical debug
      data debug/.false./

      if(debug)then
        write(6,*)' casinfoset :'
        write(6,*)' ------------'
        write(6,*)' iorocc_c  :',iorocc_c
        write(6,*)' iorclos_c :',iorclos_c
        write(6,*)' iorcore_c :',iorcore_c
        write(6,*)' nstsym_c  :',nstsym_c
        write(6,*)' weight_c  :',weight_c
        write(6,*)' istnel_c  :',istnel_c
        write(6,*)' istsy_c   :',istsy_c
        write(6,*)' istms2_c  :',istms2_c
        write(6,*)' nstats_c  :',nstats_c
        write(6,*)' strtint_c :',strtint_c
        write(6,*)' strtci_c  :',strtci_c
        write(6,*)' strtmo_c  :',strtmo_c
        write(6,*)' iorocc_d  :',iorocc_d
        write(6,*)' iorclos_d :',iorclos_d
        write(6,*)' iorcore_d :',iorcore_d
        write(6,*)' nstsym_d  :',nstsym_d
        write(6,*)' weight_d  :',weight_d
        write(6,*)' istnel_d  :',istnel_d
        write(6,*)' istsy_d   :',istsy_d
        write(6,*)' istms2_d  :',istms2_d
        write(6,*)' nstats_d  :',nstats_d
        write(6,*)' strtint_d :',strtint_d
        write(6,*)' strtci_d  :',strtci_d
        write(6,*)' strtmo_d  :',strtmo_d
      endif
      hadinput=.false.
      do 100 i=1,mxstsy_ci
      if(iorocc_d(i).ne.-1)hadinput=.true.
      if(iorclos_d(i).ne.-1)hadinput=.true.
      if(iorcore_d(i).ne.-1)hadinput=.true.
100   continue
      if(hadinput)then
        do 200 i=1,mxstsy_ci
        if(iorocc_d(i).eq.-1)iorocc_d(i)=0
        if(iorclos_d(i).eq.-1)iorclos_d(i)=0
        if(iorcore_d(i).eq.-1)iorcore_d(i)=0
200     continue
      else
        call imove_cvb(iorcore_c,iorcore_d,mxstsy_ci)
        call imove_cvb(iorclos_c,iorclos_d,mxstsy_ci)
        call imove_cvb(iorocc_c,iorocc_d,mxstsy_ci)
      endif

      mcore_d=0
c  Ensure no negative number of orbitals :
      do 300 i=1,mxirrep
      iorclos_d(i)=iorclos_d(i)+iorcore_d(i)
      iorocc_d(i)=iorocc_d(i)+iorclos_d(i)
      mcore_d=mcore_d+iorclos_d(i)
300   continue

      if(nstsym_d.eq.0)then
        nstsym_d=nstsym_c
        call imove_cvb(nstats_c,nstats_d,mxstsy_ci)
        call imove_cvb(istnel_c,istnel_d,mxstsy_ci)
        call imove_cvb(istsy_c,istsy_d,mxstsy_ci)
        call imove_cvb(istms2_c,istms2_d,mxstsy_ci)
        call fmove_cvb(weight_c,weight_d,mxstt_ci*mxstsy_ci)
        if(mcore_d.ne.mcore_c)then
c  Different number of core orbitals input -> assume NELTOT the same :
          do 400 i=1,mxstsy_ci
          if(istnel_d(i).ne.0)istnel_d(i)=istnel_d(i)
     >      +2*(mcore_c-mcore_d)
400       continue
        endif
      endif
      strtint_d=strtint
      strtmo_d=strtmo
      strtci_d=strtci
      if(.not.valid_cvb(strtint_d))strtint_d=strtint_c
      if(.not.valid_cvb(strtmo_d))strtmo_d=strtmo_c
      if(.not.valid_cvb(strtci_d))strtci_d=strtci_c
      strtint=strtint_d
      strtmo=strtmo_d
      strtci=strtci_d

c  Set active space information
      sum=zero
      do 600 i=1,nstsym_d
      do 601 j=1,nstats_d(i)
      if(weight_d(j,i).lt.zero)then
        write(6,'(a,f10.4,i3,a,i1)')
     >    ' Fatal error: WEIGHT factor negative :',weight_d(j,i),j,'.',i
        call abend_cvb()
      endif
      sum=sum+weight_d(j,i)
601   continue
600   continue
      sum=one/sum
      call dscal_(mxstt_ci*mxstsy_ci,sum,weight_d,1)
      nel_d=-1
      i2s_d=-1
      call izero(isymv,mxirrep)
      do 700 i=1,nstsym_d
      do 800 j=1,nstats_d(i)
      if(weight_d(j,i).gt.1.d-20)goto 900
800   continue
      goto 700
900   continue
      if(nel_d.ne.-1.and.nel_d.ne.istnel_d(i))then
        write(6,*)' Fatal error: ELEC varies in WF cards!'
        call abend_cvb()
      endif
      if(i2s_d.ne.-1.and.i2s_d.ne.istms2_d(i))then
        write(6,*)' Fatal error: SPIN varies in WF cards!'
        call abend_cvb()
      endif
      nel_d=istnel_d(i)
      i2s_d=istms2_d(i)
      isym_d=istsy_d(i)
      isymv(isym_d)=1
700   continue
      nsym=0
      do 1000 is=1,mxirrep
      if(isymv(is).eq.1)nsym=nsym+1
1000  continue

      nel=nel_d
      isym=isym_d

      call izero(ityp,mxorb)
      mcore=0
      norb=0
      incr=0
      do 1100 i=1,mxirrep
      ioc=iorocc_d(i)-iorclos_d(i)
      norb=norb+ioc
      mcore=mcore+iorclos_d(i)
      do 1200 j=1,ioc
      ityp(j+incr)=i
1200  continue
      incr=incr+ioc
1100  continue
c  Set NIRREP :
      nirrep=1
      do 1350 irrep=1,mxirrep
      if(iorcore_d(irrep).gt.0.or.iorclos_d(irrep).gt.0.or.
     >   iorocc_d(irrep).gt.0)nirrep=irrep
1350  continue
      if(nirrep.eq.3)nirrep=4
      if(nirrep.gt.4)nirrep=8

      noe=max(norb,nel)
      nbet=(nel-i2s_d)/2
      nalf=nel-nbet
c  Basic checks
      if(nel.lt.0.or.norb.lt.0.or.i2s_d.lt.0.or.nel.gt.2*norb.or.
     >  mod(nel,2).ne.mod(i2s_d,2))then
        write(6,*)' Impossible numbers: active electrons :',nel
        write(6,*)'                     active orbitals  :',norb
        write(6,*)'                     total spin       :',
     >    DBLE(nalf-nbet)/two
        call abend_cvb()
      endif
      if(isym.eq.0)then
        write(6,*)' WARNING: State symmetry not found - assuming A1.'
        isym=1
        nsym=1
        call izero(isymv,mxirrep)
        isymv(1)=1
      endif
      return
      end
