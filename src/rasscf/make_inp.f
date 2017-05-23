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
* Copyright (C) 2015, Giovanni Li Manni                                *
************************************************************************
      subroutine make_inp(iter,rotmax)
************************************************************************
*                                                                      *
*     Objective: generate a standardized Input for NECI                *
*                                                                      *
*     Author:    G. Li Manni, Max-Planck Institute January 2015        *
*                                                                      *
************************************************************************
#include "trafo_fciqmc.fh"
#include "fciqmc_global.fh"
#include "WrkSpc.fh"
#include "fciqmc.fh"
      Logical :: dbg
      integer, parameter :: nWalkRef=3500 ! My initial choice
!      integer, parameter ::nTWlk=500000, nWalkRef=3500, iTime =200,nmcyc=50000 !going to Molcas Runtime choice
      integer, parameter :: iDM=3
!      integer, parameter ::IterQMC=3000,IterEnQMC=500,nmcyc=15000,iDM=3
      real, parameter :: sDamp= 0.02
!      integer :: iSym
      character(len=10) :: formt

*----------------------------------------------------------------------*
*---- Check that RunTime variables are passed correctly
*---- This is checked by Molcas verify
*----------------------------------------------------------------------*
      call add_info('Default number of total walkers',real(nTWlk),1,6)
      call add_info('Default number of cycles ',real(nmcyc),1,6)
      call add_info('Default value for iTime  ',real(iTime),1,6)
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*---- Start program and say Hello                                      *
*----------------------------------------------------------------------*
      Call qEnter('make_inp')
      Dbg=.true.
*----------------------------------------------------------------------*
*----  Open LuFCINP for NECI interface
*----------------------------------------------------------------------*
      LuFCINP= isFreeUnit(39)
      call Molcas_Open(LuFCINP,'FCINP')
*----------------------------------------------------------------------*
*----  Fill it up!
*----------------------------------------------------------------------*
c
      Call get_iscalar('iSpin',iSpin)
      write(LuFCINP,'(A5)') 'Title'
      write(LuFCINP,'(A1)') ''
      write(LuFCINP,'(A11)') 'System read'
      Call Get_iScalar('nActel',nActEl)
      write(LuFCINP,'(A10,I3)') 'electrons ', nActEl
      write(LuFCINP,'("nonuniformrandexcits 4ind-weighted-2")')
      write(LuFCINP,'(A18)') 'nobrillouintheorem'
c      if (btest(iSpin, 0)) then
c       ! iSpin is odd --> Use HPHF
c          write(LuFCINP, "('HPHF ', i1)") mod(iSpin-1, 4) / 2
c      end if
      if(iSpin.ne.1) then
       write(LuFCINP,'(A13,I5)')'spin-restrict',int((ISPIN-1.0d0))
      end if
!      write(LuFCINP,'(A3)') 'SYM'
!      Call Get_iScalar('LSym',iSym)
!      write(LuFCINP,'(4I3)') 0,0,0,iSym-1
      write(LuFCINP,'(A10)') 'freeformat'
      write(LuFCINP,'(A6)') 'endsys'
      write(LuFCINP,'(A1)') ''
      write(LuFCINP,'(A4)') 'calc'
      if(DefineDet) then
        write(formt,'(I5)') nActEl
        formt='(A,('//trim(adjustl(formt))//'I5))'
        write(LuFCINP,formt)'definedet', (iWork(ipDet+i),i=0,nActEl-1)
      end if
!     &19,20,22,27,28,29,30,31,32,39,40,41,42'
! This one for Sym2
!        write(LuFCINP,'A') 'definedet 1,2,3,4,5,6,7,8,9,15,16,17,18,
!     &19,20,22,27,28,29,30,31,32,39,40,41,42'
! This one for Sym3
!     write(LuFCINP,'(A)') 'definedet 1,2,3,4,5,6,7,8,9,15,16,17,18,
!     &19,20,27,28,29,30,31,32,34,39,40,41,42'
! This one for Sym4
!      write(LuFCINP,'(A)') 'definedet 1,2,3,4,5,6,7,8,9,15,16,17,18,
!     &19,20,27,28,29,30,31,32,39,40,41,42,44'
      write(LuFCINP,'(A1)') ''
      write(LuFCINP,'(A7)') 'methods'
      write(LuFCINP,'(A19)') 'method vertex fcimc'
      write(LuFCINP,'(A10)') 'endmethods'
      write(LuFCINP,'(A1)') ' '
      write(LuFCINP,'(A13,I10)') 'totalwalkers ', nTWlk
!      write(LuFCINP,'(A13,I10)') 'maxnoathf ', nWalkRef
      write(LuFCINP,'(A9,F5.2)') 'diagshift', diagshift
      write(LuFCINP,'(A9,F5.2)') 'shiftdamp', sDamp
      write(LuFCINP,'(A6,I10)') 'nmcyc ', nmcyc
!      if(iter.gt.1) then
!        write(LuFCINP,'(A13)') 'stepsshift 1'
!      else
        write(LuFCINP,'(A13)') 'stepsshift 10'
!      end if
      write(LuFCINP,'(A19)') 'proje-changeref 1.2'
!      write(LuFCINP,'(A19)') 'targetgrowrate 0.01 50'
!      next four entries have been commented off for testing purposes
      write(LuFCINP,'(A14)') 'truncinitiator'
      write(LuFCINP,'(A16)') 'addtoinitiator 3'
      write(LuFCINP,'(A12)') 'allrealcoeff'
      write(LuFCINP,'(A16,F5.2)') 'realspawncutoff',realspawncutoff
!      write(LuFCINP,'(A20)') 'realspawncutoff 0.3' ! Suggested by Ali
      write(LuFCINP,'(A10)') 'jump-shift'
      write(LuFCINP,'(A15)') 'tau 0.01 search'
      write(LuFCINP,'(A12)') 'max-tau 0.02'
      write(LuFCINP,'(A16)') 'maxwalkerbloom 1'
!      write(LuFCINP,'(A19)') 'memoryfacspawn 5.0' ! Suggested by Ali
      write(LuFCINP,'(A19)') 'memoryfacspawn 10.0'
      write(LuFCINP,'(A17)') 'memoryfacpart 5.0'
      write(LuFCINP,'(A5,I4)') 'time ', itime
      write(LuFCINP,'("startsinglepart 10")')
      write(LuFCINP, '("semi-stochastic 1000")')
      write(LuFCINP, '("pops-core 10000")')
      write(LuFCINP, '("rdmsamplingiters 25000")')
!      write(LuFCINP,'(A18)') 'trial-wavefunction'
!      write(LuFCINP,'(A14)') 'pops-trial 500'
c      if(abs(rotmax).le.1.0d-1.and.iter.ne.1) then
c          write(LuFCINP,'(A)') 'readpops'
c      else
!          write(LuFCINP, '("doubles-core")')
c          write(LuFCINP, '("startmp1", i10)') int(nwalkers / 2) ! Simon agreed to remove it (for now!)
c      end if
      write(LuFCINP,'(A7)') 'endcalc'
      write(LuFCINP,'(A1)') ' '
      write(LuFCINP,'(A7)') 'logging'
!      if(iter.gt.1) write(LuFCINP,'(A12)') 'fcimcdebug 6'
!      write(LuFCINP,'(A20)') 'INSTANT-S2-FULL 1000'
      write(LuFCINP,'(A17)') 'Highlypopwrite 50'
      write(LuFCINP,'(A10)') 'binarypops'
      write(LuFCINP,'(A11)') 'printonerdm'
      write(LuFCINP, '("diagflyonerdm")')
      write(LuFCINP,'(A13,3I10)')
     &     'calcrdmonfly ',iDM,IterFillRDM,IterSampleRDM
      write(LuFCINP,'(A6)') 'endlog'
      write(LuFCINP,'(A3)') 'end'

      close(LuFCINP)
      Call qExit('make_inp')

      return
c Avoid unused argument warnings
      if (.false.) then
        Call Unused_integer(iter)
        Call Unused_real(rotmax)
      end if
      End
