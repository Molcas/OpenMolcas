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
      integer, parameter :: iDM=3
      real, parameter :: sDamp= 0.02
      character(len=10) :: formt

*----------------------------------------------------------------------*
*---- Check that RunTime variables are passed correctly
*---- This is checked by Molcas verify
*----------------------------------------------------------------------*
      call add_info('Default number of total walkers',dble(nTWlk),1,6)
      call add_info('Default number of cycles ',dble(nmcyc),1,6)
      call add_info('Default value for iTime  ',dble(iTime),1,6)
*----------------------------------------------------------------------*
*---- Start program and say Hello                                      *
*----------------------------------------------------------------------*
      Call qEnter('make_inp')
      Dbg=.false.
*----------------------------------------------------------------------*
*----  Open LuFCINP for NECI interface
*----------------------------------------------------------------------*
      LuFCINP= isFreeUnit(39)
      call Molcas_Open(LuFCINP,'FCINP')
*----------------------------------------------------------------------*
*----  Fill it up!
*----------------------------------------------------------------------*
      Call get_iscalar('iSpin',iSpin)
      write(LuFCINP,'(A5)') 'Title'
      write(LuFCINP,'(A1)') ''
      write(LuFCINP,'(A11)') 'System read'
      Call Get_iScalar('nActel',nActEl)
      write(LuFCINP,'(A10,I3)') 'electrons ', nActEl
      write(LuFCINP,'("nonuniformrandexcits 4ind-weighted-2")')
      write(LuFCINP,'(A18)') 'nobrillouintheorem'
      if(iSpin.ne.1) then
       write(LuFCINP,'(A13,I5)')'spin-restrict',int((ISPIN-1.0d0))
      end if
      write(LuFCINP,'(A10)') 'freeformat'
      write(LuFCINP,'(A6)') 'endsys'
      write(LuFCINP,'(A1)') ''
      write(LuFCINP,'(A4)') 'calc'
      if(DefineDet) then
        write(formt,'(I5)') nActEl
        formt='(A,('//trim(adjustl(formt))//'I5))'
        write(LuFCINP,formt)'definedet', (iWork(ipDet+i),i=0,nActEl-1)
      end if
      write(LuFCINP,'(A1)') ''
      write(LuFCINP,'(A7)') 'methods'
      write(LuFCINP,'(A19)') 'method vertex fcimc'
      write(LuFCINP,'(A10)') 'endmethods'
      write(LuFCINP,'(A1)') ' '
      write(LuFCINP,'(A13,I10)') 'totalwalkers ', nTWlk
      write(LuFCINP,'(A9,F5.2)') 'diagshift', diagshift
      write(LuFCINP,'(A9,F5.2)') 'shiftdamp', sDamp
      write(LuFCINP,'(A6,I10)') 'nmcyc ', nmcyc
      write(LuFCINP,'(A13)') 'stepsshift 10'
      write(LuFCINP,'(A19)') 'proje-changeref 1.2'
      write(LuFCINP,'(A14)') 'truncinitiator'
      write(LuFCINP,'(A16)') 'addtoinitiator 3'
      write(LuFCINP,'(A12)') 'allrealcoeff'
      write(LuFCINP,'(A16,F5.2)') 'realspawncutoff',realspawncutoff
      write(LuFCINP,'(A10)') 'jump-shift'
      write(LuFCINP,'(A15)') 'tau 0.01 search'
      write(LuFCINP,'(A12)') 'max-tau 0.02'
      write(LuFCINP,'(A16)') 'maxwalkerbloom 1'
      write(LuFCINP,'(A19)') 'memoryfacspawn 10.0'
      write(LuFCINP,'(A17)') 'memoryfacpart 5.0'
      write(LuFCINP,'(A5,I4)') 'time ', itime
      write(LuFCINP,'("startsinglepart 10")')
      write(LuFCINP, '("semi-stochastic 1000")')
      write(LuFCINP, '("pops-core 10000")')
      write(LuFCINP, '("rdmsamplingiters 25000")')
      write(LuFCINP,'(A18)') 'trial-wavefunction'
      write(LuFCINP,'(A14)') 'pops-trial 500'
c      if(abs(rotmax).le.1.0d-1.and.iter.ne.1) then
c          write(LuFCINP,'(A)') 'readpops'
c      end if
      write(LuFCINP,'(A7)') 'endcalc'
      write(LuFCINP,'(A1)') ' '
      write(LuFCINP,'(A7)') 'logging'
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
