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
*               2019, Oskar Weser                                      *
************************************************************************

      module fciqmc_make_inp
        integer ::
     &    totalwalkers = 500000,
! Default value for time per NECI run
     &    Time = 200,
! Practically this means no trial_wavefunction by default
     &    trial_wavefunction = 1000000000,
     &    nmcyc = 50000,
     &    pops_trial = 1000,
     &    rdmsamplingiters =  25000,
     &    stepsshift = 10,
     &    addtoinitiator = 3,
     &    maxwalkerbloom = 1,
     &    semi_stochastic = 500,
     &    highlypopwrite = 50,
     &    calcrdmonfly(3) = [3, 500, 500]
        integer, allocatable ::
     &    definedet(:)
        real(kind=8) ::
     &    proje_changeref = 1.2d0,
     &    max_tau = 0.02d0,
     &    memoryfacpart = 5.0d0,
     &    memoryfacspawn = 10.0d0,
     &    realspawncutoff = 0.3d0, ! Default value for NECI RealSpawnCutOff
     &    diagshift = 0.00d0, ! Default value for NECI diagonal shift value
     &    shiftdamp = 0.02d0
        save
      contains

!>  @brief
!>    generate a standardized Input for NECI
!>
!>  @author
!>    G. Li Manni, Oskar Weser
      subroutine make_inp
      use general_data, only : nActEl, iSpin
      use stdalloc, only : mma_deallocate
      implicit none
      integer :: i, isFreeUnit, file_id
      character(len=10) :: formt
      character(*), parameter ::
     &    I_fmt = '(A, 1x, I0)',
     &    R_fmt = '(A, 1x, F0.2)'

*----------------------------------------------------------------------*
*---- Check that RunTime variables are passed correctly
*---- This is checked by Molcas verify
*----------------------------------------------------------------------*
      call add_info('Default number of total walkers',
     &  dble(totalwalkers), 1, 6)
      call add_info('Default number of cycles ',dble(nmcyc),1,6)
      call add_info('Default value for Time  ',dble(Time),1,6)

      call qEnter('make_inp')

      file_id = isFreeUnit(39)
      call Molcas_Open(file_id,'FCINP')

      write(file_id,'(A5)') 'Title'
      write(file_id,'(A1)') ''
      write(file_id,'(A11)') 'System read'
      write(file_id,'(A10,I3)') 'electrons ', nActEl
      write(file_id,'("nonuniformrandexcits 4ind-weighted-2")')
      write(file_id,'(A18)') 'nobrillouintheorem'
      if(iSpin .ne. 1) write(file_id, I_fmt) 'spin-restrict', iSpin - 1
      write(file_id,'(A10)') 'freeformat'
      write(file_id,'(A6)') 'endsys'
      write(file_id,'(A1)') ''
      write(file_id,'(A4)') 'calc'
      if (allocated(DefineDet)) then
        write(formt,'(I5)') nActEl
        formt='(A,('//trim(adjustl(formt))//'I5))'
        write(file_id,formt) 'definedet', (definedet(i), i = 1, nActEl)
        call mma_deallocate(definedet)
      end if
      write(file_id,'(A1)') ''
      write(file_id,'(A7)') 'methods'
      write(file_id,'(A19)') 'method vertex fcimc'
      write(file_id,'(A10)') 'endmethods'
      write(file_id,'(A1)') ' '
      write(file_id, I_fmt) 'totalwalkers', totalwalkers
      write(file_id, R_fmt) 'diagshift', diagshift
      write(file_id, R_fmt) 'shiftdamp', shiftdamp
      write(file_id, I_fmt) 'nmcyc', nmcyc
      write(file_id, I_fmt) 'stepsshift', stepsshift
      write(file_id, R_fmt) 'proje-changeref', proje_changeref
      write(file_id,'(A14)') 'truncinitiator'
      write(file_id, I_fmt) 'addtoinitiator ', addtoinitiator
      write(file_id,'(A12)') 'allrealcoeff'
      write(file_id, R_fmt) 'realspawncutoff', realspawncutoff
      write(file_id,'(A10)') 'jump-shift'
      write(file_id,'(A15)') 'tau 0.01 search'
      write(file_id, R_fmt) 'max-tau', max_tau
      write(file_id, I_fmt) 'maxwalkerbloom', maxwalkerbloom
      write(file_id, R_fmt) 'memoryfacspawn', memoryfacspawn
      write(file_id, R_fmt) 'memoryfacpart', memoryfacpart
      write(file_id,'(A, 1x, I4)') 'time', time
      write(file_id,'("startsinglepart 10")')
      write(file_id, I_fmt) 'semi-stochastic', semi_stochastic
      write(file_id, '("pops-core 10000")')
      write(file_id, I_fmt), 'rdmsamplingiters', rdmsamplingiters
!      if(KeyTRIA) write(file_id, I_fmt)
!     & 'trial-wavefunction', trial_wavefunction
!      write(file_id, I_fmt) 'pops-trial', pops_trial
c      if(abs(rotmax).le.1.0d-1.and.iter.ne.1) then
c          write(file_id,'(A)') 'readpops'
c      end if
      write(file_id,'(A7)') 'endcalc'
      write(file_id,'(A1)') ' '
      write(file_id,'(A7)') 'logging'
      write(file_id, I_fmt) 'Highlypopwrite', Highlypopwrite
      write(file_id,'(A24)') 'Print-Spin-Resolved-RDMS'
      write(file_id,'(A10)') 'binarypops'
      write(file_id,'(A11)') 'printonerdm'
      write(file_id, '("diagflyonerdm")')
      write(file_id,'(A,1x,I0,1x,I0,1x,I0)')
     &     'calcrdmonfly', (calcrdmonfly(i), i=1,3)
      write(file_id,'(A6)') 'endlog'
      write(file_id,'(A3)') 'end'


      close(file_id)
      call qExit('make_inp')

      return
      end subroutine make_inp
      end module fciqmc_make_inp
