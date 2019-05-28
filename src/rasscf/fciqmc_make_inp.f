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
! No default value on purpose
     &    totalwalkers,
     &    calcrdmonfly(2),
     &    rdmsamplingiters,
! Default value for time per NECI run
     &    Time = 200,
! Practically this means no trial_wavefunction by default
     &    trial_wavefunction = 1000000000,
     &    nmcyc = 50000,
     &    pops_trial = 1000,
     &    stepsshift = 10,
     &    addtoinitiator = 3,
     &    maxwalkerbloom = 1,
     &    semi_stochastic = 1000,
     &    highlypopwrite = 50,
     &    startsinglepart = 10,
     &    pops_core =  10000
        integer, allocatable ::
     &    definedet(:)
        real*8 ::
     &    proje_changeref = 1.2d0,
     &    max_tau = 0.02d0,
     &    memoryfacpart = 5.0d0,
     &    memoryfacspawn = 10.0d0,
! Default value for NECI RealSpawnCutOff
     &    realspawncutoff = 0.3d0,
! Default value for NECI diagonal shift value
     &    diagshift = 0.00d0,
     &    shiftdamp = 0.02d0
        save
      contains

!>  @brief
!>    Generate a standardized Input for NECI
!>
!>  @author
!>    G. Li Manni, Oskar Weser
!>
!>  @paramin[in] readpops  If true the readpops option for NECI is set.
      subroutine make_inp(readpops)
      use general_data, only : nActEl, iSpin
      use stdalloc, only : mma_deallocate
      use fortran_strings, only : str
      implicit none
      logical, intent(in), optional :: readpops
      logical :: readpops_
      integer :: i, isFreeUnit, file_id, indentlevel
      integer, parameter :: indentstep = 4
      character(len=10) :: formt

      readpops_ = merge(readpops, .false., present(readpops))

*----------------------------------------------------------------------*
*---- Check that RunTime variables are passed correctly
*---- This is checked by Molcas verify
*----------------------------------------------------------------------*
      call add_info('Default number of total walkers',
     &  [dble(totalwalkers)], 1, 6)
      call add_info('Default number of cycles ',[dble(nmcyc)],1,6)
      call add_info('Default value for Time  ',[dble(Time)],1,6)

      call qEnter('make_inp')

      file_id = isFreeUnit(39)
      call Molcas_Open(file_id,'FCINP')

      indentlevel = 0
      write(file_id, A_fmt()) 'Title'
      write(file_id, A_fmt()) ''
      write(file_id, A_fmt()) 'System read'
      call indent()
        write(file_id, I_fmt()) 'electrons ', nActEl
        write(file_id,A_fmt()) 'nonuniformrandexcits 4ind-weighted-2'
        write(file_id,A_fmt()) 'nobrillouintheorem'
        if(iSpin /= 1) then
          write(file_id, I_fmt()) 'spin-restrict', iSpin - 1
        end if
        write(file_id, A_fmt()) 'freeformat'
      call dedent()
      write(file_id, A_fmt()) 'endsys'
      write(file_id, A_fmt()) ''
      write(file_id, A_fmt()) 'calc'
      call indent()
        if (allocated(DefineDet)) then
          write(formt,I_fmt()) nActEl
          formt='(A,('//trim(adjustl(formt))//'I5))'
          write(file_id,formt) 'definedet', (definedet(i), i = 1,nActEl)
          call mma_deallocate(definedet)
          write(file_id, A_fmt()) ''
        end if
        write(file_id, A_fmt()) 'methods'
        call indent()
          write(file_id, A_fmt()) 'method vertex fcimc'
        call dedent()
        write(file_id, A_fmt()) 'endmethods'
        write(file_id, A_fmt()) ''
        if (readpops_) write(file_id, A_fmt()) 'readpops'
        if (readpops_) write(file_id, A_fmt()) 'walkcontgrow'
        write(file_id, I_fmt()) 'semi-stochastic', semi_stochastic
        write(file_id, A_fmt()) ''
        write(file_id, I_fmt()) 'totalwalkers', totalwalkers
        write(file_id, R_fmt()) 'diagshift', diagshift
        write(file_id, R_fmt()) 'shiftdamp', shiftdamp
        write(file_id, I_fmt()) 'nmcyc', nmcyc
        write(file_id, I_fmt()) 'stepsshift', stepsshift
        write(file_id, R_fmt()) 'proje-changeref', proje_changeref
        write(file_id, A_fmt()) 'truncinitiator'
        write(file_id, I_fmt()) 'addtoinitiator ', addtoinitiator
        write(file_id, A_fmt()) 'allrealcoeff'
        write(file_id, R_fmt()) 'realspawncutoff', realspawncutoff
        write(file_id, A_fmt()) 'jump-shift'
        write(file_id, A_fmt()) 'tau 0.01 search'
        write(file_id, R_fmt()) 'max-tau', max_tau
        write(file_id, I_fmt()) 'maxwalkerbloom', maxwalkerbloom
        write(file_id, R_fmt()) 'memoryfacspawn', memoryfacspawn
        write(file_id, R_fmt()) 'memoryfacpart', memoryfacpart
        write(file_id, I_fmt()) 'time', time
        write(file_id, I_fmt()) 'startsinglepart', startsinglepart
        write(file_id, I_fmt()) 'pops-core', pops_core
        write(file_id, I_fmt()) 'rdmsamplingiters', rdmsamplingiters
      call dedent()
      write(file_id, A_fmt()) 'endcalc'
      write(file_id, A_fmt()) ' '
      write(file_id, A_fmt()) 'logging'
      call indent()
        write(file_id, I_fmt()) 'Highlypopwrite', Highlypopwrite
        write(file_id, A_fmt()) 'Print-Spin-Resolved-RDMS'
        write(file_id, A_fmt()) 'hdf5-pops'
        write(file_id, A_fmt()) 'printonerdm'
        write(file_id, A_fmt()) '(diagflyonerdm'
        write(file_id,'('//str(indentlevel)//'x, A,1x,I0,1x,I0,1x,I0)')
     &     'calcrdmonfly', 3, (calcrdmonfly(i), i=1,2)
      call dedent()
      write(file_id, A_fmt()) 'endlog'
      write(file_id, A_fmt()) 'end'


      close(file_id)
      call qExit('make_inp')

      contains

        function I_fmt() result(res)
          implicit none
          character(:), allocatable :: res
          if (indentlevel /= 0) then
            res = '('//str(indentlevel)//'x, A, 1x, I0)'
          else
            res = '(A, 1x, I0)'
          end if
        end function

        function R_fmt() result(res)
          implicit none
          character(:), allocatable :: res
          if (indentlevel /= 0) then
            res = '('//str(indentlevel)//'x, A, 1x, F0.2)'
          else
            res = '(A, 1x, F0.2)'
          end if
        end function

        function A_fmt() result(res)
          implicit none
          character(:), allocatable :: res
          if (indentlevel /= 0) then
            res = '('//str(indentlevel)//'x, A)'
          else
            res = '(A)'
          end if
        end function

        subroutine indent()
          implicit none
          indentlevel = indentlevel + indentstep
        end subroutine

        subroutine dedent()
          implicit none
          indentlevel = indentlevel - indentstep
        end subroutine
      end subroutine make_inp

      end module fciqmc_make_inp
