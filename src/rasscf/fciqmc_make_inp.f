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
        use stdalloc, only : mma_deallocate
        implicit none
        private
        public :: make_inp, cleanup
        integer, public ::
! No default value on purpose
     &    totalwalkers,
! &    calcrdmonfly(2),
! &    rdmsamplingiters,
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
        integer, allocatable, public ::
     &    definedet(:)
        real*8, public ::
     &    proje_changeref = 1.2d0,
     &    max_tau = 0.02d0,
     &    memoryfacpart = 5.0d0,
     &    memoryfacspawn = 10.0d0,
! Default value for NECI RealSpawnCutOff
     &    realspawncutoff = 0.3d0,
! Default value for NECI diagonal shift value
     &    diagshift = 0.00d0,
     &    shiftdamp = 0.02d0

        type, public :: t_RDMsampling
          sequence
          integer :: start, n_samples, step
        end type t_RDMsampling

        type(t_RDMsampling), public :: RDMsampling
        save
      contains

!>  @brief
!>    Generate an Input for NECI
!>
!>  @author
!>    G. Li Manni, Oskar Weser
!>
!>  @paramin[in] readpops  If true the readpops option for NECI is set.
      subroutine make_inp(path, readpops, doGAS)
      use general_data, only : nActEl, iSpin
      use stdalloc, only : mma_deallocate
      use fortran_strings, only : str
      character(*), intent(in) :: path
      logical, intent(in), optional :: readpops, doGAS
      logical :: readpops_, doGAS_
      integer :: i, isFreeUnit, file_id, indentlevel
      integer, parameter :: indentstep = 4

      if (present(readpops)) then
        readpops_ = readpops
      else
        readpops_ = .false.
      end if
      if (present(doGAS)) then
        doGAS_ = doGAS
      else
        doGAS_ = .false.
      end if

      call add_info('Default number of total walkers',
     &  [dble(totalwalkers)], 1, 6)
      call add_info('Default number of cycles ',[dble(nmcyc)],1,6)
      call add_info('Default value for Time  ',[dble(Time)],1,6)

      call qEnter('make_inp')

      file_id = isFreeUnit(39)
      call Molcas_Open(file_id, path)

      indentlevel = 0
      write(file_id, A_fmt()) 'Title'
      write(file_id, *)
      write(file_id, A_fmt()) 'System read'
      call indent()
        write(file_id, I_fmt()) 'electrons ', nActEl
        write(file_id,A_fmt()) 'nonuniformrandexcits 4ind-weighted-2'
        write(file_id,A_fmt()) 'nobrillouintheorem'
        if(iSpin /= 1) then
          write(file_id, I_fmt()) 'spin-restrict', iSpin - 1
        end if
        write(file_id, A_fmt()) 'freeformat'
        if (doGas_) write(file_id, A_fmt()) 'part-conserving-gas'
      call dedent()
      write(file_id, A_fmt()) 'endsys'
      write(file_id, *)
      write(file_id, A_fmt()) 'calc'
      call indent()
        if (allocated(DefineDet)) then
          write(file_id, kw_fmt('('//str(nActEl)//'I5)'))
     &        'definedet', (definedet(i), i = 1,nActEl)
          write(file_id, A_fmt()) ''
        end if
        write(file_id, *)
        write(file_id, I_fmt()) 'totalwalkers', totalwalkers
        write(file_id, A_fmt()) merge(' ', '(', readpops_)//'readpops'
        write(file_id, A_fmt())merge(' ', '(',readpops_)//'walkcontgrow'
        write(file_id, I_fmt()) 'semi-stochastic', semi_stochastic
        write(file_id, *)
        write(file_id, A_fmt()) 'methods'
        call indent()
          write(file_id, A_fmt()) 'method vertex fcimc'
        call dedent()
        write(file_id, A_fmt()) 'endmethods'
        write(file_id, *)
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
        write(file_id, I_fmt())
      call dedent()
      write(file_id, A_fmt()) 'endcalc'
      write(file_id, *)
      write(file_id, A_fmt()) 'logging'
      call indent()
        write(file_id, I_fmt()) 'highlypopwrite', Highlypopwrite
        write(file_id, A_fmt()) 'print-spin-resolved-RDMs'
        write(file_id, A_fmt()) 'hdf5-pops'
        write(file_id, A_fmt()) 'printonerdm'
       write(file_id,'('//str(indentlevel)//'x, A,1x,I0,1x,I0,1x,I0)')
     &     'RDMlinspace',
     &      RDMsampling%start, RDMsampling%n_samples, RDMsampling%step
      call dedent()
      write(file_id, A_fmt()) 'endlog'
      write(file_id, A_fmt()) 'end'


      close(file_id)
      call qExit('make_inp')

      contains

        function indent_fmt() result(res)
          character(:), allocatable :: res
          if (indentlevel /= 0) then
            res = str(indentlevel)//'x, '
          else
            res = ''
          end if
        end function

        function kw_fmt(value_fmt) result(res)
          character(*), intent(in) :: value_fmt
          character(:), allocatable :: res
          res = '('//indent_fmt()//'A, 1x, '//value_fmt//')'
        end function

        function I_fmt() result(res)
          character(:), allocatable :: res
          res = kw_fmt('I0')
        end function

        function R_fmt() result(res)
          character(:), allocatable :: res
          res = kw_fmt('F0.2')
        end function

        function A_fmt() result(res)
          character(:), allocatable :: res
          res = kw_fmt('A')
        end function

        subroutine indent()
          indentlevel = indentlevel + indentstep
        end subroutine

        subroutine dedent()
          indentlevel = indentlevel - indentstep
        end subroutine
      end subroutine make_inp

      subroutine cleanup()
        if (allocated(definedet)) call mma_deallocate(definedet)
      end subroutine cleanup

      end module fciqmc_make_inp
