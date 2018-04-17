!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2017, Stefan Knecht                                    *
!***********************************************************************
  subroutine set_as_solver()

! module dependencies
  use active_space_solver_cfg
#ifdef _DMRG_
  use qcmaquis_interface_cfg
#endif

  implicit none
! ----------------------------------------------------------------------
  character(len=8)   :: is_solver
  character(len=2)   :: onoff
  integer            :: luspool
  character(len=180) :: line, line2, blank

  external           :: get_ln, isfreeunit
  character(len=180) :: get_ln
  integer            :: isfreeunit
! ----------------------------------------------------------------------


  luspool = 18
  luspool = isfreeunit(luspool)
  call spoolinp(luspool)

  rewind(luspool)
  call rdnlst(luspool,'DMRGSCF')
  blank=' '

  !> find line as_solver; next unblank line should be the active space solver
  do
    line = get_ln(luspool)
    call leftad(line)
    call upcase(line)

    if(line == blank) cycle

    if(line(1:4) == 'ACTI')then
      do
        line2 = get_ln(luspool)
        call leftad(line2)
        call upcase(line2)
        if(line2 /= blank)then
          is_solver(1:8) = line2(1:8)
          exit
        end if
      end do
      !exit
    end if
    if(line(1:4) == 'CIDE')then
      do
        line2 = get_ln(luspool)
        call leftad(line2)
        call upcase(line2)
        if(line2 /= blank)then
          onoff(1:2) = line2(1:2)
#ifdef _DMRG_
          if(onoff(1:2) == 'ON') dmrg_warmup%doCIDEAS  = .true.
#endif
          exit
        end if
      end do
      !exit
    end if
    if(line(1:4) == 'FIED')then
      do
        line2 = get_ln(luspool)
        call leftad(line2)
        call upcase(line2)
        if(line2 /= blank)then
          onoff(1:2) = line2(1:2)
#ifdef _DMRG_
          if(onoff(1:2) == 'ON') dmrg_warmup%doFIEDLER  = .true.
#endif
          exit
        end if
      end do
      !exit
    end if
    if(line(1:6) == 'END OF' .or. line(1:6) == 'DMRGSE' .or. line(1:6) == 'OPTIMI') exit

  end do


  !> set user-selected DMRG driver as active space solver
  call upcase(is_solver)

  select case(is_solver)

    !> QCMaquis
    case('QCMAQUIS')

      as_solver(1:8) = 'qcmaquis'

    !> BLOCK
    case('BLOCK')

      as_solver(1:8) = 'block   '

    !> CheMPS2
    case('CHEMPS2')

      as_solver(1:8) = 'chemps2 '

    case default
      write(6,*) 'unknown DMRG active space solver'
      call abend()
  end select

  write(6,'(/5x,a,a )') ' DMRGSCF: active space solver is set to ',as_solver(1:8)
  write(6,'( 5x,a,a/)') ' -------                                ','--------'

  if(as_solver(1:8) == 'qcmaquis') then
#ifdef _DMRG_
    doDMRG = .true.
#endif
  end if

  call Close_LuSpool(LuSpool)

  end subroutine set_as_solver
! ----------------------------------------------------------------------
