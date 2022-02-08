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
use active_space_solver_cfg, only: as_solver
#ifdef _DMRG_
use qcmaquis_interface_cfg, only: dmrg_warmup
#endif
use Definitions, only: iwp, u6
use rasscf_data, only: doDMRG

implicit none
! ----------------------------------------------------------------------
character(len=8)   :: is_solver
character(len=2)   :: onoff
integer(kind=iwp)  :: luspool
character(len=180) :: line, line2

character(len=180), external :: get_ln
integer(kind=iwp), external  :: isfreeunit
! ----------------------------------------------------------------------

luspool = isfreeunit(18)
call spoolinp(luspool)

rewind(luspool)
call rdnlst(luspool,'DMRGSCF')

!> find line as_solver; next unblank line should be the active space solver
do
  line = get_ln(luspool)
  call upcase(line)
  line = adjustl(line)

  if (line == '') cycle

  if (line(1:4) == 'ACTI') then
    do
      line2 = get_ln(luspool)
      call upcase(line2)
      line2 = adjustl(line2)
      if (line2 /= '') then
        is_solver(1:8) = line2(1:8)
        exit
      end if
    end do
    !exit
  end if
  if (line(1:4) == 'CIDE') then
    do
      line2 = get_ln(luspool)
      call upcase(line2)
      line2 = adjustl(line2)
      if (line2 /= '') then
        onoff(1:2) = line2(1:2)
#       ifdef _DMRG_
        if (onoff(1:2) == 'ON') dmrg_warmup%doCIDEAS = .true.
#       endif
        exit
      end if
    end do
    !exit
  end if
  if (line(1:4) == 'FIED') then
    do
      line2 = get_ln(luspool)
      call upcase(line2)
      line2 = adjustl(line2)
      if (line2 /= '') then
        onoff(1:2) = line2(1:2)
#       ifdef _DMRG_
        if (onoff(1:2) == 'ON') dmrg_warmup%doFIEDLER = .true.
#       endif
        exit
      end if
    end do
    !exit
  end if
  if ((line(1:6) == 'END OF') .or. (line(1:6) == 'DMRGSE') .or. (line(1:6) == 'OOPTIM')) exit

end do

!> set user-selected DMRG driver as active space solver
call upcase(is_solver)

select case (is_solver)

  !> QCMaquis
  case ('QCMAQUIS')

    as_solver(1:8) = 'qcmaquis'

  !> BLOCK
  case ('BLOCK')

    as_solver(1:8) = 'block   '

  !> CheMPS2
  case ('CHEMPS2')

    as_solver(1:8) = 'chemps2 '

  case default
    write(u6,*) 'unknown DMRG active space solver'
    call abend()
end select

write(u6,'(/5x,a,a )') ' DMRGSCF: active space solver is set to ',as_solver(1:8)
write(u6,'( 5x,a,a/)') ' -------                                ','--------'

#ifdef _DMRG_
if (as_solver(1:8) == 'qcmaquis') doDMRG = .true.
#endif

call Close_LuSpool(luspool)

end subroutine set_as_solver
! ----------------------------------------------------------------------
