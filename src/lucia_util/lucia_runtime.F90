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
! Copyright (C) 2026, Meng Wang                                        *
!***********************************************************************

module lucia_runtime

implicit none
private

public :: LUCIA_OPTIMIZATIONS_ENABLED

logical, save :: Initialized = .false.
logical, save :: OptimizationsEnabled = .false.

contains

logical function LUCIA_OPTIMIZATIONS_ENABLED()
  character(len=32) :: Value
  integer :: Length, Status

  if (.not. Initialized) then
    Value = ''
    call get_environment_variable('MOLCAS_RASSCF_LUCIA_OPT',Value,length=Length,status=Status)

    if ((Status == -1) .or. (Length > len(Value))) then
      call SYSABENDMSG('lucia_util/lucia_runtime','MOLCAS_RASSCF_LUCIA_OPT value is too long','')
    else if ((Status > 0) .and. (Length > 0)) then
      call SYSABENDMSG('lucia_util/lucia_runtime','Could not read MOLCAS_RASSCF_LUCIA_OPT','')
    else if (Length > 0) then
      call UpCase(Value)
      select case (trim(adjustl(Value)))
      case ('1','ON','YES','TRUE')
        OptimizationsEnabled = .true.
      case ('0','OFF','NO','FALSE','')
        OptimizationsEnabled = .false.
      case default
        call SYSABENDMSG('lucia_util/lucia_runtime','Invalid MOLCAS_RASSCF_LUCIA_OPT value',trim(Value))
      end select
    end if

    Initialized = .true.
  end if

  LUCIA_OPTIMIZATIONS_ENABLED = OptimizationsEnabled
end function LUCIA_OPTIMIZATIONS_ENABLED

end module lucia_runtime
