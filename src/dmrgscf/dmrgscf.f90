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
  subroutine dmrgscf(iReturn)

  use rasscf_data, only: doDMRG
  implicit none

  integer, intent(inout) :: iReturn
! ----------------------------------------------------------------------
  !> set DMRG driver as active space solver
  call set_as_solver()

  !> read DMRG settings (driver-specific input)
  call set_dmrg_settings()

  !> call wave function optimizer
  iReturn = 0
  call rasscf(iReturn)

#ifdef _DMRG_
  !> reset in case we call RASSCF (or RASSI or CASPT2) afterwards requesting a CI driver
  if(doDMRG) doDMRG = .false.
#endif

  end subroutine dmrgscf
! ----------------------------------------------------------------------
