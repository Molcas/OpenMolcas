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
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Apr. 24, 2022, created this file.               *
! ****************************************************************

module hybridpdft
  use definitions, only: wp

  implicit none
  private

  logical :: Do_Hybrid = .false.
  real(kind=wp) :: Ratio_WF = 0.0d0
  Real(kind=wp) :: E_nohyb ! PDFT energy if it is not hybrid

  public :: do_hybrid, ratio_wf, e_nohyb
end module hybridpdft

