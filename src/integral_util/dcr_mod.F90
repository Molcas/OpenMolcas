!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

module dcr_mod

use Definitions, only: iwp

implicit none
private

integer(kind=iwp), parameter :: MxIndx = 50, Mx = MxIndx*(MxIndx+1)/2
integer(kind=iwp) :: iDCR_all(0:7,Mx), Indx(MxIndx), Lambda_all(Mx), mDCR_all(Mx), nIndx
logical(kind=iwp) :: Done(Mx)

public :: DCR_Init, Done, iDCR_all, Indx, Lambda_all, mDCR_all, nIndx

contains

subroutine DCR_Init()

  nIndx = 0
  Done(:) = .false.

end subroutine DCR_Init

end module dcr_mod
