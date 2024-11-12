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
! Copyright (C) 1992,2020, Roland Lindh                                *
!***********************************************************************

subroutine GetInf(DoRys,nDiff)
!***********************************************************************
!                                                                      *
! Object: to read all input information on the file INFO.              *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             January 1992                                             *
!***********************************************************************

use Real_Spherical, only: lMax_Internal, Sphere
use Her_RW, only: nPrp
use External_Centers, only: nOrdEF
use Gateway_global, only: Test
use DKH_Info, only: DKroll
use Sizes_of_Seward, only: S
use rctfld_module, only: lMax, PCM_Info_Get
use Definitions, only: iwp

implicit none
logical(kind=iwp), intent(in) :: DoRys
integer(kind=iwp), intent(inout) :: nDiff

! Load the dynamic input area.

call Get_Info_Dynamic()

! Load the static input area.

call Get_Info_Static()
!                                                                      *
!***********************************************************************
!                                                                      *
! Reaction field parameters

call PCM_Info_Get()

!                                                                      *
!***********************************************************************
!                                                                      *
! Generate the transformation matrices

if (S%iAngMx-1 >= lMax) then
  call Sphere(S%iAngMx)
  lmax_internal = S%iAngMx
else
  call Sphere(lMax)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Set the highest number of differentiations which will be
! applied to the basis functions. In this case 2 + 1 ( the
! kinetic energy operator and a differentiaion with respect
! to a nuclear coordinate.

nPrp = max(lMax,3)

! Setup of tables for coefficients of the Rys roots and weights.

if (S%iAngMx == 0) nDiff = 2
if (DKroll .and. (nOrdEF > 0)) nDiff = nDiff+nOrdEF
if (.not. Test) call Setup_RW(DoRys,nDiff)
!                                                                      *
!***********************************************************************
!                                                                      *
! Set up for contracted calculation

call Flip_Flop(.false.)
!                                                                      *
!***********************************************************************
!                                                                      *
call Get_EFP()
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine GetInf
