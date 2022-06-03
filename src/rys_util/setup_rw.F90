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
! Copyright (C) 1990,1996, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

subroutine SetUp_RW(DoRys,nDiff)
!***********************************************************************
!                                                                      *
! Object: to setup tables for auxiliary functions to be used direct    *
!         in the recurrence relations of integral form or indirectly   *
!         to compute the Rys roots and weights which are used in the   *
!         recurrence relations of integrand form. For the lower order  *
!         Rys polynomials the roots and weight are computed from expa- *
!         nsion coefficients.                                          *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             January '90                                              *
!                                                                      *
!             Unified version August '96, RL.                          *
!***********************************************************************

use External_Centers, only: nOrdEF, XF
use Sizes_of_Seward, only: S
use Gateway_Info, only: GIAO
use Definitions, only: iwp

implicit none
logical(kind=iwp), intent(in) :: DoRys
integer(kind=iwp), intent(in) :: nDiff
integer(kind=iwp) :: iAng2, mRys

! Compute max sum of angular momentum index

iAng2 = 4*S%iAngMx

! Set up roots and weights for Hermite polynomials.

call SetHer(nDiff)

! Set up coefficients for Rys polynomials.

! 1) for two-electron integrals
! 2) for external field and nuclear attraction

mRys = (iAng2+2+nDiff)/2
if (allocated(XF) .or. (nOrdEF == 1) .or. GIAO) mRys = max(mRys,(2*S%iAngMx+1+2+nDiff)/2)
if (nOrdEF == 2) mRys = max(mRys,(2*S%iAngMx+2+2+nDiff)/2)
if (DoRys) call SetUpR(mRys)

return

end subroutine SetUp_RW
