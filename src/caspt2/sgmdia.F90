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
! Copyright (C) 2011, Steven Vancoillie                                *
!***********************************************************************
!***********************************************************************
! Written by Steven Vancoillie, May 2011
! A set of subroutines that can handle RHS arrays in either a serial or
! parallel environment, depending on the situation.
!***********************************************************************
! --> when running serially, the RHS arrays are stored on LUSOLV and are
! loaded into the WORK array when needed.
! --> when running in parallel, the RHS arrays are stored on disk as
! disk resident arrays (DRAs) with filename RHS_XX_XX_XX, where XX is a
! number referring to the case, symmetry, and RHS vector respectively,
! and are loaded onto a global array when needed.
!***********************************************************************

subroutine sgmdia(nRow,nCol,W,LDW,dIn,dIs)

use caspt2_global, only: imag_shift, real_shift, sigma_p_epsilon, sigma_p_exponent
use Constants, only: Zero, One
use definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nRow, nCol, LDW
real(kind=wp), intent(inout) :: W(LDW,nCol)
real(kind=wp), intent(in) :: dIn(nRow), dIs(nCol)
integer(kind=iwp) :: i, j, p
real(kind=wp) :: delta, sigma, epsilon

do j=1,nCol
  do i=1,nRow
    ! energy denominator plus real shift
    delta = dIn(i)+dIs(j)+real_shift
    ! add the imaginary shift
    delta = delta+imag_shift**2/delta
    ! multiply by sigma-p regularizer
    epsilon = sigma_p_epsilon
    p = sigma_p_exponent
    if (epsilon > Zero) then
      sigma = One/epsilon**p
      delta = delta/(One-exp(-sigma*abs(delta)**p))
    end if

    W(i,j) = delta*W(i,j)
  end do
end do

end subroutine sgmdia
