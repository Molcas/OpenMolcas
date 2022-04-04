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

! One diffuse, the other not diffuse.
subroutine ABOne(iLdiff,iLpoi,dMul,Ep,R,Rinv,Colle,lDiffA)

use Index_Functions, only: nTri_Elem1
use Constants, only: Zero, One, Two, Three, Four, Nine, OneHalf
use Definitions, only: wp, iwp, u6

! Maximum multipole implemented
#define _MxM_ 2

implicit none
integer(kind=iwp), intent(in) :: iLdiff, iLpoi
real(kind=wp), intent(in) :: dMul(nTri_Elem1(_MxM_)), Ep, R, Rinv
real(kind=wp), intent(out) :: Colle(3)
logical(kind=iwp), intent(in) :: lDiffA
real(kind=wp) :: DAMP, er, Ex, Pi1, Pi2, Sigma
real(kind=wp), parameter :: d3 = sqrt(Three)
#include "warnings.h"

! The omnipresent exponential and distance-exponent product.

er = Ep*R
Ex = exp(-Two*er)
Colle(:) = Zero

if ((iLdiff == 0) .and. (iLpoi == 0)) then
  ! s-s; see ABBoth for comments on sigma and similar below.

  Sigma = dMul(1)
  DAMP = (One+er)*Ex
  Colle(1) = Sigma*Rinv*(One-DAMP)

else if ((iLdiff == 0) .and. (iLpoi == 1)) then
  ! s-p

  Sigma = dMul(3)
  if (lDiffA) then
    Sigma = dMul(1)
  end if
  DAMP = (One+Two*er+Two*er**2)*Ex
  Colle(1) = Sigma*Rinv**2*(One-DAMP)

else if ((iLdiff == 0) .and. (iLpoi == 2)) then
  ! s-d

  Sigma = dMul(3)
  if (lDiffA) then
    Sigma = dMul(1)
  end if
  DAMP = (One+Two*er+Two*er**2+Four*er**3/Three)*Ex
  Colle(1) = Sigma*Rinv**3*(One-DAMP)

else if ((iLdiff == 1) .and. (iLpoi == 0)) then
  ! p-s

  Sigma = dMul(1)
  if (lDiffA) then
    Sigma = dMul(3)
  end if
  DAMP = (One+Two*er+Two*er**2+er**3)*Ex
  Colle(1) = Sigma*Rinv**2*(One-DAMP)

else if ((iLdiff == 1) .and. (iLpoi == 1)) then
  ! p-p

  Sigma = dMul(3)
  Pi1 = dMul(1)
  Pi2 = dMul(2)
  DAMP = (One+Two*er+Two*er**2+OneHalf*er**3+er**4)*Ex
  Colle(1) = Two*Sigma*Rinv**3*(One-DAMP)
  DAMP = (One+Two*er+Two*er**2+er**3)*Ex
  Colle(2) = Pi1*Rinv**3*(One-DAMP)
  Colle(3) = Pi2*Rinv**3*(One-DAMP)

else if ((iLdiff == 1) .and. (iLpoi == 2)) then
  ! p-d

  Sigma = dMul(3)
  Pi1 = dMul(2)
  Pi2 = dMul(4)
  if (lDiffA) then
    Pi1 = dMul(1)
    Pi2 = dMul(2)
  end if
  DAMP = (One+Two*er+Two*er**2+Four*er**3/Three+Two*er**4/Three+Four*er**5/Nine)*Ex
  Colle(1) = Three*Sigma*Rinv**4*(One-DAMP)
  DAMP = (One+Two*er+Two*er**2+Four*er**3/Three+Two*er**4/Three)*Ex
  Colle(2) = d3*Pi1*Rinv**4*(One-DAMP)
  Colle(3) = d3*Pi2*Rinv**4*(One-DAMP)
  !Colle = Colle1+Colle2+Colle3

else
  ! Higher moments.

  write(u6,*)
  write(u6,*) 'Too high momentum!'
  call Quit(_RC_IO_ERROR_READ_)
end if

return

end subroutine ABOne
