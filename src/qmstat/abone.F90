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

implicit real*8(a-h,o-z)
parameter(MxMltp=2)
dimension dMul((MxMltp+1)*(MxMltp+2)/2), Colle(3)
logical lDiffA
#include "warnings.h"

! The omnipresent exponential and distance-exponent product.

er = Ep*R
Ex = exp((-2.0d0)*er)
d3 = sqrt(3.0d0)
do i=1,3
  Colle(i) = 0.0d0
end do

if ((iLdiff == 0) .and. (iLpoi == 0)) then
  ! s-s; see ABBoth for comments on sigma and similar below.

  Sigma = dMul(1)
  DAMP = (1.0d0+er)*Ex
  Colle(1) = Sigma*Rinv*(1.0d0-DAMP)

else if ((iLdiff == 0) .and. (iLpoi == 1)) then
  ! s-p

  Sigma = dMul(3)
  if (lDiffA) then
    Sigma = dMul(1)
  end if
  DAMP = (1.0d0+2.0d0*er+2.0d0*er**2)*Ex
  Colle(1) = Sigma*Rinv**2*(1.0d0-DAMP)

else if ((iLdiff == 0) .and. (iLpoi == 2)) then
  ! s-d

  Sigma = dMul(3)
  if (lDiffA) then
    Sigma = dMul(1)
  end if
  DAMP = (1.0d0+2.0d0*er+2.0d0*er**2+4.0d0*er**3/3.0d0)*Ex
  Colle(1) = Sigma*Rinv**3*(1.0d0-DAMP)

else if ((iLdiff == 1) .and. (iLpoi == 0)) then
  ! p-s

  Sigma = dMul(1)
  if (lDiffA) then
    Sigma = dMul(3)
  end if
  DAMP = (1.0d0+2.0d0*er+2.0d0*er**2+er**3)*Ex
  Colle(1) = Sigma*Rinv**2*(1.0d0-DAMP)

else if ((iLdiff == 1) .and. (iLpoi == 1)) then
  ! p-p

  Sigma = dMul(3)
  Pi1 = dMul(1)
  Pi2 = dMul(2)
  DAMP = (1.0d0+2.0d0*er+2.0d0*er**2+3.0d0*er**3/2.0d0+er**4)*Ex
  Colle(1) = 2.0d0*Sigma*Rinv**3*(1.0d0-DAMP)
  DAMP = (1.0d0+2.0d0*er+2.0d0*er**2+er**3)*Ex
  Colle(2) = Pi1*Rinv**3*(1.0d0-DAMP)
  Colle(3) = Pi2*Rinv**3*(1.0d0-DAMP)

else if ((iLdiff == 1) .and. (iLpoi == 2)) then
  ! p-d

  Sigma = dMul(3)
  Pi1 = dMul(2)
  Pi2 = dMul(4)
  if (lDiffA) then
    Pi1 = dMul(1)
    Pi2 = dMul(2)
  end if
  DAMP = (1.0d0+2.0d0*er+2.0d0*er**2+4.0d0*er**3/3.0d0+2.0d0*er**4/3.0d0+4.0d0*er**5/9.0d0)*Ex
  Colle(1) = 3.0d0*Sigma*Rinv**4*(1.0d0-DAMP)
  DAMP = (1.0d0+2.0d0*er+2.0d0*er**2+4.0d0*er**3/3.0d0+2.0d0*er**4/3.0d0)*Ex
  Colle(2) = d3*Pi1*Rinv**4*(1.0d0-DAMP)
  Colle(3) = d3*Pi2*Rinv**4*(1.0d0-DAMP)
  !Colle = Colle1+Colle2+Colle3

else
  ! Higher moments.

  write(6,*)
  write(6,*) 'Too high momentum!'
  call Quit(_RC_IO_ERROR_READ_)
end if

return

end subroutine ABOne
