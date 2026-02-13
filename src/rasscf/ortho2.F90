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

subroutine ORTHO2(S,U,V,N)
! Purpose: normalize vector U and calculate V=S*U.
!
! Called from: ORTHO1.
!
!      ****** IBM 3090 MOLCAS Release: 90 02 22 ******

use Constants, only: Zero, One
use Definitions, only: wp, u6

implicit none
integer N
real*8 S(*), U(*), V(*)
real*8 THR, SUM, X
real*8, external :: DDot_
integer I
#include "warnings.h"

THR = 1.0e-10_wp
if (N == 0) return
call DGEMM_('N','N',N,1,N,One,S,N,U,N,Zero,V,N)
SUM = DDOT_(N,U,1,V,1)
if (SUM < THR) then
  write(u6,*) ' TEST IN ORTHO2: N=',N
  write(u6,'(1X,5G16.8)') (U(I),I=1,N)
  write(u6,'(1X,5G16.8)') (V(I),I=1,N)
  write(u6,*) ' Error in ORTHO2. Norm=',SUM
  write(u6,*) ' RASSCF tried to orthonormalize orbitals, but'
  write(u6,*) ' failed due to a condition that should not be'
  write(u6,*) ' possible in a low-level subroutine. Either'
  write(u6,*) ' some extremely strange orbitals have been'
  write(u6,*) ' produced, or something is seriously wrong'
  write(u6,*) ' with the program. Please check, and consider'
  write(u6,*) ' issuing a bug report.'
  call QUIT(_RC_GENERAL_ERROR_)
end if
X = One/sqrt(SUM)
do I=1,N
  U(I) = X*U(I)
  V(I) = X*V(I)
end do

end subroutine ORTHO2
