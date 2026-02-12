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

use output_ras, only: LF

implicit none
integer N
real*8 S(*), U(*), V(*)
real*8 THR, SUM, X
real*8, external :: DDot_
integer I
#include "warnings.h"

THR = 1.D-10
if (N == 0) return
call DGEMM_('N','N',N,1,N,1.0d0,S,N,U,N,0.0d0,V,N)
SUM = DDOT_(N,U,1,V,1)
if (SUM < THR) then
  write(LF,*) ' TEST IN ORTHO2: N=',N
  write(LF,'(1X,5G16.8)') (U(I),I=1,N)
  write(LF,'(1X,5G16.8)') (V(I),I=1,N)
  write(LF,*) ' Error in ORTHO2. Norm=',SUM
  write(LF,*) ' RASSCF tried to orthonormalize orbitals, but'
  write(LF,*) ' failed due to a condition that should not be'
  write(LF,*) ' possible in a low-level subroutine. Either'
  write(LF,*) ' some extremely strange orbitals have been'
  write(LF,*) ' produced, or something is seriously wrong'
  write(LF,*) ' with the program. Please check, and consider'
  write(LF,*) ' issuing a bug report.'
  call QUIT(_RC_GENERAL_ERROR_)
end if
X = 1.0d0/sqrt(SUM)
do I=1,N
  U(I) = X*U(I)
  V(I) = X*V(I)
end do

end subroutine ORTHO2
