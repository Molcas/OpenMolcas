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

subroutine DIAG_R2_RASSI(MATRIX,NBTOT,INFO,W1,Z1)
! THIS ROUTINE PERFORMS THE DIAGONALIZATION OF A REAL SQUARE
! MATRIX WITH THE DIMENSION NBTOT. THE EIGENVALUES OF THE DIAGONALIZATION
! ARE DIRECTED INTO W1 AND THE REAL EIGENVECTORS ARE WRITTEN TO Z1.

use Constants, only: Zero

implicit none
integer INFO, I, J, NBTOT
real*8 AP(Nbtot*(Nbtot+1)/2), WORK(3*Nbtot), W1(NBTOT), W(Nbtot), Z(Nbtot,Nbtot), Z1(NBTOT,NBTOT), MATRIX(NBTOT,NBTOT)

! initializations
INFO = 0
AP(:) = Zero
WORK(:) = Zero
W1(:) = Zero
W(:) = Zero
Z1(:,:) = Zero
Z(:,:) = Zero

do j=1,Nbtot
  do i=1,j
    AP(i+(j-1)*j/2) = MATRIX(i,j)
  end do
end do

call dspev_('V','U',Nbtot,AP,W,Z,Nbtot,WORK,INFO)
do I=1,Nbtot
  W1(I) = W(I)
  do J=1,Nbtot
    Z1(J,I) = Z(J,I)
  end do
end do

end subroutine DIAG_R2_RASSI
