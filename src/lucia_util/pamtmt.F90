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
! Copyright (C) 1988, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine PAMTMT(X,T,NORB)
! GENERATE PER AKE'S T MATRIX FROM AN
! ORBITAL ROTATION MATRIX X
!
! T IS OBTAINED AS A STRICTLY LOWER TRIANGULAR
! MATRIX TL AND AN UPPER TRIANGULAR MATRIX TU
!
!         TL = 1 - L
!         TU = U ** -1
!
! WHERE L AND U ARISES FROM A LU DECOMPOSITION OF
! X :
!         X = L * U
! WITH L BEING A LOWER TRIANGULAR MATRIX WITH UNIT ON THE
! DIAGONAL AND U IS AN UPPER TRIANGULAR MATRIX
!
! JEPPE OLSEN OCTOBER 1988

use stdalloc, only: mma_allocate, mma_deallocate
use Index_Functions, only: nTri_Elem
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nOrb
real(kind=wp), intent(in) :: X(NORB,NORB)
real(kind=wp), intent(out) :: T(NORB,NORB)

integer(kind=iwp) :: I, ISING, J, KLL, KLU
real(kind=wp), allocatable:: SCR(:)

! Allocate local memory
Call mma_allocate(SCR,nOrb**2+nTri_Elem(nOrb),Label='SCR')
KLL = 1
KLU = KLL+nTri_Elem(NORB)
! LU factorize X
call LULU(X,SCR(KLL),SCR(KLU),NORB)
! Expand U to full matrix
T(:,:) = Zero
do I=1,NORB
  do J=I,NORB
    T(I,J) = SCR(KLU-1+nTri_Elem(J-1)+I)
  end do
end do
! Invert U
call INVMAT(T,SCR(KLU),NORB,ISING)
! Subtract L
do I=1,NORB
  do J=1,I-1
    T(I,J) = -SCR(KLL-1+nTri_Elem(I-1)+J)
  end do
end do
Call mma_deallocate(SCR)

end subroutine PAMTMT
