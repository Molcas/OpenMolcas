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
subroutine PAMTMT(X,T,SCR,NORB)
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

use Index_Functions, only: nTri_Elem
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nOrb
real(kind=wp), intent(in) :: X(NORB,NORB)
real(kind=wp), intent(out) :: T(NORB,NORB), SCR(nOrb**2+nTri_Elem(nOrb))
integer(kind=iwp) :: I, ISING, J, KLFREE, KLL, KLU

#ifdef _DEBUGPRINT_
write(u6,*) ' Welcome to PAMTMT'
write(u6,*) ' ================='
write(u6,*)
#endif
! Allocate local memory
KLFREE = 1
!KLL = KFLREE
KLL = KLFREE
KLFREE = KLL+nTri_Elem(NORB)
KLU = KLFREE
KLFREE = KLU+NORB**2
! LU factorize X
call LULU(X,SCR(KLL),SCR(KLU),NORB)
! Expand U to full matrix
T(:,:) = Zero
do I=1,NORB
  do J=I,NORB
    T(I,J) = SCR(KLU-1+nTri_Elem(J-1)+I)
  end do
end do
#ifdef _DEBUGPRINT_
write(u6,*) ' MATRIX TO BE INVERTED'
call WRTMAT(T,NORB,NORB,NORB,NORB)
#endif
! Invert U
call INVMAT(T,SCR(KLU),NORB,ISING)
#ifdef _DEBUGPRINT_
write(u6,*) ' INVERTED MATRIX'
call WRTMAT(T,NORB,NORB,NORB,NORB)
#endif
! Subtract L
do I=1,NORB
  do J=1,I-1
    T(I,J) = -SCR(KLL-1+nTri_Elem(I-1)+J)
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' INPUT X MATRIX'
call WRTMAT(X,NORB,NORB,NORB,NORB)
write(u6,*) ' T MATRIX'
call WRTMAT(T,NORB,NORB,NORB,NORB)
#endif

end subroutine PAMTMT
