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

implicit real*8(A-H,O-Z)
integer nOrb
real*8 X(NORB,NORB), T(NORB,NORB)
real*8 SCR(nOrb**2+nOrb*(nOrb+1)/2)

NTEST = 0
if (NTEST >= 2) then
  write(6,*) ' Wellcome to PAMTMT'
  write(6,*) ' =================='
  write(6,*)
end if
! Allocate local memory
KLFREE = 1
!KLL = KFLREE
KLL = KLFREE
KLFREE = KLL+NORB*(NORB+1)/2
KLU = KLFREE
KLFREE = KLU+NORB**2
! LU factorize X
call LULU(X,SCR(KLL),SCR(KLU),NORB)
! Expand U to full matrix
call SETVEC(T,0.0d0,NORB**2)
do I=1,NORB
  do J=I,NORB
    T(I,J) = SCR(KLU-1+J*(J-1)/2+I)
  end do
end do
if (NTEST >= 100) then
  write(6,*) ' MATRIX TO BE INVERTED'
  call WRTMAT(T,NORB,NORB,NORB,NORB)
end if
! Invert U
call INVMAT(T,SCR(KLU),NORB,NORB,ISING)
if (NTEST >= 100) then
  write(6,*) ' INVERTED MATRIX'
  call WRTMAT(T,NORB,NORB,NORB,NORB)
end if
! Subtract L
do I=1,NORB
  do J=1,I-1
    T(I,J) = -SCR(KLL-1+I*(I-1)/2+J)
  end do
end do

if (NTEST >= 2) then
  write(6,*) ' INPUT X MATRIX'
  call WRTMAT(X,NORB,NORB,NORB,NORB)
  write(6,*) ' T MATRIX'
  call WRTMAT(T,NORB,NORB,NORB,NORB)
end if

end subroutine PAMTMT
