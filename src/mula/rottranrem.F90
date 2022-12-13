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
! Copyright (C) 1995, Niclas Forsberg                                  *
!***********************************************************************

subroutine RotTranRem(Sinv,S,Mass,AtCoord,NumOfAt,NumInt)
!  Purpose:
!    Project total translation and total rotation out of inverted S matrix.
!
!  Input:
!    S        : Real three dimensional array.
!    Mass     : Real array - masses of the atoms.
!    AtCoord  : Real two dimensional array - coordinates of the atoms.
!
!  Output:
!    Sinv     : Real three dimensional array - Inverted S matrix with rotation and translation projected out.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, uToau
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NumOfAt, NumInt
real(kind=wp), intent(out) :: Sinv(3,NumOfAt,NumInt)
real(kind=wp), intent(in) :: S(3,NumOfAt,NumInt), Mass(NumOfAt), AtCoord(3,NumOfAt)
integer(kind=iwp) :: i, iAtom, j, k, n
real(kind=wp) :: Det, X
real(kind=wp), allocatable :: Ainv(:,:), Amat(:,:), AmatMass(:,:), Stemp(:,:,:), Temp1(:,:), Temp2(:,:)
integer(kind=iwp), parameter :: nFree = 6

! Initialize.
n = 3*NumOfAt
call mma_allocate(Amat,n,nFree,label='Amat')
Amat(:,:) = Zero

! Invert S matrix.
call mma_allocate(Temp2,NumInt,NumInt,label='Temp2')

call DGEMM_('T','N',NumInt,NumInt,3*NumOfAt,One,S,3*NumOfAt,S,3*NumOfAt,Zero,Temp2,NumInt)
! Invert, by solving eq Temp2*X=Temp1. Solution computed in-place.
! Thus, Temp1=Unit matrix(in) and contains solution (out).
call mma_allocate(Temp1,NumInt,NumInt,label='Temp1')

call unitmat(Temp1,NumInt)
call Dool_MULA(Temp2,NumInt,NumInt,Temp1,NumInt,NumInt,det)
!PAM01 Replacement for Dool_MULA, if superstable solution is wanted:
!Eps = 1.0e-8_wp
!call SymSolve(Temp2,Temp1,Temp3,Temp4,Eps)
!Temp1(:,:) = Temp3
!PAM01 End of replacement code.

call DGEMM_('N','N',3*NumOfAt,NumInt,NumInt,One,S,3*NumOfAt,Temp1,NumInt,Zero,Sinv,3*NumOfAt)
call mma_deallocate(Temp1)
call mma_deallocate(Temp2)

! Pure translation.
k = 0
do iAtom=1,NumOfAt
  Amat(k+1,1) = One
  Amat(k+2,2) = One
  Amat(k+3,3) = One
  k = k+3
end do

! Pure rotation.
k = 0
do iAtom=1,NumOfAt
  Amat(k+1,4) = -AtCoord(2,iAtom)
  Amat(k+2,4) = AtCoord(1,iAtom)
  Amat(k+3,4) = Zero
  Amat(k+1,5) = AtCoord(3,iAtom)
  Amat(k+2,5) = Zero
  Amat(k+3,5) = -AtCoord(1,iAtom)
  Amat(k+1,6) = Zero
  Amat(k+2,6) = -AtCoord(3,iAtom)
  Amat(k+3,6) = AtCoord(2,iAtom)
  k = k+3
end do

! Scale with the mass of the atom.
n = 3*NumOfAt
call mma_allocate(AmatMass,n,nFree,label='AmatMass')
AmatMass(:,:) = Amat
!PAM04: Replace the following code section...
!do i=1,nFree
!  Acol => AmatMass(:,i)
!  k = 1
!  do j=1,NumOfAt
!    jMass = (k+2)/3
!    Acol(k) = uToau*Mass(jMass)*Acol(k)
!    Acol(k+1) = uToau*Mass(jMass)*Acol(k+1)
!    Acol(k+2) = uToau*Mass(jMass)*Acol(k+2)
!    k = k+3
!  end do
!end do
!PAM04: ... with the following...
do i=1,nFree
  do j=1,NumOfAt
    X = uToau*Mass(j)
    AmatMass(1+3*(j-1),i) = X*AmatMass(1+3*(j-1),i)
    AmatMass(2+3*(j-1),i) = X*AmatMass(2+3*(j-1),i)
    AmatMass(3+3*(j-1),i) = X*AmatMass(3+3*(j-1),i)
  end do
end do
!PAM04: ... until here.

! Project rotation and translation out of S matrix.
call mma_allocate(Temp1,nFree,nFree,label='Temp1')
call DGEMM_('T','N',nFree,nFree,3*NumOfAt,One,AmatMass,3*NumOfAt,Amat,3*NumOfAt,Zero,Temp1,nFree)
call mma_allocate(Ainv,nFree,nFree,label='Ainv')
call unitmat(Ainv,nFree)
call Dool_MULA(Temp1,nFree,nFree,Ainv,nFree,nFree,det)
call mma_deallocate(Temp1)

call mma_allocate(Temp2,nFree,NumInt,label='Temp2')
call DGEMM_('T','N',nFree,NumInt,3*NumOfAt,One,AmatMass,3*NumOfAt,Sinv,3*NumOfAt,Zero,Temp2,nFree)
call mma_deallocate(AmatMass)

call mma_allocate(Temp1,nFree,NumInt,label='Temp1')

call DGEMM_('N','N',nFree,NumInt,nFree,One,Ainv,nFree,Temp2,nFree,Zero,Temp1,nFree)
call mma_deallocate(Ainv)
call mma_deallocate(Temp2)

call mma_allocate(Stemp,3,NumOfAt,NumInt,label='Stemp')
call DGEMM_('N','N',3*NumOfAt,NumInt,nFree,One,Amat,3*NumOfAt,Temp1,nFree,Zero,Stemp,3*NumOfAt)
call mma_deallocate(Amat)
call mma_deallocate(Temp1)

Sinv(:,:,:) = Sinv-Stemp

call mma_deallocate(Stemp)

end subroutine RotTranRem
