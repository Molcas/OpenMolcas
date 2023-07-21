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
! Copyright (C) 2016, Giovanni Li Manni                                *
!***********************************************************************

subroutine DecoMat(MAT,dimens,eigenvec,NumV,rc)
!************ by G. Li Manni Stuttgart March 2016 *************
!
! MAT:      copy of the D1A matrix (1RDM) in AO basis.
!           It will be destroyed in eigen_molcas
!
! DIMENS:   dimension of the square MAT = nbas(isym)
!
! EIGENVEC: First it is used as scratch inside eigen_molcas. Next it will contains eigenvectors
!           of the eigen-decomposed matrix in Cholesky form. Instead of eigenvectors X such that MAT = XDX^t,
!           eigenvectors are presented as Y = X(D**0.5) such that MAT = YY^t.
!           Negative eigenvalues are set to zero and eigenvalue larger than 2.0d0 set to 2.0d0
!
! NUMV    : Number of non negative eigenvalues

implicit none
integer dimens, NumV, rc, i, j
real*8 MAT(dimens,dimens), eigenvec(dimens,dimens)
real*8 eigenval(dimens)
character*12 routine
parameter(routine='DecoNegatMat')

rc = 0
NumV = 0

if (dimens < 1) then
  rc = -1
  write(6,*) 'matrix size < 1'
  Go To 10
end if

! Step 1: diagonalize MAT. MAT is destroyed and replaced by the eigenvectors values.
! IMPORTANT: At this stage eigenvec is used as scratch.
call eigen_molcas(dimens,MAT,eigenval,eigenvec)
! Move MAT to eigenvec. where it belongs. I could not wait!
call dcopy_(dimens**2,MAT,1,eigenvec,1)
! Set to zero negative eigenvalue and to TWO values larger than 2.0d0.
! Count only the positive ones (counter NumV)
do j=1,dimens
  if (eigenval(j) > 1.0d-12) then
    NumV = NumV+1
    if (eigenval(j) > 2.0d0) eigenval(j) = 2.0d0
  else
    eigenval(j) = 0.0d0
  end if
end do
! Sort eigenvalues in decreasing order of occupation number
call IncrSort(eigenval,eigenvec,dimens)
! Compute sqrt(eigenvalues)
do i=1,dimens
  eigenval(i) = sqrt(eigenval(i))
end do
! Generate Y = X(D**0.5)
do j=1,dimens
  do i=1,dimens
    eigenvec(i,j) = eigenvec(i,j)*eigenval(j)
  end do
end do
!***************** Exit ****************
10 continue
return

end subroutine DecoMat
