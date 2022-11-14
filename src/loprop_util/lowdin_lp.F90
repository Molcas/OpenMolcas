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

subroutine Lowdin_LP(S,C,nDim)
! S: full-storage overlap (S) matrix
! C: on exit, the S^-1/2 matrix

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nDim
real(kind=wp), intent(in) :: S(nDim,nDim)
real(kind=wp), intent(out) :: C(nDim,nDim)
integer(kind=iwp) :: i, j, k
real(kind=wp) :: eigenv, sij, toosml
real(kind=wp), allocatable :: Blk(:,:), Eval(:)
real(kind=wp), parameter :: DIAGTH = 1.0e-12_wp, DANGER = 1.0e3_wp

! diagth  threshold for matrix diagonalization used in
!         subroutine jacobi.  in jacobi, this  constant
!         is called doneth.
! danger  criterion for deciding that the job should be
!         aborted due to numerical problems caused by near
!         linear dependencies in the basis set.
!         all eigenvalues of the weighted overlap
!         matrix must be greater than DIAGTH*DANGER.

toosml = diagth*danger

call mma_allocate(Eval,nDim*(nDim+1)/2,label='Eval')
call mma_allocate(Blk,nDim,nDim,label='Blk')
call unitmat(Blk,nDim)
do i=1,nDim
  do j=1,i
    Eval(i*(i-1)/2+j) = S(i,j)
  end do
end do

! diagonalize overlap matrix

call Jacob(Eval,Blk,nDim,nDim)

!call RecPrt('Blk',' ',Blk,nDim,nDim)
!call TriPrt('Eval',' ',Eval,nDim)

! form the inverse sqrt of the overlap matrix of the vectors:
! (avoid numerical problems of linear dependence (too small eigenvalues)
! by prescreening
! the  eigenvalues)

do i=1,nDim
  eigenv = Eval(i*(i+1)/2)
  if (eigenv < toosml) then
    write(u6,910) eigenv,toosmL
    return
  end if
  Eval(i*(i+1)/2) = One/sqrt(eigenv)
end do

! Compute S^-1/2

do i=1,nDim
  do j=1,i
    sij = Zero
    do k=1,nDim
      sij = sij+Eval(k*(k+1)/2)*Blk(i,k)*Blk(j,k)
    end do
    C(i,j) = sij
    C(j,i) = sij
  end do
end do

call mma_deallocate(Eval)
call mma_deallocate(Blk)

!call RecPrt('C',' ',C,nDim,nDim)

return

910 format(/1X,'An eigenvalue of the overlap matrix of the symmetrized Jacobi transf. matrix of ',E13.5,' has been found.'/1X, &
            'This is lower than the allowed threshold of ',E13.5)

end subroutine Lowdin_LP
