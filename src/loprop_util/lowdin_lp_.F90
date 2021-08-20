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

subroutine Lowdin_LP_(S,Eval,C,nDim,nDim2,Blk)
! S: full-storage overlap matrix (it will be destroyed!)
! C: on exit, the S^-1/2 matrix

implicit real*8(A-H,O-Z)
#include "real.fh"
real*8 C(nDim,nDim), S(nDim,nDim), Blk(nDim,nDim),Eval(nDim*(nDim+1)/2)
data DIAGTH,DANGER/1.0D-12,1.0d3/

! diagth  threshold for matrix diagonalization used in
!         subroutine jacobi.  in jacobi, this  constant
!         is called doneth.
! danger  criterion for deciding that the job should be
!         aborted due to numerical problems caused by near
!         linear dependencies in the basis set.
!         all eigenvalues of the weighted overlap
!         matrix must be greater than DIAGTH*DANGER.

toosml = diagth*danger

! diagonalize overlap matrix

call Jacob(Eval,Blk,nDim,nDim)

!call RecPrt('Blk',' ',Blk,nDim,nDim)
!call TriPrt('Eval',' ',Eval,nDim)

! form the inverse sqrt of the overlap matrix of the vectors:
! (avoid numerical problems of linear dependence (too small eigenvals)
! by prescreening
! the  eigenvalues)

do i=1,nDim
  eigenv = eval(i*(i+1)/2)
  if (eigenv < toosml) go to 900
  eval(i*(i+1)/2) = ONE/sqrt(eigenv)
end do

! Compute S^-1/2

do i=1,nDim
  do j=1,i
    sij = zero
    do k=1,nDim
      sij = sij+eval(k*(k+1)/2)*blk(i,k)*blk(j,k)
    end do
    C(i,j) = sij
    C(j,i) = sij
  end do
end do

!call RecPrt('C',' ',C,nDim,nDim)

return

900 write(6,910) eigenv,toosmL
910 format(/1X,'An eigenvalue of the overlap matrix of the symmetrized Jacobi transf. matrix of ',E13.5,' has been found.'/1X, &
            'This is lower than the allowed threshold of ',E13.5)
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(S)
  call Unused_integer(nDim2)
end if

end subroutine Lowdin_LP_
