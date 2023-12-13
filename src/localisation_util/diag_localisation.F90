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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Diag_Localisation(A,EVR,EVI,n,iGetVecs)
! Thomas Bondo Pedersen, Dec. 2005.
!
! Purpose: diagonalize a real general matrix A and return
!          the real and imaginary part of eigenvalues in EVR and
!          EVI. If iGetVecs=0 no eigenvectors are computed; else,
!          eigenvectors are returned in A (such that column one
!          corresponds to eigenvalue 1 in EVR and EVI). Uses XEIGEN
!          from the linalg_util directory (which uses LAPACK).
!          See LAPACK for details about the storage of real and
!          complex eigenvectors.

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n, iGetVecs
real(kind=wp), intent(inout) :: A(n,n)
real(kind=wp), intent(out) :: EVR(n), EVI(n)
integer(kind=iwp) :: iErr
real(kind=wp), allocatable :: Vecs(:,:)
character(len=*), parameter :: SecNam = 'Diag_Localisation'

call mma_allocate(Vecs,n,n,label='Vecs')

iErr = 0
call xEigen(iGetVecs,n,n,A,EVR,EVI,Vecs,iErr)
if (iErr /= 0) then
  write(u6,*) SecNam,': xEigen returned ',iErr
  call SysAbendMsg(SecNam,'Error in xEigen',' ')
end if

if (iGetVecs /= 0) A(:,:) = Vecs

call mma_deallocate(Vecs)

end subroutine Diag_Localisation
