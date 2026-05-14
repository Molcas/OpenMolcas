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
!********************************************************
!** Public-domain library routines used by casvb only. **
!********************************************************
!**********************
!** EISPACK ROUTINES **
!**********************

subroutine rs(nm,n,a,w,matz,z,ierr)
! this subroutine calls the recommended sequence of
! subroutines from the eigensystem subroutine package (eispack)
! to find the eigenvalues and eigenvectors (if desired)
! of a real symmetric matrix.
!
! on input
!
!    nm  must be set to the row dimension of the two-dimensional
!    array parameters as declared in the calling program
!    dimension statement.
!
!    n  is the order of the matrix  a.
!
!    a  contains the real symmetric matrix.
!
!    matz  is an integer variable set equal to zero if
!    only eigenvalues are desired.  otherwise it is set to
!    any non-zero integer for both eigenvalues and eigenvectors.
!
! on output
!
!    w  contains the eigenvalues in ascending order.
!
!    z  contains the eigenvectors if matz is not zero.
!
!    ierr  is an integer output variable set equal to an error
!       completion code described in the documentation for tqlrat
!       and tql2.  the normal completion code is zero.
!
! questions and comments should be directed to burton s. garbow,
! mathematics and computer science div, argonne national laboratory
!
! this version dated august 1983.
!
! Updated to Fortran 90+ (Sep. 2023)
! ----------------------------------------------------------------------

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nm, n, matz
real(kind=wp), intent(inout) :: a(nm,n)
real(kind=wp), intent(out) :: w(n), z(nm,n)
integer(kind=iwp), intent(out) :: ierr
real(kind=wp), allocatable :: fv1(:), fv2(:)

if (n > nm) then
  ierr = 10*n
  return
end if

call mma_allocate(fv1,n,label='fv1')
if (matz /= 0) then
  ! .......... find both eigenvalues and eigenvectors ..........
  call tred2(nm,n,a,w,fv1,z)
  call tql2(nm,n,w,fv1,z,ierr)
else
  ! .......... find eigenvalues only ..........
  call mma_allocate(fv2,n,label='fv2')
  call tred1(nm,n,a,w,fv1,fv2)
  ! tqlrat encounters catastrophic underflow on the Vax
  !call tqlrat(n,w,fv2,ierr)
  call tql1(n,w,fv1,ierr)
  call mma_deallocate(fv2)
end if
call mma_deallocate(fv1)

return

end subroutine rs
