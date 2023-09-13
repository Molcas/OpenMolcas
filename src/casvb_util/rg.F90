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

subroutine rg(nm,n,a,wr,wi,matz,z,iv1,fv1,ierr)
! this subroutine calls the recommended sequence of
! subroutines from the eigensystem subroutine package (eispack)
! to find the eigenvalues and eigenvectors (if desired)
! of a real general matrix.
!
! on input
!
!    nm  must be set to the row dimension of the two-dimensional
!    array parameters as declared in the calling program
!    dimension statement.
!
!    n  is the order of the matrix  a.
!
!    a  contains the real general matrix.
!
!    matz  is an integer variable set equal to zero if
!    only eigenvalues are desired.  otherwise it is set to
!    any non-zero integer for both eigenvalues and eigenvectors.
!
! on output
!
!    wr  and  wi  contain the real and imaginary parts,
!    respectively, of the eigenvalues.  complex conjugate
!    pairs of eigenvalues appear consecutively with the
!    eigenvalue having the positive imaginary part first.
!
!    z  contains the real and imaginary parts of the eigenvectors
!    if matz is not zero.  if the j-th eigenvalue is real, the
!    j-th column of  z  contains its eigenvector.  if the j-th
!    eigenvalue is complex with positive imaginary part, the
!    j-th and (j+1)-th columns of  z  contain the real and
!    imaginary parts of its eigenvector.  the conjugate of this
!    vector is the eigenvector for the conjugate eigenvalue.
!
!    ierr  is an integer output variable set equal to an error
!       completion code described in the documentation for hqr
!       and hqr2.  the normal completion code is zero.
!
!    iv1  and  fv1  are temporary storage arrays.
!
! questions and comments should be directed to burton s. garbow,
! mathematics and computer science div, argonne national laboratory
!
! this version dated august 1983.
!
! ----------------------------------------------------------------------

integer n, nm, is1, is2, ierr, matz
real*8 a(nm,n), wr(n), wi(n), z(nm,n), fv1(n)
integer iv1(n)

if (n <= nm) go to 10
ierr = 10*n
go to 50

10 call balanc(nm,n,a,is1,is2,fv1)
call elmhes(nm,n,is1,is2,a,iv1)
if (matz /= 0) go to 20
! .......... find eigenvalues only ..........
call hqr(nm,n,is1,is2,a,wr,wi,ierr)
go to 50
! .......... find both eigenvalues and eigenvectors ..........
20 call eltran(nm,n,is1,is2,a,iv1,z)
call hqr2(nm,n,is1,is2,a,wr,wi,z,ierr)
if (ierr /= 0) go to 50
call balbak(nm,n,is1,is2,fv1,n,z)

50 return

end subroutine rg
