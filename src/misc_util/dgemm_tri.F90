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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************
!  dGeMM_Tri
!
!> @brief
!>   Compute matrix product, store result in triangular storage format.
!>   May be used exactly as level 3 BLAS library routine ::DGEMM with \p C stored triangularly.
!> @author Thomas Bondo Pedersen
!>
!> @details
!> This subroutine computes the matrix product of matrices \p A and \p B
!> in a manner similar to the level 3 BLAS routine ::DGEMM. Unlike
!> ::DGEMM, however, the result matrix \p C is computed in triangular
!> storage (i.e., only upper triangle [identical to lower triangle
!> if \p C is symmetric] including diagonal elements). The argument
!> list is identical to that of ::DGEMM. Note, however, that \p m must
!> be equal to \p n. Moreover, \p ldC is not used here (it is included
!> only to have the same interface as ::DGEMM), but must be at least
!> ``1``.
!>
!> @param[in]     TransA Transposition of \p A
!> @param[in]     TransB Transposition of \p B
!> @param[in]     m      Row dimension of \p C
!> @param[in]     n      Column dimension of \p C
!> @param[in]     k      Dimension of summation index
!> @param[in]     alpha  Scale factor
!> @param[in]     A      Factor matrix \p A
!> @param[in]     ldA    Leading dimension of \p A
!> @param[in]     B      Factor matrix \p B
!> @param[in]     ldB    Leading dimension of\p  B
!> @param[in]     beta   Scale factor
!> @param[in,out] C      Result matrix
!> @param[in]     ldC    Leading dimension of \p C
!***********************************************************************

subroutine dGeMM_Tri(TransA,TransB,m,n,k,alpha,A,ldA,B,ldB,beta,C,ldC)

use Index_Functions, only: nTri_Elem
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
character, intent(in) :: TransA, TransB
integer(kind=iwp), intent(in) :: m, n, k, ldA, ldB, ldC
real(kind=wp), intent(in) :: alpha, A(ldA,*), B(ldB,*), beta
real(kind=wp), intent(inout) :: C(*)
integer(kind=iwp) :: iArg, j, kOff, ldARef, ldBRef
logical(kind=iwp) :: NoAB, NoC
character(len=2) :: ArgNum
character(len=*), parameter :: ArgErr = ' Illegal argument number ', SecNam = 'dGeMM_Tri'

! Test input parameters.
! ----------------------

if ((TransA == 'N') .or. (TransA == 'n')) then
  ldARef = m
else if ((TransA == 'T') .or. (TransA == 't')) then
  ldARef = k
else
  ArgNum = ' 1'
  call SysAbendMsg(SecNam,ArgErr,ArgNum)
  ldARef = 1
end if
ldARef = max(ldARef,1)

if ((TransB == 'N') .or. (TransB == 'n')) then
  ldBRef = k
else if ((TransB == 'T') .or. (TransB == 't')) then
  ldBRef = n
else
  ArgNum = ' 2'
  call SysAbendMsg(SecNam,ArgErr,ArgNum)
  ldBRef = 1
end if
ldBRef = max(ldBRef,1)

iArg = 0
if (m < 0) then
  iArg = 3
else if (n /= m) then
  iArg = 4
else if (k < 0) then
  iArg = 5
else if (ldA < ldARef) then
  iArg = 8
else if (ldB < ldBRef) then
  iArg = 10
else if (ldC < 1) then
  iArg = 13
end if
if (iArg /= 0) then
  write(ArgNum,'(I2)') iArg
  call SysAbendMsg(SecNam,ArgErr,ArgNum)
end if

! Scale C and return if possible.
! -------------------------------

NoC = n == 0
NoAB = (alpha == Zero) .or. (k == 0)
if (NoC .or. (NoAB .and. (beta == one))) return
if (beta == Zero) then
  C(1:nTri_Elem(n)) = Zero
else if (beta /= One) then
  C(1:nTri_Elem(n)) = beta*C(1:nTri_Elem(n))
end if
if (NoAB) return

! Compute using level 2 BLAS library (matrix-vector products).
! ------------------------------------------------------------

kOff = 1
if ((TransB == 'N') .or. (TransB == 'n')) then
  if ((TransA == 'N') .or. (TransA == 'n')) then
    do j=1,n
      call dGeMV_('N',j,k,alpha,A,ldA,B(:,j),1,One,C(kOff),1)
      kOff = kOff+j
    end do
  else
    do j=1,n
      call dGeMV_('T',k,j,alpha,A,ldA,B(:,j),1,One,C(kOff),1)
      kOff = kOff+j
    end do
  end if
else
  if ((TransA == 'N') .or. (TransA == 'n')) then
    do j=1,n
      call dGeMV_('N',j,k,alpha,A,ldA,B(j,1:k),1,One,C(kOff),1)
      kOff = kOff+j
    end do
  else
    do j=1,n
      call dGeMV_('T',k,j,alpha,A,ldA,B(j,1:k),1,One,C(kOff),1)
      kOff = kOff+j
    end do
  end if
end if

end subroutine dGeMM_Tri
