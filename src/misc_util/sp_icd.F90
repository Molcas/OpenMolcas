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
! Copyright (C) 2014, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  Sp_ICD
!
!> @ingroup Sparse
!> @brief
!>   Perform an incomplete Cholesky decomposition of sparse matrix \p A
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Perform an incomplete Cholesky decomposition of a symmetric matrix \p A, in sparse format.
!> On output, matrix \p B contains a lower triangular matrix such that:
!> \f$ A \simeq B B^\text{T} \f$
!> "Incomplete" means that only non-zero elements in \p A will be computed in \p B.
!> This is useful as a preconditioner.
!>
!> @param[in]  n   Size of the square matrices
!> @param[in]  A   Input matrix in sparse format
!> @param[in]  ija Index vector of matrix \p A
!> @param[out] B   Output matrix in sparse format
!> @param[out] ijb Index vector of matrix \p B
!***********************************************************************

subroutine Sp_ICD(n,A,ija,B,ijb)

use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: n, ija(*)
real(kind=wp), intent(in) :: A(*)
real(kind=wp), intent(_OUT_) :: B(*)
integer(kind=iwp), intent(_OUT_) :: ijb(*)
integer(kind=iwp) :: nijb, i, j, k, kk, kkb, l
real(kind=wp) :: Ljk
logical(kind=iwp) :: GoOn, Sym
real(kind=wp), parameter :: Thr = 1.0e-12_wp
integer(kind=iwp), external :: ip_of_Work

Sym = (A(n+1) > Zero)
if (ip_of_Work(A(1)) == ip_of_Work(B(1))) then
  if (.not. Sym) call SysAbendMsg('Sp_ICD','In-place decomposition only allowed with symmetric-stored matrix.','')
end if
nijb = n+1
ijb(1) = n+2
do i=1,n
  B(i) = A(i)

  ! Loop all elements in row i
  do k=ija(i),ija(i+1)-1
    j = ija(k)
    if (j < i) then
      nijb = nijb+1
      B(nijb) = A(k)
      ijb(nijb) = ija(k)

      ! Loop all previous elements of row i
      do kk=ija(i),k-1
        Ljk = Zero
        GoOn = .true.
        l = ijb(j)

        ! This loop to find an element in row j that belongs
        ! to the same column as each of the parent loop
        do while (GoOn)
          if (ijb(l) >= j) GoOn = .false.
          if (ijb(l) == ija(kk)) then
            Ljk = B(l)
            GoOn = .false.
          end if
          l = l+1
          if (l >= ijb(j+1)) GoOn = .false.
        end do

        ! As the number of elements per row in A and B can be
        ! different, the proper offset must be calculated
        kkb = kk-ija(i)+ijb(i)
        B(nijb) = B(nijb)-B(kkb)*Ljk
      end do
      if (B(j) > Thr) then
        B(nijb) = B(nijb)/B(j)
      else
        B(nijb) = Zero
      end if
      B(i) = B(i)-B(nijb)**2
    end if
  end do
  B(i) = sqrt(abs(B(i)))
  ijb(i+1) = nijb+1
end do

! The lower triangular matrix is not symmetric, obviously
B(n+1) = Zero

end subroutine Sp_ICD
