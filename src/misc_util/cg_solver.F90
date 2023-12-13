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
!  CG_Solver
!
!> @brief
!>   Solve a system of linear equations with a preconditioned conjugate gradients method
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Given a system of linear equations \f$ A x = b \f$, this routine attempts to solve
!> it with a preconditioned conjugate gradients method.
!> The matrix \p A should be symmetric positive definite.
!> The method is useful with a sparse matrix, where \f$ A x \f$ can be computed efficiently.
!> The preconditioner is an incomplete Cholesky decomposition, which becomes complete
!> with a dense matrix, and then it should converge at the first iteration.
!> If a dense matrix is used, the \p ija(1) value must be `0`.
!> On output, \p info contains the number of iterations needed for convergence,
!> or `-1` if it did not converge.
!>
!> @param[in]     n    Size of the system
!> @param[in]     nij  Size of \p A and \p ija
!> @param[in]     A    Input matrix in dense or sparse format
!> @param[in]     ija  Index vector of matrix \p A
!> @param[in]     b    Vector of independent terms
!> @param[in,out] x    Solution vector
!> @param[in,out] info Status info
!***********************************************************************

subroutine CG_Solver(n,nij,A,ija,b,x,info)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Ten
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, nij, ija(*)
real(kind=wp), intent(in) :: A(nij), b(n)
real(kind=wp), intent(inout) :: x(n)
integer(kind=iwp), intent(inout) :: info
integer(kind=iwp) :: k, maxk, recomp
real(kind=wp) :: alpha, beta, rr, Thr, RelThr
logical(kind=iwp) :: Sparse
integer(kind=iwp), allocatable :: ijLo(:), ijUp(:)
real(kind=wp), allocatable :: Ap(:), Lo(:), p(:), r(:), Up(:), y(:), z(:)
parameter(Thr=1.0e-20_wp)
real(kind=wp), external :: ddot_

call mma_allocate(r,n)
call mma_allocate(p,n)
call mma_allocate(Ap,n)
call mma_allocate(z,n)
call mma_allocate(y,n)

! Same algorithm, with different calls for sparse or dense matrices
Sparse = ija(1) > 0

maxk = max(10,n*n)
recomp = max(50,int(n/Ten))
r(:) = b

if (Sparse) then
  call mma_allocate(Lo,nij,label='Lo')
  call mma_allocate(ijLo,nij,label='ijLo')
  call mma_allocate(Up,nij,label='Up')
  call mma_allocate(ijUp,nij,label='ijUp')

  ! Compute the preconditioner
  call Sp_ICD(n,A,ija,Lo,ijLo)
  call Sp_Transpose(n,Lo,ijLo,Up,ijUp,nij)

  ! Initial guess: r = A*x-b
  call Sp_MV(n,-One,A,ija,x,One,r)
  call Sp_TriSolve(n,'L',Lo,ijLo,r,y)
  call Sp_TriSolve(n,'U',Up,ijUp,y,z)
  p(:) = z
  rr = DDot_(n,z,1,r,1)
  RelThr = Thr*max(rr,One)
  k = 1
  do while ((abs(rr) >= RelThr) .and. (k <= maxk))
    call Sp_MV(n,One,A,ija,p,Zero,Ap)
    alpha = rr/DDot_(n,p,1,Ap,1)
    x(:) = x+alpha*p
    beta = rr

    ! Recompute or update the residual
    if (mod(k,recomp) == 0) then
      r(:) = b
      call Sp_MV(n,-One,A,ija,x,One,r)
    else
      r(:) = r-alpha*Ap
    end if
    call Sp_TriSolve(n,'L',Lo,ijLo,r,y)
    call Sp_TriSolve(n,'U',Up,ijUp,y,z)
    rr = DDot_(n,z,1,r,1)
    p(:) = p*rr/beta+z
    k = k+1
  end do
  call mma_deallocate(Lo)
  call mma_deallocate(ijLo)
  call mma_deallocate(Up)
  call mma_deallocate(ijUp)
else
  call mma_allocate(Lo,nij,label='Lo')

  ! With a dense matrix, the preconditioner could be replaced with
  ! something else, otherwise this is just solving the system with
  ! a direct method
  Lo(1:n*n) = A(1:n*n)
  call dpotrf_('L',n,Lo,n,info)

  call DSyMV('L',n,-One,A,n,x,1,One,r,1)
  z(:) = r
  call DPoTrS('L',n,1,Lo,n,z,n,info)
  p(:) = z
  rr = DDot_(n,z,1,r,1)
  RelThr = Thr*max(rr,One)
  k = 1
  do while ((abs(rr) >= RelThr) .and. (k <= maxk))
    call dGeMV_('N',n,n,One,A,n,p,1,Zero,Ap,1)
    alpha = rr/DDot_(n,p,1,Ap,1)
    x(:) = x+alpha*p
    beta = rr
    if (mod(k,recomp) == 0) then
      r(:) = b
      call dGeMV_('N',n,n,-One,A,n,x,1,One,r,1)
    else
      r(:) = r-alpha*Ap
    end if
    z(:) = r
    call DPoTrS('L',n,1,Lo,n,z,n,info)
    rr = DDot_(n,z,1,r,1)
    p(:) = p*rr/beta+z
    k = k+1
  end do
  call mma_deallocate(Lo)
end if

! Set the return value
if (k <= maxk) then
  info = k
else
  info = -1
end if

call mma_deallocate(r)
call mma_deallocate(p)
call mma_deallocate(Ap)
call mma_deallocate(z)
call mma_deallocate(y)

end subroutine CG_Solver
