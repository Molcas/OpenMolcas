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
integer(kind=iwp) :: n, nij, ija(*), info
real(kind=wp) :: A(nij), b(n), x(n)
#include "WrkSpc.fh"
integer(kind=iwp) :: ipLo, ipUp, ipijLo, ipijUp, k, maxk, recomp
real(kind=wp) :: alpha, beta, rr, Thr, RelThr
logical(kind=iwp) :: Sparse
real(kind=wp), allocatable :: Ap(:), p(:), r(:), y(:), z(:)
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
call dcopy_(n,b,1,r,1)

if (Sparse) then
  call Allocate_Work(ipLo,nij)
  call Allocate_Work(ipUp,nij)
  call Allocate_iWork(ipijLo,nij)
  call Allocate_iWork(ipijUp,nij)

  ! Compute the preconditioner
  call Sp_ICD(n,A,ija,Work(ipLo),iWork(ipijLo))
  call Sp_Transpose(n,Work(ipLo),iWork(ipijLo),Work(ipUp),iWork(ipijUp),nij)

  ! Initial guess: r = A*x-b
  call Sp_MV(n,-One,A,ija,x,One,r)
  call Sp_TriSolve(n,'L',Work(ipLo),iWork(ipijLo),r,y)
  call Sp_TriSolve(n,'U',Work(ipUp),iWork(ipijUp),y,z)
  call dcopy_(n,z,1,p,1)
  rr = DDot_(n,z,1,r,1)
  RelThr = Thr*max(rr,One)
  k = 1
  do while ((abs(rr) >= RelThr) .and. (k <= maxk))
    call Sp_MV(n,One,A,ija,p,Zero,Ap)
    alpha = rr/DDot_(n,p,1,Ap,1)
    call daxpy_(n,alpha,p,1,x,1)
    beta = rr

    ! Recompute or update the residual
    if (mod(k,recomp) == 0) then
      call dcopy_(n,b,1,r,1)
      call Sp_MV(n,-One,A,ija,x,One,r)
    else
      call daxpy_(n,-alpha,Ap,1,r,1)
    end if
    call Sp_TriSolve(n,'L',Work(ipLo),iWork(ipijLo),r,y)
    call Sp_TriSolve(n,'U',Work(ipUp),iWork(ipijUp),y,z)
    rr = DDot_(n,z,1,r,1)
    call dscal_(n,rr/beta,p,1)
    call daxpy_(n,One,z,1,p,1)
    k = k+1
  end do
  call Free_Work(ipLo)
  call Free_Work(ipUp)
  call Free_iWork(ipijLo)
  call Free_iWork(ipijUp)
else
  call Allocate_Work(ipLo,nij)

  ! With a dense matrix, the preconditioner could be replaced with
  ! something else, otherwise this is just solving the system with
  ! a direct method
  call dcopy_(n*n,A,1,Work(ipLo),1)
  call dpotrf_('L',n,Work(ipLo),n,info)

  call DSyMV('L',n,-One,A,n,x,1,One,r,1)
  call dcopy_(n,r,1,z,1)
  call DPoTrS('L',n,1,Work(ipLo),n,z,n,info)
  call dcopy_(n,z,1,p,1)
  rr = DDot_(n,z,1,r,1)
  RelThr = Thr*max(rr,One)
  k = 1
  do while ((abs(rr) >= RelThr) .and. (k <= maxk))
    call dGeMV_('N',n,n,One,A,n,p,1,Zero,Ap,1)
    alpha = rr/DDot_(n,p,1,Ap,1)
    call daxpy_(n,alpha,p,1,x,1)
    beta = rr
    if (mod(k,recomp) == 0) then
      call dcopy_(n,b,1,r,1)
      call dGeMV_('N',n,n,-One,A,n,x,1,One,r,1)
    else
      call daxpy_(n,-alpha,Ap,1,r,1)
    end if
    call dcopy_(n,r,1,z,1)
    call DPoTrS('L',n,1,Work(ipLo),n,z,n,info)
    rr = DDot_(n,z,1,r,1)
    call dscal_(n,rr/beta,p,1)
    call daxpy_(n,One,z,1,p,1)
    k = k+1
  end do
  call Free_Work(ipLo)
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
