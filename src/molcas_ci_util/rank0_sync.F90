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
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

!> Copy rank 0's values to all processes (Rank0_Bcast).
!>
!> RASSCF keeps the orbitals, the CI vectors and the Davidson subspace
!> replicated on every process, and every process takes the iteration
!> decisions (convergence, new trial vectors) itself. Rank-local linear
!> algebra (diagonalizations, threaded BLAS) need not agree bitwise
!> between processes, and in a degenerate subspace the difference is not
!> round-off but a different eigenvector. The processes then hold
!> different data, combine it in the global sums, and take different
!> branches (wrong results, or a deadlock in the next global sum).
!> Rank0_Bcast makes rank 0's copy the one all processes use. It does
!> nothing in a serial run; all processes must call it at the same point.
!> Rank0_Active tells whether there is anything to do (to skip work that
!> only prepares a broadcast).

module Rank0_Sync

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
use MPI_Wrapper, only: MPI_Bcast, MPI_COMM_WORLD, MPI_REAL8
#ifdef _I8_
use MPI_Wrapper, only: MPI_INTEGER8
#else
use MPI_Wrapper, only: MPI_INTEGER4
#endif
use Definitions, only: MPIInt
#endif
use Definitions, only: wp, iwp

implicit none
private

public :: Rank0_Active, Rank0_Bcast

interface Rank0_Bcast
  module procedure :: Bcast_r0, Bcast_r1, Bcast_r2, Bcast_i0, Bcast_i1
end interface Rank0_Bcast

#ifdef _MOLCAS_MPP_
! the MPI count is a default integer: long arrays go in pieces
integer(kind=iwp), parameter :: MxPiece = 2_iwp**30
#endif

#include "macros.fh"

contains

function Rank0_Active()

  logical(kind=iwp) :: Rank0_Active

#ifdef _MOLCAS_MPP_
  Rank0_Active = Is_Real_Par()
#else
  Rank0_Active = .false.
#endif

end function Rank0_Active

subroutine Bcast_rbuf(A,n)

  integer(kind=iwp), intent(in) :: n
  real(kind=wp), intent(inout) :: A(n)
#ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: i0, m
  integer(kind=MPIInt) :: iErr

  if (.not. Is_Real_Par()) return
  do i0=1,n,MxPiece
    m = min(MxPiece,n-i0+1)
    call MPI_Bcast(A(i0),int(m,kind=MPIInt),MPI_REAL8,0_MPIInt,MPI_COMM_WORLD,iErr)
  end do
#else
  unused_var(A)
#endif

end subroutine Bcast_rbuf

subroutine Bcast_ibuf(iA,n)

  integer(kind=iwp), intent(in) :: n
  integer(kind=iwp), intent(inout) :: iA(n)
#ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: i0, m
  integer(kind=MPIInt) :: iErr

  if (.not. Is_Real_Par()) return
  do i0=1,n,MxPiece
    m = min(MxPiece,n-i0+1)
#   ifdef _I8_
    call MPI_Bcast(iA(i0),int(m,kind=MPIInt),MPI_INTEGER8,0_MPIInt,MPI_COMM_WORLD,iErr)
#   else
    call MPI_Bcast(iA(i0),int(m,kind=MPIInt),MPI_INTEGER4,0_MPIInt,MPI_COMM_WORLD,iErr)
#   endif
  end do
#else
  unused_var(iA)
#endif

end subroutine Bcast_ibuf

subroutine Bcast_r0(x)

  real(kind=wp), intent(inout) :: x
  real(kind=wp) :: Buf(1)

  Buf(1) = x
  call Bcast_rbuf(Buf,1_iwp)
  x = Buf(1)

end subroutine Bcast_r0

subroutine Bcast_r1(A)

  real(kind=wp), contiguous, intent(inout) :: A(:)

  call Bcast_rbuf(A,size(A,kind=iwp))

end subroutine Bcast_r1

subroutine Bcast_r2(A)

  real(kind=wp), contiguous, target, intent(inout) :: A(:,:)
  real(kind=wp), pointer, contiguous :: A1(:)

  A1(1:size(A)) => A
  call Bcast_rbuf(A1,size(A,kind=iwp))
  nullify(A1)

end subroutine Bcast_r2

subroutine Bcast_i0(i)

  integer(kind=iwp), intent(inout) :: i
  integer(kind=iwp) :: Buf(1)

  Buf(1) = i
  call Bcast_ibuf(Buf,1_iwp)
  i = Buf(1)

end subroutine Bcast_i0

subroutine Bcast_i1(iA)

  integer(kind=iwp), contiguous, intent(inout) :: iA(:)

  call Bcast_ibuf(iA,size(iA,kind=iwp))

end subroutine Bcast_i1

end module Rank0_Sync
