
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Francesco Aquilante                                    *
!               2002, Thomas Bondo Pedersen                            *
!***********************************************************************
!  CHO_GET_ORD_bj
!
!> @brief
!>   Computes the decomposition pattern of the 2nd-rank orbital energy denominators according to Eq. (6) of \cite Koc2000-JCP-113-508
!> @author F. Aquilante
!> @modified_by Thomas Bondo Pedersen, December 2012: code restructured, fix of minor bug
!>
!> @details
!> Computes the decomposition pattern of the 2nd-rank orbital energy denominators
!> according to Eq. (6) of \cite Koc2000-JCP-113-508
!>
!> The decomposition is considered converged when
!> either of the two criteria (\p MaxNVec or \p thr) is fulfilled.
!> Therefore these two arguments should have appropriate
!> values depending on what the user expects:
!>
!> - \p MaxNVec = \p nOV if the user specifies a meaningful \p thr
!> - \p thr = ``0.0``    if the user requires a given number of vectors (or percentage of \p nOV)
!>
!> @param[in]  nOV     Number of (occ,vir) pairs matching a given compound symmetry
!> @param[in]  MaxNVec Max number of Cholesky vectors
!> @param[in]  thr     Threshold for the Cholesky decomposition
!> @param[in]  W       Array (\p nOV) of the orbital energy differences \p W(bj) = ``e(b)-e(j)``
!> @param[out] ID_bj   Index array (\p MaxNVec) of the decomposition pattern
!> @param[out] nVec    Number of resulting Cholesky vectors
!> @param[out] Dmax    Max value of the remainder after the decomposition
!***********************************************************************

subroutine CHO_GET_ORD_bj(nOV,MaxNVec,thr,W,ID_bj,nVec,Dmax)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nOV, MaxNVec
real(kind=wp), intent(in) :: thr, W(*)
integer(kind=iwp), intent(_OUT_) :: ID_bj(*)
integer(kind=iwp), intent(out) :: NVec
real(kind=wp), intent(out) :: Dmax
integer(kind=iwp) :: ip, Jm
real(kind=wp), allocatable :: Diag(:)

! Initialize
! ----------
nVec = 0
if (nOV < 1) then
  Dmax = -9.987654321_wp ! dummy value
  return
end if

! Allocate diagonal
! -----------------
call mma_allocate(Diag,NOV,Label='Diag')

! Compute diagonals
! -----------------
do ip=1,nOV
  if (W(ip) > Zero) then
    Diag(ip) = Half/W(ip)
  else ! tbp: perhaps we should stop it here (matrix not PSD)
    Diag(ip) = Zero
  end if
end do

! Find ID (Jm) of max diagonal
! ----------------------------
Jm = 1
do ip=2,NOV
  if (Diag(ip) > Diag(Jm)) Jm = ip
end do

! Find CD pattern using the diagonal update:
!   Diag(p)[k] = Diag(p)[k-1] * ((W(p) - W(J[k-1]))/(W(p) + W(J[k-1])))^2
! -----------------------------------------------------------------
do while ((nVec < MaxNVec) .and. (Diag(Jm) > thr))
  ! update vector counter
  nVec = nVec+1
  ! save ID of vector = ID of current max diagonal
  ID_bj(nVec) = Jm
  ! update diagonals
  Diag(1:nOV) = Diag(1:nOV)*((W(1:nOV)-W(Jm))/(W(1:nOV)+W(Jm)))**2
  ! find ID (Jm) of max updated diagonal
  Jm = 1
  do ip=2,nOV
    if (Diag(ip) > Diag(Jm)) Jm = ip
  end do
end do

! Set max diagonal (Dmax) to return to caller
! -------------------------------------------
Dmax = Diag(Jm)

! Deallocate diagonal
! -------------------
call mma_deallocate(Diag)

return

end subroutine CHO_GET_ORD_bj
