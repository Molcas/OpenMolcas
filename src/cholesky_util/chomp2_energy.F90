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
! Copyright (C) 2004, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_Energy(irc,EMP2,EOcc,EVir,Sorted,DelOrig)
!
! Thomas Bondo Pedersen, Dec. 2004.
!
! Purpose: compute MP2 energy correction from MO Cholesky vectors,
!          constructing (ai|bj) integrals on the fly. Flag Sorted
!          refers to whether or not the MO vectors have been sorted
!          into the sizes of the batches over occupied orbitals.

use ChoMP2, only: nBatch
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(out) :: EMP2
real(kind=wp), intent(in) :: EOcc(*), EVir(*)
logical(kind=iwp), intent(in) :: Sorted, DelOrig
integer(kind=iwp) :: lWrk
real(kind=wp), allocatable :: Wrk(:)
character(len=*), parameter :: SecNam = 'ChoMP2_Energy'

irc = 0

call mma_maxDBLE(lWrk)
call mma_allocate(Wrk,lWrk,Label='Wrk')

if (Sorted) then
  call ChoMP2_Energy_Srt(irc,DelOrig,EMP2,EOcc,EVir,Wrk,lWrk)
  if (irc /= 0) write(u6,*) SecNam,': ChoMP2_Energy_Srt returned ',irc
else if (nBatch == 1) then
  call ChoMP2_Energy_Fll(irc,DelOrig,EMP2,EOcc,EVir,Wrk,lWrk)
  if (irc /= 0) write(u6,*) SecNam,': ChoMP2_Energy_Fll returned ',irc
else
  call ChoMP2_Energy_Org(irc,DelOrig,EMP2,EOcc,EVir,Wrk,lWrk)
  if (irc /= 0) write(u6,*) SecNam,': ChoMP2_Energy_Org returned ',irc
end if

call mma_deallocate(Wrk)

end subroutine ChoMP2_Energy
