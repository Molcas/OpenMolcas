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
! Copyright (C) 2008, Francesco Aquilante                              *
!***********************************************************************

subroutine ChoMP2_FNO(irc,D_ab,D_ii,EOcc,EVir,Sorted,DelOrig)
!
! F. Aquilante, Geneva May 2008  (snick in Pedersen's code)

use ChoMP2, only: nBatch
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(inout) :: D_ab(*), D_ii(*)
real(kind=wp), intent(in) :: EOcc(*), EVir(*)
logical(kind=iwp), intent(in) :: Sorted, DelOrig
integer(kind=iwp) :: lWrk
real(kind=wp), allocatable :: Wrk(:)
character(len=*), parameter :: SecNam = 'ChoMP2_FNO'

irc = 0

call mma_maxDBLE(lWrk)
call mma_allocate(Wrk,lWrk,Label='Wrk')

if (Sorted) then
  call ChoMP2_fno_Srt(irc,DelOrig,D_ab,D_ii,EOcc,EVir,Wrk,lWrk)
  if (irc /= 0) write(u6,*) SecNam,': ChoMP2_fno_Srt returned ',irc
else if (nBatch == 1) then
  call ChoMP2_fno_Fll(irc,DelOrig,D_ab,D_ii,EOcc,EVir,Wrk,lWrk)
  if (irc /= 0) write(u6,*) SecNam,': ChoMP2_fno_Fll returned ',irc
else
  call ChoMP2_fno_Org(irc,DelOrig,D_ab,D_ii,EOcc,EVir,Wrk,lWrk)
  if (irc /= 0) write(u6,*) SecNam,': ChoMP2_fno_Org returned ',irc
end if

call mma_deallocate(Wrk)

end subroutine ChoMP2_FNO
