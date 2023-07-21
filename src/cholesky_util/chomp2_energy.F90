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

use stdalloc

implicit none
real*8 EMP2
real*8 EOcc(*), EVir(*)
integer irc
logical Sorted, DelOrig
#include "chomp2.fh"
#include "chomp2_cfg.fh"
character(len=6), parameter :: ThisNm = 'Energy'
character(len=13), parameter :: SecNam = 'ChoMP2_Energy'
integer lWrk
real*8, allocatable :: Wrk(:)

irc = 0

call mma_maxDBLE(lWrk)
call mma_allocate(Wrk,lWrk,Label='Wrk')

if (Sorted) then
  call ChoMP2_Energy_Srt(irc,DelOrig,EMP2,EOcc,EVir,Wrk,lWrk)
  if (irc /= 0) then
    write(6,*) SecNam,': ChoMP2_Energy_Srt returned ',irc
    Go To 1 ! exit
  end if
else
  if (nBatch == 1) then
    call ChoMP2_Energy_Fll(irc,DelOrig,EMP2,EOcc,EVir,Wrk,lWrk)
    if (irc /= 0) then
      write(6,*) SecNam,': ChoMP2_Energy_Fll returned ',irc
      Go To 1 ! exit
    end if
  else
    call ChoMP2_Energy_Org(irc,DelOrig,EMP2,EOcc,EVir,Wrk,lWrk)
    if (irc /= 0) then
      write(6,*) SecNam,': ChoMP2_Energy_Org returned ',irc
      Go To 1 ! exit
    end if
  end if
end if

1 call mma_deallocate(Wrk)

end subroutine ChoMP2_Energy
