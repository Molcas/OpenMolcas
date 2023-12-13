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
!               Francesco Aquilante                                    *
!***********************************************************************
!  Cho_X_CheckDiag
!
!> @brief
!>   Check diagonal
!> @author Thomas Bondo Pedersen
!> @modified_by F. Aquilante (add ::OneCenter_ChkDiag)
!> @modified_by T.B. Pedersen (If ``(Cho_1Center)``: \p Err only contains 1-center errors on exit)
!>
!> @details
!> This routine reads and analyzes (histogram and statistics)
!> the exact integral diagonal, computes and analyzes the
!> diagonal from Cholesky vectors,
!> and the difference between the two (exact minus Cholesky).
!>
!> The statistics printed are: minimum value, maximum
!> value, mean value, mean absolute value, variance (wrt mean
!> value), and standard deviation (wrt mean value).
!>
!> On exit:
!>
!> - \p Err(1) = min error
!> - \p Err(2) = max error
!> - \p Err(3) = average error
!> - \p Err(4) = RMS error
!>
!> Return code is ``0`` if successful execution. If \p irc is non-zero,
!> the contents or \p Err are ill-defined.
!> Results will only be printed to output if \c iPrint is ``-5`` or
!> greater.
!>
!> @param[out] irc Return code
!> @param[out] Err min, max, average, and RMS error
!***********************************************************************

subroutine Cho_X_CheckDiag(irc,Err)

use Cholesky, only: Cho_1Center, IPRINT, nnBstRT, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(out) :: Err(4)
integer(kind=iwp) :: i
real(kind=wp), allocatable :: Bin(:), CD(:), Stat(:), XD(:)
integer(kind=iwp), parameter :: iPrThr = -5
character(len=*), parameter :: SecNam = 'Cho_X_CheckDiag'
real(kind=wp), external :: ddot_

! Set return code.
! ----------------

irc = 0
if (nnBstRT(1) < 1) then
  Err(:) = Zero
  return
end if

! Allocations.
! ------------

call mma_allocate(XD,nnBstRT(1),Label='XD')
call mma_allocate(CD,nnBstRT(1),Label='CD')
call mma_allocate(Bin,16,Label='Bin')
call mma_allocate(Stat,7,Label='Stat')

! Set bins for histograms.
! ------------------------

Bin(1) = One
do i=1,size(Bin)-1
  Bin(1+i) = Bin(i)*1.0e-1_wp
end do

! Read exact diagonal.
! --------------------

call Cho_IODiag(XD,2)

! Print histogram of exact diagonal and get statistics.
! -----------------------------------------------------

if (iPrint >= iPrThr) then
  call Cho_Head('Analysis of Exact Integral Diagonal','=',80,u6)
  call Cho_AnaSize(XD,size(XD),Bin,size(Bin),u6)
  call Statistics(XD,size(XD),Stat,1,2,3,4,5,6,7)
  call Cho_PrtSt(XD,size(XD),Stat)
end if

! Calculate Cholesky diagonal.
! ----------------------------

call Cho_X_CalcChoDiag(irc,CD)
if (irc /= 0) then
  write(u6,*) SecNam,': Cho_X_CalcChoDiag returned ',irc
else

  ! Print histogram of Cholesky diagonal and get statistics.
  ! --------------------------------------------------------

  if (iPrint >= iPrThr) then
    call Cho_Head('Analysis of Cholesky Integral Diagonal','=',80,u6)
    call Cho_AnaSize(CD,size(CD),Bin,size(Bin),u6)
    call Statistics(CD,size(CD),Stat,1,2,3,4,5,6,7)
    call Cho_PrtSt(CD,size(CD),Stat)
  end if

  ! Subtract Cholesky diagonal from exact diagonal.
  ! -----------------------------------------------

  XD(:) = XD(:)-CD(:)

  ! Print histogram of difference array and get statistics.
  ! -------------------------------------------------------

  if (iPrint >= iPrThr) then
    call Cho_Head('Analysis of Difference (Exact-Cholesky)','=',80,u6)
    call Cho_AnaSize(XD,size(XD),Bin,size(Bin),u6)
  end if
  call Statistics(XD,size(XD),Stat,1,2,3,4,5,6,7)
  if (iPrint >= iPrThr) call Cho_PrtSt(XD,size(XD),Stat)

  ! Set Err array.
  ! --------------

  Err(1) = Stat(3)
  Err(2) = Stat(4)
  Err(3) = Stat(1)
  Err(4) = sqrt(dDot_(nnBstRT(1),XD,1,XD,1)/real(nnBstRT(1),kind=wp))

  if (iPrint >= iPrThr) then
    write(u6,'(/,1X,A,ES15.6)') 'Minimum error   : ',Err(1)
    write(u6,'(1X,A,ES15.6)') 'Maximum error   : ',Err(2)
    write(u6,'(1X,A,ES15.6)') 'Average error   : ',Err(3)
    write(u6,'(1X,A,ES15.6)') 'RMS error       : ',Err(4)
  end if

  ! Error analysis for the 1-center diagonals only.
  ! If this is a one-center calculation, use statistics from 1-center
  ! diagonals only as elements of Err array.
  ! -----------------------------------------------------------------

  if (nSym == 1) then
    call OneCenter_ChkDiag(XD,size(XD),Stat,iPrint >= iPrThr)
    if (Cho_1Center) then
      Err(1) = Stat(3)
      Err(2) = Stat(4)
      Err(3) = Stat(1)
      Err(4) = sqrt(dDot_(nnBstRT(1),XD,1,XD,1)/real(nnBstRT(1),kind=wp))
    end if
  end if

end if

! Deallocations.
! --------------

call mma_deallocate(Stat)
call mma_deallocate(Bin)
call mma_deallocate(CD)
call mma_deallocate(XD)

end subroutine Cho_X_CheckDiag
