!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine TRPMT3(XIN,NROW,NCOL,XOUT)
! XOUT(I,J) = XIN(J,I)
!
! With a few considerations for large scale cases with cache minimization

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NROW, NCOL
real(kind=wp), intent(in) :: XIN(NROW,NCOL)
real(kind=wp), intent(out) :: XOUT(NCOL,NROW)
integer(kind=iwp) :: ICBLK, ICEND, ICOFF, IRBLK, IREND, IROFF, IROW, LCBLK, LRBLK, NCBLK, NRBLK
integer(kind=iwp), parameter :: IMET = 2

if (IMET == 1) then
  ! Straightforward, no blocking
  call TRNSPS(NROW,NCOL,XIN,XOUT)
else if (IMET == 2) then
  ! Simple blocking of matrix
  LRBLK = 40
  LCBLK = 40
  NRBLK = NROW/LRBLK
  NCBLK = NCOL/LCBLK
  if (LRBLK*NRBLK /= NROW) NRBLK = NRBLK+1
  if (LCBLK*NCBLK /= NCOL) NCBLK = NCBLK+1

  do IRBLK=1,NRBLK
    if (IRBLK == 1) then
      IROFF = 1
    else
      IROFF = IROFF+LRBLK
    end if
    IREND = min(NROW,IROFF+LRBLK-1)
    do ICBLK=1,NCBLK
      if (ICBLK == 1) then
        ICOFF = 1
      else
        ICOFF = ICOFF+LCBLK
      end if
      ICEND = min(NCOL,ICOFF+LCBLK-1)

      do IROW=IROFF,IREND
        XOUT(ICOFF:ICEND,IROW) = XIN(IROW,ICOFF:ICEND)
      end do

    end do
  end do
end if

end subroutine TRPMT3
