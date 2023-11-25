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

subroutine CHO_ANADIA(DIAG,BIN1,STEP,NUMBIN,FULL)
!
! Purpose: analyze diagonal (histogram).

use Cholesky, only: DIAMNZ, IABMNZ, LuPri, nnBstRT, nnZTot, ThrCom
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Diag(*), BIN1, STEP
integer(kind=iwp), intent(in) :: NUMBIN
logical(kind=iwp), intent(in) :: FULL
integer(kind=iwp), parameter :: MAXBIN = 50, NUMSTA = 7
integer(kind=iwp) :: IAB, IBIN, MBIN, NBIN, NCONV
real(kind=wp) :: BIN(MAXBIN), BINLOC, STAT(NUMSTA), STPLOC
logical(kind=iwp) :: FOUND

! Print header.
! -------------

call CHO_HEAD('Histogram of Diagonal Elements','=',80,LUPRI)

! Set up size bins for analysis of diagonal.
! ------------------------------------------

if (NUMBIN < 1) then
  MBIN = min(10,MAXBIN)
  BINLOC = 1.0e2_wp
  STPLOC = 1.0e-2_wp
else
  MBIN = min(NUMBIN,MAXBIN)
  BINLOC = BIN1
  STPLOC = STEP
end if

BIN(1) = BINLOC
do IBIN=2,MBIN
  BIN(IBIN) = BIN(IBIN-1)*STPLOC
end do

! Set smallest BIN according to full.
! -----------------------------------

if (FULL) then
  NBIN = MBIN
else
  NBIN = MBIN
  IBIN = MBIN
  FOUND = .false.
  do while ((IBIN > 1) .and. (.not. FOUND))
    IBIN = IBIN-1
    if (THRCOM >= BIN(IBIN)) then
      NBIN = IBIN+1
    else
      FOUND = .true.
    end if
  end do
end if

! Histogram.
! ----------

call CHO_ANASIZE(DIAG,NNBSTRT(1),BIN,NBIN,LUPRI)

! Count converged.
! ----------------

NCONV = 0
do IAB=1,NNBSTRT(1)
  if (DIAG(IAB) <= THRCOM) NCONV = NCONV+1
end do
write(LUPRI,'(/,1X,A,I10,/,1X,A,I10)') 'Converged  : ',NCONV,'Unconverged: ',NNBSTRT(1)-NCONV

! Print total number of negative zeroed diagonal as well as the most
! negative one.
! ------------------------------------------------------------------

write(LUPRI,'(/,1X,A,5X,I10)') 'Total number of zeroed negative diagonals: ',NNZTOT
if (NNZTOT > 0) then
  if (IABMNZ < 1) then
    write(LUPRI,'(1X,A)') 'WARNING: most negative zeroed diagonal has not been stored!'
  else
    write(LUPRI,'(1X,A,ES15.6)') '- most negative zeroed diagonal          : ',DIAMNZ
  end if
end if

! Print statistics.
! -----------------

call STATISTICS(DIAG,NNBSTRT(1),STAT,1,2,3,4,5,6,7)
write(LUPRI,'(/,1X,A,ES15.6)') 'Minimum diagonal: ',STAT(3)
write(LUPRI,'(1X,A,ES15.6)') 'Maximum diagonal: ',STAT(4)
write(LUPRI,'(1X,A,ES15.6)') 'Mean value      : ',STAT(1)
write(LUPRI,'(1X,A,ES15.6)') 'Mean abs. value : ',STAT(2)
write(LUPRI,'(1X,A,ES15.6)') 'Biased variance : ',STAT(6)
write(LUPRI,'(1X,A,ES15.6,A)') 'Standard dev.   : ',STAT(7),' (unbiased variance)'

end subroutine CHO_ANADIA
