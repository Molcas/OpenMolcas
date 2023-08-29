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

subroutine CHO_ANASIZE(VEC,LVEC,BIN,LBIN,LUPRI)
!
! Purpose: analyse vector (histogram).

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LVEC, LBIN, LUPRI
real(kind=wp), intent(in) :: VEC(LVEC)
real(kind=wp), intent(inout) :: BIN(LBIN)
integer(kind=iwp), parameter :: MBIN = 20
integer(kind=iwp) :: I, IBIN, ICOUNT(MBIN), IJOB, JCOUNT, NBIN, NLOW, NNEG, NZER
real(kind=wp) :: TEST, TOPCT, XNEG
logical(kind=iwp) :: FOUND

! Return if nothing to do.
! ------------------------

if ((LVEC < 1) .or. (LBIN < 1)) return

! Ensure that BIN is in descending order.
! ---------------------------------------

IJOB = -1
call CHO_ORDER(BIN,LBIN,IJOB)

! Test that BIN is positive.
! --------------------------

if (BIN(1) <= Zero) return

! Analysis.
! ---------

NBIN = min(LBIN,MBIN)
ICOUNT(1:NBIN) = 0
NLOW = 0
NZER = 0
NNEG = 0
XNEG = Zero

do I=1,LVEC

  TEST = VEC(I)

  if (TEST < Zero) then
    NNEG = NNEG+1
    XNEG = min(XNEG,TEST)
  else if (TEST == Zero) then
    NZER = NZER+1
  end if

  IBIN = 0
  FOUND = .false.
  do while ((.not. FOUND) .and. (IBIN < NBIN))
    IBIN = IBIN+1
    if (TEST >= BIN(IBIN)) then
      ICOUNT(IBIN) = ICOUNT(IBIN)+1
      FOUND = .true.
    end if
  end do
  if (.not. FOUND) NLOW = NLOW+1

end do

! Print.
! ------

TOPCT = 1.0e2_wp/real(LVEC,kind=wp)

JCOUNT = ICOUNT(1)
write(LUPRI,'(/,1X,A,11X,D11.4,A,I12,1X,F7.2,A,3X,A,F7.2,A)') 'Larger than ',BIN(1),':',ICOUNT(1),real(ICOUNT(1),kind=wp)*TOPCT, &
                                                              '%','Accumulated: ',real(JCOUNT,kind=wp)*TOPCT,'%'
do IBIN=2,NBIN
  JCOUNT = JCOUNT+ICOUNT(IBIN)
  write(LUPRI,'(1X,A,D11.4,A,D11.4,A,I12,1X,F7.2,A,3X,A,F7.2,A)') 'Between ',BIN(IBIN-1),' and ',BIN(IBIN),':',ICOUNT(IBIN), &
                                                                  real(ICOUNT(IBIN),kind=wp)*TOPCT,'%','Accumulated: ', &
                                                                  real(JCOUNT,kind=wp)*TOPCT,'%'
end do
JCOUNT = JCOUNT+NLOW
write(LUPRI,'(1X,A,10X,D11.4,A,I12,1X,F7.2,A,3X,A,F7.2,A)') 'Smaller than ',BIN(NBIN),':',NLOW,real(NLOW,kind=wp)*TOPCT,'%', &
                                                            'Accumulated: ',real(JCOUNT,kind=wp)*TOPCT,'%'

write(LUPRI,'(/,1X,A,I12,1X,F7.2,A)') 'Number of elements exactly 0.0   :',NZER,real(NZER,kind=wp)*TOPCT,'%'
write(LUPRI,'(1X,A,I12,1X,F7.2,A)') 'Number of negative elements      :',NNEG,real(NNEG,kind=wp)*TOPCT,'%'
if (NNEG > 0) write(LUPRI,'(1X,A,D12.4)') ' - numerically largest           :',XNEG

end subroutine CHO_ANASIZE
