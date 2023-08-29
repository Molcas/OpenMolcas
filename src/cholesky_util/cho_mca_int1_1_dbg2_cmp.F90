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

subroutine CHO_MCA_INT1_1_DBG2_CMP(XINT1,XINT2,NI,NJ,ERRMIN,IMN,JMN,ERRMAX,IMX,JMX,ITST,IERR,THR,PRTERR,LUPRI)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NI, NJ, LUPRI
real(kind=wp), intent(in) :: XINT1(NI,NJ), XINT2(NJ,NI), THR
real(kind=wp), intent(out) :: ERRMIN, ERRMAX
integer(kind=iwp), intent(out) :: IMN, JMN, IMX, JMX
integer(kind=iwp), intent(inout) :: ITST, IERR
logical(kind=iwp), intent(in) :: PRTERR
integer(kind=iwp) :: I, J, JERR
real(kind=wp) :: DIFF

if ((NI < 1) .or. (NJ < 1)) then
  ERRMAX = Zero
  ERRMIN = Zero
  IMN = 0
  JMN = 0
  IMX = 0
  JMX = 0
  return
end if

ERRMAX = XINT1(1,1)-XINT2(1,1)
ERRMIN = XINT1(1,1)-XINT2(1,1)
IMN = 1
JMN = 1
IMX = 1
JMX = 1

JERR = 0
do J=1,NJ
  do I=1,NI
    ITST = ITST+1
    DIFF = XINT1(I,J)-XINT2(J,I)
    if (abs(DIFF) > THR) then
      JERR = JERR+1
      if (PRTERR) write(LUPRI,*) '      Error: ',I,J,DIFF
    end if
    if (DIFF < ERRMIN) then
      ERRMIN = DIFF
      IMN = I
      JMN = J
    end if
    if (DIFF > ERRMAX) then
      ERRMAX = DIFF
      IMX = I
      JMX = J
    end if
  end do
end do
IERR = IERR+JERR

if ((JERR /= 0) .and. (NI == NJ)) then
  JERR = 0
  write(LUPRI,*) '         Checking for identity...'
  do J=1,NJ
    do I=1,NI
      DIFF = XINT1(I,J)-XINT2(I,J)
      if (abs(DIFF) > 1.0e-14_wp) JERR = JERR+1
    end do
  end do
  if (JERR /= 0) then
    write(LUPRI,*) '      ...not identical!!'
  else
    write(LUPRI,*) '      ...identical!!'
  end if
end if

end subroutine CHO_MCA_INT1_1_DBG2_CMP
