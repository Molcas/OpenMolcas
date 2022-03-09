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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine THETSET(ICASE,THE,NII)
!PAM06 This routine is called from SDCI if this is an MCPF calculation

use cpf_global, only: IPRINT, IRC, IREF0, LN, LWSP
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ICASE(*), NII
real(kind=wp), intent(out) :: THE(NII,NII)
integer(kind=iwp) :: I, II, II1, IINT, IIOR, IJ, IK, IL, IOCR(100), IP, IQ, JJ, JOJ, NI, NJ
integer(kind=iwp), external :: ICUNP

IIOR = 0
II1 = (IREF0-1)*LN
do I=1,LN
  JOJ = ICUNP(ICASE,II1+I)
  IIOR = IIOR+1
  IOCR(IIOR) = JOJ
end do
if (IPRINT > 5) write(u6,888) IREF0,(IOCR(I),I=1,LN)

IINT = IRC(4)
THE(1:IINT,1:IINT) = Zero
do IP=1,IINT
  II = 0
  IJ = 0
  do I=1,LN
    JJ = (IP-1)*LN+I
    if ((ICUNP(ICASE,JJ) == IOCR(I)) .or. (ICUNP(ICASE,JJ) == 3)) cycle
    if (LWSP .and. (ICUNP(ICASE,JJ)*IOCR(I) == 2)) cycle
    if (II == 0) II = I
    IJ = I
  end do
  !PAM06 BUG: What if we come down here with II == 0 still??
  ! the IOCR will be accessed below first element. Provisional fix:
  if (II /= 0) then
    NI = IOCR(II)
    if (NI > 1) NI = NI-1
    NJ = IOCR(IJ)
    if (NJ > 1) NJ = NJ-1
    do IQ=1,IINT
      IK = 0
      IL = 0
      do I=1,LN
        JJ = (IQ-1)*LN+I
        if ((ICUNP(ICASE,JJ) == IOCR(I)) .or. (ICUNP(ICASE,JJ) == 3)) cycle
        if (LWSP .and. (ICUNP(ICASE,JJ)*IOCR(I) == 2)) cycle
        if (IK == 0) IK = I
        IL = I
      end do
      if ((II == IK) .and. (IJ == IL)) THE(IQ,IP) = One
    end do
  end if
end do

return

888 format(5X,'IREF0=',I3/5X,'IOCR=',10I5)

end subroutine THETSET
