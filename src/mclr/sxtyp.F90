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

!#define _DEBUGPRINT_
subroutine SXTYP(NSXTP,ITP,JTP,LEL1,LEL3,REL1,REL3)
! Types of creation and annihilation  operators so
! <L!a+ a!R> is nonvanishing
!
! L is defined by LEL1,LEL3
! R is defined by REL1,REL3

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(out) :: NSXTP, ITP(3), JTP(3)
integer(kind=iwp), intent(in) :: LEL1, LEL3, REL1, REL3
integer(kind=iwp) :: I1, I123, I3, IJ1, IJ3, J123
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: I
#endif

NSXTP = 0

! To get rid of annoying and incorrect compiler warnings
I1 = 0
I3 = 0
IJ1 = 0
IJ3 = 0

do I123=1,3
  if (I123 == 1) then
    I1 = 1
    I3 = 0
  else if (I123 == 2) then
    I1 = 0
    I3 = 0
  else if (I123 == 3) then
    I1 = 0
    I3 = 1
  end if
  if (LEL1-I1 < 0) cycle
  if (LEL3-I3 < 0) cycle
  do J123=1,3
    if (J123 == 1) then
      IJ1 = I1-1
      IJ3 = I3
    else if (J123 == 2) then
      IJ1 = I1
      IJ3 = I3
    else if (J123 == 3) then
      IJ1 = I1
      IJ3 = I3-1
    end if
    if ((REL1+IJ1 == LEL1) .and. (REL3+IJ3 == LEL3)) then
      NSXTP = NSXTP+1
      ITP(NSXTP) = I123
      JTP(NSXTP) = J123
    end if
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,'(A,4I4)') ' SX  connecting LEL1,LEL3,REL1,REL3 ',LEL1,LEL3,REL1,REL3
write(u6,*) ' Number of connections obtained ',NSXTP
write(u6,*) ' ITYPE JTYPE'
write(u6,*) ' ==========='
do I=1,NSXTP
  write(u6,'(2I5)') ITP(I),JTP(I)
end do
#endif

return

end subroutine SXTYP
