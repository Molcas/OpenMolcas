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
subroutine DXTYP(NDXTP,ITYP,JTYP,KTYP,LTYP,LEL1,LEL3,REL1,REL3)
! Obtain types of I,J,K,l so
! <L!a+I a+K a L a J!R> is nonvanishing
! only combinations with type(I) >= type(K) and type(L) >= type(J)
! are included

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: NDXTP, ITYP(36), JTYP(36), KTYP(36), LTYP(36), LEL1, LEL3, REL1, REL3
integer(kind=iwp) :: I1, I3, IK1, IK3, IKL1, IKL3, IKLJ1, IKLJ3, ITP, JTP, KTP, LTP

! To get rid of annoying and incorrect compiler warnings
I1 = 0
I3 = 0
IK1 = 0
IK3 = 0
IKL1 = 0
IKL3 = 0
IKLJ1 = 0
IKLJ3 = 0

NDXTP = 0
do ITP=1,3
  if (ITP == 1) then
    I1 = 1
    I3 = 0
  else if (ITP == 2) then
    I1 = 0
    I3 = 0
  else if (ITP == 3) then
    I1 = 0
    I3 = 1
  end if
  do KTP=1,ITP
    if (KTP == 1) then
      IK1 = I1+1
      IK3 = I3
    else if (KTP == 2) then
      IK1 = I1
      IK3 = I3
    else if (KTP == 3) then
      IK1 = I1
      IK3 = I3+1
    end if
    if (LEL1-IK1 < 0) cycle
    if (LEL3-IK3 < 0) cycle
    do LTP=1,3
      if (LTP == 1) then
        IKL1 = IK1-1
        IKL3 = IK3
      else if (LTP == 2) then
        IKL1 = IK1
        IKL3 = IK3
      else if (LTP == 3) then
        IKL1 = IK1
        IKL3 = IK3-1
      end if
      do JTP=1,3
        if (JTP == 1) then
          IKLJ1 = IKL1-1
          IKLJ3 = IKL3
        else if (JTP == 2) then
          IKLJ1 = IKL1
          IKLJ3 = IKL3
        else if (JTP == 3) then
          IKLJ1 = IKL1
          IKLJ3 = IKL3-1
        end if
        if ((IKLJ1+REL1 == LEL1) .and. (IKLJ3+REL3 == LEL3)) then
          NDXTP = NDXTP+1
          ITYP(NDXTP) = ITP
          KTYP(NDXTP) = KTP
          LTYP(NDXTP) = LTP
          JTYP(NDXTP) = JTP
        end if
      end do
    end do
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,'(A,4I4)') ' Double excitations connecting LEL1,LEL3,LEL1,LEL3 ',LEL1,LEL3,REL1,REL3
write(u6,*) '  ITYP KTYP LTYP JTYP'
write(u6,*) '  ===================='
do IDX=1,NDXTP
  write(u6,'(1X,5I5)') ITYP(IDX),KTYP(IDX),LTYP(IDX),JTYP(IDX)
end do
#endif

return

end subroutine DXTYP
