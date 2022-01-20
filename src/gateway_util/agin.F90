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

subroutine AGIN()

use AMatrix, only: DFAC, KOSUU, lp1, lp12, NYU, RCA
use Constants, only: Zero, One, Two, Three, Four, Five, Seven, Nine, Ten
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: I, ICOL, IJ, IROW, IVAL, J, L

DFAC(1) = One
DFAC(2) = One
do I=3,lp12
  DFAC(I) = DFAC(I-2)*real(I-1,kind=wp)
end do

! RCA(i,j) = c(k) (la,0;lb,0) * sqrt((2*la+1)(2*lb+1))
!   j identifies la,lb:  ss,ps,pp,ds,dp,dd,... 1,2,3,4,5,6,...
!   i identifies k: i=1 for the lowest k whose c(k) is non-zero,
!                   i=2 for the next   k  "     "    "   ", etc.

RCA(:,:) = Zero
RCA(1,1) = One
RCA(1,2) = One/Three
RCA(1,3) = One/Three
RCA(2,3) = Two/15.0_wp
RCA(1,4) = One/Five
RCA(1,5) = Two/15.0_wp
RCA(2,5) = Three/35.0_wp
RCA(1,6) = One/Five
RCA(2,6) = Two/35.0_wp
RCA(3,6) = Two/35.0_wp
RCA(1,7) = One/Seven
RCA(1,8) = Three/35.0_wp
RCA(2,8) = Four/63.0_wp
RCA(1,9) = Three/35.0_wp
RCA(2,9) = Four/105.0_wp
RCA(3,9) = Ten/231.0_wp
RCA(1,10) = One/Seven
RCA(2,10) = Four/105.0_wp
RCA(3,10) = Two/77.0_wp
RCA(4,10) = 100.0_wp/3003.0_wp
RCA(1,11) = One/Nine
RCA(1,12) = Four/63.0_wp
RCA(2,12) = Ten/198.0_wp
RCA(1,13) = Two/35.0_wp
RCA(2,13) = 20.0_wp/693.0_wp
RCA(3,13) = Ten/286.0_wp
RCA(1,14) = Four/63.0_wp
RCA(2,14) = Two/77.0_wp
RCA(3,14) = 20.0_wp/1001.0_wp
RCA(4,14) = 70.0_wp/2574.0_wp
RCA(1,15) = One/Nine
RCA(2,15) = 20.0_wp/693.0_wp
RCA(3,15) = 162.0_wp/9009.0_wp
RCA(4,15) = 20.0_wp/1287.0_wp
RCA(5,15) = 490.0_wp/21879.0_wp

IJ = 0
do I=1,lp1
  do J=1,I
    IJ = IJ+1
    KOSUU(IJ) = J
  end do
end do

ICOL = 0
do L=1,lp1
  do I=1,L
    ICOL = ICOL+1
    IVAL = L-I-2
    do IROW=1,I
      IVAL = IVAL+2
      NYU(IROW,ICOL) = IVAL
    end do
  end do
end do

return

end subroutine AGIN
