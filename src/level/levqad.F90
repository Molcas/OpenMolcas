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

!***********************************************************************
! Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
!***********************************************************************
subroutine LEVQAD(Y1,Y2,Y3,H,RT,ANS1,ANS2)
!** Subroutine "LEVQAD" fits quadratic  Y = A + B*X + C*X**2  through
!  function values  Y1, Y2, Y3  at equally spaced points separated by
!  distance H, where  Y1 < 0  and (Y2,Y3 >= 0), locates the function
!  zero (at RT, relative to  X1 < X2 = 0) between points X1 & X2, and
!  evaluates the integral from RT to R3 of   1/sqrt(Y)  , called
!  ANS1, and the integral (same range) of  sqrt(Y) , which is ANS2
!** Alternately, if Y1 & Y3 both  < 0  and only the middle point
!  Y2 >= 0 ,   fit the points to:  Y = A - B*(X-X0)**2 , locate the
!  turning points between which  Y(X) > 0  and evaluate these integrals
!  on this interval.  *************************************************
!----------------------------------------------------------------------

use Constants, only: Zero, One, Two, Three, Four, Half, Pi
use Definitions, only: wp, u6

implicit none
real(kind=wp), intent(in) :: Y1, Y2, Y3, H
real(kind=wp), intent(out) :: RT, ANS1, ANS2
real(kind=wp) :: A, B, C, CQ, R1, R2, RCQ, RR, SL3, SLT, X0, ZT

if ((Y1 >= 0) .or. (Y2 < 0)) then
  write(u6,602) Y1,Y2
  ANS1 = Zero
  ANS2 = Zero
else if (Y3 < Zero) then
  ! Here treat case when only 'Y2' is non-negative
  RR = (Y2-Y1)/(Y2-Y3)
  X0 = H*(RR-ONe)/((RR+ONe)*Two)
  B = (Y2-Y1)/(H*(Two*X0+H))
  A = Y2+B*X0**2
  ZT = sqrt(A/B)
  RT = X0-ZT
  ANS1 = PI/sqrt(B)
  ANS2 = ANS1*A*Half
else
  ! Here treat case where both 'Y2' & 'Y3' are positive
  if (abs((Y2-Y1)/(Y3-Y2)-One) < 1.0e-10_wp) then
    ! ... special case of true (to 1/10^10) linearity ...
    RT = -H*Y2/(Y2-Y1)
    ANS1 = Two*(H-RT)/sqrt(Y3)
    ANS2 = ANS1*Y3/Three
    return
  end if
  C = (Y3-Two*Y2+Y1)/(Two*H*H)
  B = (Y3-Y2)/H-C*H
  A = Y2
  CQ = B**2-Four*A*C
  RCQ = sqrt(CQ)
  R1 = (-B-RCQ)/(Two*C)
  R2 = R1+RCQ/C
  if ((R2 <= Zero) .and. (R2 >= -H)) RT = R2
  if ((R1 <= Zero) .and. (R1 >= -H)) RT = R1
  SL3 = Two*C*H+B
  SLT = Two*C*RT+B
  if (C < Zero) then
    ANS1 = -(asin(SL3/RCQ)-sign(Half*Pi,SLT))/sqrt(-C)
  else
    ANS1 = log((Two*sqrt(C*Y3)+SL3)/SLT)/sqrt(C)
  end if
  ANS2 = (SL3*sqrt(Y3)-CQ*ANS1*Half)/(Four*C)
  if (RT >= H) write(u6,601) H,R1,R2
end if

return

601 format(' *** CAUTION *** in LEVQAD, turning point not between points 1 & 2.   H =',F9.6,'   R1 =',F9.6,'   R2 =',F9.6)
602 format(' *** ERROR in LEVQAD *** No turning point between 1-st two points as   Y1=',ES10.3,'   Y2=',ES10.3)

end subroutine LEVQAD
