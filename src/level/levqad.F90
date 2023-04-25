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

real*8 A, ANS1, ANS2, B, C, CQ, H, HPI, R1, R2, RCQ, RR, RT, SL3, SLT, X0, Y1, Y2, Y3, ZT
data HPI/1.570796326794896d0/

if ((Y1 >= 0) .or. (Y2 < 0)) go to 99
if (Y3 < 0.d0) go to 50
! Here treat case where both 'Y2' & 'Y3' are positive
if (dabs((Y2-Y1)/(Y3-Y2)-1.d0) < 1.d-10) then
  ! ... special case of true (to 1/10^10) linearity ...
  RT = -H*Y2/(Y2-Y1)
  ANS1 = 2.d0*(H-RT)/dsqrt(Y3)
  ANS2 = ANS1*Y3/3.d0
  return
end if
C = (Y3-2.d0*Y2+Y1)/(2.d0*H*H)
B = (Y3-Y2)/H-C*H
A = Y2
CQ = B**2-4.d0*A*C
RCQ = dsqrt(CQ)
R1 = (-B-RCQ)/(2.d0*C)
R2 = R1+RCQ/C
if ((R2 <= 0.d0) .and. (R2 >= -H)) RT = R2
if ((R1 <= 0.d0) .and. (R1 >= -H)) RT = R1
SL3 = 2.d0*C*H+B
SLT = 2.d0*C*RT+B
if (C < 0.d0) go to 10
ANS1 = dlog((2.d0*dsqrt(C*Y3)+SL3)/SLT)/dsqrt(C)
go to 20
10 continue
ANS1 = -(dasin(SL3/RCQ)-dsign(HPI,SLT))/dsqrt(-C)
20 continue
ANS2 = (SL3*dsqrt(Y3)-CQ*ANS1/2.d0)/(4.d0*C)
if (RT >= H) write(6,601) H,R1,R2
601 format(' *** CAUTION *** in LEVQAD, turning point not between points 1 & 2.   H =',F9.6,'   R1 =',F9.6,'   R2 =',F9.6)

return

! Here treat case when only 'Y2' is non-negative
50 continue
RR = (Y2-Y1)/(Y2-Y3)
X0 = H*(RR-1.d0)/((RR+1.d0)*2.d0)
B = (Y2-Y1)/(H*(2.d0*X0+H))
A = Y2+B*X0**2
ZT = dsqrt(A/B)
RT = X0-ZT
ANS1 = 2.d0*HPI/dsqrt(B)
ANS2 = ANS1*A*0.5d0

return

99 continue
write(6,602) Y1,Y2
602 format(' *** ERROR in LEVQAD *** No turning point between 1-st two points as   Y1=',D10.3,'   Y2=',D10.3)

ANS1 = 0.d0
ANS2 = 0.d0

return

end subroutine LEVQAD
