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
subroutine PLYINTRP(XI,YI,NPT,RR,C,NCFT,IDER)
!* From the NPT known mesh points (XI,YI), given in order of increasing
!  or decreasing XI(I), select the NCFT points (XJ,YJ) surrounding the
!  given point RR, and by fitting an (NCFT-1)-th degree polynomial through
!  them, interpolate to find the function CC(1) and its first IDER
!  derivatives (CC(I+1),I=1,IDER) evaluated at RR.
!* Adapted by  R.J. Le Roy  from algorithm #416,Comm.A.C.M.;  27/02/1988
!=======================================================================

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NPT, NCFT
real(kind=wp), intent(in) :: XI(NPT), YI(NPT), RR
real(kind=wp), intent(out) :: C(NCFT)
integer(kind=iwp), intent(out) :: IDER
integer(kind=iwp) :: I, I1, I2, IFC, IM, J, K, NH
real(kind=wp) :: XJ(20), XX, YJ(20)

IM = 0
if ((NCFT > 20) .or. (NCFT > NPT)) then
  write(u6,601) NCFT,NCFT,NPT
  !stop
  call ABEND()
end if
NH = NCFT/2
! First locate the known mesh points (XJ,YJ) bracketing RR
I1 = 1
I2 = NCFT
if (NCFT /= NPT) then
  if (XI(NPT) <= XI(1)) then
    do I=1,NPT
      IM = I
      if (XI(I) < RR) exit
    end do
  else
    do I=1,NPT
      IM = I
      if (XI(I) > RR) exit
    end do
  end if
  I1 = IM-NH
  if (I1 <= 0) I1 = 1
  I2 = I1+NCFT-1
  if (I2 > NPT) then
    I2 = NPT
    I1 = I2-NCFT+1
  end if
end if
XJ(1:I2-I1+1) = XI(I1:I2)-RR
YJ(1:I2-I1+1) = YI(I1:I2)
! Now determine polynomial coefficients C(I).
do I=2,NCFT
  I1 = I-1
  K = I1+1
  do J=1,I1
    K = K-1
    YJ(K) = (YJ(K+1)-YJ(K))/(XJ(I)-XJ(K))
  end do
end do
C(1) = YJ(1)
do I=2,NCFT
  XX = XJ(I)
  C(I) = C(I-1)
  if (I /= 2) then
    I1 = I-1
    K = I1+1
    do J=2,I1
      K = K-1
      C(K) = -XX*C(K)+C(K-1)
    end do
  end if
  C(1) = YJ(I)-XX*C(1)
end do
! Finally, convert polynomial coefficients to derivatives at RR.
IFC = 1
if (IDER >= NCFT) IDER = NCFT-1
if (IDER > 1) then
  do I=2,IDER
    J = I+1
    IFC = IFC*I
    C(J) = C(J)*IFC
  end do
  C(J+1:NCFT) = Zero
end if

return

601 format(/' *** Dimensioning ERROR in PLYINTRP :  either   (NCFT=',I2,' > 20)   or   (NCFT=',I2,' > NPT=',I3,')')

end subroutine PLYINTRP
