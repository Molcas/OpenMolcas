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
subroutine LEVEL_SPLINE(X,Y,N,IOPT,C,N4,IER)
!** Subroutine for generating cubic spline coefficients
!  C(J), (J=1,N4=4*N) through the N points X(I), Y(I).
!** C(I+M*N), M=0-3  are the coefficients of order  0-3  of cubic
!  polynomial expanded about X(I) so as to describe the interval:
!             -  X(I) to X(I+1)  , if  X(I)  in increasing order
!             -  X(I-1) to X(I)  , if  X(I)  in decreasing order.
!** IOPT indicates boundary conditions used in creating the  spline .
!*  If (IOPT=0)  second derivatives = zero at both ends of range.
!*  If (IOPT=1)  1st derivative at first point X(1) fixed at C(1),
!                and 2nd derivative at X(N) = zero.
!*  If (IOPT=2)  1st derivative at last point X(N) fixed at C(2),
!                and 2nd derivative at X(1) = zero.
!*  If (IOPT=3)  constrain first derivatives at end points to have
!                (read in) values  C(1)  at  X(1)  &  C(2)  at  X(N)
!** IER is the error flag.  IER=0  on return if routine successful.
!-----------------------------------------------------------------------

integer I, II, IER, IOH, IOL, IOPT, J, J1, J2, J3, NER, N, N4, JMP
real*8 A, H, R, DY2, DYA, DYB, XB, XC, YA, YB, X(N), Y(N), C(N4)

J1 = 0
II = 0
! Make sure II is "referenced":
write(6,*) II
JMP = 1
NER = 1000
if (N <= 1) go to 250
! Initialization
XC = X(1)
YB = Y(1)
H = 0.d0
A = 0.d0
R = 0.d0
DYB = 0.d0
NER = 2000

! IOL=0 - given derivative at firstpoint
! IOH=0 - given derivative at last point

IOL = IOPT-1
IOH = IOPT-2
if (IOH == 1) then
  IOL = 0
  IOH = 0
end if
DY2 = C(2)

! Form the system of linear equations
! and eliminate subsequentially

J = 1
do I=1,N
  J2 = N+I
  J3 = J2+N
  A = H*(2.d0-A)
  DYA = DYB+H*R
  if (I >= N) then

    ! set derivative dy2 at last point

    DYB = DY2
    H = 0.d0
    if (IOH == 0) goto 200
    DYB = DYA
    goto 220
  end if
  J = J+JMP
  XB = XC
  XC = X(J)
  H = XC-XB

  ! II=0 - increasing abscissae
  ! II=1 - decreasing abscissae

  II = 0
  if (H < 0) II = 1
  if (H == 0) go to 250
  YA = YB
  YB = Y(J)
  DYB = (YB-YA)/H
  if (I <= 1) then
    J1 = II
    if (IOL /= 0) go to 220
    DYA = C(1)
  end if
  200 continue
  if (J1 /= II) go to 250
  A = 1.d0/(H+H+A)
  220 continue
  R = A*(DYB-DYA)
  C(J3) = R
  A = H*A
  C(J2) = A
  C(I) = DYB
end do

! back substitution of the system of linear equations
! and computation of the other coefficients

A = 1.d0
J1 = J3+N+II-II*N
I = N
do IOL=1,N
  XB = X(J)
  H = XC-XB
  XC = XB
  A = A+H
  YB = R
  R = C(J3)-R*C(J2)
  YA = R+R
  C(J3) = YA+R
  C(J2) = C(I)-H*(YA+YB)
  C(J1) = (YB-R)/A
  C(I) = Y(J)
  A = 0.d0
  J = J-JMP
  I = I-1
  J2 = J2-1
  J3 = J3-1
  J1 = J3+N+II
end do
IER = 0

return

250 continue
IER = NER

return

end subroutine LEVEL_SPLINE
