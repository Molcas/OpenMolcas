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
! Copyright (C) 1994, Niclas Forsberg                                  *
!***********************************************************************

subroutine BondStr(R,i1,i2,j,S,NumOfAt,NumInt)
!  Purpose:
!    Calculate contribution to the S-matrix due to bond stretching.
!
!  Input:
!    R        : Array of Real*8 real -  contains the
!               cartesian coordinates of the bond.
!    i1,i2    : Integer - the number of the atom.
!    j        : Integer - the number of the internal coordinate.
!
!  Output:
!    S        : Real*8 three dimensional array - the
!               contributions to S for the parameters specified
!               in the input.
!
!  Uses:
!    Linalg
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1994.

implicit real*8(a-h,o-z)
#include "Constants_mula.fh"
real*8 S(3,NumOfAt,NumInt)
real*8 R(3)

! Contributions to S.
SR = sqrt(R(1)**2+R(2)**2+R(3)**2)
do k=1,3
  S(k,i1,j) = -R(k)/SR
  S(k,i2,j) = -S(k,i1,j)
  call NaNChk(S(k,i1,j),k,i1,j)
  call NaNChk(S(k,i2,j),k,i1,j)
end do

end subroutine BondStr
!####
subroutine NaNChk(X,k,i,j)

use Definitions, only: u6

implicit real*8(a-h,o-z)
real*8 X
character*16 str1, str2

write(str2,'(G16.8)') X
call Normalize(str2,str1)
if (str1(1:3) == 'NAN') then
  write(u6,*) ' CalcS subroutine produced Not-a-Number!'
  write(u6,*) ' Internal coordinate nr j=',j
  write(u6,*) ' Atom nr.               i=',i
  write(u6,*) ' Component              k=',k
  write(u6,*) ' S(k,i,j)=',X
end if

return

end subroutine NaNChk
!####
subroutine AngBend(R1,R2,i1,i2,i3,j,S,NumOfAt,NumInt)
!  Purpose:
!    Calculate contribution to the S-matrix due to valence
!    angle bending.
!
!  Input:
!    R1,R2    : Array of Real*8 real - contains the
!               cartesian coordinates of the bond.
!    i1,i2,i3 : Integer - the number of the atom.
!    j        : Integer - the number of the internal coordinate.
!
!  Output:
!    S        : Real*8 three dimensional array - the
!               contributions to S for the parameters specified
!               in the input.
!
!  Uses:
!    Linalg
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1994.

use Constants, only: One

implicit real*8(a-h,o-z)
#include "Constants_mula.fh"
real*8 S(3,NumOfAt,NumInt)
real*8 R1(3), R2(3)
real*8 NR1(3), NR2(3)
real*8 CosTheta, SinTheta, F1, F2, F3
real*8 SR1, SR2

! Unit vectors.
SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
do k=1,3
  NR1(k) = R1(k)/SR1
  NR2(k) = R2(k)/SR2
end do

! Angle Theta.
Theta = acos(dDot_(3,NR1,1,NR2,1))
CosTheta = cos(Theta)
SinTheta = sin(Theta)

! Contributions to S.
F1 = One/(SR1*SinTheta)
S(1,i1,j) = (CosTheta*NR1(1)-NR2(1))*F1
S(2,i1,j) = (CosTheta*NR1(2)-NR2(2))*F1
S(3,i1,j) = (CosTheta*NR1(3)-NR2(3))*F1

F2 = SR1-SR2*CosTheta
F3 = SR2-SR1*CosTheta
F4 = F1/SR2
S(1,i2,j) = (F2*NR1(1)+F3*NR2(1))*F4
S(2,i2,j) = (F2*NR1(2)+F3*NR2(2))*F4
S(3,i2,j) = (F2*NR1(3)+F3*NR2(3))*F4

F5 = One/(SR2*SinTheta)
S(1,i3,j) = (CosTheta*NR2(1)-NR1(1))*F5
S(2,i3,j) = (CosTheta*NR2(2)-NR1(2))*F5
S(3,i3,j) = (CosTheta*NR2(3)-NR1(3))*F5

call NaNChk(S(1,i1,j),1,i1,j)
call NaNChk(S(2,i1,j),2,i1,j)
call NaNChk(S(3,i1,j),3,i1,j)
call NaNChk(S(1,i2,j),1,i2,j)
call NaNChk(S(2,i2,j),2,i2,j)
call NaNChk(S(3,i2,j),3,i2,j)
call NaNChk(S(1,i3,j),1,i3,j)
call NaNChk(S(2,i3,j),2,i3,j)
call NaNChk(S(3,i3,j),3,i3,j)

end subroutine AngBend
!####
subroutine LinBend(R1,R2,i1,i2,i3,j,S,NumOfAt,NumInt)
!  Purpose:
!    Calculate contribution to the S-matrix due to valence
!    angle bending for a linear molecule.
!
!  Input:
!    R1,R2    : Array of Real*8 real -  contains the
!               cartesian coordinates of the bond.
!    i1,i2,i3 : Integer - the number of the atom.
!    j        : Integer - the number of the internal coordinate.
!
!  Output:
!    S        : Real*8 three dimensional array - the
!               contributions to S for the parameters specified
!               in the input.
!
!  Uses:
!    Linalg
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1994.

use Constants, only: Zero, One

implicit real*8(a-h,o-z)
#include "Constants_mula.fh"
real*8 S(3,NumOfAt,NumInt)
real*8 R1(3), R2(3)

! Length of vectors.
SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)

! Contributions to S.
F1 = One/SR1
F2 = One/SR2

S(1,i1,j) = F1
S(2,i1,j) = Zero
S(3,i1,j) = Zero

S(1,i2,j) = -F1-F2
S(2,i2,j) = Zero
S(3,i2,j) = Zero

S(1,i3,j) = F2
S(2,i3,j) = Zero
S(3,i3,j) = Zero

S(1,i1,j+1) = Zero
S(2,i1,j+1) = F1
S(3,i1,j+1) = Zero

S(1,i2,j+1) = Zero
S(2,i2,j+1) = -F1-F2
S(3,i2,j+1) = Zero

S(1,i3,j+1) = Zero
S(2,i3,j+1) = F2
S(3,i3,j+1) = Zero

call NaNChk(S(1,i1,j),1,i1,j)
call NaNChk(S(2,i1,j),2,i1,j)
call NaNChk(S(3,i1,j),3,i1,j)
call NaNChk(S(1,i2,j),1,i2,j)
call NaNChk(S(2,i2,j),2,i2,j)
call NaNChk(S(3,i2,j),3,i2,j)
call NaNChk(S(1,i3,j),1,i3,j)
call NaNChk(S(2,i3,j),2,i3,j)
call NaNChk(S(3,i3,j),3,i3,j)

end subroutine LinBend
!####
subroutine Torsion(R1,R2,R3,i1,i2,i3,i4,j,S,NumOfAt,NumInt)
!  Purpose:
!    Calculate contribution to the S-matrix due to torsion.
!
!  Input:
!    R1,R2,R3 : Array of Real*8 real -  contains the
!               cartesian coordinates of the bond.
!    i1,i2,
!    i3,i4    : Integer - the number of the atom.
!    j        : Integer - the number of the internal coordinate.
!
!  Output:
!    S        : Real*8 three dimensional array - the
!               contributions to S for the parameters specified
!               in the input.
!
!  Uses:
!    Linalg
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1994.

use Constants, only: One

implicit real*8(a-h,o-z)
#include "Constants_mula.fh"
real*8 S(3,NumOfAt,NumInt)
real*8 R1(3), R2(3), R3(3)
real*8 NR1(3), NR2(3), NR3(3)

! Unit vectors.
SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
SR3 = sqrt(R3(1)**2+R3(2)**2+R3(3)**2)
do k=1,3
  NR1(k) = R1(k)/SR1
  NR2(k) = R2(k)/SR2
  NR3(k) = R3(k)/SR3
end do

! Angles Theta1 and Theta2.
Theta1 = acos(-dDot_(3,NR1,1,NR2,1))
Theta2 = acos(-dDot_(3,NR2,1,NR3,1))
CosTheta1 = cos(Theta1)
SinTheta1 = sin(Theta1)
CosTheta2 = cos(Theta2)
SinTheta2 = sin(Theta2)

! Contributions to S.
F1 = One/(SR1*(SinTheta1)**2)
F2 = NR1(2)*NR2(3)-NR1(3)*NR2(2)
F3 = NR1(3)*NR2(1)-NR1(1)*NR2(3)
F4 = NR1(1)*NR2(2)-NR1(2)*NR2(1)
S(1,i1,j) = -F2*F1
S(2,i1,j) = -F3*F1
S(3,i1,j) = -F4*F1

F5 = (SR2-SR1*CosTheta1)/(SR1*SR2*(SinTheta1)**2)
F6 = CosTheta2/(SR2*(SinTheta2)**2)
F7 = NR3(2)*NR2(3)-NR3(3)*NR2(2)
F8 = NR3(3)*NR2(1)-NR3(1)*NR2(3)
F9 = NR3(1)*NR2(2)-NR3(2)*NR2(1)
S(1,i2,j) = F5*F2+F6*F7
S(2,i2,j) = F5*F3+F6*F8
S(3,i2,j) = F5*F4+F6*F9

F10 = (SR2-SR3*CosTheta2)/(SR2*SR3*SinTheta2**2)
F11 = CosTheta1/(SR2*(SinTheta1)**2)
S(1,i3,j) = F10*F7+F11*F2
S(2,i3,j) = F10*F8+F11*F3
S(3,i3,j) = F10*F9+F11*F4

F12 = One/(SR3*(SinTheta2)**2)
S(1,i4,j) = -F7*F12
S(2,i4,j) = -F8*F12
S(3,i4,j) = -F9*F12

call NaNChk(S(1,i1,j),1,i1,j)
call NaNChk(S(2,i1,j),2,i1,j)
call NaNChk(S(3,i1,j),3,i1,j)
call NaNChk(S(1,i2,j),1,i2,j)
call NaNChk(S(2,i2,j),2,i2,j)
call NaNChk(S(3,i2,j),3,i2,j)
call NaNChk(S(1,i3,j),1,i3,j)
call NaNChk(S(2,i3,j),2,i3,j)
call NaNChk(S(3,i3,j),3,i3,j)
call NaNChk(S(1,i1,j),1,i4,j)
call NaNChk(S(2,i1,j),2,i4,j)
call NaNChk(S(3,i1,j),3,i4,j)

end subroutine Torsion
!####
subroutine OutOfPl(R1,R2,R3,i1,i2,i3,i4,j,S,NumOfAt,NumInt)
!
!  Purpose:
!    Calculate contribution to the S-matrix due to the angle between
!    a bond and a plane defined by two bonds.
!
!  Input:
!    R1,R2,R3 : Array of Real*8 real -  contains the
!               cartesian coordinates of the bond.
!    i1,i2,
!    i3,i4    : Integer - the number of the atom.
!    j        : Integer - the number of the internal coordinate.
!
!  Output:
!    S        : Real*8 three dimensional array - the
!               contributions to S for the parameters specified
!               in the input.
!  Uses:
!    Linalg
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1994.

use Constants, only: One

implicit real*8(a-h,o-z)
#include "Constants_mula.fh"
real*8 S(3,NumOfAt,NumInt)
real*8 R1(3), R2(3), R3(3)
real*8 NR1(3), NR2(3), NR3(3), F0(3)

! Unit vectors.
SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
SR3 = sqrt(R3(1)**2+R3(2)**2+R3(3)**2)
do k=1,3
  NR1(k) = R1(k)/SR1
  NR2(k) = R2(k)/SR2
  NR3(k) = R3(k)/SR3
end do

! Calculation of the angle Theta.
Dot = Ddot_(3,NR2,1,NR3,1)
Theta = acos(Dot)
CosTheta = cos(Theta)
SinTheta = sin(Theta)

! Calculation of the out of plane angle Phi.
F1 = NR2(2)*NR3(3)-NR2(3)*NR3(2)
F2 = NR2(3)*NR3(1)-NR2(1)*NR3(3)
F3 = NR2(1)*NR3(2)-NR2(2)*NR3(1)
F0(1) = F1/SinTheta
F0(2) = F2/SinTheta
F0(3) = F3/SinTheta
Phi = asin(Ddot_(3,F0,1,NR1,1))
CosPhi = cos(Phi)
TanPhi = tan(Phi)

! Contributions to S.
F4 = One/(SR1*CosPhi*SinTheta)
F5 = TanPhi/SR1
S(1,i1,j) = F1*F4-NR1(1)*F5
S(2,i1,j) = F2*F4-NR1(2)*F5
S(3,i1,j) = F3*F4-NR1(3)*F5

F6 = NR3(2)*NR1(3)-NR3(3)*NR1(2)
F7 = NR3(3)*NR1(1)-NR3(1)*NR1(3)
F8 = NR3(1)*NR1(2)-NR3(2)*NR1(1)
F9 = One/(SR2*CosPhi*SinTheta)
F10 = TanPhi/(SR2*(SinTheta)**2)
S(1,i2,j) = F6*F9-NR2(1)*F10+CosTheta*F10*NR3(1)
S(2,i2,j) = F7*F9-NR2(2)*F10+CosTheta*F10*NR3(2)
S(3,i2,j) = F8*F9-NR2(3)*F10+CosTheta*F10*NR3(3)

F11 = NR1(2)*NR2(3)-NR1(3)*NR2(2)
F12 = NR1(3)*NR2(1)-NR1(1)*NR2(3)
F13 = NR1(1)*NR2(2)-NR1(2)*NR2(1)
F14 = One/(SR3*CosPhi*SinTheta)
F15 = TanPhi/(SR3*(SinTheta)**2)
S(1,i3,j) = F11*F14-NR3(1)*F15+CosTheta*F15*NR2(1)
S(2,i3,j) = F12*F14-NR3(2)*F15+CosTheta*F15*NR2(2)
S(3,i3,j) = F13*F14-NR3(3)*F15+CosTheta*F15*NR2(3)

S(1,i4,j) = -S(1,i1,j)-S(1,i2,j)-S(1,i3,j)
S(2,i4,j) = -S(2,i1,j)-S(2,i2,j)-S(2,i3,j)
S(3,i4,j) = -S(3,i1,j)-S(3,i2,j)-S(3,i3,j)

call NaNChk(S(1,i1,j),1,i1,j)
call NaNChk(S(2,i1,j),2,i1,j)
call NaNChk(S(3,i1,j),3,i1,j)
call NaNChk(S(1,i2,j),1,i2,j)
call NaNChk(S(2,i2,j),2,i2,j)
call NaNChk(S(3,i2,j),3,i2,j)
call NaNChk(S(1,i3,j),1,i3,j)
call NaNChk(S(2,i3,j),2,i3,j)
call NaNChk(S(3,i3,j),3,i3,j)
call NaNChk(S(1,i1,j),1,i4,j)
call NaNChk(S(2,i1,j),2,i4,j)
call NaNChk(S(3,i1,j),3,i4,j)

end subroutine OutOfPl
!####
subroutine CalcS(AtCoord,InterVec,S,NumInt,NumOfAt)
!
!  Purpose:
!    Calculate the contribution to the S-matrix for each of the
!    internal coordinates.
!
!  Input:
!    AtCoord  : Two dimensional Real*8 array - contains
!               the cartesian coordinates of the atoms.
!    InterVec : Integer array - contains the atoms that are used
!               in the calculations of each internal coordinate.
!    S        : Real*8 array filled with zeros.
!
!  Output:
!    S        : Real*8 array - contains the contributions
!               of each of the internal coordinates.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1994.

implicit real*8(a-h,o-z)
#include "Constants_mula.fh"
real*8 AtCoord(3,NumOfAt)
integer InterVec(*)
real*8 S(3,NumOfAt,NumInt)
real*8 R(3), R1(3), R2(3), R3(3)

! Calculate the contributions to the S-matrix for each internal
! coordinate specified in the vector InterVec.
k = 1
IntType = InterVec(k)
nInt = 1
do while (nInt <= NumInt)
  if (IntType == 1) then
    ! Bond Stretching.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    R(:) = AtCoord(:,i2)-AtCoord(:,i1)
    call BondStr(R,i1,i2,nInt,S,NumOfAt,NumInt)
    k = k+3
  end if
  if (IntType == 2) then
    ! Valence Angle Bending.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i2)
    R2(:) = AtCoord(:,i3)-AtCoord(:,i2)
    call AngBend(R1,R2,i1,i2,i3,nInt,S,NumOfAt,NumInt)
    k = k+4
  end if
  if (IntType == 3) then
    ! Linear Valence Angle.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i2)
    R2(:) = AtCoord(:,i3)-AtCoord(:,i2)
    call LinBend(R1,R2,i1,i2,i3,nInt-1,S,NumOfAt,NumInt)
    k = k+4
  end if
  if (IntType == 4) then
    ! Torsion.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    i4 = InterVec(k+4)
    R1(:) = AtCoord(:,i2)-AtCoord(:,i1)
    R2(:) = AtCoord(:,i3)-AtCoord(:,i2)
    R3(:) = AtCoord(:,i4)-AtCoord(:,i3)
    call Torsion(R1,R2,R3,i1,i2,i3,i4,nInt,S,NumOfAt,NumInt)
    k = k+5
  end if
  if (IntType == 5) then
    ! Out of Plane Angle Bending.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    i4 = InterVec(k+4)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i4)
    R2(:) = AtCoord(:,i2)-AtCoord(:,i4)
    R3(:) = AtCoord(:,i3)-AtCoord(:,i4)
    call OutOfPl(R1,R2,R3,i1,i2,i3,i4,nInt,S,NumOfAt,NumInt)
    k = k+5
  end if
  IntType = InterVec(k)
  nInt = nInt+1
end do

end subroutine CalcS
!####
subroutine CalcG(G,Mass,S,NumInt,NumOfAt)
!  Purpose:
!    Calculate G-matrix.
!
!  Input:
!    Mass     : Real*8 array - the mass of the atoms.
!    S        : Real*8 three dimensional array - the
!               contributions from all the internal coordinates.
!
!  Output:
!    G        : Real*8 two dimensional array.
!
!  Uses:
!    Constants
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1994.

use Constants, only: Zero, One

implicit real*8(a-h,o-z)
#include "Constants_mula.fh"
real*8 G(NumInt,NumInt)
real*8 Mass(NumOfAt)
real*8 S(3,NumOfAt,NumInt)

do i=1,NumInt
  do j=1,NumInt
    GSum = Zero
    do k=1,NumOfAt
      GSum = GSum+(One/(uToAu*Mass(k)))*(S(1,k,i)*S(1,k,j)+S(2,k,i)*S(2,k,j)+S(3,k,i)*S(3,k,j))
    end do
    G(i,j) = GSum
  end do
end do

end subroutine CalcG
