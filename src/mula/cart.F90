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
! Copyright (C) 1997, Niclas Forsberg                                  *
!               1997,2000, Per Ake Malmqvist                           *
!***********************************************************************

subroutine Cart_to_Int1(InterVec,AtCoord,xvec,BMatrix,NumOfAt,NumInt)
!  Purpose:
!    Calculate internal coordinates, and their gradients.
!
!  Input:
!    InterVec : Integer array - contains the atoms that are used
!               in the calculations of each internal coordinate.
!    AtCoord  : Two dimensional Real*8 array - contains
!               the cartesian coordinates of the atoms.
!
!  Output:
!    xvec     : Real*8 array - internal coordinates.
!   BMatrix   : Real*8 array - Their gradients.
!
!  Written by:
!    Niclas Forsberg, Per Ake Malmqvist,
!    Dept. of Theoretical Chemistry, Lund University, 1997.
!
!  Uses:
!    Linalg, Constants

implicit none
#include "Constants_mula.fh"
integer inttype, numint, numofat
real*8 AtCoord(3,NumOfAt)
integer InterVec(*)
real*8 xvec(NumInt)
real*8 BMatrix(3,NumOfAt,NumInt)
real*8 R(3), R1(3), R2(3), R3(3)
real*8 NR1(3), NR2(3), NR3(3)
real*8 Normal123(3), Normal234(3)
real*8 G1(3), G2(3), G3(3)
integer i1, i2, k, nInt, i3, i4
integer iv
real*8 SR1, SR2, SR3, SR, R2R3, Area
real*8 cosv, sinv, cosv0, sinv0, cosdv, sindv
real*8 cos123, sin123, cot123, cos234, sin234, cot234

k = 1
IntType = InterVec(k)
nInt = 1
do while (nInt <= NumInt)
  if (IntType == 1) then
    ! Bond.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    do iv=1,3
      R(iv) = (AtCoord(iv,i1)-AtCoord(iv,i2))
    end do
    SR = sqrt(R(1)**2+R(2)**2+R(3)**2)
    xvec(nInt) = SR
    do iv=1,3
      BMatrix(iv,i1,nInt) = R(iv)/SR
      BMatrix(iv,i2,nInt) = -R(iv)/SR
    end do
    !write(6,*) ' Bond. i1,i2=',i1,i2
    !write(6,'(1x,a,3F16.8)') ' R=',R
    !write(6,'(1x,a,3F16.8)') 'SR=',SR
    !write(6,'(1x,a,3F16.8)') 'BMatrix(:,i1,nInt)=',BMatrix(:,i1,nInt)
    !write(6,'(1x,a,3F16.8)') 'BMatrix(:,i2,nInt)=',BMatrix(:,i2,nInt)
    k = k+3
    nInt = nInt+1
  end if
  if (IntType == 2) then
    ! Valence Angle.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    do iv=1,3
      R1(iv) = (AtCoord(iv,i1)-AtCoord(iv,i2))
      R2(iv) = (AtCoord(iv,i2)-AtCoord(iv,i3))
    end do
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    cosv = -(R1(1)*R2(1)+R1(2)*R2(2)+R1(3)*R2(3))/(SR1*SR2)
    if (cosv > 1.0d0) cosv = 1.0d0
    if (cosv < -1.0d0) cosv = -1.0d0
    sinv = sqrt(1.0d0-cosv**2)
    xvec(nInt) = atan2(sinv,cosv)
    do iv=1,3
      BMatrix(iv,i1,nInt) = (SR1*R2(iv)+SR2*cosv*R1(iv))/(SR1**2*SR2*sinv)
      BMatrix(iv,i3,nInt) = -(SR2*R1(iv)+SR1*cosv*R2(iv))/(SR2**2*SR1*sinv)
      BMatrix(iv,i2,nInt) = -(BMatrix(iv,i1,nInt)+BMatrix(iv,i3,nInt))
    end do
    !write(6,*) ' Angle. i1,i2,i3=',i1,i2,i3
    !write(6,'(1x,a,3F16.8)') ' R1=',R1
    !write(6,'(1x,a,3F16.8)') ' R2=',R2
    !write(6,'(1x,a,3F16.8)') 'BMatrix(:,i1,nInt)=',BMatrix(:,i1,nInt)
    !write(6,'(1x,a,3F16.8)') 'BMatrix(:,i2,nInt)=',BMatrix(:,i2,nInt)
    !write(6,'(1x,a,3F16.8)') 'BMatrix(:,i3,nInt)=',BMatrix(:,i3,nInt)
    k = k+4
    nInt = nInt+1
  end if
  if (IntType == 3) then
    ! Linear Valence Angle.
    !vv?? stop
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    do iv=1,3
      R1(iv) = (AtCoord(iv,i1)-AtCoord(iv,i2))
      R2(iv) = (AtCoord(iv,i2)-AtCoord(iv,i3))
    end do
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    cosv = -(R1(1)*R2(1)+R1(2)*R2(2)+R1(3)*R2(3))/(SR1*SR2)
    if (cosv > 1.0d0) cosv = 1.0d0
    if (cosv < -1.0d0) cosv = -1.0d0
    sinv = sqrt(1.0d0-cosv**2)
    xvec(nInt) = atan2(sinv,cosv)
    BMatrix(1,i1,nInt) = 1.0d0/SR1
    BMatrix(1,i3,nInt) = 1.0d0/SR2
    BMatrix(1,i2,nInt) = -(BMatrix(1,i1,nInt)+BMatrix(1,i3,nInt))
    BMatrix(2,i1,nInt+1) = 1.0d0/SR1
    BMatrix(2,i3,nInt+1) = 1.0d0/SR2
    BMatrix(2,i2,nInt+1) = -(BMatrix(1,i1,nInt)+BMatrix(1,i3,nInt))
    nInt = nInt+2
    k = k+4
  end if
  if (IntType == 4) then
    ! Torsion.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    i4 = InterVec(k+4)
    do iv=1,3
      R1(iv) = (AtCoord(iv,i1)-AtCoord(iv,i2))
      R2(iv) = (AtCoord(iv,i2)-AtCoord(iv,i3))
      R3(iv) = (AtCoord(iv,i3)-AtCoord(iv,i4))
    end do
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    SR3 = sqrt(R3(1)**2+R3(2)**2+R3(3)**2)
    do iv=1,3
      NR1(iv) = R1(iv)/SR1
      NR2(iv) = R2(iv)/SR2
      NR3(iv) = R3(iv)/SR3
    end do
    cos123 = -(NR1(1)*NR2(1)+NR1(2)*NR2(2)+NR1(3)*NR2(3))
    if (cos123 > 1.0d0) cos123 = 1.0d0
    if (cos123 < -1.0d0) cos123 = -1.0d0
    sin123 = sqrt(1.0d0-cos123**2)
    cot123 = cos123/sin123
    Normal123(1) = (NR1(2)*NR2(3)-NR1(3)*NR2(2))/sin123
    Normal123(2) = (NR1(3)*NR2(1)-NR1(1)*NR2(3))/sin123
    Normal123(3) = (NR1(1)*NR2(2)-NR1(2)*NR2(1))/sin123
    cos234 = -(NR2(1)*NR3(1)+NR2(2)*NR3(2)+NR2(3)*NR3(3))
    if (cos234 > 1.0d0) cos234 = 1.0d0
    if (cos234 < -1.0d0) cos234 = -1.0d0
    sin234 = sqrt(1.0d0-cos234**2)
    cot234 = cos234/sin234
    Normal234(1) = (NR2(2)*NR3(3)-NR2(3)*NR3(2))/sin234
    Normal234(2) = (NR2(3)*NR3(1)-NR2(1)*NR3(3))/sin234
    Normal234(3) = (NR2(1)*NR3(2)-NR2(2)*NR3(1))/sin234
    sinv = NR2(1)*(Normal123(2)*Normal234(3)-Normal123(3)*Normal234(2))+ &
           NR2(2)*(Normal123(3)*Normal234(1)-Normal123(1)*Normal234(3))+ &
           NR2(3)*(Normal123(1)*Normal234(2)-Normal123(2)*Normal234(1))
    if (sinv > 1.0d0) sinv = 1.0d0
    if (sinv < -1.0d0) sinv = -1.0d0
    cosv = Normal123(1)*Normal234(1)+Normal123(2)*Normal234(2)+Normal123(3)*Normal234(3)
    if (cosv > 1.0d0) cosv = 1.0d0
    if (cosv < -1.0d0) cosv = -1.0d0
    cosv0 = cos(xvec(nInt))
    sinv0 = sin(xvec(nInt))
    sindv = sinv*cosv0-cosv*sinv0
    cosdv = cosv*cosv0+sinv*sinv0
    xvec(nInt) = xvec(nInt)+atan2(sindv,cosdv)
    do iv=1,3
      G1(iv) = Normal123(iv)/(SR1*sin123)
      G3(iv) = Normal234(iv)/(SR3*sin234)
      G2(iv) = (cot123*Normal123(iv)+cot234*Normal234(iv))/SR2
    end do
    do iv=1,3
      BMatrix(iv,i1,nInt) = G1(iv)
      BMatrix(iv,i2,nInt) = -G1(iv)+G2(iv)
      BMatrix(iv,i3,nInt) = -G2(iv)+G3(iv)
      BMatrix(iv,i4,nInt) = -G3(iv)
    end do
    !write(6,*) ' Torsion. i1,i2,i3,i4=',i1,i2,i3,i4
    !write(6,'(1x,a,3F16.8)') ' R1=',R1
    !write(6,'(1x,a,3F16.8)') ' R2=',R2
    !write(6,'(1x,a,3F16.8)') ' R3=',R3
    !write(6,'(1x,a,3F16.8)') 'BMatrix(:,i1,nInt)=',BMatrix(:,i1,nInt)
    !write(6,'(1x,a,3F16.8)') 'BMatrix(:,i2,nInt)=',BMatrix(:,i2,nInt)
    !write(6,'(1x,a,3F16.8)') 'BMatrix(:,i3,nInt)=',BMatrix(:,i3,nInt)
    !write(6,'(1x,a,3F16.8)') 'BMatrix(:,i4,nInt)=',BMatrix(:,i4,nInt)
    k = k+5
    nInt = nInt+1
  end if
  if (IntType == 5) then
    ! Out of Plane Angle.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    i4 = InterVec(k+4)
    do iv=1,3
      R1(iv) = (AtCoord(iv,i1)-AtCoord(iv,i2))
      R2(iv) = (AtCoord(iv,i2)-AtCoord(iv,i3))
      R3(iv) = (AtCoord(iv,i3)-AtCoord(iv,i4))
    end do
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    SR3 = sqrt(R3(1)**2+R3(2)**2+R3(3)**2)
    R2R3 = (R2(1)*R3(1)+R2(2)*R3(2)+R2(3)*R3(3))
    cos234 = -R2R3/(SR2*SR3)
    if (cos234 > 1.0d0) cos234 = 1.0d0
    if (cos234 < -1.0d0) cos234 = -1.0d0
    sin234 = sqrt(1.0d0-cos234**2)
    Area = SR2*SR3*sin234
    Normal234(1) = (R2(2)*R3(3)-R2(3)*R3(2))/Area
    Normal234(2) = (R2(3)*R3(1)-R2(1)*R3(3))/Area
    Normal234(3) = (R2(1)*R3(2)-R2(2)*R3(1))/Area
    sinv = (R1(1)*Normal234(1)+R1(2)*Normal234(2)+R1(3)*Normal234(3))/SR1
    if (sinv > 1.0d0) sinv = 1.0d0
    if (sinv < -1.0d0) sinv = -1.0d0
    cosv = sqrt(1.0d0-sinv**2)
    xvec(nInt) = atan2(sinv,cosv)
    G1(1) = (R2(2)*R3(3)-R3(2)*R2(3)-(Area*sinv/SR1)*R1(1))/(SR1*Area*cosv)
    G1(2) = (R2(3)*R3(1)-R3(3)*R2(1)-(Area*sinv/SR1)*R1(2))/(SR1*Area*cosv)
    G1(3) = (R2(1)*R3(2)-R3(1)*R2(2)-(Area*sinv/SR1)*R1(3))/(SR1*Area*cosv)
    G2(1) = (R3(2)*R1(3)-R1(2)*R3(3)-(SR1*sinv/Area)*(SR3*R2(1)-R2R3*R3(1)))/(SR1*Area*cosv)
    G2(2) = (R3(3)*R1(1)-R1(3)*R3(1)-(SR1*sinv/Area)*(SR3*R2(2)-R2R3*R3(2)))/(SR1*Area*cosv)
    G2(3) = (R3(1)*R1(2)-R1(1)*R3(2)-(SR1*sinv/Area)*(SR3*R2(3)-R2R3*R3(3)))/(SR1*Area*cosv)
    G3(1) = (R1(2)*R2(3)-R2(2)*R1(3)-(SR1*sinv/Area)*(SR2*R3(1)-R2R3*R2(1)))/(SR1*Area*cosv)
    G3(2) = (R1(3)*R2(1)-R2(3)*R1(1)-(SR1*sinv/Area)*(SR2*R3(2)-R2R3*R2(2)))/(SR1*Area*cosv)
    G3(3) = (R1(1)*R2(2)-R2(1)*R1(2)-(SR1*sinv/Area)*(SR2*R3(3)-R2R3*R2(3)))/(SR1*Area*cosv)
    do iv=1,3
      BMatrix(iv,i1,nInt) = G1(iv)
      BMatrix(iv,i2,nInt) = -G1(iv)+G2(iv)
      BMatrix(iv,i3,nInt) = -G2(iv)+G3(iv)
      BMatrix(iv,i4,nInt) = -G3(iv)
    end do
    k = k+5
    nInt = nInt+1
  end if
  IntType = InterVec(k)
end do

end subroutine Cart_to_Int1
!####
subroutine Cart_to_Int0(InterVec,AtCoord,xvec,NumOfAt,NumInt)
!
!  Purpose:
!    Calculate internal coordinates only (not their gradients).
!
!  Input:
!    InterVec : Integer array - contains the atoms that are used
!               in the calculations of each internal coordinate.
!    AtCoord  : Two dimensional Real*8 array - contains
!               the cartesian coordinates of the atoms.
!
!  Output:
!    xvec     : Real*8 array - internal coordinates.
!
!  Written by:
!    Niclas Forsberg, Per Ake Malmqvist,
!    Dept. of Theoretical Chemistry, Lund University, 1997.
!
!  Uses:
!    Linalg, Constants

implicit none
#include "Constants_mula.fh"
integer inttype, numint, numofat
real*8 AtCoord(3,NumOfAt)
integer InterVec(*)
real*8 xvec(NumInt)
real*8 R(3), R1(3), R2(3), R3(3)
real*8 NR1(3), NR2(3), NR3(3)
real*8 Normal123(3), Normal234(3)
!real*8 G1(3), G2(3), G3(3)
integer i1, i2, k, nInt, i3, i4, iv
real*8 SR1, SR2, SR3, SR, R2R3
real*8 cosv, sinv, cosv0, sinv0, cosdv, sindv
real*8 cos123, sin123, cos234, sin234

k = 1
IntType = InterVec(k)
nInt = 1
do while (nInt <= NumInt)
  if (IntType == 1) then
    ! Bond.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    do iv=1,3
      R(iv) = (AtCoord(iv,i1)-AtCoord(iv,i2))
    end do
    SR = sqrt(R(1)**2+R(2)**2+R(3)**2)
    xvec(nInt) = SR
    k = k+3
    nInt = nInt+1
  end if
  if (IntType == 2) then
    ! Valence Angle.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    do iv=1,3
      R1(iv) = (AtCoord(iv,i1)-AtCoord(iv,i2))
      R2(iv) = (AtCoord(iv,i2)-AtCoord(iv,i3))
    end do
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    cosv = -(R1(1)*R2(1)+R1(2)*R2(2)+R1(3)*R2(3))/(SR1*SR2)
    if (cosv > 1.0d0) cosv = 1.0d0
    if (cosv < -1.0d0) cosv = -1.0d0
    sinv = sqrt(1.0d0-cosv**2)
    xvec(nInt) = atan2(sinv,cosv)
    k = k+4
    nInt = nInt+1
  end if
  if (IntType == 3) then
    ! Linear Valence Angle.
    !vv? stop
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    do iv=1,3
      R1(iv) = (AtCoord(iv,i1)-AtCoord(iv,i2))
      R2(iv) = (AtCoord(iv,i3)-AtCoord(iv,i2))
    end do
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    do iv=1,3
      NR1(iv) = R1(iv)/SR1
      NR2(iv) = R2(iv)/SR2
    end do
    cosv = NR1(1)*NR2(1)+NR1(2)*NR2(2)+NR1(3)*NR2(3)
    if (cosv > 1.0d0) cosv = 1.0d0
    if (cosv < -1.0d0) cosv = -1.0d0
    xvec(nInt) = acos(cosv)
    k = k+4
    nInt = nInt+2
  end if
  if (IntType == 4) then
    ! Torsion.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    i4 = InterVec(k+4)
    do iv=1,3
      R1(iv) = (AtCoord(iv,i1)-AtCoord(iv,i2))
      R2(iv) = (AtCoord(iv,i2)-AtCoord(iv,i3))
      R3(iv) = (AtCoord(iv,i3)-AtCoord(iv,i4))
    end do
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    SR3 = sqrt(R3(1)**2+R3(2)**2+R3(3)**2)
    do iv=1,3
      NR1(iv) = R1(iv)/SR1
      NR2(iv) = R2(iv)/SR2
      NR3(iv) = R3(iv)/SR3
    end do
    cos123 = -(NR1(1)*NR2(1)+NR1(2)*NR2(2)+NR1(3)*NR2(3))
    if (cos123 > 1.0d0) cos123 = 1.0d0
    if (cos123 < -1.0d0) cos123 = -1.0d0
    sin123 = sqrt(1.0d0-cos123**2)
    Normal123(1) = (NR1(2)*NR2(3)-NR1(3)*NR2(2))/sin123
    Normal123(2) = (NR1(3)*NR2(1)-NR1(1)*NR2(3))/sin123
    Normal123(3) = (NR1(1)*NR2(2)-NR1(2)*NR2(1))/sin123
    cos234 = -(NR2(1)*NR3(1)+NR2(2)*NR3(2)+NR2(3)*NR3(3))
    if (cos234 > 1.0d0) cos234 = 1.0d0
    if (cos234 < -1.0d0) cos234 = -1.0d0
    sin234 = sqrt(1.0d0-cos234**2)
    Normal234(1) = (NR2(2)*NR3(3)-NR2(3)*NR3(2))/sin234
    Normal234(2) = (NR2(3)*NR3(1)-NR2(1)*NR3(3))/sin234
    Normal234(3) = (NR2(1)*NR3(2)-NR2(2)*NR3(1))/sin234
    sinv = NR2(1)*(Normal123(2)*Normal234(3)-Normal123(3)*Normal234(2))+ &
           NR2(2)*(Normal123(3)*Normal234(1)-Normal123(1)*Normal234(3))+ &
           NR2(3)*(Normal123(1)*Normal234(2)-Normal123(2)*Normal234(1))
    if (sinv > 1.0d0) sinv = 1.0d0
    if (sinv < -1.0d0) sinv = -1.0d0
    cosv = Normal123(1)*Normal234(1)+Normal123(2)*Normal234(2)+Normal123(3)*Normal234(3)
    if (cosv > 1.0d0) cosv = 1.0d0
    if (cosv < -1.0d0) cosv = -1.0d0
    cosv0 = cos(xvec(nInt))
    sinv0 = sin(xvec(nInt))
    sindv = sinv*cosv0-cosv*sinv0
    cosdv = cosv*cosv0+sinv*sinv0
    xvec(nInt) = xvec(nInt)+atan2(sindv,cosdv)
    k = k+5
    nInt = nInt+1
  end if
  if (IntType == 5) then
    ! Out of Plane Angle.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    i4 = InterVec(k+4)
    do iv=1,3
      R1(iv) = (AtCoord(iv,i1)-AtCoord(iv,i2))
      R2(iv) = (AtCoord(iv,i2)-AtCoord(iv,i3))
      R3(iv) = (AtCoord(iv,i3)-AtCoord(iv,i4))
    end do
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    SR3 = sqrt(R3(1)**2+R3(2)**2+R3(3)**2)
    R2R3 = (R2(1)*R3(1)+R2(2)*R3(2)+R2(3)*R3(3))
    cos234 = -R2R3/(SR2*SR3)
    if (cos234 > 1.0d0) cos234 = 1.0d0
    if (cos234 < -1.0d0) cos234 = -1.0d0
    sin234 = sqrt(1.0d0-cos234**2)
    Normal234(1) = (R2(2)*R3(3)-R2(3)*R3(2))/(SR2*SR3*sin234)
    Normal234(2) = (R2(3)*R3(1)-R2(1)*R3(3))/(SR2*SR3*sin234)
    Normal234(3) = (R2(1)*R3(2)-R2(2)*R3(1))/(SR2*SR3*sin234)
    sinv = (R1(1)*Normal234(1)+R1(2)*Normal234(2)+R1(3)*Normal234(3))/SR1
    if (sinv > 1.0d0) sinv = 1.0d0
    if (sinv < -1.0d0) sinv = -1.0d0
    xvec(nInt) = asin(sinv)
    k = k+5
    nInt = nInt+1
  end if
  IntType = InterVec(k)
end do

end subroutine Cart_to_Int0
!####
subroutine Int_to_Cart1(InterVec,xvec,AtCoord,NumOfAt,NumInt)
!  Purpose:
!    Transform a set of internal coordinates to cartesian coordinates.
!
!  Written by:
!    Per Ake Malmqvist
!    Dept. of Theoretical Chemistry, Lund University, 2000.

!implicit none
#include "Constants_mula.fh"
integer NumInt, NumOfAt
!integer ncart, iter
!real*8 sum, RMSErr
integer InterVec(*)
real*8 xvec(NumInt)
real*8 AtCoord(3,NumOfAt)
#include "WrkSpc.fh"

call GetMem('Bmatrix','Allo','Real',ipBmatrix,3*NumOfAt*NumInt)
call Int_to_Cart1_a(InterVec,xvec,AtCoord,NumOfAt,NumInt,Work(ipBmatrix))
call GetMem('Bmatrix','Free','Real',ipBmatrix,3*NumOfAt*NumInt)

end subroutine Int_to_Cart1
!####
subroutine Int_to_Cart1_a(InterVec,xvec,AtCoord,NumOfAt,NumInt,Bmatrix)
!  Purpose:
!    Transform a set of internal coordinates to cartesian coordinates.
!
!  Written by:
!    Per Ake Malmqvist
!    Dept. of Theoretical Chemistry, Lund University, 2000.

!implicit none
#include "Constants_mula.fh"
integer NumInt, NumOfAt
integer ncart, iter, j, k, ii, i1, i2, jj, j1, j2
real*8 sum, RMSErr, det
integer InterVec(*)
real*8 xvec(NumInt)
real*8 AtCoord(3,NumOfAt)
real*8 BMatrix(3,NumOfAt,NumInt)
#include "WrkSpc.fh"

! Initialize.

! Purpose: For a given array InterVec, which describes the set of
! internal coordinates, and a given array xvec with values of these
! internal coordinates, find the structure which has the required
! internal coordinates and which is as close as possible to the
! input structure.
! Method: Newton-Raphson iteration, modified by Gauss' approximation
! in each iteration.
ncart = 3*NumOfAt
!vv call GetMem('Bmatrix','Allo','Real',ipBmatrix,3*NumOfAt*NumInt)

!BMatrix = 0.0D0
call dcopy_(3*NumOfAt*NumInt,[0.0d0],0,Bmatrix,1)
call GetMem('EqMat','Allo','Real',ipEqMat,ncart*ncart)
call GetMem('EqRHS','Allo','Real',ipEqRHS,ncart)
call GetMem('xvectmp','Allo','Real',ipxvectmp,NumInt)
!write(6,*) ' Int_to_Cart1 has been called with target internal coordinates:'
!write(6,'(1x,5F16.8)') xvec
iter = 0
10 continue
iter = iter+1
if (iter > 50) then
  write(6,*) ' Int_to_Cart1 fails to converge.'
  call abend()
end if
! Find the current internal coordinates and B matrix:
!xvectmp = xvec
call dcopy_(NumInt,xvec,1,Work(ipxvectmp),1)
call Cart_to_Int1(InterVec,AtCoord,Work(ipxvectmp),BMatrix,NumOfAt,NumInt)

! Present error:
sum = 0.0d0
do k=1,NumInt
  sum = sum+(Work(ipxvectmp+k-1)-xvec(k))**2
end do
RMSErr = sqrt(sum)
! The B matrix is the matrix of partial derivatives of xvectmp with
! respect to AtCoord.
! Form the matrix B(transpose)*B -- This will be the matrix of the
! linearized equation system.
jj = 0
do j1=1,NumOfAt
  do j2=1,3
    jj = jj+1
    ii = 0
    do i1=1,NumOfAt
      do i2=1,3
        ii = ii+1
        sum = 0.0d0
        do k=1,NumInt
          sum = sum+BMatrix(i2,i1,k)*BMatrix(j2,j1,k)
        end do
        Work(ipEqMat+ii+ncart*(jj-1)-1) = sum
      end do
    end do
  end do
end do
! The right-hand side of the equation system:
ii = 0
do i1=1,NumOfAt
  do i2=1,3
    ii = ii+1
    sum = 0.0d0
    do k=1,NumInt
      sum = sum+BMatrix(i2,i1,k)*(xvec(k)-Work(ipxvectmp+k-1))
    end do
    Work(ipEqRHS+ii-1) = sum
  end do
end do
!D write(6,*) ' EqRHS:'
!D write(6,'(1x,5f16.8)') EqRHS
! Automatic step limitation, and in fact ensuring that we go in the
! direction of minimization always, is achieved indirectly:
do j=1,ncart
  Work(ipEqMat+j+ncart*(j-1)-1) = Work(ipEqMat+j+ncart*(j-1)-1)+1.0D-12
end do
! Solve the linear equation system:
call Dool_MULA(Work(ipEqMat),ncart,ncart,Work(ipEqRHS),ncart,1,det)
!D write(6,*) ' Solution vector:'
!D write(6,'(1x,5f16.8)') EqRHS
! After return from dool, EqMat is destroyed and EqRHS is solution.
! Update the cartesian coordinates:
ii = 0
do i1=1,NumOfAt
  do i2=1,3
    ii = ii+1
    AtCoord(i2,i1) = AtCoord(i2,i1)+Work(ipEqRHS+ii-1)
  end do
end do
!write(6,*) ' Iteration iter=',iter
!write(6,*) ' AtCoord array:'
!do i1=1,NumOfAt
!  write(6,'(1x,3f16.8)') (AtCoord(i2,i1),i2=1,3)
!end do
!write(6,*) 'RMSErr:',RMSErr
if (RMSErr > 1.0D-12) goto 10
call GetMem('EqMat','Free','Real',ipEqMat,ncart*ncart)
call GetMem('EqRHS','Free','Real',ipEqRHS,ncart)
call GetMem('xvectmp','Free','Real',ipxvectmp,NumInt)

return

end subroutine Int_to_Cart1_a
