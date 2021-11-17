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
!    InterVec : Integer array - contains the atoms that are used in the calculations of each internal coordinate.
!    AtCoord  : Two dimensional Real array - contains the cartesian coordinates of the atoms.
!
!  Output:
!    xvec     : Real array - internal coordinates.
!    BMatrix  : Real array - Their gradients.
!
!  Written by:
!    Niclas Forsberg, Per Ake Malmqvist,
!    Dept. of Theoretical Chemistry, Lund University, 1997.
!
!  Uses:
!    Linalg, Constants

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: InterVec(*), NumOfAt, NumInt
real(kind=wp), intent(in) :: AtCoord(3,NumOfAt)
real(kind=wp), intent(inout) :: xvec(NumInt)
real(kind=wp), intent(inout) :: BMatrix(3,NumOfAt,NumInt)
integer(kind=iwp) :: i1, i2, i3, i4, inttype, k, n_Int
real(kind=wp) :: Area, cos123, cos234, cosdv, cosv, cosv0, cot123, cot234, G1(3), G2(3), G3(3), Normal123(3), Normal234(3), &
                 NR1(3), NR2(3), NR3(3), R(3), R1(3), R2(3), R2R3, R3(3), sin123, sin234, sindv, sinv, sinv0, SR, SR1, SR2, SR3

k = 1
IntType = InterVec(k)
n_Int = 1
do while (n_Int <= NumInt)
  if (IntType == 1) then
    ! Bond.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    R(:) = AtCoord(:,i1)-AtCoord(:,i2)
    SR = sqrt(R(1)**2+R(2)**2+R(3)**2)
    xvec(n_Int) = SR
    BMatrix(:,i1,n_Int) = R/SR
    BMatrix(:,i2,n_Int) = -R/SR
    !write(u6,*) ' Bond. i1,i2=',i1,i2
    !write(u6,'(1x,a,3F16.8)') ' R=',R
    !write(u6,'(1x,a,3F16.8)') 'SR=',SR
    !write(u6,'(1x,a,3F16.8)') 'BMatrix(:,i1,n_Int)=',BMatrix(:,i1,n_Int)
    !write(u6,'(1x,a,3F16.8)') 'BMatrix(:,i2,n_Int)=',BMatrix(:,i2,n_Int)
    k = k+3
    n_Int = n_Int+1
  end if
  if (IntType == 2) then
    ! Valence Angle.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i2)
    R2(:) = AtCoord(:,i2)-AtCoord(:,i3)
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    cosv = -(R1(1)*R2(1)+R1(2)*R2(2)+R1(3)*R2(3))/(SR1*SR2)
    if (cosv > One) cosv = One
    if (cosv < -One) cosv = -One
    sinv = sqrt(One-cosv**2)
    xvec(n_Int) = atan2(sinv,cosv)
    BMatrix(:,i1,n_Int) = (SR1*R2+SR2*cosv*R1)/(SR1**2*SR2*sinv)
    BMatrix(:,i3,n_Int) = -(SR2*R1+SR1*cosv*R2)/(SR2**2*SR1*sinv)
    BMatrix(:,i2,n_Int) = -(BMatrix(:,i1,n_Int)+BMatrix(:,i3,n_Int))
    !write(u6,*) ' Angle. i1,i2,i3=',i1,i2,i3
    !write(u6,'(1x,a,3F16.8)') ' R1=',R1
    !write(u6,'(1x,a,3F16.8)') ' R2=',R2
    !write(u6,'(1x,a,3F16.8)') 'BMatrix(:,i1,n_Int)=',BMatrix(:,i1,n_Int)
    !write(u6,'(1x,a,3F16.8)') 'BMatrix(:,i2,n_Int)=',BMatrix(:,i2,n_Int)
    !write(u6,'(1x,a,3F16.8)') 'BMatrix(:,i3,n_Int)=',BMatrix(:,i3,n_Int)
    k = k+4
    n_Int = n_Int+1
  end if
  if (IntType == 3) then
    ! Linear Valence Angle.
    !vv?? stop
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i2)
    R2(:) = AtCoord(:,i2)-AtCoord(:,i3)
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    cosv = -(R1(1)*R2(1)+R1(2)*R2(2)+R1(3)*R2(3))/(SR1*SR2)
    if (cosv > One) cosv = One
    if (cosv < -One) cosv = -One
    sinv = sqrt(One-cosv**2)
    xvec(n_Int) = atan2(sinv,cosv)
    BMatrix(1,i1,n_Int) = One/SR1
    BMatrix(1,i3,n_Int) = One/SR2
    BMatrix(1,i2,n_Int) = -(BMatrix(1,i1,n_Int)+BMatrix(1,i3,n_Int))
    BMatrix(2,i1,n_Int+1) = One/SR1
    BMatrix(2,i3,n_Int+1) = One/SR2
    BMatrix(2,i2,n_Int+1) = -(BMatrix(1,i1,n_Int)+BMatrix(1,i3,n_Int))
    n_Int = n_Int+2
    k = k+4
  end if
  if (IntType == 4) then
    ! Torsion.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    i4 = InterVec(k+4)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i2)
    R2(:) = AtCoord(:,i2)-AtCoord(:,i3)
    R3(:) = AtCoord(:,i3)-AtCoord(:,i4)
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    SR3 = sqrt(R3(1)**2+R3(2)**2+R3(3)**2)
    NR1(:) = R1/SR1
    NR2(:) = R2/SR2
    NR3(:) = R3/SR3
    cos123 = -(NR1(1)*NR2(1)+NR1(2)*NR2(2)+NR1(3)*NR2(3))
    if (cos123 > One) cos123 = One
    if (cos123 < -One) cos123 = -One
    sin123 = sqrt(One-cos123**2)
    cot123 = cos123/sin123
    Normal123(1) = (NR1(2)*NR2(3)-NR1(3)*NR2(2))/sin123
    Normal123(2) = (NR1(3)*NR2(1)-NR1(1)*NR2(3))/sin123
    Normal123(3) = (NR1(1)*NR2(2)-NR1(2)*NR2(1))/sin123
    cos234 = -(NR2(1)*NR3(1)+NR2(2)*NR3(2)+NR2(3)*NR3(3))
    if (cos234 > One) cos234 = One
    if (cos234 < -One) cos234 = -One
    sin234 = sqrt(One-cos234**2)
    cot234 = cos234/sin234
    Normal234(1) = (NR2(2)*NR3(3)-NR2(3)*NR3(2))/sin234
    Normal234(2) = (NR2(3)*NR3(1)-NR2(1)*NR3(3))/sin234
    Normal234(3) = (NR2(1)*NR3(2)-NR2(2)*NR3(1))/sin234
    sinv = NR2(1)*(Normal123(2)*Normal234(3)-Normal123(3)*Normal234(2))+ &
           NR2(2)*(Normal123(3)*Normal234(1)-Normal123(1)*Normal234(3))+ &
           NR2(3)*(Normal123(1)*Normal234(2)-Normal123(2)*Normal234(1))
    if (sinv > One) sinv = One
    if (sinv < -One) sinv = -One
    cosv = Normal123(1)*Normal234(1)+Normal123(2)*Normal234(2)+Normal123(3)*Normal234(3)
    if (cosv > One) cosv = One
    if (cosv < -One) cosv = -One
    cosv0 = cos(xvec(n_Int))
    sinv0 = sin(xvec(n_Int))
    sindv = sinv*cosv0-cosv*sinv0
    cosdv = cosv*cosv0+sinv*sinv0
    xvec(n_Int) = xvec(n_Int)+atan2(sindv,cosdv)
    G1(:) = Normal123/(SR1*sin123)
    G3(:) = Normal234/(SR3*sin234)
    G2(:) = (cot123*Normal123+cot234*Normal234)/SR2
    BMatrix(:,i1,n_Int) = G1
    BMatrix(:,i2,n_Int) = -G1+G2
    BMatrix(:,i3,n_Int) = -G2+G3
    BMatrix(:,i4,n_Int) = -G3
    !write(u6,*) ' Torsion. i1,i2,i3,i4=',i1,i2,i3,i4
    !write(u6,'(1x,a,3F16.8)') ' R1=',R1
    !write(u6,'(1x,a,3F16.8)') ' R2=',R2
    !write(u6,'(1x,a,3F16.8)') ' R3=',R3
    !write(u6,'(1x,a,3F16.8)') 'BMatrix(:,i1,n_Int)=',BMatrix(:,i1,n_Int)
    !write(u6,'(1x,a,3F16.8)') 'BMatrix(:,i2,n_Int)=',BMatrix(:,i2,n_Int)
    !write(u6,'(1x,a,3F16.8)') 'BMatrix(:,i3,n_Int)=',BMatrix(:,i3,n_Int)
    !write(u6,'(1x,a,3F16.8)') 'BMatrix(:,i4,n_Int)=',BMatrix(:,i4,n_Int)
    k = k+5
    n_Int = n_Int+1
  end if
  if (IntType == 5) then
    ! Out of Plane Angle.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    i4 = InterVec(k+4)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i2)
    R2(:) = AtCoord(:,i2)-AtCoord(:,i3)
    R3(:) = AtCoord(:,i3)-AtCoord(:,i4)
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    SR3 = sqrt(R3(1)**2+R3(2)**2+R3(3)**2)
    R2R3 = (R2(1)*R3(1)+R2(2)*R3(2)+R2(3)*R3(3))
    cos234 = -R2R3/(SR2*SR3)
    if (cos234 > One) cos234 = One
    if (cos234 < -One) cos234 = -One
    sin234 = sqrt(One-cos234**2)
    Area = SR2*SR3*sin234
    Normal234(1) = (R2(2)*R3(3)-R2(3)*R3(2))/Area
    Normal234(2) = (R2(3)*R3(1)-R2(1)*R3(3))/Area
    Normal234(3) = (R2(1)*R3(2)-R2(2)*R3(1))/Area
    sinv = (R1(1)*Normal234(1)+R1(2)*Normal234(2)+R1(3)*Normal234(3))/SR1
    if (sinv > One) sinv = One
    if (sinv < -One) sinv = -One
    cosv = sqrt(One-sinv**2)
    xvec(n_Int) = atan2(sinv,cosv)
    G1(1) = (R2(2)*R3(3)-R3(2)*R2(3)-(Area*sinv/SR1)*R1(1))/(SR1*Area*cosv)
    G1(2) = (R2(3)*R3(1)-R3(3)*R2(1)-(Area*sinv/SR1)*R1(2))/(SR1*Area*cosv)
    G1(3) = (R2(1)*R3(2)-R3(1)*R2(2)-(Area*sinv/SR1)*R1(3))/(SR1*Area*cosv)
    G2(1) = (R3(2)*R1(3)-R1(2)*R3(3)-(SR1*sinv/Area)*(SR3*R2(1)-R2R3*R3(1)))/(SR1*Area*cosv)
    G2(2) = (R3(3)*R1(1)-R1(3)*R3(1)-(SR1*sinv/Area)*(SR3*R2(2)-R2R3*R3(2)))/(SR1*Area*cosv)
    G2(3) = (R3(1)*R1(2)-R1(1)*R3(2)-(SR1*sinv/Area)*(SR3*R2(3)-R2R3*R3(3)))/(SR1*Area*cosv)
    G3(1) = (R1(2)*R2(3)-R2(2)*R1(3)-(SR1*sinv/Area)*(SR2*R3(1)-R2R3*R2(1)))/(SR1*Area*cosv)
    G3(2) = (R1(3)*R2(1)-R2(3)*R1(1)-(SR1*sinv/Area)*(SR2*R3(2)-R2R3*R2(2)))/(SR1*Area*cosv)
    G3(3) = (R1(1)*R2(2)-R2(1)*R1(2)-(SR1*sinv/Area)*(SR2*R3(3)-R2R3*R2(3)))/(SR1*Area*cosv)
    BMatrix(:,i1,n_Int) = G1
    BMatrix(:,i2,n_Int) = -G1+G2
    BMatrix(:,i3,n_Int) = -G2+G3
    BMatrix(:,i4,n_Int) = -G3
    k = k+5
    n_Int = n_Int+1
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
!    InterVec : Integer array - contains the atoms that are used in the calculations of each internal coordinate.
!    AtCoord  : Two dimensional Real array - contains the cartesian coordinates of the atoms.
!
!  Output:
!    xvec     : Real array - internal coordinates.
!
!  Written by:
!    Niclas Forsberg, Per Ake Malmqvist,
!    Dept. of Theoretical Chemistry, Lund University, 1997.
!
!  Uses:
!    Linalg, Constants

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: InterVec(*), NumOfAt, NumInt
real(kind=wp), intent(in) :: AtCoord(3,NumOfAt)
real(kind=wp), intent(inout) :: xvec(NumInt)
integer(kind=iwp) :: i1, i2, i3, i4, inttype, k, n_Int
real(kind=wp) :: cos123, cos234, cosdv, cosv, cosv0, Normal123(3), Normal234(3), NR1(3), NR2(3), NR3(3), R(3), R1(3), R2(3), R2R3, &
                 R3(3), sin123, sin234, sindv, sinv, sinv0, SR, SR1, SR2, SR3

k = 1
IntType = InterVec(k)
n_Int = 1
do while (n_Int <= NumInt)
  if (IntType == 1) then
    ! Bond.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    R(:) = AtCoord(:,i1)-AtCoord(:,i2)
    SR = sqrt(R(1)**2+R(2)**2+R(3)**2)
    xvec(n_Int) = SR
    k = k+3
    n_Int = n_Int+1
  end if
  if (IntType == 2) then
    ! Valence Angle.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i2)
    R2(:) = AtCoord(:,i2)-AtCoord(:,i3)
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    cosv = -(R1(1)*R2(1)+R1(2)*R2(2)+R1(3)*R2(3))/(SR1*SR2)
    if (cosv > One) cosv = One
    if (cosv < -One) cosv = -One
    sinv = sqrt(One-cosv**2)
    xvec(n_Int) = atan2(sinv,cosv)
    k = k+4
    n_Int = n_Int+1
  end if
  if (IntType == 3) then
    ! Linear Valence Angle.
    !vv? stop
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i2)
    R2(:) = AtCoord(:,i3)-AtCoord(:,i2)
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    NR1(:) = R1/SR1
    NR2(:) = R2/SR2
    cosv = NR1(1)*NR2(1)+NR1(2)*NR2(2)+NR1(3)*NR2(3)
    if (cosv > One) cosv = One
    if (cosv < -One) cosv = -One
    xvec(n_Int) = acos(cosv)
    k = k+4
    n_Int = n_Int+2
  end if
  if (IntType == 4) then
    ! Torsion.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    i4 = InterVec(k+4)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i2)
    R2(:) = AtCoord(:,i2)-AtCoord(:,i3)
    R3(:) = AtCoord(:,i3)-AtCoord(:,i4)
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    SR3 = sqrt(R3(1)**2+R3(2)**2+R3(3)**2)
    NR1(:) = R1/SR1
    NR2(:) = R2/SR2
    NR3(:) = R3/SR3
    cos123 = -(NR1(1)*NR2(1)+NR1(2)*NR2(2)+NR1(3)*NR2(3))
    if (cos123 > One) cos123 = One
    if (cos123 < -One) cos123 = -One
    sin123 = sqrt(One-cos123**2)
    Normal123(1) = (NR1(2)*NR2(3)-NR1(3)*NR2(2))/sin123
    Normal123(2) = (NR1(3)*NR2(1)-NR1(1)*NR2(3))/sin123
    Normal123(3) = (NR1(1)*NR2(2)-NR1(2)*NR2(1))/sin123
    cos234 = -(NR2(1)*NR3(1)+NR2(2)*NR3(2)+NR2(3)*NR3(3))
    if (cos234 > One) cos234 = One
    if (cos234 < -One) cos234 = -One
    sin234 = sqrt(One-cos234**2)
    Normal234(1) = (NR2(2)*NR3(3)-NR2(3)*NR3(2))/sin234
    Normal234(2) = (NR2(3)*NR3(1)-NR2(1)*NR3(3))/sin234
    Normal234(3) = (NR2(1)*NR3(2)-NR2(2)*NR3(1))/sin234
    sinv = NR2(1)*(Normal123(2)*Normal234(3)-Normal123(3)*Normal234(2))+ &
           NR2(2)*(Normal123(3)*Normal234(1)-Normal123(1)*Normal234(3))+ &
           NR2(3)*(Normal123(1)*Normal234(2)-Normal123(2)*Normal234(1))
    if (sinv > One) sinv = One
    if (sinv < -One) sinv = -One
    cosv = Normal123(1)*Normal234(1)+Normal123(2)*Normal234(2)+Normal123(3)*Normal234(3)
    if (cosv > One) cosv = One
    if (cosv < -One) cosv = -One
    cosv0 = cos(xvec(n_Int))
    sinv0 = sin(xvec(n_Int))
    sindv = sinv*cosv0-cosv*sinv0
    cosdv = cosv*cosv0+sinv*sinv0
    xvec(n_Int) = xvec(n_Int)+atan2(sindv,cosdv)
    k = k+5
    n_Int = n_Int+1
  end if
  if (IntType == 5) then
    ! Out of Plane Angle.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    i4 = InterVec(k+4)
    R1(:) = AtCoord(:,i1)-AtCoord(:,i2)
    R2(:) = AtCoord(:,i2)-AtCoord(:,i3)
    R3(:) = AtCoord(:,i3)-AtCoord(:,i4)
    SR1 = sqrt(R1(1)**2+R1(2)**2+R1(3)**2)
    SR2 = sqrt(R2(1)**2+R2(2)**2+R2(3)**2)
    SR3 = sqrt(R3(1)**2+R3(2)**2+R3(3)**2)
    R2R3 = (R2(1)*R3(1)+R2(2)*R3(2)+R2(3)*R3(3))
    cos234 = -R2R3/(SR2*SR3)
    if (cos234 > One) cos234 = One
    if (cos234 < -One) cos234 = -One
    sin234 = sqrt(One-cos234**2)
    Normal234(1) = (R2(2)*R3(3)-R2(3)*R3(2))/(SR2*SR3*sin234)
    Normal234(2) = (R2(3)*R3(1)-R2(1)*R3(3))/(SR2*SR3*sin234)
    Normal234(3) = (R2(1)*R3(2)-R2(2)*R3(1))/(SR2*SR3*sin234)
    sinv = (R1(1)*Normal234(1)+R1(2)*Normal234(2)+R1(3)*Normal234(3))/SR1
    if (sinv > One) sinv = One
    if (sinv < -One) sinv = -One
    xvec(n_Int) = asin(sinv)
    k = k+5
    n_Int = n_Int+1
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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: InterVec(*), NumOfAt, NumInt
real(kind=wp), intent(in) :: xvec(NumInt)
real(kind=wp), intent(inout) :: AtCoord(3,NumOfAt)
integer(kind=iwp) :: i1, i2, ii, iter, j, j1, j2, jj, k, ncart
real(kind=wp) :: det, RMSErr, rsum
real(kind=wp), allocatable :: BMatrix(:,:,:), EqMat(:,:), EqRHS(:), xvectmp(:)

! Initialize.

call mma_allocate(BMatrix,3,NumOfAt,NumInt,label='BMatrix')

! Purpose: For a given array InterVec, which describes the set of
! internal coordinates, and a given array xvec with values of these
! internal coordinates, find the structure which has the required
! internal coordinates and which is as close as possible to the
! input structure.
! Method: Newton-Raphson iteration, modified by Gauss' approximation
! in each iteration.
ncart = 3*NumOfAt

BMatrix(:,:,:) = Zero
call mma_allocate(EqMat,ncart,ncart,label='EqMat')
call mma_allocate(EqRHS,ncart,label='EqRHS')
call mma_allocate(xvectmp,NumInt,label='xvectmp')
!write(u6,*) ' Int_to_Cart1 has been called with target internal coordinates:'
!write(u6,'(1x,5F16.8)') xvec
iter = 0
RMSErr = huge(RMSErr)
do while (RMSErr > 1.0e-12_wp)
  iter = iter+1
  if (iter > 50) then
    write(u6,*) ' Int_to_Cart1 fails to converge.'
    call abend()
  end if
  ! Find the current internal coordinates and B matrix:
  !xvectmp = xvec
  xvectmp(:) = xvec
  call Cart_to_Int1(InterVec,AtCoord,xvectmp,BMatrix,NumOfAt,NumInt)

  ! Present error:
  rsum = Zero
  do k=1,NumInt
    rsum = rsum+(xvectmp(k)-xvec(k))**2
  end do
  RMSErr = sqrt(rsum)
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
          rsum = Zero
          do k=1,NumInt
            rsum = rsum+BMatrix(i2,i1,k)*BMatrix(j2,j1,k)
          end do
          EqMat(ii,jj) = rsum
        end do
      end do
    end do
  end do
  ! The right-hand side of the equation system:
  ii = 0
  do i1=1,NumOfAt
    do i2=1,3
      ii = ii+1
      rsum = Zero
      do k=1,NumInt
        rsum = rsum+BMatrix(i2,i1,k)*(xvec(k)-xvectmp(k))
      end do
      EqRHS(ii) = rsum
    end do
  end do
  !D write(u6,*) ' EqRHS:'
  !D write(u6,'(1x,5f16.8)') EqRHS
  ! Automatic step limitation, and in fact ensuring that we go in the
  ! direction of minimization always, is achieved indirectly:
  do j=1,ncart
    EqMat(j,j) = EqMat(j,j)+1.0e-12_wp
  end do
  ! Solve the linear equation system:
  call Dool_MULA(EqMat,ncart,ncart,EqRHS,ncart,1,det)
  !D write(u6,*) ' Solution vector:'
  !D write(u6,'(1x,5f16.8)') EqRHS
  ! After return from dool, EqMat is destroyed and EqRHS is solution.
  ! Update the cartesian coordinates:
  ii = 0
  do i1=1,NumOfAt
    do i2=1,3
      ii = ii+1
      AtCoord(i2,i1) = AtCoord(i2,i1)+EqRHS(ii)
    end do
  end do
  !write(u6,*) ' Iteration iter=',iter
  !write(u6,*) ' AtCoord array:'
  !do i1=1,NumOfAt
  !  write(u6,'(1x,3f16.8)') (AtCoord(i2,i1),i2=1,3)
  !end do
  !write(u6,*) 'RMSErr:',RMSErr
end do
call mma_deallocate(EqMat)
call mma_deallocate(EqRHS)
call mma_deallocate(xvectmp)
call mma_deallocate(BMatrix)

return

end subroutine Int_to_Cart1
