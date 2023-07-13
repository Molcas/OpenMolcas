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

subroutine OutofP(xyz,nCent,Teta,Bt,lWrite,lWarn,Label,dBt,ldB)

use Constants, only: Zero, One, Two, Nine, Pi, deg2rad
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nCent
real(kind=wp), intent(inout) :: xyz(3,nCent)
real(kind=wp), intent(out) :: Teta, Bt(3,nCent)
logical(kind=iwp), intent(in) :: lWrite, lWarn, ldB
character(len=8), intent(in) :: Label
real(kind=wp), intent(inout) :: dBt(3,nCent,3,nCent)
integer(kind=iwp) :: ix, iy, iz, jx, jy, jz, mCent
real(kind=wp) :: BR14X(3,3), C14X(3,3), CosFi1, CosFi2, CosFi3, dBR14X(3,3,3,3), dFi1, dFi2, dFi3, dTeta, e41x, e41y, e41z, e42x, &
                 e42y, e42z, e43x, e43y, e43z, Fi1, Fi2, Fi3, Q41, Q42, Q43, R41KV, R42(3), R42KV, R43(3), R43KV, RX1, RX2, RX3, &
                 RY1, RY2, RY3, RZ1, RZ2, RZ3
!#define _Test_Numerical_
#ifdef _Test_Numerical_
real(kind=wp) :: Bt_temp_m(3,4); Bt_temp_p(3,4)
#endif
real(kind=wp), external :: ArCos

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
! First some diagnostics

! 4-->1 (Bond)
RX1 = xyz(1,1)-xyz(1,4)
RY1 = xyz(2,1)-xyz(2,4)
RZ1 = xyz(3,1)-xyz(3,4)
R41KV = RX1*RX1+RY1*RY1+RZ1*RZ1
Q41 = sqrt(R41KV)
e41x = RX1/Q41
e41y = RY1/Q41
e41z = RZ1/Q41
! 4-->2 (Bond in plane)
RX2 = xyz(1,2)-xyz(1,4)
RY2 = xyz(2,2)-xyz(2,4)
RZ2 = xyz(3,2)-xyz(3,4)
R42KV = RX2*RX2+RY2*RY2+RZ2*RZ2
Q42 = sqrt(R42KV)
e42x = RX2/Q42
e42y = RY2/Q42
e42z = RZ2/Q42
! 4-->3 (Bond in plane)
RX3 = xyz(1,3)-xyz(1,4)
RY3 = xyz(2,3)-xyz(2,4)
RZ3 = xyz(3,3)-xyz(3,4)
R43KV = RX3*RX3+RY3*RY3+RZ3*RZ3
Q43 = sqrt(R43KV)
e43x = RX3/Q43
e43y = RY3/Q43
e43z = RZ3/Q43

! Get the angle between e43 and e42

CosFi1 = e43x*e42x+e43y*e42y+e43z*e42z

Fi1 = ArCos(CosFi1)
if (abs(CosFi1) > One) call RecPrt('xyz(1)',' ',xyz,3,4)
dFi1 = Fi1/deg2rad
if (lWarn .and. ((dFi1 > 177.5_wp) .or. (dFi1 < 2.5_wp))) write(u6,*) 'Warning: auxiliary Angle close to end of range'

! Dirty exit! This happens when an earlier structure is ill defined.

if (abs(Fi1-Pi) < 1.0e-13_wp) then
  Teta = Zero
  Bt(:,:) = Zero
  return
end if

! Get the angle between e41 and e43

CosFi2 = e41x*e43x+e41y*e43y+e41z*e43z

Fi2 = ArCos(CosFi2)
if (abs(CosFi2) > One) call RecPrt('xyz(2)',' ',xyz,3,4)
dFi2 = Fi2/deg2rad
if (lWarn .and. ((dFi2 > 177.5_wp) .or. (dFi2 < 2.5_wp))) write(u6,*) 'Warning: auxiliary Angle close to end of range'

! Get the angle between e41 and e42

CosFi3 = e41x*e42x+e41y*e42y+e41z*e42z

Fi3 = ArCos(CosFi3)
if (abs(CosFi3) > One) call RecPrt('xyz(3)',' ',xyz,3,4)
dFi3 = Fi3/deg2rad
if (lWarn .and. ((dFi3 > 177.5_wp) .or. (dFi3 < 2.5_wp))) write(u6,*) 'Warning: auxiliary Angle close to end of range'
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt('xyz',' ',xyz,3,nCent)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! The first two centers are trivially

C14X(:,1) = xyz(:,1)
C14X(:,2) = xyz(:,4)

! The 3rd is

R42(:) = xyz(:,2)-xyz(:,4)
R43(:) = xyz(:,3)-xyz(:,4)
C14X(1,3) = R42(2)*R43(3)-R42(3)*R43(2)
C14X(2,3) = R42(3)*R43(1)-R42(1)*R43(3)
C14X(3,3) = R42(1)*R43(2)-R42(2)*R43(1)

! Exit if 2-3-4 are collinear
! (equivalent to the above check, but this is more concrete)

if ((C14X(1,3)**2+C14X(2,3)**2+C14X(3,3)**2) < 1.0e-10_wp) then
  Teta = Zero
  Bt(:,:) = Zero
  return
end if
C14X(:,3) = C14X(:,3)+xyz(:,4)

mCent = 3
call Bend(C14X,mCent,Teta,BR14X,.false.,.false.,Label,dBR14X,ldB)

Teta = Teta-Pi/Two
dTeta = Teta/deg2rad
if (lWarn .and. ((dTeta > 87.5_wp) .or. (dTeta < -87.5_wp))) write(u6,*) 'Warning: Out of plane angle close to end of range'
if (lWrite) write(u6,'(1X,A,A,F10.4,A,F10.4,A)') Label,' : Out of plane angle=',dTeta,'/degree, ',Teta,'/rad'

! Compute the WDC matrix

do ix=1,3
  iy = mod(ix+1,4)+(ix+1)/4
  iz = mod(iy+1,4)+(iy+1)/4

  Bt(ix,1) = -BR14X(ix,1)
  Bt(ix,2) = R43(iz)*BR14X(iy,3)-R43(iy)*BR14X(iz,3)
  Bt(ix,3) = -R42(iz)*BR14X(iy,3)+R42(iy)*BR14X(iz,3)

  Bt(ix,4) = -(Bt(ix,1)+Bt(ix,2)+Bt(ix,3))

end do
#ifdef _DEBUGPRINT_
call RecPrt('Outofp: R43',' ',R43,1,3)
call RecPrt('Outofp: R43',' ',R42,1,3)
call RecPrt('Outofp: BR14X',' ',BR14X,3,3)
call RecPrt('Outofp: B matrix',' ',Bt,3,nCent)
#endif

if (ldB) then

  ! Compute the derivative of the WDC matrix.

  dBt(:,:,:,:) = Nine
  do ix=1,3
    iy = mod(ix+1,4)+(ix+1)/4
    iz = mod(iy+1,4)+(iy+1)/4
    do jx=1,ix
      jy = mod(jx+1,4)+(jx+1)/4
      jz = mod(jy+1,4)+(jy+1)/4
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Do Block (1,1), (2,1), (3,1). Construct
      ! block (1,2) and (1,3) by symmetry.

      dBt(ix,1,jx,1) = -dBR14X(ix,1,jx,1)

      dBt(ix,2,jx,1) = +R43(iz)*dBR14X(iy,3,jx,1)-R43(iy)*dBR14X(iz,3,jx,1)
      dBt(jx,1,ix,2) = dBt(ix,2,jx,1)

      dBt(ix,3,jx,1) = -R42(iz)*dBR14X(iy,3,jx,1)+R42(iy)*dBR14X(iz,3,jx,1)
      dBt(jx,1,ix,3) = dBt(ix,3,jx,1)

      ! Do block (4,1) by translational invariance and
      ! (1,4) by symmetry

      dBt(ix,4,jx,1) = -(dBt(ix,1,jx,1)+dBt(ix,2,jx,1)+dBt(ix,3,jx,1))
      dBt(jx,1,ix,4) = dBt(ix,4,jx,1)
      if (ix /= jx) then
        dBt(jx,1,ix,1) = dBt(ix,1,jx,1)
        dBt(jx,2,ix,1) = +R43(jz)*dBR14X(jy,3,ix,1)-R43(jy)*dBR14X(jz,3,ix,1)
        dBt(ix,1,jx,2) = dBt(jx,2,ix,1)
        dBt(jx,3,ix,1) = -R42(jz)*dBR14X(jy,3,ix,1)+R42(jy)*dBR14X(jz,3,ix,1)
        dBt(ix,1,jx,3) = dBt(jx,3,ix,1)
        dBt(jx,4,ix,1) = -(dBt(jx,1,ix,1)+dBt(jx,2,ix,1)+dBt(jx,3,ix,1))
        dBt(ix,1,jx,4) = dBt(jx,4,ix,1)
      end if
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Do block (2,2), and (3,2). Construct block
      ! (2,3) by symmetry

      dBt(ix,2,jx,2) = -R43(iz)*(R43(jz)*dBR14X(iy,3,jy,3)-R43(jy)*dBR14X(iy,3,jz,3))+ &
                       R43(iy)*(R43(jz)*dBR14X(iz,3,jy,3)-R43(jy)*dBR14X(iz,3,jz,3))
      dBt(ix,3,jx,2) = R42(iz)*(R43(jz)*dBR14X(iy,3,jy,3)-R43(jy)*dBR14X(iy,3,jz,3))- &
                       R42(iy)*(R43(jz)*dBR14X(iz,3,jy,3)-R43(jy)*dBR14X(iz,3,jz,3))
      if (ix == jz) dBt(ix,3,jx,2) = dBt(ix,3,jx,2)+BR14X(jy,3)
      if (ix == jy) dBt(ix,3,jx,2) = dBt(ix,3,jx,2)-BR14X(jz,3)
      dBt(jx,2,ix,3) = dBt(ix,3,jx,2)

      ! Do block (4,2) by translational invariance and (2,4) by symmetry

      dBt(ix,4,jx,2) = -(dBt(ix,1,jx,2)+dBt(ix,2,jx,2)+dBt(ix,3,jx,2))
      dBt(jx,2,ix,4) = dBt(ix,4,jx,2)
      if (ix /= iy) then
        dBt(jx,2,ix,2) = dBt(ix,2,jx,2)

        dBt(jx,3,ix,2) = R42(jz)*(R43(iz)*dBR14X(jy,3,iy,3)-R43(iy)*dBR14X(jy,3,iz,3))- &
                         R42(jy)*(R43(iz)*dBR14X(jz,3,iy,3)-R43(iy)*dBR14X(jz,3,iz,3))
        if (jx == iz) dBt(jx,3,ix,2) = dBt(jx,3,ix,2)+BR14X(iy,3)
        if (jx == iy) dBt(jx,3,ix,2) = dBt(jx,3,ix,2)-BR14X(iz,3)
        dBt(ix,2,jx,3) = dBt(jx,3,ix,2)
        dBt(jx,4,ix,2) = -(dBt(jx,1,ix,2)+dBt(jx,2,ix,2)+dBt(jx,3,ix,2))
        dBt(ix,2,jx,4) = dBt(jx,4,ix,2)
      end if
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Do block (3,3)

      dBt(ix,3,jx,3) = -R42(iz)*(R42(jz)*dBR14X(iy,3,jy,3)-R42(jy)*dBR14X(iy,3,jz,3))+ &
                       R42(iy)*(R42(jz)*dBR14X(iz,3,jy,3)-R42(jy)*dBR14X(iz,3,jz,3))

      ! Do (4,3) byh translational invariance and (3,4) by
      ! symmetry

      dBt(ix,4,jx,3) = -(dBt(ix,1,jx,3)+dBt(ix,2,jx,3)+dBt(ix,3,jx,3))
      dBt(jx,3,ix,4) = dBt(ix,4,jx,3)
      if (ix /= iy) then
        dBt(jx,3,ix,3) = dBt(ix,3,jx,3)
        dBt(jx,4,ix,3) = -(dBt(jx,1,ix,3)+dBt(jx,2,ix,3)+dBt(jx,3,ix,3))
        dBt(ix,3,jx,4) = dBt(jx,4,ix,3)
      end if
      !                                                                *
      !*****************************************************************
      !                                                                *
      !        Finally do (4,4) by translational invariance and symmetry

      dBt(ix,4,jx,4) = -(dBt(ix,1,jx,4)+dBt(ix,2,jx,4)+dBt(ix,3,jx,4))
      if (ix /= jx) dBt(jx,4,ix,4) = dBt(ix,4,jx,4)
      !                                                                *
      !*****************************************************************
      !                                                                *
    end do
  end do
# ifdef _DEBUGPRINT_
  call RecPrt('dBt','(4(3F7.2,2X))',dBt,12,12)
# endif
end if
Bt(:,:) = -Bt(:,:)
!dBt(:,:,:,:) = -dBt(:,:,:,:)
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _Test_Numerical_

! test of Bt

delta = 1.0e-5_wp
do iAtom=1,4
  do iCar=1,3
    tmp = xyz(iCar,iAtom)
    xyz(iCar,iAtom) = tmp+delta
    call OutofP0(xyz,4,Tetap,ldB)
    xyz(iCar,iAtom) = tmp-delta
    call OutofP0(xyz,4,Tetam,ldB)
    xyz(iCar,iAtom) = tmp

    dbdx = (Tetap-Tetam)/(Two*delta)
    if (abs(dbdx-Bt(iCar,iAtom)) > delta) then
      write(u6,*) dbdx,Bt(iCar,iAtom)
      call Abend()
    end if
  end do
end do

! test of dBt

if (ldb) then
  do iAtom=1,4
    do iCar=1,3
      tmp = xyz(iCar,iAtom)
      xyz(iCar,iAtom) = tmp+delta
      call OutofP1(xyz,4,Tetap,Bt_temp_p,ldB)
      xyz(iCar,iAtom) = tmp-delta
      call OutofP1(xyz,4,Tetam,Bt_temp_m,ldB)
      xyz(iCar,iAtom) = tmp

      do jAtom=1,4
        do jCar=1,3
          ddbddx = (Bt_temp_p(jCar,jAtom)-Bt_temp_m(jCar,jAtom))/(Two*delta)
          if (abs(ddbddx-dBt(iCar,iAtom,jCar,jAtom)) > delta) then
            write(u6,*) ddbddx,dBt(iCar,iAtom,jCar,jAtom)
            call Abend()
          end if
        end do
      end do

    end do ! iCar
  end do   ! iAtom
end if
#endif

return

end subroutine OutofP
