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
! Copyright (C) 1996, Roland Lindh                                     *
!***********************************************************************

subroutine Trsn(xyz,nCent,Tau,Bt,lWrite,lWarn,Label,dBt,ldB)
!***********************************************************************
!                                                                      *
! Reference: Molecular Vibrations, E. Bright Wilson, Jr, J. C. Decicius*
!            and Paul C. Cross, Sec. 4-1, Eq. 20-24                    *
!                                                                      *
! R.Lindh May-June '96                                                 *
!***********************************************************************

use Constants, only: Zero, One, Two, Pi, deg2rad
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nCent
real(kind=wp), intent(in) :: xyz(3,nCent)
real(kind=wp), intent(out) :: Tau, Bt(3,nCent)
logical(kind=iwp), intent(in) :: lWrite, lWarn, ldB
character(len=8), intent(in) :: Label
real(kind=wp), intent(inout) :: dBt(3,nCent,3,nCent)
integer(kind=iwp) :: ix, iy, iz, jx, jy, jz, mCent
real(kind=wp) :: Bf2(3,3), Bf3(3,3), BRij(3,2), BRjk(3,2), BRkl(3,2), CosFi2, CosFi3, CosTau, dBRij(3,2,3,2), dBRjk(3,2,3,2), &
                 dBRkl(3,2,3,2), dFi2, dFi3, dTau, Dum(1), Fac, Fi2, Fi3, Rij1, Rjk1, Rkl1, SinFi2, SinFi3, SinTau

mCent = 2
call Strtch(xyz(1,1),mCent,Rij1,BRij,.false.,Label,dBRij,ldB)
call Strtch(xyz(1,2),mCent,Rjk1,BRjk,.false.,Label,dBRjk,ldB)
call Strtch(xyz(1,3),mCent,Rkl1,BRkl,.false.,Label,dBRkl,ldB)
mCent = 3
call Bend(xyz(1,1),mCent,Fi2,Bf2,.false.,.false.,Label,Dum,.false.)
SinFi2 = sin(Fi2)
CosFi2 = cos(Fi2)
call Bend(xyz(1,2),mCent,Fi3,Bf3,.false.,.false.,Label,Dum,.false.)
SinFi3 = sin(Fi3)
CosFi3 = cos(Fi3)

if (SinFi2*SinFi3 < 1.0e-13_wp) then
  Tau = Zero
  dTau = Zero
  if (lWrite) write(u6,1) Label,-dTau,-Tau
  return
end if

! Get the angle between the two planes, i.e. the
! angle between the normal vectors.
!
! r123 * r234 = CosTau

CosTau = ((BRij(2,1)*BRjk(3,2)-BRij(3,1)*BRjk(2,2))*(BRjk(2,1)*BRkl(3,2)-BRjk(3,1)*BRkl(2,2))+ &
          (BRij(3,1)*BRjk(1,2)-BRij(1,1)*BRjk(3,2))*(BRjk(3,1)*BRkl(1,2)-BRjk(1,1)*BRkl(3,2))+ &
          (BRij(1,1)*BRjk(2,2)-BRij(2,1)*BRjk(1,2))*(BRjk(1,1)*BRkl(2,2)-BRjk(2,1)*BRkl(1,2)))/(SinFi2*SinFi3)

! For the vector product of the two vectors. This
! will give a vector parallel to e23. The direction
! relative to e23 defines the sign.
!
! e123 X e234 = SinTau * e23

SinTau = (BRij(1,2)*(BRjk(2,1)*BRkl(3,2)-BRjk(3,1)*BRkl(2,2))+BRij(2,2)*(BRjk(3,1)*BRkl(1,2)-BRjk(1,1)*BRkl(3,2))+ &
          BRij(3,2)*(BRjk(1,1)*BRkl(2,2)-BRjk(2,1)*BRkl(1,2)))/(SinFi2*SinFi3)

! (-Pi < Tau <= Pi)

Tau = atan2(SinTau,CosTau)
if (abs(Tau) == Pi) Tau = Pi

dTau = Tau/deg2rad
dFi2 = Fi2/deg2rad
dFi3 = Fi3/deg2rad
if (lWarn) then
  if ((dTau > 177.5_wp) .or. (dTau < -177.5_wp)) call WarningMessage(1,' Warning: dihedral angle close to end of range')
  if ((dFi2 > 177.5_wp) .or. (dFi2 < 2.5_wp)) call WarningMessage(1,' Warning: bond angle 2 close to end of range')
  if ((dFi3 > 177.5_wp) .or. (dFi3 < 2.5_wp)) call WarningMessage(1,' Warning: bond angle 3 close to end of range')
end if
if (lWrite) write(u6,1) Label,-dTau,-Tau

! Compute the WDC matrix.

do ix=1,3
  iy = ix+1
  if (iy > 3) iy = iy-3
  iz = iy+1
  if (iz > 3) iz = iz-3
  Bt(ix,1) = (BRij(iy,2)*BRjk(iz,2)-BRij(iz,2)*BRjk(iy,2))/(Rij1*SinFi2**2)
  Bt(ix,4) = (BRkl(iy,1)*BRjk(iz,1)-BRkl(iz,1)*BRjk(iy,1))/(Rkl1*SinFi3**2)
  Bt(ix,2) = -((Rjk1-Rij1*CosFi2)*Bt(ix,1)+Rkl1*CosFi3*Bt(ix,4))/Rjk1
  Bt(ix,3) = -(Bt(ix,1)+Bt(ix,2)+Bt(ix,4))
end do

if (ldB) then

  ! Compute the derivative of the WDC matrix.

  do ix=1,3
    iy = ix+1
    if (iy > 3) iy = iy-3
    iz = iy+1
    if (iz > 3) iz = iz-3
    do jx=1,ix
      jy = jx+1
      if (jy > 3) jy = jy-3
      jz = jy+1
      if (jz > 3) jz = jz-3

      dBt(ix,1,jx,1) = (dBRij(ix,1,jy,2)*BRjk(jz,2)-dBRij(ix,1,jz,2)*BRjk(jy,2)- &
                        Bt(jx,1)*(BRij(ix,1)*SinFi2**2+Rij1*Two*SinFi2*CosFi2*Bf2(ix,1)))/(Rij1*SinFi2**2)
      dBt(ix,1,jx,2) = -((-BRij(ix,1)*CosFi2+Rij1*SinFi2*Bf2(ix,1))*Bt(jx,1)+(Rjk1-Rij1*CosFi2)*dBt(ix,1,jx,1))/Rjk1
      dBt(jx,2,ix,1) = dBt(ix,1,jx,2)
      dBt(ix,1,jx,4) = Zero
      dBt(jx,4,ix,1) = dBt(ix,1,jx,4)
      dBt(ix,1,jx,3) = -(dBt(ix,1,jx,1)+dBt(ix,1,jx,2))
      dBt(jx,3,ix,1) = dBt(ix,1,jx,3)
      dBt(ix,4,jx,4) = (dBRkl(ix,2,jy,1)*BRjk(jz,1)-dBRkl(ix,2,jz,1)*BRjk(jy,1)- &
                        Bt(jx,4)*(BRkl(ix,2)*SinFi3**2+Rkl1*Two*SinFi3*CosFi3*Bf3(ix,3)))/(Rkl1*SinFi3**2)
      dBt(ix,4,jx,3) = -((-BRkl(ix,2)*CosFi3+Rkl1*SinFi3*Bf3(ix,3))*Bt(jx,4)+(Rjk1-Rkl1*CosFi3)*dBt(ix,4,jx,4))/Rjk1
      dBt(jx,3,ix,4) = dBt(ix,4,jx,3)
      dBt(ix,4,jx,2) = -(dBt(ix,4,jx,4)+dBt(ix,4,jx,3))
      dBt(jx,2,ix,4) = dBt(ix,4,jx,2)
      if (ix /= jx) then
        dBt(jx,1,ix,1) = dBt(ix,1,jx,1)
        dBt(ix,4,jx,1) = Zero
        dBt(jx,4,ix,4) = dBt(ix,4,jx,4)
        dBt(jx,1,ix,4) = dBt(ix,4,jx,1)
        dBt(jx,1,ix,2) = -((-BRij(jx,1)*CosFi2+Rij1*SinFi2*Bf2(jx,1))*Bt(ix,1)+(Rjk1-Rij1*CosFi2)*dBt(jx,1,ix,1))/Rjk1
        dBt(ix,2,jx,1) = dBt(jx,1,ix,2)
        dBt(ix,3,jx,1) = -(dBt(ix,1,jx,1)+dBt(ix,2,jx,1)+dBt(ix,4,jx,1))
        dBt(jx,1,ix,3) = dBt(ix,3,jx,1)
        dBt(jx,4,ix,3) = -((-BRkl(jx,2)*CosFi3+Rkl1*SinFi3*Bf3(jx,3))*Bt(ix,4)+(Rjk1-Rkl1*CosFi3)*dBt(jx,4,ix,4))/Rjk1
        dBt(ix,3,jx,4) = dBt(jx,4,ix,3)
        dBt(ix,2,jx,4) = -(dBt(ix,4,jx,4)+dBt(ix,3,jx,4))
        dBt(jx,4,ix,2) = dBt(ix,2,jx,4)
      end if
      dBt(ix,2,jx,3) = -((BRjk(ix,1)+Rkl1*SinFi3*Bf3(ix,1))*Bt(jx,4)+(Rjk1-Rkl1*CosFi3)*dBt(ix,2,jx,4)+ &
                         (BRij(ix,2)*CosFi2-Rij1*SinFi2*Bf2(ix,2))*Bt(jx,1)+Rij1*CosFi2*dBt(ix,2,jx,1)+Bt(jx,3)*BRjk(ix,1))/Rjk1
      dBt(jx,3,ix,2) = dBt(ix,2,jx,3)
      dBt(ix,2,jx,2) = -(dBt(ix,2,jx,1)+dBt(ix,2,jx,4)+dBt(ix,2,jx,3))
      dBt(ix,3,jx,3) = -(dBt(ix,2,jx,3)+dBt(ix,1,jx,3)+dBt(ix,4,jx,3))
      if (ix /= jx) then
        dBt(ix,3,jx,2) = -(dBt(ix,2,jx,2)+dBt(ix,1,jx,2)+dBt(ix,4,jx,2))
        dBt(jx,2,ix,3) = dBt(ix,3,jx,2)
        dBt(jx,2,ix,2) = dBt(ix,2,jx,2)
        dBt(jx,3,ix,3) = dBt(ix,3,jx,3)
      end if

    end do
  end do

end if

Fac = -One
Tau = Fac*Tau
Bt(:,:) = Fac*Bt(:,:)
if (ldB) dBt(:,:,:,:) = Fac*dBt(:,:,:,:)

return

1 format(1X,A,' : Dihedral= ',F10.4,'   / Degree  ',F10.6,' / rad')

end subroutine Trsn
