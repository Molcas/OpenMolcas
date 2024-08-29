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

subroutine hmod(gx,gy,gz,V,EFx,EFy,EFz,Cavxyz,lmax_)
!***********************************************************************
!                                                                      *
!     Object: to compute the potential and electric field in a point   *
!             due to the multipole moment expansion at origin.         *
!                                                                      *
!***********************************************************************

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: lmax_
real(kind=wp), intent(in) :: gx, gy, gz, Cavxyz((lMax_+1)*(lMax_+2)*(lMax_+3)/6)
real(kind=wp), intent(out) :: V, EFx, EFy, EFz
integer(kind=iwp) :: ix, iy, iz, lMax
real(kind=wp) :: aX, aY, aZ, xEff, xyEff, xyzEff
!---- Statement functions
integer(kind=iwp) :: iOff, Indx
iOff(ix,iy,iz) = (ix+iy+iz)*(ix+iy+iz+1)*(ix+iy+iz+2)/6
Indx(ix,iy,iz) = iOff(ix,iy,iz)+(iy+iz)*(iy+iz+1)/2+iz+1

V = Zero
EFx = Zero
EFy = Zero
EFz = Zero

lmax = lmax_-1
do ix=0,lmax
  if (ix == 0) then
    xeff = One
  else
    xeff = gx**ix
  end if
  ax = real(ix+1,kind=wp)
  do iy=0,lmax-ix
    if (iy == 0) then
      xyeff = xeff
    else
      xyeff = xeff*gy**iy
    end if
    ay = real(iy+1,kind=wp)
    do iz=0,lmax-ix-iy
      if (iz == 0) then
        xyzeff = xyeff
      else
        xyzeff = xyeff*gz**iz
      end if
      az = real(iz+1,kind=wp)
      !write(u6,*) ix,iy,iz,Indx(ix,iy,iz),Indx(ix+1,iy,iz),Indx(ix,iy+1,iz),Indx(ix,iy,iz+1)

      ! Charge term

      V = V+xyzeff*Cavxyz(Indx(ix,iy,iz))

      ! Dipole terms

      EFx = EFx+ax*xyzeff*Cavxyz(Indx(ix+1,iy,iz))
      EFy = EFy+ay*xyzeff*Cavxyz(Indx(ix,iy+1,iz))
      EFz = EFz+az*xyzeff*Cavxyz(Indx(ix,iy,iz+1))

    end do
  end do
end do

return

end subroutine hmod
