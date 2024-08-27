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

implicit none
integer lmax_
real*8 gx, gy, gz, V, EFx, EFy, EFz
real*8 Cavxyz((lMax_+1)*(lMax_+2)*(lMax_+3)/6)
integer ix, iy, iz, iOff, Index, lMax
real*8 xEff, xyEff, xyzEff, aX, aY, aZ
!---- Statement functions
iOff(ix,iy,iz) = (ix+iy+iz)*(ix+iy+iz+1)*(ix+iy+iz+2)/6
index(ix,iy,iz) = iOff(ix,iy,iz)+(iy+iz)*(iy+iz+1)/2+iz+1

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
  ax = dble(ix)+One
  do iy=0,lmax-ix
    if (iy == 0) then
      xyeff = xeff
    else
      xyeff = xeff*gy**iy
    end if
    ay = dble(iy)+One
    do iz=0,lmax-ix-iy
      if (iz == 0) then
        xyzeff = xyeff
      else
        xyzeff = xyeff*gz**iz
      end if
      az = dble(iz)+One
      !write(6,*) ix,iy,iz,Index(ix,iy,iz),Index(ix+1,iy,iz),Index(ix,iy+1,iz),Index(ix,iy,iz+1)

      ! Charge term

      V = V+xyzeff*Cavxyz(index(ix,iy,iz))

      ! Dipole terms

      EFx = EFx+ax*xyzeff*Cavxyz(index(ix+1,iy,iz))
      EFy = EFy+ay*xyzeff*Cavxyz(index(ix,iy+1,iz))
      EFz = EFz+az*xyzeff*Cavxyz(index(ix,iy,iz+1))

    end do
  end do
end do

return

end subroutine hmod
