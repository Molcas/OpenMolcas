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
! Copyright (C) 2000, Gunnar Karlstrom                                 *
!               2000, Roland Lindh                                     *
!***********************************************************************

subroutine qlm(gx,gy,gz,qa,dax,day,daz,lmax_,Cavxyz)
!***********************************************************************
!                                                                      *
!     Object: to reexpand the charge and the dipole moment at a given  *
!             point as a multipole moment expansion at origin.         *
!                                                                      *
!     Authors: G. Karlstroem                                           *
!              Dept. of Theor. Chem., Univ. of Lund, Sweden.           *
!                                                                      *
!              and                                                     *
!                                                                      *
!              R. Lindh                                                *
!              Dept. of Chem. Phys., Univ. of Lund, Sweden.            *
!                                                                      *
!              March 2000                                              *
!***********************************************************************

use Index_Functions, only: C3_Ind3, nTri3_Elem1
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: gx, gy, gz, qa, dax, day, daz
integer(kind=iwp), intent(in) :: lmax_
real(kind=wp), intent(inout) :: Cavxyz(nTri3_Elem1(lMax_))
integer(kind=iwp) :: ix, iy, iz, lMax
real(kind=wp) :: ax, ay, az, xeff, xyeff, xyzeff

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
        xyzeff = xyeff*gz**iz
      else
        xyzeff = xyeff*gz**iz
      end if
      az = real(iz+1,kind=wp)

      ! Charge term

      Cavxyz(C3_Ind3(ix,iy,iz)) = Cavxyz(C3_Ind3(ix,iy,iz))+xyzeff*qa

      ! Dipole terms

      Cavxyz(C3_Ind3(ix+1,iy,iz)) = Cavxyz(C3_Ind3(ix+1,iy,iz))+xyzeff*dax*ax
      Cavxyz(C3_Ind3(ix,iy+1,iz)) = Cavxyz(C3_Ind3(ix,iy+1,iz))+xyzeff*day*ay
      Cavxyz(C3_Ind3(ix,iy,iz+1)) = Cavxyz(C3_Ind3(ix,iy,iz+1))+xyzeff*daz*az

    end do
  end do
end do

return

end subroutine qlm
