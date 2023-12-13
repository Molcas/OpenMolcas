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

subroutine Box_On_Sphere(x_Min_,x_Max_,y_Min_,y_Max_,z_Min_,z_Max_,xMin_,xMax_,yMin_,yMax_,zMin_,zMax_)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: x_Min_, x_Max_, y_Min_, y_Max_, z_Min_, z_Max_
real(kind=wp), intent(out) :: xMin_, xMax_, yMin_, yMax_, zMin_, zMax_
integer(kind=iwp) :: i, ix, iy, iz, j, ny_Roots, nz_Roots
real(kind=wp) :: r, Roots(3,3), x, x_r, xMax, xMin, xyz(3,2), xyz0(3,2), y, z
real(kind=wp), parameter :: Delta = 1.0e-15_wp

xyz(1,1) = x_min_
xyz(1,2) = x_max_
xyz(2,1) = y_min_
xyz(2,2) = y_max_
xyz(3,1) = z_min_
xyz(3,2) = z_max_
!write(u6,*)
!write(u6,*) 'Box limits'
!write(u6,*) 'x:',xyz(1,1),xyz(1,2)
!write(u6,*) 'y:',xyz(2,1),xyz(2,2)
!write(u6,*) 'z:',xyz(3,1),xyz(3,2)

! Set extremal values

xyz0(1,1) = One
xyz0(1,2) = -One
xyz0(2,1) = One
xyz0(2,2) = -One
xyz0(3,1) = One
xyz0(3,2) = -One

do ix=1,3
  iy = ix+1
  if (iy > 3) iy = 1
  iz = iy+1
  if (iz > 3) iz = 1

  xMax = xyz(ix,2)
  xMin = xyz(ix,1)
  !write(u6,*)
  Roots(1,iy) = xyz(iy,1)
  Roots(2,iy) = xyz(iy,2)
  if (xyz(iy,1)*xyz(iy,2) < Zero) then
    ny_Roots = 3
    Roots(3,iy) = Zero
  else
    ny_Roots = 2
  end if
  Roots(1,iz) = xyz(iz,1)
  Roots(2,iz) = xyz(iz,2)
  if (xyz(iz,1)*xyz(iz,2) < Zero) then
    nz_Roots = 3
    Roots(3,iz) = Zero
  else
    nz_Roots = 2
  end if
  !call RecPrt('Roots','(3G25.12)',Roots,3,3)

  do i=1,ny_Roots
    !write(u6,*) 'i=',i,ny_Roots
    !write(u6,*)
    y = Roots(i,iy)
    do j=1,nz_Roots
      !write(u6,*) 'j=',j,nz_Roots
      z = Roots(j,iz)

      x = xMin
      r = sqrt(x**2+y**2+z**2)
      !write(u6,*) x/r
      if (r == Zero) then
        x_r = Zero
      else
        x_r = x/r
      end if
      xyz0(ix,1) = min(xyz0(ix,1),x_r)
      xyz0(ix,2) = max(xyz0(ix,2),x_r)

      x = xMax
      r = sqrt(x**2+y**2+z**2)
      !write(u6,*) x/r
      if (r == Zero) then
        x_r = Zero
      else
        x_r = x/r
      end if
      xyz0(ix,1) = min(xyz0(ix,1),x_r)
      xyz0(ix,2) = max(xyz0(ix,2),x_r)

    end do
  end do
end do

!write(u6,*) 'xMin=',xyz0(1,1)
!write(u6,*) 'xMax=',xyz0(1,2)
!write(u6,*) 'yMin=',xyz0(2,1)
!write(u6,*) 'yMax=',xyz0(2,2)
!write(u6,*) 'zMin=',xyz0(3,1)
!write(u6,*) 'zMax=',xyz0(3,2)
xMin_ = xyz0(1,1)-Delta
xMax_ = xyz0(1,2)+Delta
yMin_ = xyz0(2,1)-Delta
yMax_ = xyz0(2,2)+Delta
zMin_ = xyz0(3,1)-Delta
zMax_ = xyz0(3,2)+Delta

return

end subroutine Box_On_Sphere
