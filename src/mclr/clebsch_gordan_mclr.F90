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
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************

function Clebsch_Gordan_mclr(j1,m1,j2,m2,j,m)
! Calculates the Clebsch-Gordan coefficients

use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: Clebsch_Gordan_mclr
real(kind=wp), intent(in) :: j1, m1, j2, m2, j, m
integer(kind=iwp) :: z, zmax, zmin
real(kind=wp) :: dz, Fct, Fct1, Fct2, N, r, rsum, T
real(kind=wp), external :: Fact

if ((j1 < Zero) .or. (j2 < Zero) .or. (j < Zero)) then
  write(u6,*) 'Error J is lower than 0'
  call Abend()
end if
r = abs(Two*j1-int(Two*j1))+abs(Two*j2-int(Two*j2))+abs(Two*j-int(Two*j))+abs(Two*m1-int(Two*m1))+abs(Two*m2-int(Two*m2))+ &
    abs(Two*m-int(Two*m))
if (r > 1.0e-6_wp) then
  write(u6,*) 'CG provided with not half integer'
  call Abend()
end if
if (m1+m2 == m) then
  Fct1 = (Two*j+One)*Fact(j1+j2-j)*Fact(j1-j2+j)*Fact(-j1+j2+j)
  Fct2 = Fact(j1+j2+j+One)
  Fct = sqrt(Fct1/fct2)*sqrt(Fact(j1+m1)*Fact(j1-m1)*Fact(j2+m2)*Fact(j2-m2)*Fact(j+m)*Fact(j-m))
  zmax = nint(min(j1+j2-j,j1-m1,j2+m2))
  zmin = -nint(min(j-j2+m1,j-j1-m2))
  rsum = Zero
  zmin = max(0,zmin)
  do z=zmin,zmax
    dz = real(z,kind=wp)
    T = (-One)**z
    N = Fact(dz)*Fact(j1+j2-j-dz)*Fact(j1-m1-dz)*Fact(j2+m2-dz)*Fact(j-j2+m1+dz)*Fact(j-j1-m2+dz)
    rsum = rsum+T/N
  end do
  Clebsch_Gordan_mclr = rsum*Fct
else
  Clebsch_Gordan_mclr = Zero
end if

return

end function Clebsch_Gordan_mclr
