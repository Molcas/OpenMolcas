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

real*8 function Clebsch_Gordan_mclr(j1,m1,j2,m2,j,m)
! Calculates the Clebsch-Gordan coefficients

use Constants, only: Zero, One, Two
use Definitions, only: wp, u6

implicit real*8(a-z)
integer z, zmax, zmin

if ((j1 < Zero) .or. (j2 < Zero) .or. (j < Zero)) then
  write(u6,*) 'Error J is lower than 0'
  call Abend()
end if
r = abs(Two*j1-aint(Two*j1))+abs(Two*j2-aint(Two*j2))+abs(Two*j-aint(Two*j))+abs(Two*m1-aint(Two*m1))+abs(Two*m2-aint(Two*m2))+ &
    abs(Two*m-aint(Two*m))
if (r > 1.0e-6_wp) then
  write(u6,*) 'CG provided with not half integer'
  call Abend()
end if
if (m1+m2 == m) then
  Fct1 = (Two*j+One)*Fact(j1+j2-j)*Fact(j1-j2+j)*Fact(-j1+j2+j)
  Fct2 = Fact(j1+j2+j+One)
  Fct = sqrt(Fct1/fct2)
  Fct = Fct*sqrt(Fact(j1+m1)*Fact(j1-m1)*Fact(j2+m2)*Fact(j2-m2)*Fact(j+m)*FacT(j-m))
  zmax = nint(min(j1+j2-j,j1-m1,j2+m2))
  zmin = -nint(min(j-j2+m1,j-j1-m2))
  sum = Zero
  zmin = max(0,zmin)
  do z=zmin,zmax
    dz = real(z,kind=wp)
    T = (-One)**z
    N = Facti(z)*FacT(j1+j2-j-dz)*fact(j1-m1-dz)*fact(j2+m2-dz)*fact(j-j2+m1+dz)*fact(j-j1-m2+dz)
    sum = sum+T/N
  end do
  Clebsch_Gordan_mclr = sum*fct
else
  Clebsch_Gordan_mclr = Zero
end if

return

end function Clebsch_Gordan_mclr
