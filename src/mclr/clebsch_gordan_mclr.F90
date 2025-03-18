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

implicit real*8(a-z)
integer z, zmax, zmin

if ((j1 < 0.0d0) .or. (j2 < 0.0d0) .or. (j < 0.0d0)) then
  write(6,*) 'Error J is lower than 0'
  call Abend()
end if
r = abs(2.0d0*j1-dint(2.0d0*j1))+abs(2.0d0*j2-dint(2.0d0*j2))+abs(2.0d0*j-dint(2.0d0*j))+abs(2.0d0*m1-dint(2.0d0*m1))+ &
    abs(2.0d0*m2-dint(2.0d0*m2))+abs(2.0d0*m-dint(2.0d0*m))
if (r > 1.0d-6) then
  write(6,*) 'CG provided with not half integer'
  call Abend()
end if
if (m1+m2 == m) then
  Fct1 = (2.0d0*j+1.0d0)*Fact(j1+j2-j)*Fact(j1-j2+j)*Fact(-j1+j2+j)
  Fct2 = Fact(j1+j2+j+1.0d0)
  Fct = sqrt(Fct1/fct2)
  Fct = Fct*sqrt(Fact(j1+m1)*Fact(j1-m1)*Fact(j2+m2)*Fact(j2-m2)*Fact(j+m)*FacT(j-m))
  zmax = nint(min(j1+j2-j,j1-m1,j2+m2))
  zmin = -nint(min(j-j2+m1,j-j1-m2))
  sum = 0.0d0
  zmin = max(0,zmin)
  do z=zmin,zmax
    T = dble((-1)**z)
    N = Facti(z)*FacT(j1+j2-j-dble(z))*fact(j1-m1-dble(z))*fact(j2+m2-dble(z))*fact(j-j2+m1+dble(z))*fact(j-j1-m2+dble(z))
    sum = sum+T/N
  end do
  Clebsch_Gordan_mclr = sum*fct
else
  Clebsch_Gordan_mclr = 0.0d0
end if

return

end function Clebsch_Gordan_mclr
