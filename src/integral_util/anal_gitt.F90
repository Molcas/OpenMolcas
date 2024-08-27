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

function Anal_Gitt(cordsi,latato)
!***********************************************************************
!                                                                      *
!     Object:                                                          *
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

use Constants, only: Zero, One

implicit none
real*8 Anal_Gitt
integer latato
real*8 cordsi(3,latato)
real*8 gAtom
integer i, j
real*8 faktor, x1, y1, z1, r2

Anal_Gitt = Zero

! analyze the lattice

gatom = Zero
do i=1,latato
  faktor = One
  x1 = cordsi(1,i)+One
  y1 = cordsi(2,i)
  z1 = cordsi(3,i)
  do j=1,latato
    r2 = (x1-cordsi(1,j))**2+(y1-cordsi(2,j))**2+(z1-cordsi(3,j))**2
    if (r2 < 0.01d0) faktor = faktor+One
  end do

  x1 = cordsi(1,i)-One
  y1 = cordsi(2,i)
  z1 = cordsi(3,i)
  do j=1,latato
    r2 = (x1-cordsi(1,j))**2+(y1-cordsi(2,j))**2+(z1-cordsi(3,j))**2
    if (r2 < 0.01d0) faktor = faktor+One
  end do

  x1 = cordsi(1,i)
  y1 = cordsi(2,i)+One
  z1 = cordsi(3,i)
  do j=1,latato
    r2 = (x1-cordsi(1,j))**2+(y1-cordsi(2,j))**2+(z1-cordsi(3,j))**2
    if (r2 < 0.01d0) faktor = faktor+One
  end do

  x1 = cordsi(1,i)
  y1 = cordsi(2,i)-One
  z1 = cordsi(3,i)
  do j=1,latato
    r2 = (x1-cordsi(1,j))**2+(y1-cordsi(2,j))**2+(z1-cordsi(3,j))**2
    if (r2 < 0.01d0) faktor = faktor+One
  end do

  x1 = cordsi(1,i)
  y1 = cordsi(2,i)
  z1 = cordsi(3,i)+One
  do j=1,latato
    r2 = (x1-cordsi(1,j))**2+(y1-cordsi(2,j))**2+(z1-cordsi(3,j))**2
    if (r2 < 0.01d0) faktor = faktor+One
  end do

  x1 = cordsi(1,i)
  y1 = cordsi(2,i)
  z1 = cordsi(3,i)-One
  do j=1,latato
    r2 = (x1-cordsi(1,j))**2+(y1-cordsi(2,j))**2+(z1-cordsi(3,j))**2
    if (r2 < 0.01d0) faktor = faktor+One
  end do
  gatom = gatom+One/faktor
end do
Anal_Gitt = gatom

return

end function Anal_Gitt
