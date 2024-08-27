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

subroutine pgamma()
!***********************************************************************
!                                                                      *
! Object: to compute the arrays gammath and gammaph for the angular    *
!         integration within a R-matrix run                            *
!                                                                      *
!***********************************************************************

use Constants, only: Pi, Two, Zero
use rmat, only: gammath, gammaph, lgamma

implicit none
integer m, n

! initialize arrays
gammath(:,:) = Zero
gammaph(:,:) = Zero

! Set intitial values for recursion of gammath
gammath(0,0) = Two
gammath(1,0) = Pi/two

m = 0
!m_gam = m
do n=0,2*lgamma+2,2
  gammath(0,n+2) = dble(n+1)/dble(m+n+3)*gammath(0,n)
  !n_gam = n
  !hgt = gammat(x)
  !write(6,*) ' m,n,gammath,hgt',m,n,gammath(m,n),hgt
end do
do n=1,2*lgamma+2,2
  gammath(0,n) = 0.0d0
  !n_gam = n
  !hgt = gammat(x)
  !write(6,*) ' m,n,gammath,hgt',m,n,gammath(m,n),hgt
end do

! m > 0
do m=1,2*lgamma+2
  !m_gam = m
  do n=0,2*lgamma+2,2
    gammath(m,n+2) = dble(n+1)/dble(m+n+3)*gammath(m,n)
    !n_gam = n
    !hgt = gammat(x)
    !write(6,*) ' m,n,gammath,hgt',m,n,gammath(m,n),hgt
  end do
  do n=1,2*lgamma+2,2
    gammath(m,n) = 0.0d0
    !n_gam = n
    !hgt = gammat(x)
    !write(6,*) ' m,n,gammath,hgt',m,n,gammath(m,n),hgt
  end do
  gammath(m+1,0) = dble(m+1)/dble(m+2)*gammath(m-1,0)
end do

! Set intitial values for recursion of gammaph
gammaph(0,0) = Two*pi
gammaph(1,0) = Zero
gammaph(0,1) = Zero

m = 0
!m_gam = m
do n=0,2*lgamma+2
  gammaph(0,n+2) = dble(n+1)/dble(m+n+2)*gammaph(0,n)
  !n_gam = n
  !hgp = gammaf(x)
  !write(6,*) ' m,n,gammaph,hgp',m,n,gammaph(m,n),hgp
end do

! m > 0
do m=1,2*lgamma+2
  !m_gam = m
  do n=0,2*lgamma+2
    gammaph(m,n+2) = dble(n+1)/dble(m+n+2)*gammaph(m,n)
    !n_gam = n
    !hgp = gammaf(x)
    !write(6,*) ' m,n,gammaph,hgp',m,n,gammaph(m,n),hgp
  end do
  gammaph(m+1,0) = dble(m)/dble(m+1)*gammaph(m-1,0)
end do
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine pgamma
