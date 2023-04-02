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

subroutine diish2(rdiis1,ndiis,cdiis)
! this rouine calculates diis coefficients by solving
!
! B -1   c  =  0
! -1 0   l    -1
!
! and consequent normalization of coefficients
!
! r1diis  - matrix of amp. overlap of ndiid+1 iterations (I)
! ndiis   - size of diis (2-4) (I)
! cdiis   - final diis coefficients (O)

real*8 rdiis1(1:4,1:4)
real*8 cdiis(1:4)
integer ndiis
! help variables
integer p, q
real*8 scalar
real*8 rdiis2(1:5,1:5)
real*8 bb(1:5)
real*8 ci(1:5)

!1 vanish rdiis2 file
call mv0zero(25,25,rdiis2)

!2.1 make rdiis2 matrix

do p=1,ndiis
  do q=1,ndiis
    rdiis2(p,q) = rdiis1(p,q)
  end do
end do

do p=1,ndiis
  rdiis2(p,ndiis+1) = -1.0d0
  rdiis2(ndiis+1,p) = -1.0d0
  bb(p) = 0.0d0
end do

bb(ndiis+1) = -1.0d0

!2.2 modify rdiis2 matrix

! scale matrix
scalar = sqrt(rdiis2(1,1)*rdiis2(ndiis,ndiis))
do q=1,ndiis
  do p=1,ndiis
    rdiis2(p,q) = rdiis2(p,q)/scalar
  end do
end do

! add penalty function
scalar = 0.01d0*rdiis2(ndiis,ndiis)
!rdiis2(ndiis,ndiis) = rdiis2(ndiis,ndiis)+scalar
!bb(ndiis) = scalar

!3 solve SLE

!3.1 vanish ci
do p=1,ndiis+1
  ci(p) = 0.0d0
end do

!3.2 get cdiis
call gauss(ndiis+1,5,rdiis2,ci,bb)

!4 final modification of cdiis coefficients

!FUE if (rc == 1) then
! matrix R2 was singular, no extrapolation
!FUE write(6,*) ' SINGULAR DIIS MATRIX, NO EXTRAPOLATION'
!FUE cdiis(1) = 1.0d0
!FUE do p=2,ndiis
!FUE   cdiis(p) = 0.0d0
!FUE end do

!FUE else
! DIIS procedure was successful, renormalize coef.

scalar = 0.0d0
do p=1,ndiis
  scalar = scalar+ci(p)
end do

do p=1,ndiis
  cdiis(p) = ci(p)/scalar
end do

scalar = 0.0d0
do p=1,ndiis
  scalar = scalar+cdiis(p)
end do

!FUE end if

!FUE write(6,*) cdiis(1),cdiis(2),cdiis(3),cdiis(4),scalar
!51 format (5(i2,d12.7))

return

end subroutine diish2
