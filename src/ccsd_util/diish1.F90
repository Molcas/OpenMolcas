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

subroutine diish1(wrk,wrksize,nind,rdiis1,mapd1,mapd2,mapd3,mapd4,mapi1,mapi2,mapi3,mapi4,ndiis,szkey)
! this routine upgrades rdiis1(p,q) matrix
! rdiis1(p,q) = szkey*rdiis1(p,q) + (xp|xq)
!
! nind   - number of indices in vectors (i)
! rdiis1 - overlap matrix of 1-5 vectors (O)
! mapd1  - direct map of vector 1 (i)
! mapd2  - direct map of vector 2 (i)
! mapd3  - direct map of vector 3 (i)
! mapd4  - direct map of vector 4 (i)
! mapi1  - inverse map of vector 1 (i)
! mapi2  - inverse map of vector 2 (i)
! mapi3  - inverse map of vector 3 (i)
! mapi4  - inverse map of vector 4 (i)
! if there is less than 5 vectors, use any map
! ndiis  - dimension of DIIS (2-4) (I)
! szkey  - 0 - no vanishing rdiis1 at the beginning
!          1 - vanishing rdiis1 at the beginning

#include "wrk.fh"
#include "ccsd1.fh"
real*8 rdiis1(1:4,1:4)
integer mapd1(0:512,1:6)
integer mapd2(0:512,1:6)
integer mapd3(0:512,1:6)
integer mapd4(0:512,1:6)
integer mapi1(1:8,1:8,1:8)
integer mapi2(1:8,1:8,1:8)
integer mapi3(1:8,1:8,1:8)
integer mapi4(1:8,1:8,1:8)
integer nind, ndiis, szkey
! help variables
integer nhelp
real*8 scalar
integer rc, num

num = ndiis+1

if (szkey == 1) then
  nhelp = 4*4
  call mv0zero(nhelp,nhelp,rdiis1)
end if

if (num > 0) then
  ! calc X11
  call multdot(wrk,wrksize,nind,mapd1,mapi1,1,mapd1,mapi1,1,scalar,rc)
  rdiis1(1,1) = rdiis1(1,1)+scalar
end if

if (num > 1) then
  ! X21
  call multdot(wrk,wrksize,nind,mapd2,mapi2,1,mapd1,mapi1,1,scalar,rc)
  rdiis1(2,1) = rdiis1(2,1)+scalar
  rdiis1(1,2) = rdiis1(1,2)+scalar
  ! X22
  call multdot(wrk,wrksize,nind,mapd2,mapi2,1,mapd2,mapi2,1,scalar,rc)
  rdiis1(2,2) = rdiis1(2,2)+scalar
end if

if (num > 2) then
  ! X31
  call multdot(wrk,wrksize,nind,mapd3,mapi3,1,mapd1,mapi1,1,scalar,rc)
  rdiis1(3,1) = rdiis1(3,1)+scalar
  rdiis1(1,3) = rdiis1(1,3)+scalar
  ! X32
  call multdot(wrk,wrksize,nind,mapd3,mapi3,1,mapd2,mapi2,1,scalar,rc)
  rdiis1(3,2) = rdiis1(3,2)+scalar
  rdiis1(2,3) = rdiis1(2,3)+scalar
  ! X33
  call multdot(wrk,wrksize,nind,mapd3,mapi3,1,mapd3,mapi3,1,scalar,rc)
  rdiis1(3,3) = rdiis1(3,3)+scalar
end if

if (num > 3) then
  ! X41
  call multdot(wrk,wrksize,nind,mapd4,mapi4,1,mapd1,mapi1,1,scalar,rc)
  rdiis1(4,1) = rdiis1(4,1)+scalar
  rdiis1(1,4) = rdiis1(1,4)+scalar
  ! X42
  call multdot(wrk,wrksize,nind,mapd4,mapi4,1,mapd2,mapi2,1,scalar,rc)
  rdiis1(4,2) = rdiis1(4,2)+scalar
  rdiis1(2,4) = rdiis1(2,4)+scalar
  ! X43
  call multdot(wrk,wrksize,nind,mapd4,mapi4,1,mapd3,mapi3,1,scalar,rc)
  rdiis1(4,3) = rdiis1(4,3)+scalar
  rdiis1(3,4) = rdiis1(3,4)+scalar
  ! X44
  call multdot(wrk,wrksize,nind,mapd4,mapi4,1,mapd4,mapi4,1,scalar,rc)
  rdiis1(4,4) = rdiis1(4,4)+scalar
end if

return

end subroutine diish1
