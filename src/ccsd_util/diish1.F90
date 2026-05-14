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

subroutine diish1(wrk,wrksize,nind,rdiis1,v1,v2,v3,v4,ndiis,szkey)
! this routine upgrades rdiis1(p,q) matrix
! rdiis1(p,q) = szkey*rdiis1(p,q) + (xp|xq)
!
! nind   - number of indices in vectors (i)
! rdiis1 - overlap matrix of 1-4 vectors (O)
! v1     - map type of vector 1 (i)
! v2     - map type of vector 1 (i)
! v3     - map type of vector 1 (i)
! v4     - map type of vector 1 (i)
! if there is less than 4 vectors, use any map
! ndiis  - dimension of DIIS (2-4) (I)
! szkey  - 0 - no vanishing rdiis1 at the beginning
!          1 - vanishing rdiis1 at the beginning

use ccsd_global, only: Map_Type
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, nind, ndiis, szkey
real(kind=wp), intent(in) :: wrk(wrksize)
real(kind=wp), intent(inout) :: rdiis1(4,4)
type(Map_Type), intent(in) :: v1, v2, v3, v4
integer(kind=iwp) :: num, rc
real(kind=wp) :: scalar

num = ndiis+1

if (szkey == 1) rdiis1(:,:) = Zero

if (num > 0) then
  ! calc X11
  call multdot(wrk,wrksize,nind,v1,1,v1,1,scalar,rc)
  rdiis1(1,1) = rdiis1(1,1)+scalar
end if

if (num > 1) then
  ! X21
  call multdot(wrk,wrksize,nind,v2,1,v1,1,scalar,rc)
  rdiis1(2,1) = rdiis1(2,1)+scalar
  rdiis1(1,2) = rdiis1(1,2)+scalar
  ! X22
  call multdot(wrk,wrksize,nind,v2,1,v2,1,scalar,rc)
  rdiis1(2,2) = rdiis1(2,2)+scalar
end if

if (num > 2) then
  ! X31
  call multdot(wrk,wrksize,nind,v3,1,v1,1,scalar,rc)
  rdiis1(3,1) = rdiis1(3,1)+scalar
  rdiis1(1,3) = rdiis1(1,3)+scalar
  ! X32
  call multdot(wrk,wrksize,nind,v3,1,v2,1,scalar,rc)
  rdiis1(3,2) = rdiis1(3,2)+scalar
  rdiis1(2,3) = rdiis1(2,3)+scalar
  ! X33
  call multdot(wrk,wrksize,nind,v3,1,v3,1,scalar,rc)
  rdiis1(3,3) = rdiis1(3,3)+scalar
end if

if (num > 3) then
  ! X41
  call multdot(wrk,wrksize,nind,v4,1,v1,1,scalar,rc)
  rdiis1(4,1) = rdiis1(4,1)+scalar
  rdiis1(1,4) = rdiis1(1,4)+scalar
  ! X42
  call multdot(wrk,wrksize,nind,v4,1,v2,1,scalar,rc)
  rdiis1(4,2) = rdiis1(4,2)+scalar
  rdiis1(2,4) = rdiis1(2,4)+scalar
  ! X43
  call multdot(wrk,wrksize,nind,v4,1,v3,1,scalar,rc)
  rdiis1(4,3) = rdiis1(4,3)+scalar
  rdiis1(3,4) = rdiis1(3,4)+scalar
  ! X44
  call multdot(wrk,wrksize,nind,v4,1,v4,1,scalar,rc)
  rdiis1(4,4) = rdiis1(4,4)+scalar
end if

return

end subroutine diish1
