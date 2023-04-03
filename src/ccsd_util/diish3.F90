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

subroutine diish3(wrk,wrksize,mapd0,mapd1,mapd2,mapd3,mapd4,cdiis,ndiis)
! this routine produces new vector
! v(0) = sum(i) [cdiis(i) . V(i)]
!
! mapd0 - direct map of V0 vector (I)
! mapd1 - direct map of V1 vector (I)
! mapd2 - direct map of V2 vector (I)
! mapd3 - direct map of V3 vector (I)
! mapd4 - direct map of V4 vector (I)
! cdiis - vector of diis coefficients (I)
! ndiis - size of diis (2-4) (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize, mapd0(0:512,6), mapd1(0:512,6), mapd2(0:512,6), mapd3(0:512,6), mapd4(0:512,6), ndiis
real(kind=wp) :: wrk(wrksize), cdiis(4)
integer(kind=iwp) :: length, nhelp, poss0, poss1, poss2, poss3, poss4

if (ndiis == 2) then
  ! 2 dimensional DIIS

  poss0 = mapd0(1,1)
  poss1 = mapd1(1,1)
  poss2 = mapd2(1,1)

  nhelp = mapd1(0,5)
  length = mapd1(nhelp,1)+mapd1(nhelp,2)-mapd1(1,1)

  if (length > 0) then
    do nhelp=0,length-1
      wrk(poss0+nhelp) = cdiis(1)*wrk(poss1+nhelp)+cdiis(2)*wrk(poss2+nhelp)
    end do
  end if

else if (ndiis == 3) then
  ! 3 dimensional DIIS

  poss0 = mapd0(1,1)
  poss1 = mapd1(1,1)
  poss2 = mapd2(1,1)
  poss3 = mapd3(1,1)

  nhelp = mapd1(0,5)
  length = mapd1(nhelp,1)+mapd1(nhelp,2)-mapd1(1,1)

  if (length > 0) then
    do nhelp=0,length-1
      wrk(poss0+nhelp) = cdiis(1)*wrk(poss1+nhelp)+cdiis(2)*wrk(poss2+nhelp)+cdiis(3)*wrk(poss3+nhelp)
    end do
  end if

else if (ndiis == 4) then
  ! 4 dimensional DIIS

  poss0 = mapd0(1,1)
  poss1 = mapd1(1,1)
  poss2 = mapd2(1,1)
  poss3 = mapd3(1,1)
  poss4 = mapd4(1,1)

  nhelp = mapd1(0,5)
  length = mapd1(nhelp,1)+mapd1(nhelp,2)-mapd1(1,1)

  if (length > 0) then
    do nhelp=0,length-1
      wrk(poss0+nhelp) = cdiis(1)*wrk(poss1+nhelp)+cdiis(2)*wrk(poss2+nhelp)+cdiis(3)*wrk(poss3+nhelp)+cdiis(4)*wrk(poss4+nhelp)
    end do
  end if

end if

return

end subroutine diish3
