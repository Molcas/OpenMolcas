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

subroutine diish3(wrk,wrksize,v0,v1,v2,v3,v4,cdiis,ndiis)
! this routine produces new vector
! v(0) = sum(i) [cdiis(i) . V(i)]
!
! v0    - map type of V0 vector (I)
! v1    - map type of V1 vector (I)
! v2    - map type of V2 vector (I)
! v3    - map type of V3 vector (I)
! v4    - map type of V4 vector (I)
! cdiis - vector of diis coefficients (I)
! ndiis - size of diis (2-4) (I)

use ccsd_global, only: Map_Type
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, ndiis
real(kind=wp), intent(inout) :: wrk(wrksize)
type(Map_Type), intent(in) :: v0, v1, v2, v3, v4
real(kind=wp), intent(in) :: cdiis(4)
integer(kind=iwp) :: length, ll, nhelp, pos0, pos1, pos2, pos3, pos4

if (ndiis == 2) then
  ! 2 dimensional DIIS

  pos0 = v0%d(1,1)
  pos1 = v1%d(1,1)
  pos2 = v2%d(1,1)

  nhelp = v1%d(0,5)
  length = v1%d(nhelp,1)+v1%d(nhelp,2)-v1%d(1,1)

  if (length > 0) then
    ll = length-1
    wrk(pos0:pos0+ll) = cdiis(1)*wrk(pos1:pos1+ll)+cdiis(2)*wrk(pos2:pos2+ll)
  end if

else if (ndiis == 3) then
  ! 3 dimensional DIIS

  pos0 = v0%d(1,1)
  pos1 = v1%d(1,1)
  pos2 = v2%d(1,1)
  pos3 = v3%d(1,1)

  nhelp = v1%d(0,5)
  length = v1%d(nhelp,1)+v1%d(nhelp,2)-v1%d(1,1)

  if (length > 0) then
    ll = length-1
    wrk(pos0:pos0+ll) = cdiis(1)*wrk(pos1:pos1+ll)+cdiis(2)*wrk(pos2:pos2+ll)+cdiis(3)*wrk(pos3:pos3+ll)
  end if

else if (ndiis == 4) then
  ! 4 dimensional DIIS

  pos0 = v0%d(1,1)
  pos1 = v1%d(1,1)
  pos2 = v2%d(1,1)
  pos3 = v3%d(1,1)
  pos4 = v4%d(1,1)

  nhelp = v1%d(0,5)
  length = v1%d(nhelp,1)+v1%d(nhelp,2)-v1%d(1,1)

  if (length > 0) then
    ll = length-1
    wrk(pos0:pos0+ll) = cdiis(1)*wrk(pos1:pos1+ll)+cdiis(2)*wrk(pos2:pos2+ll)+cdiis(3)*wrk(pos3:pos3+ll)+cdiis(4)*wrk(pos4:pos4+ll)
  end if

end if

return

end subroutine diish3
