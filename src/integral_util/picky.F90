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

subroutine Picky(iBasi,iBsInc,iPrimi,iBasAO,iBasn,jBasj,jBsInc,jPrimj,jBasAO,jBasn,iCmpi,jCmpj,iShell,jShell,mDCRij,ipDij,ipDDij, &
                 mDij,DeDe,nDeDe)

use Symmetry_Info, only: nIrrep
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iBasi, iBsInc, iPrimi, iBasAO, iBasn, jBasj, jBsInc, jPrimj, jBasAO, jBasn, iCmpi, jCmpj, iShell, &
                                 jShell, mDCRij, ipDij, nDeDe
integer(kind=iwp), intent(inout) :: ipDDij
integer(kind=iwp), intent(out) :: mDij
real(kind=wp), intent(inout) :: DeDe(nDeDe)
integer(kind=iwp) :: ii1, ii2, ii3, jj1, jj2, jj3, i1, i2, i3, j1, j2, j3

if (nIrrep == 1) then
  ii1 = 0
  ii2 = 1
  ii3 = 0
  jj1 = 0
  jj2 = 1
  jj3 = 0
else
  ii1 = iBasi
  ii2 = iBasAO
  ii3 = iBasn
  jj1 = jBasj
  jj2 = jBasAO
  jj3 = jBasn
end if
if (mDCRij /= 0) then
  if (iShell >= jShell) then
    i1 = ii1
    i2 = ii2
    i3 = ii3
    j1 = jj1
    j2 = jj2
    j3 = jj3
  else
    i1 = jj1
    i2 = jj2
    i3 = jj3
    j1 = ii1
    j2 = ii2
    j3 = ii3
  end if
  if ((iBasi == iBsInc) .and. (jBasj == jBsInc)) then
    ipDDij = ipDij
  else
    call Picky_inner(DeDe(ipDij),i1,j1,iPrimi*jPrimj,iCmpi*jCmpj,mDCRij,i2,i2+i3-1,j2,j2+j3-1,DeDe(ipDDij))
  end if
end if
mDij = (ii3*jj3+1)*iCmpi*jCmpj+iPrimi*jPrimj+1

return

end subroutine Picky
