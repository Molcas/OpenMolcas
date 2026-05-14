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

subroutine Picky_Mck(nSD,iSD4,i,j,nMethod)

use Symmetry_Info, only: nIrrep
use Dens_Stuff, only: ipDDij, ipDDij2, ipDDik, ipDDik2, ipDDil, ipDDil2, ipDDjk, ipDDjk2, ipDDjl, ipDDjl2, ipDDkl, ipDDkl2, ipDij, &
                      ipDij2, ipDik, ipDik2, ipDil, ipDil2, ipDjk, ipDjk2, ipDjl, ipDjl2, ipDkl, ipDkl2, mDCRij, mDCRik, mDCRil, &
                      mDCRjk, mDCRjl, mDCRkl, mDij, mDik, mDil, mDjk, mDjl, mDkl
use k2_arrays, only: DeDe, DeDe2
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSD, iSD4(0:nSD,4), i, j, nMethod
integer(kind=iwp) :: i1, i2, i3, iBasAO, iBasi, iBasn, iCmpi, ii1, ii2, ii3, iPrimi, iShell, j1, j2, j3, jBasAO, jBasj, jBasn, &
                     jCmpj, jj1, jj2, jj3, jPrimj, jShell
integer(kind=iwp), pointer :: ipD2_ij, ipD_ij, ipDD2_ij, ipDD_ij, mD_ij, mDCR_ij
integer(kind=iwp), parameter :: SCF = 1, RASSCF = 2

iCmpi = iSD4(2,i)
iBasi = iSD4(3,i)
iPrimi = iSD4(5,i)
iShell = iSD4(11,i)
iBasAO = iSD4(8,i)+1
iBasn = iSD4(19,i)

jCmpj = iSD4(2,j)
jBasj = iSD4(3,j)
jPrimj = iSD4(5,j)
jShell = iSD4(11,j)
jBasAO = iSD4(8,j)+1
jBasn = iSD4(19,j)

if ((i == 1) .and. (j == 2)) then
  mDCR_ij => mDCRij
  ipD_ij => ipDij
  ipDD_ij => ipDDij
  mD_ij => mDij
  ipD2_ij => ipDij2
  ipDD2_ij => ipDDij2
else if ((i == 1) .and. (j == 3)) then
  mDCR_ij => mDCRik
  ipD_ij => ipDik
  ipDD_ij => ipDDik
  mD_ij => mDik
  ipD2_ij => ipDik2
  ipDD2_ij => ipDDik2
else if ((i == 1) .and. (j == 4)) then
  mDCR_ij => mDCRil
  ipD_ij => ipDil
  ipDD_ij => ipDDil
  mD_ij => mDil
  ipD2_ij => ipDil2
  ipDD2_ij => ipDDil2
else if ((i == 2) .and. (j == 3)) then
  mDCR_ij => mDCRjk
  ipD_ij => ipDjk
  ipDD_ij => ipDDjk
  mD_ij => mDjk
  ipD2_ij => ipDjk2
  ipDD2_ij => ipDDjk2
else if ((i == 2) .and. (j == 4)) then
  mDCR_ij => mDCRjl
  ipD_ij => ipDjl
  ipDD_ij => ipDDjl
  mD_ij => mDjl
  ipD2_ij => ipDjl2
  ipDD2_ij => ipDDjl2
else if ((i == 3) .and. (j == 4)) then
  mDCR_ij => mDCRkl
  ipD_ij => ipDkl
  ipDD_ij => ipDDkl
  mD_ij => mDkl
  ipD2_ij => ipDkl2
  ipDD2_ij => ipDDkl2
else
  write(u6,*) 'Picky: illegal i and j combination'
  write(u6,*) 'i,j=',i,j
  call Abend()
  nullify(mDCR_ij)
  nullify(ipD_ij)
  nullify(ipDD_ij)
  nullify(mD_ij)
  nullify(ipD2_ij)
  nullify(ipDD2_ij)
end if

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
if (mDCR_ij /= 0) then
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
  call Picky_inner(DeDe(ipD_ij),i1,j1,iPrimi*jPrimj,iCmpi*jCmpj,mDCR_ij,i2,i2+i3-1,j2,j2+j3-1,DeDe(ipDD_ij))
  if (nMethod == RASSCF) &
    call Picky_inner(DeDe2(ipD2_ij),i1,j1,iPrimi*jPrimj,iCmpi*jCmpj,mDCR_ij,i2,i2+i3-1,j2,j2+j3-1,DeDe2(ipDD2_ij))
end if
mD_ij = (ii3*jj3+1)*iCmpi*jCmpj+iPrimi*jPrimj+1

end subroutine Picky_Mck
