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

subroutine Gen_GeoList()

use GeoList, only: Centr, Chrg, Mass
use Basis_Info, only: dbsc, nCnttp
use Center_Info, only: dc
use Symmetry_Info, only: nIrrep, iChCar
use Sizes_of_Seward, only: S
use Gateway_Info, only: TMass, qNuc, CoM, CoC
use stdalloc, only: mma_allocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: i, jCnt, jCnttp, mCnt, nc, nchr, ndc
real(kind=wp) :: Z

!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(Centr,3,S%mCentr,label='Centr')
call mma_allocate(Mass,S%mCentr,label='Mass')
call mma_allocate(Chrg,S%mCentr,label='Chrg')
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate the center list.

S%kCentr = 0

nc = 1
do jCnttp=1,nCnttp
  mCnt = dbsc(jCnttp)%nCntr

  ! Do not include Auxiliary basis sets, or fragment basis sets

  if (dbsc(jCnttp)%Aux .or. dbsc(jCnttp)%Frag) cycle

  ! Do not include ECP basis sets which does not have any valence basis set.

  if (dbsc(jCnttp)%ECP .and. (dbsc(jCnttp)%nVal == 0)) cycle

  do jCnt=1,mCnt
    ndc = jCnt+dbsc(jCnttp)%mdci
    do i=0,nIrrep/dc(ndc)%nStab-1
      call OA(dc(ndc)%iCoSet(i,0),dbsc(jCnttp)%Coor(1:3,jCnt),Centr(1:3,nc))
      nchr = dbsc(jCnttp)%AtmNr
      if (nchr >= 0) then
        Mass(nc) = dbsc(jCnttp)%CntMass
      else
        Mass(nc) = Zero
      end if
      nchr = dbsc(jCnttp)%AtmNr
      if (nchr >= 0) then
        Chrg(nc) = real(nchr,kind=wp)
      else
        Chrg(nc) = Zero
      end if
      nc = nc+1
    end do
    S%kCentr = S%kCentr+nIrrep/dc(ndc)%nStab
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute Total Charge and Center of Charge centroid

call CoW(Centr,CoC,Chrg,S%kCentr,qNuc)
if (iChCar(1) /= 0) CoC(1) = Zero
if (iChCar(2) /= 0) CoC(2) = Zero
if (iChCar(3) /= 0) CoC(3) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
! Modify charges to effective charges.

nc = 1
do jCnttp=1,nCnttp
  Z = dbsc(jCnttp)%Charge
  mCnt = dbsc(jCnttp)%nCntr
  if (dbsc(jCnttp)%Aux .or. dbsc(jCnttp)%Frag) cycle
  if (dbsc(jCnttp)%ECP .and. (dbsc(jCnttp)%nVal == 0)) cycle
  do jCnt=1,mCnt
    ndc = jCnt+dbsc(jCnttp)%mdci
    do i=0,nIrrep/dc(ndc)%nStab-1
      nchr = dbsc(jCnttp)%AtmNr
      Chrg(nc) = Z
      nc = nc+1
    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute Total Mass and Center of Mass

call CoW(Centr,CoM,Mass,S%kCentr,TMass)
if (iChCar(1) /= 0) CoM(1) = Zero
if (iChCar(2) /= 0) CoM(2) = Zero
if (iChCar(3) /= 0) CoM(3) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Gen_GeoList
