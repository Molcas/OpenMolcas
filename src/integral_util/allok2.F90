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
! Copyright (C) 1992, Roland Lindh                                     *
!               1995, Martin Schuetz                                   *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine AlloK2()
!***********************************************************************
!                                                                      *
!  Object: Allocate space for K2 entities.                             *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, Sweden. November '92                 *
!             Martin Schuetz, Dept. of Theoretical Chemistry,          *
!             University of Lund, Sweden. Jun '95                      *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem, nTri_Elem1, nTri3_Elem1
use k2_arrays, only: DoGrad_, DoHess_, MaxDe, nDeDe
use iSD_data, only: iSD
use Basis_Info, only: Shells
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep
use k2_structure, only: Allocate_k2data, Allocate_k2data_in, IndK2, k2_Processed, k2data, ZZZ_i, ZZZ_r
use stdalloc, only: mma_allocate
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: iAng, iAO, iBas, iCmp, iDeSiz, iIrrep, ijCmp, ijS, ik2, iPrim, iS, iShell, iShll, iSMLbl, jAng, jAO, jBas, &
                     jCmp, jPrim, jS, jShell, jShll, nData, nHm, nIndK2, nk2_integer, nk2_real, Nr_of_Densities, nSkal, nSO, nZeta
integer(kind=iwp), external :: MemSO1

#ifdef _DEBUGPRINT_
if (allocated(k2Data) .and. k2_processed) then
  write(u6,*) 'Enter Allok2, k2_Status=Processed'
else if (allocated(k2Data)) then
  write(u6,*) 'Enter Allok2, k2_Status=Active'
else
  write(u6,*) 'Enter Allok2, k2_Status=InActive'
end if
#endif
if (allocated(k2Data) .or. k2_processed) return

call Nr_Shells(nSkal)

ik2 = 0
nk2_real = 0
nk2_integer = 0
do iS=1,nSkal
  iShll = iSD(0,iS)
  if (Shells(iShll)%Aux .and. (iS /= nSkal)) cycle
  iPrim = iSD(5,iS)
  iAng = iSD(1,iS)
  iCmp = iSD(2,iS)
  do jS=1,iS
    jShll = iSD(0,jS)
    if (Shells(iShll)%Aux .and. (.not. Shells(jShll)%Aux)) cycle
    if (Shells(jShll)%Aux .and. (jS == nSkal)) cycle
    jPrim = iSD(5,jS)
    jAng = iSD(1,jS)
    jCmp = iSD(2,jS)

    nZeta = iPrim*jPrim
    ijCmp = 0
    if (DoGrad_) ijCmp = nTri_Elem1(iAng)*nTri_Elem1(jAng)
    nHm = iCmp*jCmp*(nTri3_Elem1(iAng+jAng)-nTri3_Elem1(max(iAng,jAng)-1))
    if (DoHess_) nHm = 0

    ik2 = ik2+1

    nData = nZeta*(10+ijCmp*2)+nHm*nIrrep
    nk2_real = nk2_real+nData*nIrrep
    nk2_integer = nk2_integer+(nZeta+1)*nIrrep

  end do
end do

call mma_allocate(ZZZ_r,nk2_real,Label='ZZZ_r')
call mma_allocate(ZZZ_i,nk2_integer,Label='ZZZ_i')
call Allocate_k2Data(nIrrep,ik2)

! determine memory size for K2 entities
! for this, run dummy K2 loop...
ik2 = 0
nDeDe = 0
MaxDe = 0
nr_of_Densities = 1 ! Hardwired option

nIndk2 = nTri_Elem(S%nShlls)
call mma_allocate(Indk2,3,nIndk2,Label='Indk2')
Indk2(:,:) = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! Double loop over shells. These loops decide the integral type

do iS=1,nSkal
  iShll = iSD(0,iS)
  if (Shells(iShll)%Aux .and. (iS /= nSkal)) cycle
  iAng = iSD(1,iS)
  iCmp = iSD(2,iS)
  iBas = iSD(3,iS)
  iPrim = iSD(5,iS)
  iAO = iSD(7,iS)
  iShell = iSD(11,iS)

  do jS=1,iS
    jShll = iSD(0,jS)
    if (Shells(iShll)%Aux .and. (.not. Shells(jShll)%Aux)) cycle
    if (Shells(jShll)%Aux .and. (jS == nSkal)) cycle
    jAng = iSD(1,jS)
    jCmp = iSD(2,jS)
    jBas = iSD(3,jS)
    jPrim = iSD(5,jS)
    jAO = iSD(7,jS)
    jShell = iSD(11,jS)

    if (nIrrep == 1) then
      iDeSiz = 1+iPrim*jPrim+iCmp*jCmp
    else
      iDeSiz = 1+iPrim*jPrim+(iBas*jBas+1)*iCmp*jCmp
    end if

    MaxDe = max(MaxDe,iDeSiz)
    iSmLbl = 1
    nSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
    if (nSO > 0) nDeDe = nDeDe+nr_of_Densities*iDeSiz*nIrrep

    nZeta = iPrim*jPrim
    ijCmp = 0
    if (DoGrad_) ijCmp = nTri_Elem1(iAng)*nTri_Elem1(jAng)
    nHm = iCmp*jCmp*(nTri3_Elem1(iAng+jAng)-nTri3_Elem1(max(iAng,jAng)-1))
    if (DoHess_) nHm = 0

    ijS = iTri(iShell,jShell)
    ik2 = ik2+1
    Indk2(3,ijS) = ik2
    do iIrrep=1,nIrrep
      call Allocate_k2data_in(k2data(iIrrep,ik2),nZeta,ijCmp,nHm)
    end do

  end do
end do

return

end subroutine AlloK2
