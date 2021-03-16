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
! Copyright (C) Ben Swerts                                             *
!***********************************************************************

subroutine AddFragDens(Array,nDens,nDens_Valence,nBas_Valence)
!***********************************************************************
!                                                                      *
! Input: Array(Size) filled with valence density at proper positions   *
!        for a density including fragments.                            *
! Output: Updated with fragment densities at their proper positions.   *
!                                                                      *
!***********************************************************************

use Basis_Info, only: dbsc, nBas, nCnttp
use Center_Info, only: dc
use Symmetry_Info, only: nIrrep, iOper
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nDens, nDens_Valence, nBas_Valence(0:7)
real(kind=wp), intent(inout) :: Array(nDens)
#include "WrkSpc.fh"
logical(kind=iwp) :: EnergyWeight
integer(kind=iwp) :: iPrint, maxDens, iCnttp, ipFragDensAO, iDpos, ipFragDensSO, i, j, iCnt, iFpos, iFD, mdc, iIrrep, nBasC
real(kind=wp) :: rDummy(1)

if (nIrrep /= 1) then
  write(u6,*) 'AddFragDens: Symmetry not implemented yet'
  call Abend()
end if

iPrint = 0

! Each fragment needs its (symmetrized) density matrix added along the diagonal
! This density matrix first has to be constructed from the MO coefficients
! so allocate space for the largest possible density matrix
maxDens = 0
do iCnttp=1,nCnttp
  if (dbsc(iCnttp)%nFragType > 0) maxDens = max(maxDens,dbsc(iCnttp)%nFragDens*(dbsc(iCnttp)%nFragDens+1)/2)
end do
call GetMem('FragDSO','Allo','Real',ipFragDensSO,maxDens)
!if (nIrrep /= 1) Then
!  call GetMem('FragDAO','Allo','Real',ipFragDensAO,maxDens)
!else
ipFragDensAO = ipFragDensSO
!end If

iDpos = 1 ! position in the total density matrix
do iIrrep=0,nIrrep-1
  nBasC = nBas_Valence(iIrrep)
  iDpos = iDpos+nBasC*(nBasC+1)/2
  mdc = 0
  do iCnttp=1,nCnttp
    if (dbsc(iCnttp)%nFragType <= 0) then
      mdc = mdc+dbsc(iCnttp)%nCntr
      Go To 1000
    end if

    ! construct the density matrix
    EnergyWeight = .false.
    call MakeDens(dbsc(iCnttp)%nFragDens,dbsc(iCnttp)%nFragEner,dbsc(iCnttp)%FragCoef,rDummy,EnergyWeight,Work(ipFragDensAO))
    ! create the symmetry adapted version if necessary
    ! (fragment densities are always calculated without symmetry)
    !if (nIrrep /= 1) call SymmDens(Work(ipFragDensAO),Work(ipFragDensSO))
    if (iPrint >= 99) call TriPrt('Fragment density',' ',Work(ipFragDensSO),dbsc(iCnttp)%nFragDens)
    do iCnt=1,dbsc(iCnttp)%nCntr
      mdc = mdc+1
      ! only add fragment densities that are active in this irrep
      ! => the following procedure still has to be verified thoroughly
      !    but appears to be working
      if (iand(dc(mdc)%iChCnt,iIrrep) == iOper(iIrrep)) then
        ! add it at the correct location in the large custom density matrix
        iFpos = 1
        ! position in fragment density matrix
        do i=1,dbsc(iCnttp)%nFragDens
          iDpos = iDpos+nBasC
          do j=0,i-1
            Array(iDpos+j) = Work(ipFragDensSO+iFpos+j-1)
          end do
          iDpos = iDpos+i
          iFpos = iFpos+i
        end do
        nBasC = nBasC+dbsc(iCnttp)%nFragDens
      end if
    end do
1000 continue
  end do
end do
if (iPrint >= 19) then
  iFD = 1
  do iIrrep=0,nIrrep-1
    call TriPrt('Combined density',' ',Array(iFD),nBas(iIrrep))
    iFD = iFD+nBas(iIrrep)*(nBas(iIrrep)+1)/2
  end do
end if
call GetMem('FragDSO','Free','Real',ipFragDensSO,maxDens)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nDens_Valence)

end subroutine AddFragDens
