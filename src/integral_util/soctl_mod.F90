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

subroutine SOCtl_mod(Mamn,nMamn,Cnt_ico,Phase_ico)

use Basis_Info
use Center_Info
use Symmetry_Info, only: nIrrep, iChTbl, iChBas
use Real_Spherical, only: iSphCr, LblCBs, LblSBs
use Constants

implicit none
#include "Molcas.fh"
integer nMamn
character(len=LENIN8) Mamn(nMamn)
integer cnt_ico(0:7,*), phase_ico(0:7,*)
integer, external :: iPrmt
character(len=8) ChTemp
logical, external :: TstFnc
integer iSO, iIrrep, mdc, mc, iCnttp, kComp, iSh, iAng, nExpi, nBasisi, jComp, iCnt, iComp, iChBs, iCntrc, NrOpr, lComp, iCo

! Generate list of symmetry adapted or petite list basis functions

! Loop over Irreps
iSO = 0

! Loop over irreducible representations and symmetry operations,
! respectively, for SO and Petite list, respectively.

do iIrrep=0,nIrrep-1

  ! Loop over distinct shell types

  mdc = 0
  mc = 1
  do iCnttp=1,nCnttp
    if (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag) Go To 201

    ! Loop over distinct centers

    do iCnt=1,dbsc(iCnttp)%nCntr
      mdc = mdc+1

      ! Loop over shells associated with this center
      ! Start with s type shells

      kComp = 0
      iSh = dbsc(iCnttp)%iVal-1
      do iAng=0,dbsc(iCnttp)%nVal-1
        iSh = iSh+1
        nExpi = Shells(iSh)%nExp
        if (nExpi == 0) Go To 2033
        nBasisi = Shells(iSh)%nBasis
        if (nBasisi == 0) Go To 2033
        if (Shells(iSh)%Prjct) then
          jComp = 2*iAng+1
        else
          jComp = (iAng+1)*(iAng+2)/2
        end if
        do iComp=1,jComp
          lComp = kComp+iComp
          ! Get character of basis function
          iChBs = iChBas(lComp)
          if (Shells(iSh)%Transf) iChBs = iChBas(iSphCr(lComp))

          ! Skip if function not a basis of irreps.

          if (.not. TstFnc(dc(mdc)%iCoSet,iIrrep,iChBs,dc(mdc)%nStab)) Go To 204

          do iCntrc=1,nBasisi
            iSO = iSO+1
            if (iSO > nMamn) then
              call WarningMessage(2,'SOout: iSO > nMamn')
              call Abend()
            end if
            ChTemp = LblCBs(lComp)
            if (Shells(iSh)%Transf) ChTemp = LblSbs(lComp)
            do ico=0,nIrrep/dc(mdc)%nStab-1
              Cnt_ico(ico,iso) = mc+ico
              Phase_ico(ico,iso) = iPrmt(NrOpr(dc(mdc)%iCoSet(iCo,0)),iChbs)*iChTbl(iIrrep,NrOpr(dc(mdc)%iCoSet(iCo,0)))
            end do
            Mamn(iSO) = dc(mdc)%LblCnt(1:LENIN)//ChTemp(1:8)
          end do

204       continue
        end do
2033    continue
        kComp = kComp+(iAng+1)*(iAng+2)/2
      end do
      mc = mc+nIrrep/dc(mdc)%nStab
    end do

201 continue
  end do
end do

end subroutine SOCtl_mod
