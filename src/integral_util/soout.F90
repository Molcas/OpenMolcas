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

subroutine SOOUT(label,cnt_ico,phase_ico)

use Index_Functions, only: nTri_Elem1
use Basis_Info, only: dbsc, nCnttp, Shells
use Center_Info, only: dc
use Symmetry_Info, only: iChBas, iChTbl, nIrrep
use Real_Spherical, only: iSphCr, LblCBs, LblSBs
use Molcas, only: LenIn, MaxBfn, MaxBfn_aux
use Definitions, only: iwp

#include "intent.fh"

implicit none
character(len=LenIn+8), intent(out) :: Label(MaxBfn+MaxBfn_Aux)
integer(kind=iwp), intent(_OUT_) :: cnt_ico(0:7,*), phase_ico(0:7,*)
integer(kind=iwp) :: iAng, iChBs, iCnt, iCntrc, iCnttp, iCo, iComp, iIrrep, iSh, iSO, jComp, kComp, lComp, mc, mdc, nBasisi, &
                     nExpi, NrOpr
character(len=8) :: ChTemp
integer(kind=iwp), external :: iPrmt
logical(kind=iwp), external :: TstFnc

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
    if (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag) cycle

    ! Loop over distinct centers

    do iCnt=1,dbsc(iCnttp)%nCntr
      mdc = mdc+1

      ! Loop over shells associated with this center
      ! Start with s type shells

      kComp = 0
      iSh = dbsc(iCnttp)%iVal-1
      do iAng=0,dbsc(iCnttp)%nVal-1
        kComp = kComp+nTri_Elem1(iAng-1)
        iSh = iSh+1
        nExpi = Shells(iSh)%nExp
        if (nExpi == 0) cycle
        nBasisi = Shells(iSh)%nBasis
        if (nBasisi == 0) cycle
        if (Shells(iSh)%Prjct) then
          jComp = 2*iAng+1
        else
          jComp = nTri_Elem1(iAng)
        end if
        do iComp=1,jComp
          lComp = kComp+iComp
          ! Get character of basis function
          iChBs = iChBas(lComp)
          if (Shells(iSh)%Transf) iChBs = iChBas(iSphCr(lComp))

          ! Skip if function not a basis of irreps.

          if (.not. TstFnc(dc(mdc)%iCoSet,iIrrep,iChBs,dc(mdc)%nStab)) cycle

          do iCntrc=1,nBasisi
            iSO = iSO+1
            if (iSO > size(Label)) then
              call WarningMessage(2,'SOout: iSO > size(Label)')
              call Abend()
            end if
            ChTemp = LblCBs(lComp)
            if (Shells(iSh)%Transf) ChTemp = LblSbs(lComp)
            do ico=0,nIrrep/dc(mdc)%nStab-1
              Cnt_ico(ico,iso) = mc+ico
              Phase_ico(ico,iso) = iPrmt(NrOpr(dc(mdc)%iCoSet(iCo,0)),iChbs)*iChTbl(iIrrep,NrOpr(dc(mdc)%iCoSet(iCo,0)))
            end do
            Label(iSO) = dc(mdc)%LblCnt(1:LenIn)//ChTemp(1:8)
          end do

        end do
      end do
      mc = mc+nIrrep/dc(mdc)%nStab
    end do

  end do
end do

end subroutine SOOUT
