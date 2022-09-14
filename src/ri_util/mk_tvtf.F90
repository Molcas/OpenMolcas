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

subroutine Mk_tVtF(TInt,nTheta_All,tVtF,nTheta,List2,mData,iPrm,nPrm,iAng,jAng,nk,Indkl,nkl,nTheta_Full)

use Index_Functions, only: nTri_Elem
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nTheta_All, nTheta, mData, List2(mData,nTheta_All), nPrm, iPrm(nPrm), iAng, jAng, nk, nkl, &
                                 Indkl(nkl), nTheta_Full
real(kind=wp), intent(in) :: TInt(nTheta_All,nTheta_All)
real(kind=wp), intent(out) :: tVtF(nTheta,nTheta_Full)
integer(kind=iwp) :: iA, ik, il, iTheta, iTheta_All, iTheta_Full, jA, jk, jl, jTheta_All, jTheta_Full, kComp, lComp, mComp, nComp

#ifdef _DEBUGPRINT_
call RecPrt('TInt',' ',TInt,nTheta_all,nTheta_all)
call iVcPrt('Indkl',' ',Indkl,nkl)
#endif
tVtF(:,:) = Zero
iA = iAng+1
jA = jAng+1
do iTheta_All=1,nTheta_All
  kComp = List2(3,iTheta_All)
  lComp = List2(4,iTheta_All)
  ik = List2(5,iTheta_All)
  il = List2(6,iTheta_All)
  if (iAng == jAng) then
    iTheta_Full = nTri_Elem(ik-1)+il
  else
    iTheta_Full = (il-1)*nk+ik
  end if
  !if ((iPrm(iTheta_Full) == 1) .and. (iAL(kComp,lComp) == 1) .and. (kComp == iA) .and. (lComp == jA)) then
  if ((iPrm(iTheta_Full) == 1) .and. (kComp == iA) .and. (lComp == jA)) then
    iTheta = Indkl(iTheta_Full)

    do jTheta_All=1,nTheta_All
      mComp = List2(3,jTheta_All)
      nComp = List2(4,jTheta_All)
      !if (iAL(mComp,nComp) == 1) then
      if ((mComp == iA) .and. (nComp == jA)) then
        jk = List2(5,jTheta_All)
        jl = List2(6,jTheta_All)
        if (iAng == jAng) then
          jTheta_Full = nTri_Elem(jk-1)+jl
        else
          jTheta_Full = (jl-1)*nk+jk
        end if

        tVtF(iTheta,jTheta_Full) = tVtF(iTheta,jTheta_Full)+TInt(iTheta_All,jTheta_All)
        !tVtF(iTheta,jTheta_Full) = tVtF(iTheta,jTheta_Full)+abs(TInt(iTheta_All,jTheta_All))

      end if
    end do

  end if
end do
#ifdef _DEBUGPRINT_
call RecPrt('tVtF',' ',tVtF,nTheta,nTheta_Full)
#endif

return

end subroutine Mk_tVtF
