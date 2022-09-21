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

subroutine Mk_Coeffs(CoeffA,nPrimA,nConA,CoeffB,nPrimB,nConB,Coeff,nTheta_Full,nPhi,iD,NumCho,List2,mData,nPhi_All,Indkl,nkl,nk, &
                     iAng,jAng,CoeffAP,CoeffBP)

use Index_Functions, only: iTri
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nPrimA, nConA, nPrimB, nConB, nTheta_Full, nPhi, NumCho, iD(NumCho), mData, nPhi_All, &
                                 List2(mData,nPhi_All), nkl, Indkl(nkl), nk, iAng, jAng
real(kind=wp), intent(in) :: CoeffA(nPrimA,nConA), CoeffB(nPrimB,nConB), CoeffAP(nPrimA,nPrimA), CoeffBP(nPrimB,nPrimB)
real(kind=wp), intent(out) :: Coeff(nTheta_Full,nPhi)
integer(kind=iwp) :: iCho, ik, il, iPhi, iPhi_All, iPhi_Full, iPrimA, iPrimB, iTheta_Full
real(kind=wp) :: Cff

!                                                                      *
!***********************************************************************
!                                                                      *

!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call RecPrt('CoeffA',' ',CoeffA,nPrimA,nConA)
call RecPrt('CoeffB',' ',CoeffB,nPrimB,nConB)
call iVcPrt('Indkl',' ',Indkl,nkl)
call iVcPrt('Mk_Coeffs: iD',' ',iD,NumCho)
write(u6,*) 'iAng,jAng=',iAng,jAng
#endif
do iCho=1,NumCho
  iPhi_All = iD(iCho)
  if ((List2(1,iPhi_All) == iAng) .and. (List2(2,iPhi_All) == jAng)) then
    ik = List2(5,iPhi_All)
    il = List2(6,iPhi_All)

    if (iAng == jAng) then
      iPhi_Full = iTri(ik,il)
      iPhi = Indkl(iPhi_Full)
      if (iPhi == 0) cycle
      do iPrimA=1,nPrimA
        do iPrimB=1,iPrimA
          Cff = (CoeffA(iPrimA,ik)*CoeffB(iPrimB,il)+CoeffA(iPrimA,il)*CoeffB(iPrimB,ik))/ &
                (CoeffAP(iPrimA,iPrimA)*CoeffBP(iPrimB,iPrimB))
          if (iPrimA == iPrimB) Cff = Half*Cff
          iTheta_Full = iTri(iPrimA,iPrimB)
          Coeff(iTheta_Full,iPhi) = Cff
        end do
      end do
    else
      iPhi_Full = (il-1)*nk+ik
      iPhi = Indkl(iPhi_Full)
      if (iPhi == 0) cycle
      do iPrimA=1,nPrimA
        do iPrimB=1,nPrimB
          iTheta_Full = (iPrimB-1)*nPrimA+iPrimA
          Coeff(iTheta_Full,iPhi) = CoeffA(iPrimA,ik)*CoeffB(iPrimB,il)/(CoeffAP(iPrimA,iPrimA)*CoeffBP(iPrimB,iPrimB))
        end do
      end do
    end if

  end if
end do
#ifdef _DEBUGPRINT_
call RecPrt('Coeff',' ',Coeff,nTheta_Full,nPhi)
#endif

return

end subroutine Mk_Coeffs
