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
                     nl,iAng,jAng,CoeffAP,CoeffBP)

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 CoeffA(nPrimA,nConA), CoeffB(nPrimB,nConB), Coeff(nTheta_Full,nPhi), CoeffAP(nPrimA,nPrimA), CoeffBP(nPrimB,nPrimB)
integer List2(mData,nPhi_All), iD(NumCho), Indkl(nkl)
! Statement function
iTri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

!                                                                      *
!***********************************************************************
!                                                                      *

!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call RecPrt('CoeffA',' ',CoeffA,nPrimA,nConA)
call RecPrt('CoeffB',' ',CoeffB,nPrimB,nConB)
call iVcPrt('Indkl',' ',Indkl,nkl)
call iVcPrt('Mk_Coeffs: iD',' ',iD,NumCho)
write(6,*) 'iAng,jAng=',iAng,jAng
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
          iTheta_Full = iPrimA*(iPrimA-1)/2+iPrimB
          Coeff(iTheta_Full,iPhi) = Cff
        end do
      end do
    else
      iPhi_Full = (il-1)*nk+ik
      iPhi = Indkl(iPhi_Full)
      if (iPhi == 0) cycle
      do iPrimA=1,nPrimA
        do iPrimB=1,nPrimB
          Cff = CoeffA(iPrimA,ik)*CoeffB(iPrimB,il)/(CoeffAP(iPrimA,iPrimA)*CoeffBP(iPrimB,iPrimB))
          iTheta_Full = (iPrimB-1)*nPrimA+iPrimA
          Coeff(iTheta_Full,iPhi) = Cff
        end do
      end do
    end if

  end if
end do
#ifdef _DEBUGPRINT_
call RecPrt('Coeff',' ',Coeff,nTheta_Full,nPhi)
#endif

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nl)

end subroutine Mk_Coeffs
