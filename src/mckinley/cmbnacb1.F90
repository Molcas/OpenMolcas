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
! Copyright (C) 2004, Anders Bernhardsson                              *
!***********************************************************************

subroutine CmbnACB1(FA1,FB1,rFinal,Fact,nAlpha,nBeta,C,nC,la,lb,iang,ifgrad,tmp,lsro,indx,mvec,idcar)
!***********************************************************************
!
! Merges the first derivatives of ECP projection/SRO integrals
! for derivatives of components
!
!***********************************************************************
!
! @parameter FA1    The first derivative of Left side. Includes no deriavtive (input)
! @parameter FB1    The first derivative of Right side. Includes no derivative (input)
! @parameter rFinal Result added up to (out)
! @parameter Fact   Factor the reult is multiplied with bef. added up (input)
! @parameter C      Coefficients for SRO (input)
! @parameter nAlpha Number of exponents LS (input)
! @parameter nBeta  Number of exponents RS (input)
! @parameter nC     Number of exponents in SRO (input)
! @parameter la     Angular monenta LS (input)
! @parameter lb     Angular monenta RS (input)
! @parameter iAng   Angular monenta SRO (input)
! @parameter Tmp    Working Area nAlpha*nC (SRO case) (scratch)
! @parameter lSRO   true for SRO false projection operator (input)
! @parameter indx   Array storing index for derivatives in rFinal (out)
! @parameter mvec   Number of derivatives calculated (out)
! @parameter idcar  Cartesiam index for current derivative (input)
!
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAlpha, nBeta, nC, la, lb, iang, idcar
logical(kind=iwp), intent(in) :: ifgrad(3,4), lsro
real(kind=wp), intent(in) :: FA1(nAlpha,nC,nTri_Elem1(la),2*iang+1,2), FB1(nC,nBeta,2*iang+1,nTri_Elem1(lb),2), Fact, &
                             C(nC,merge(nC,0,lsro))
real(kind=wp), intent(out) :: rFinal(nAlpha*nBeta,nTri_Elem1(la),nTri_Elem1(lb),6), Tmp(nAlpha,merge(nC,0,lsro))
integer(kind=iwp), intent(out) :: indx(3,4), mvec
integer(kind=iwp) :: ia, ib, iC, iCent, iFa, iFb

rFinal(:,:,:,1:6) = Zero
indx(:,:) = 0

mVec = 0
do iCent=1,2
  if (ifGrad(iDCar,iCent)) then
    mVec = mVec+1
    indx(iDcar,icent) = mvec

    if (iCent == 1) then
      iFa = 2
      iFb = 1
    else
      iFa = 1
      iFb = 2
    end if

    do ib=1,nTri_Elem1(lb)
      do ia=1,nTri_Elem1(la)

        do iC=1,2*iAng+1
          if (lsro) then
            call mult_sro(FA1(:,:,ia,ic,iFa),nAlpha,C,nC,FB1(:,:,ic,ib,iFb),nBeta,Fact,rFinal(:,ia,ib,mvec),Tmp)
          else
            call DGEMM_('N','N',nAlpha,nBeta,nC,Fact,FA1(:,:,ia,ic,iFa),nAlpha,FB1(:,:,ic,ib,iFb),nC,One,rFinal(:,ia,ib,mvec), &
                        nAlpha)
          end if

        end do ! iC
      end do   ! iA
    end do     ! ib

  end if
end do ! icent

return

end subroutine CmbnACB1
