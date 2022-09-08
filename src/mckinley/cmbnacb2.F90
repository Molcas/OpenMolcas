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

subroutine CmbnACB2(FA1,FA2,FB1,FB2,rFinal,Fact,nAlpha,nBeta,C,nC,la,lb,iang,jfhess,Tmp,lsro)
!******************************************************************************
!
! Merges the second derivatives of ECP projection/SRO integrals
! for derivatives of components
!
!******************************************************************************
!
! @parameter FA1    The first derivative of Left side. Includes no derivative
! @parameter FA2    The second derivative of Left side
! @parameter FB1    The first derivative of Right side. Includes no derivative
! @parameter FB2    The second derivative of Right side
! @parameter rFinal Result added up to (out)
! @parameter Fact   Factor the result is multiplied with before added up
! @parameter C      Coefficients for SRO
! @parameter nAlpha Number of exponents LS
! @parameter nBeta  Number of exponents RS
! @parameter nC     Number of exponents in SRO
! @parameter la     Angular monenta LS
! @parameter lb     Angular monenta RS
! @parameter iAng   Angular monenta SRO
! @parameter nBeta  Number of exponents RS
! @parameter nC     Number of exponents in SRO
! @parameter Tmp    Working Area nAlpha*nC (SRO case)
! @parameter lSRO   true for SRO false projection operator
!
!******************************************************************************

use Index_Functions, only: iTri, nTri_Elem1
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAlpha, nBeta, nC, la, lb, iang
logical(kind=iwp), intent(in) :: jfHess(4,3,4,3), lsro
real(kind=wp), intent(in) :: FA1(nAlpha,nC,nTri_Elem1(la),2*iang+1,4), FA2(nAlpha,nC,nTri_Elem1(la),2*iang+1,6), &
                             FB1(nC,nBeta,2*iang+1,nTri_Elem1(lb),4), FB2(nC,nBeta,2*iang+1,nTri_Elem1(lb),6), Fact, &
                             C(nC,merge(nC,0,lsro))
real(kind=wp), intent(inout) :: rFinal(nAlpha*nBeta,nTri_Elem1(la),nTri_Elem1(lb),21)
real(kind=wp), intent(out) :: Tmp(nAlpha,merge(nC,0,lsro))
integer(kind=iwp) :: ia, ib, iC, iCar, jCar, mVec, mVecB

!                                                                      *
!***********************************************************************
!                                                                      *
! Merge integrals with one derivative on each center

do iCar=1,3
  do jCar=1,3
    mVec = itri(iCar+3,jCar)
    if (jfHess(2,iCar,1,jCar)) then
      do ib=1,nTri_Elem1(lb)
        do ia=1,nTri_Elem1(la)
          do iC=1,2*iAng+1

            if (lsro) then
              call mult_sro(FA1(:,:,ia,ic,iCar+1),nAlpha,C,nC,FB1(:,:,ic,ib,jCar+1),nBeta,Fact,rFinal(:,ia,ib,mVec),Tmp)
            else
              call DGEMM_('N','N',nAlpha,nBeta,nC,Fact,FA1(:,:,ia,ic,iCar+1),nAlpha,FB1(:,:,ic,ib,jCar+1),nC,One, &
                          rFinal(:,ia,ib,mVec),nAlpha)
            end if
          end do
        end do
      end do
    end if
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Merge integrals with both derivative on center A
do iCar=1,3
  do jCar=1,iCar
    mVec = itri(iCar,jCar)
    if (jfHess(1,iCar,1,jCar)) then
      do ib=1,nTri_Elem1(lb)
        do ia=1,nTri_Elem1(la)
          do iC=1,2*iAng+1
            if (lsro) then
              call mult_sro(FA2(:,:,ia,ic,mVec),nAlpha,C,nC,FB1(:,:,ic,ib,1),nBeta,Fact,rFinal(:,ia,ib,mVec),Tmp)
            else
              call DGEMM_('N','N',nAlpha,nBeta,nC,Fact,FA2(:,:,ia,ic,mVec),nAlpha,FB1(:,:,ic,ib,1),nC,One,rFinal(:,ia,ib,mVec), &
                          nAlpha)
            end if
          end do
        end do
      end do
    end if
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Merge integrals with both derivative on center B
do iCar=1,3
  do jCar=1,iCar
    mVec = itri(3+iCar,3+jCar)
    mVecB = itri(iCar,jCar)
    if (jfHess(2,iCar,2,jCar)) then

      do ib=1,nTri_Elem1(lb)
        do ia=1,nTri_Elem1(la)

          do iC=1,2*iAng+1
            if (lsro) then
              call mult_sro(FA1(:,:,ia,ic,1),nAlpha,C,nC,FB2(:,:,ic,ib,mVecB),nBeta,Fact,rFinal(:,ia,ib,mVec),Tmp)
            else
              call DGEMM_('N','N',nAlpha,nBeta,nC,Fact,FA1(:,:,ia,ic,1),nAlpha,FB2(:,:,ic,ib,mVecB),nC,One,rFinal(:,ia,ib,mVec), &
                          nAlpha)
            end if
          end do
        end do
      end do
    end if
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine CmbnACB2
