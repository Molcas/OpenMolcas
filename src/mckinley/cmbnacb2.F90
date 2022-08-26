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

subroutine CmbnACB2(Fa1,Fa2,Fb1,Fb2,rFinal,Fact,nalpha,nbeta,C,nC,la,lb,iang,jfhess,Tmp,lsro)
!******************************************************************************
!
! Merges the second derivatives of ECP projection/SRO  integrals
! for derivatives of components
!
!******************************************************************************
!
! @parameter FA1    The first derivative of Left side , Includes no deriavtive
! @parameter FA2    The second derivative of Left side
! @parameter FB1    The first derivative of Right side . Includes no derivative
! @parameter FB2    The second derivative of Right side
! @parameter rFinal Result added up to (out)
! @parameter Fact   Factor the reult is multiplied with bef. added up
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
integer(kind=iwp) :: nalpha, nbeta, nC, la, lb, iang
real(kind=wp) :: FA1(nAlpha,nC,nTri_Elem1(la),2*iang+1,*), FA2(nAlpha,nC,nTri_Elem1(la),2*iang+1,*), &
                 FB1(nC,nBeta,2*iang+1,nTri_Elem1(lb),*), FB2(nC,nBeta,2*iang+1,nTri_Elem1(lb),*), &
                 rFinal(nAlpha*nbeta,nTri_Elem1(la),nTri_Elem1(lb),21), Fact, c(*), Tmp(*)
logical(kind=iwp) :: jfhess(4,3,4,3), lsro
integer(kind=iwp) :: ia, ib, iC, iCar, jCar, mvec, mvecB

!                                                                      *
!***********************************************************************
!                                                                      *
! Merge integrals with one deriavtive on each center

do iCar=1,3
  do jCar=1,3
    mvec = itri(iCar+3,jcar)
    if (jfHess(2,iCar,1,jcar)) then
      do ib=1,nTri_Elem1(lb)
        do ia=1,nTri_Elem1(la)
          do iC=1,(2*iAng+1)

            if (lsro) then
              call mult_sro(FA1(1,1,ia,ic,icar+1),nAlpha,C,nC,FB1(1,1,ic,ib,jcar+1),nBeta,Fact,rFinal(1,ia,ib,mVec),Tmp)
            else
              call DGEMM_('N','N',nAlpha,nBeta,nC,Fact,FA1(1,1,ia,ic,icar+1),nAlpha,FB1(1,1,ic,ib,jcar+1),nC,One, &
                          rFinal(1,ia,ib,mVec),nAlpha)
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
  do jCar=1,icar
    mvec = itri(iCar,jcar)
    if (jfHess(1,iCar,1,jcar)) then
      do ib=1,nTri_Elem1(lb)
        do ia=1,nTri_Elem1(la)
          do iC=1,(2*iAng+1)
            if (lsro) then
              call mult_sro(FA2(1,1,ia,ic,mvec),nAlpha,C,nC,FB1(1,1,ic,ib,1),nBeta,Fact,rFinal(1,ia,ib,mVec),Tmp)
            else
              call DGEMM_('N','N',nAlpha,nBeta,nC,Fact,FA2(1,1,ia,ic,mvec),nAlpha,FB1(1,1,ic,ib,1),nC,One,rFinal(1,ia,ib,mVec), &
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
  do jCar=1,icar
    mvec = itri(3+icar,3+jcar)
    mvecB = itri(icar,jcar)
    if (jfHess(2,iCar,2,jcar)) then

      do ib=1,nTri_Elem1(lb)
        do ia=1,nTri_Elem1(la)

          do iC=1,(2*iAng+1)
            if (lsro) then
              call mult_sro(FA1(1,1,ia,ic,1),nAlpha,C,nC,FB2(1,1,ic,ib,mvecb),nBeta,Fact,rFinal(1,ia,ib,mVec),Tmp)
            else
              call DGEMM_('N','N',nAlpha,nBeta,nC,Fact,FA1(1,1,ia,ic,1),nAlpha,FB2(1,1,ic,ib,mvecb),nC,One,rFinal(1,ia,ib,mVec), &
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
