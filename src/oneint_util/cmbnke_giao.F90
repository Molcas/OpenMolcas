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
! Copyright (C) 1991,2002, Roland Lindh                                *
!***********************************************************************

subroutine CmbnKE_GIAO(Rxyz,nZeta,la,lb,lr,Zeta,rKappa,rFinal,nComp,nB,Txyz,Wxyz,A,RB,C)
!***********************************************************************
!                                                                      *
! Object: to compute the first derivative of the kinetic energy        *
!         integrals with respect to the magnetic field.                *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!                                                                      *
!     Modified for GIAO's by RL June 2002, Tokyo, Japan.               *
!***********************************************************************

use Index_Functions, only: C_Ind, nTri_Elem1
use Constants, only: Half, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, lr, nComp, nB
real(kind=wp), intent(in) :: Rxyz(nZeta,3,0:la+1,0:lb+1,0:lr+1), Zeta(nZeta), rKappa(nZeta), Txyz(nZeta,3,0:la,0:lb,0:lr+1), &
                             Wxyz(nZeta,3,0:la,0:lb,2), A(3), RB(3), C(3)
real(kind=wp), intent(out) :: rFinal(nZeta,nComp,nTri_Elem1(la),nTri_Elem1(lb),nB)
integer(kind=iwp) :: iBx, iBy, iBz, iComp, indx(3,2), ipa, ipb, ixa, ixb, iya, iyaMax, iyb, iybMax, iza, izb, iZeta, jxa, jxb, &
                     jya, jyb, jza, jzb
real(kind=wp) :: Fact, RAB(3), temp1, temp2a, temp2b, temp2c

!iRout = 134
!iPrint = nPrint(iRout)

RAB(:) = A-RB

iComp = 1
do ixa=0,la
  iyaMax = la-ixa
  indx(1,1) = ixa
  do ixb=0,lb
    iybMax = lb-ixb
    indx(1,2) = ixb
    do iya=0,iyaMax
      iza = la-ixa-iya
      indx(2,1) = iya
      indx(3,1) = iza
      ipa = C_Ind(la,ixa,iza)
      do iyb=0,iybMax
        izb = lb-ixb-iyb
        indx(2,2) = iyb
        indx(3,2) = izb
        ipb = C_Ind(lb,ixb,izb)

        !if (iPrint >= 99) then
        !  write(u6,*)
        !  write(u6,*) ixa,iya,iza
        !  write(u6,*) ixb,iyb,izb
        !  write(u6,*) ipa,ipb
        !end if
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Loop over components of B

        do iBx=1,3
          iBy = iBx+1-((iBx+1)/4)*3
          iBz = iBy+1-((iBy+1)/4)*3
          jxa = indx(iBx,1)
          jxb = indx(iBx,2)
          jya = indx(iBy,1)
          jyb = indx(iBy,2)
          jza = indx(iBz,1)
          jzb = indx(iBz,2)
          !write(u6,*) 'iBx,iBy,iBz=',iBx,iBy,iBz
          !write(u6,*) 'nZeta=',nZeta
          !write(u6,*) jxa,jya,jza
          !write(u6,*) jxb,jyb,jzb
          !                                                            *
          !*************************************************************
          !                                                            *
          ! Combine integrals

          do iZeta=1,nZeta

            Fact = rKappa(iZeta)*Zeta(iZeta)**(-OneHalf)

            temp1 = Rxyz(iZeta,iBx,jxa,jxb,0)*(Wxyz(iZeta,iBy,jya,jyb,2)*Rxyz(iZeta,iBz,jza+1,jzb,0)- &
                                               Rxyz(iZeta,iBy,jya+1,jyb,0)*Wxyz(iZeta,iBz,jza,jzb,2)- &
                                               Wxyz(iZeta,iBy,jya,jyb,1)*Rxyz(iZeta,iBz,jza,jzb+1,0)+ &
                                               Rxyz(iZeta,iBy,jya,jyb+1,0)*Wxyz(iZeta,iBz,jza,jzb,1))

            temp2a = Txyz(iZeta,iBx,jxa,jxb,0)*(RAB(iBy)*Rxyz(iZeta,iBy,jya,jyb,0)* &
                                                (Rxyz(iZeta,iBz,jza,jzb,1)+Rxyz(iZeta,iBz,jza,jzb,0)*C(iBz))- &
                                                (Rxyz(iZeta,iBy,jya,jyb,1)+Rxyz(iZeta,iBy,jya,jyb,0)*C(iBy))* &
                                                RAB(iBz)*Rxyz(iZeta,iBz,jza,jzb,0))

            temp2b = Rxyz(iZeta,iBx,jxa,jxb,0)*(RAB(iBy)*Txyz(iZeta,iBy,jya,jyb,0)* &
                                                (Rxyz(iZeta,iBz,jza,jzb,1)+Rxyz(iZeta,iBz,jza,jzb,0)*C(iBz))- &
                                                (Txyz(iZeta,iBy,jya,jyb,1)+Txyz(iZeta,iBy,jya,jyb,0)*C(iBy))* &
                                                RAB(iBz)*Rxyz(iZeta,iBz,jza,jzb,0))

            temp2c = Rxyz(iZeta,iBx,jxa,jxb,0)*(RAB(iBy)*Rxyz(iZeta,iBy,jya,jyb,0)* &
                                                (Txyz(iZeta,iBz,jza,jzb,1)+Txyz(iZeta,iBz,jza,jzb,0)*C(iBz))- &
                                                (Rxyz(iZeta,iBy,jya,jyb,1)+Rxyz(iZeta,iBy,jya,jyb,0)*C(iBy))* &
                                                RAB(iBz)*Txyz(iZeta,iBz,jza,jzb,0))

            rFinal(iZeta,iComp,ipa,ipb,iBx) = Half*Fact*(temp1+Half*(temp2a+temp2b+temp2c))
          end do
          !write(u6,*)
          !                                                            *
          !*************************************************************
          !                                                            *
        end do ! iBx
        !                                                              *
        !***************************************************************
        !                                                              *
      end do
    end do
  end do
end do

return

end subroutine CmbnKE_GIAO
