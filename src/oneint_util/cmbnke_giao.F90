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

subroutine CmbnKE_GIAO(Rxyz,nZeta,la,lb,lr,Zeta,rKappa,final,nComp,nB,Txyz,Wxyz,A,RB,C)
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

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
real*8 final(nZeta,nComp,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nB), Zeta(nZeta), rKappa(nZeta), Rxyz(nZeta,3,0:la+1,0:lb+1,0:lr+1), &
       Txyz(nZeta,3,0:la,0:lb,0:lr+1), Wxyz(nZeta,3,0:la,0:lb,2), A(3), RB(3), RAB(3), C(3)
integer index(3,2)
! Statement function for Cartesian index
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1

!iRout = 134
!iPrint = nPrint(iRout)

RAB(1) = A(1)-RB(1)
RAB(2) = A(2)-RB(2)
RAB(3) = A(3)-RB(3)

iComp = 1
do ixa=0,la
  iyaMax = la-ixa
  index(1,1) = ixa
  do ixb=0,lb
    iybMax = lb-ixb
    index(1,2) = ixb
    do iya=0,iyaMax
      iza = la-ixa-iya
      index(2,1) = iya
      index(3,1) = iza
      ipa = Ind(la,ixa,iza)
      do iyb=0,iybMax
        izb = lb-ixb-iyb
        index(2,2) = iyb
        index(3,2) = izb
        ipb = Ind(lb,ixb,izb)

        !if (iPrint >= 99) then
        !  write(6,*)
        !  write(6,*) ixa,iya,iza
        !  write(6,*) ixb,iyb,izb
        !  write(6,*) ipa,ipb
        !end if
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Loop over components of B

        do iBx=1,3
          iBy = iBx+1-((iBx+1)/4)*3
          iBz = iBy+1-((iBy+1)/4)*3
          jxa = index(iBx,1)
          jxb = index(iBx,2)
          jya = index(iBy,1)
          jyb = index(iBy,2)
          jza = index(iBz,1)
          jzb = index(iBz,2)
          !write(6,*) 'iBx,iBy,iBz=',iBx,iBy,iBz
          !write(6,*) 'nZeta=',nZeta
          !write(6,*) jxa,jya,jza
          !write(6,*) jxb,jyb,jzb
          !                                                            *
          !*************************************************************
          !                                                            *
          ! Combine integrals

          do iZeta=1,nZeta

            Fact = rKappa(iZeta)*Zeta(iZeta)**(-Three/Two)

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

            final(iZeta,iComp,ipa,ipb,iBx) = Half*Fact*(temp1+Half*(temp2a+temp2b+temp2c))
          end do
          !write(6,*)
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
