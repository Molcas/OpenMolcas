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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine Assemble_dTdmu(nZeta,final,la,lb,Elalbp,Elalbm,Beta)
!***********************************************************************
!                                                                      *
! Object: to assemble the diamagnetic shielding integrals from         *
!         electric field integrals.                                    *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             February '91                                             *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,3), Elalbp(nZeta,(la+1)*(la+2)/2,(lb+2)*(lb+3)/2,3), &
       Elalbm(nZeta,(la+1)*(la+2)/2,(lb)*(lb+1)/2,3), Beta(nZeta)
character*80 Label
! Statement function for cartesian index
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1
nElem(ix) = (ix+1)*(ix+2)/2

iRout = 231
iPrint = nPrint(iRout)

!Fact = -1.D6*One2C2
if (iPrint >= 99) then
  write(6,*) ' In Assemble_dTdmu la,lb=',la,lb
  do ia=1,nElem(la)
    do ib=1,nElem(lb+1)
      write(Label,'(A,I2,A,I2,A)') ' Elalbp(',ia,',',ib,',x)'
      call RecPrt(Label,' ',Elalbp(1,ia,ib,1),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Elalbp(',ia,',',ib,',y)'
      call RecPrt(Label,' ',Elalbp(1,ia,ib,2),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Elalbp(',ia,',',ib,',z)'
      call RecPrt(Label,' ',Elalbp(1,ia,ib,3),nZeta,1)
    end do
  end do
  do ia=1,nElem(la)
    ib_max = nElem(lb-1)
    if (lb == 0) ib_max = 0
    do ib=1,ib_max
      write(Label,'(A,I2,A,I2,A)') ' Elalbm(',ia,',',ib,',x)'
      call RecPrt(Label,' ',Elalbm(1,ia,ib,1),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Elalbm(',ia,',',ib,',y)'
      call RecPrt(Label,' ',Elalbm(1,ia,ib,2),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Elalbm(',ia,',',ib,',z)'
      call RecPrt(Label,' ',Elalbm(1,ia,ib,3),nZeta,1)
    end do
  end do
end if

do ixa=la,0,-1
  do iya=la-ixa,0,-1
    iza = la-ixa-iya
    ipa = Ind(la,ixa,iza)

    do ixb=lb,0,-1
      do iyb=lb-ixb,0,-1
        izb = lb-ixb-iyb
        ipb = Ind(lb,ixb,izb)

        do iZeta=1,nZeta
          xyTmp = -Two*Beta(nzeta)*Elalbp(iZeta,ipa,Ind(lb+1,ixb,izb),1)
          yxTmp = -Two*Beta(nzeta)*Elalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),2)
          yzTmp = -Two*Beta(nzeta)*Elalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),2)
          zyTmp = -Two*Beta(nzeta)*Elalbp(iZeta,ipa,Ind(lb+1,ixb,izb),3)
          zxTmp = -Two*Beta(nzeta)*Elalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),3)
          xzTmp = -Two*Beta(nzeta)*Elalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),1)
          if (ixb >= 1) then
            yxTmp = yxTmp+dble(ixb)*Elalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),2)
            zxTmp = zxTmp+dble(ixb)*Elalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),3)
          end if
          if (iyb >= 1) then
            xyTmp = xyTmp+dble(iyb)*Elalbm(iZeta,ipa,Ind(lb-1,ixb,izb),1)
            zyTmp = xyTmp+dble(iyb)*Elalbm(iZeta,ipa,Ind(lb-1,ixb,izb),3)
          end if
          if (izb >= 1) then
            xzTmp = xzTmp+dble(izb)*Elalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),1)
            yzTmp = yzTmp+dble(izb)*Elalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),2)
          end if
          final(iZeta,ipa,ipb,1) = -(xyTmp-yxTmp)
          final(iZeta,ipa,ipb,2) = -(yzTmp-zyTmp)
          final(iZeta,ipa,ipb,3) = -(zxTmp-xzTmp)
        end do

      end do
    end do

  end do
end do

if (iPrint >= 49) then
  do iComp=1,3
    write(Label,'(A,I2,A)') ' Final (',iComp,') '
    call RecPrt(Label,' ',final(1,1,1,iComp),nZeta,nElem(la)*nELem(lb))
  end do
end if

return

end subroutine Assemble_dTdmu
