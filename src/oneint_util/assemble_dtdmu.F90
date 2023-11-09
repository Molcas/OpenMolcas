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

subroutine Assemble_dTdmu(nZeta,rFinal,la,lb,Elalbp,Elalbm,Beta)
!***********************************************************************
!                                                                      *
! Object: to assemble the diamagnetic shielding integrals from         *
!         electric field integrals.                                    *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             February '91                                             *
!***********************************************************************

use Index_Functions, only: C_Ind, nTri_Elem1
use Constants, only: Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb
real(kind=wp), intent(out) :: rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),3)
real(kind=wp), intent(in) :: Elalbp(nZeta,nTri_Elem1(la),nTri_Elem1(lb+1),3), Elalbm(nZeta,nTri_Elem1(la),nTri_Elem1(lb-1),3), &
                             Beta(nZeta)
#include "print.fh"
integer(kind=iwp) :: ia, ib, ib_max, iComp, ipa, ipb, iPrint, iRout, ixa, ixb, iya, iyb, iza, izb, iZeta
real(kind=wp) :: xyTmp, xzTmp, yxTmp, yzTmp, zxTmp, zyTmp
character(len=80) :: Label

iRout = 231
iPrint = nPrint(iRout)

!Fact = -1.0e6_wp*Half/c_in_au**2
if (iPrint >= 99) then
  write(u6,*) ' In Assemble_dTdmu la,lb=',la,lb
  do ia=1,nTri_Elem1(la)
    do ib=1,nTri_Elem1(lb+1)
      write(Label,'(A,I2,A,I2,A)') ' Elalbp(',ia,',',ib,',x)'
      call RecPrt(Label,' ',Elalbp(:,ia,ib,1),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Elalbp(',ia,',',ib,',y)'
      call RecPrt(Label,' ',Elalbp(:,ia,ib,2),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Elalbp(',ia,',',ib,',z)'
      call RecPrt(Label,' ',Elalbp(:,ia,ib,3),nZeta,1)
    end do
  end do
  do ia=1,nTri_Elem1(la)
    ib_max = nTri_Elem1(lb-1)
    if (lb == 0) ib_max = 0
    do ib=1,ib_max
      write(Label,'(A,I2,A,I2,A)') ' Elalbm(',ia,',',ib,',x)'
      call RecPrt(Label,' ',Elalbm(:,ia,ib,1),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Elalbm(',ia,',',ib,',y)'
      call RecPrt(Label,' ',Elalbm(:,ia,ib,2),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Elalbm(',ia,',',ib,',z)'
      call RecPrt(Label,' ',Elalbm(:,ia,ib,3),nZeta,1)
    end do
  end do
end if

do ixa=la,0,-1
  do iya=la-ixa,0,-1
    iza = la-ixa-iya
    ipa = C_Ind(la,ixa,iza)

    do ixb=lb,0,-1
      do iyb=lb-ixb,0,-1
        izb = lb-ixb-iyb
        ipb = C_Ind(lb,ixb,izb)

        do iZeta=1,nZeta
          xyTmp = -Two*Beta(nzeta)*Elalbp(iZeta,ipa,C_Ind(lb+1,ixb,izb),1)
          yxTmp = -Two*Beta(nzeta)*Elalbp(iZeta,ipa,C_Ind(lb+1,ixb+1,izb),2)
          yzTmp = -Two*Beta(nzeta)*Elalbp(iZeta,ipa,C_Ind(lb+1,ixb,izb+1),2)
          zyTmp = -Two*Beta(nzeta)*Elalbp(iZeta,ipa,C_Ind(lb+1,ixb,izb),3)
          zxTmp = -Two*Beta(nzeta)*Elalbp(iZeta,ipa,C_Ind(lb+1,ixb+1,izb),3)
          xzTmp = -Two*Beta(nzeta)*Elalbp(iZeta,ipa,C_Ind(lb+1,ixb,izb+1),1)
          if (ixb >= 1) then
            yxTmp = yxTmp+real(ixb,kind=wp)*Elalbm(iZeta,ipa,C_Ind(lb-1,ixb-1,izb),2)
            zxTmp = zxTmp+real(ixb,kind=wp)*Elalbm(iZeta,ipa,C_Ind(lb-1,ixb-1,izb),3)
          end if
          if (iyb >= 1) then
            xyTmp = xyTmp+real(iyb,kind=wp)*Elalbm(iZeta,ipa,C_Ind(lb-1,ixb,izb),1)
            zyTmp = xyTmp+real(iyb,kind=wp)*Elalbm(iZeta,ipa,C_Ind(lb-1,ixb,izb),3)
          end if
          if (izb >= 1) then
            xzTmp = xzTmp+real(izb,kind=wp)*Elalbm(iZeta,ipa,C_Ind(lb-1,ixb,izb-1),1)
            yzTmp = yzTmp+real(izb,kind=wp)*Elalbm(iZeta,ipa,C_Ind(lb-1,ixb,izb-1),2)
          end if
          rFinal(iZeta,ipa,ipb,1) = -(xyTmp-yxTmp)
          rFinal(iZeta,ipa,ipb,2) = -(yzTmp-zyTmp)
          rFinal(iZeta,ipa,ipb,3) = -(zxTmp-xzTmp)
        end do

      end do
    end do

  end do
end do

if (iPrint >= 49) then
  do iComp=1,3
    write(Label,'(A,I2,A)') ' rFinal (',iComp,') '
    call RecPrt(Label,' ',rFinal(:,:,:,iComp),nZeta,nTri_Elem1(la)*nTri_Elem1(lb))
  end do
end if

return

end subroutine Assemble_dTdmu
