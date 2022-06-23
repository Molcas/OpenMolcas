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
! Copyright (C) 1994, Bernd Artur Hess                                 *
!***********************************************************************

subroutine Ass_pXp(Beta,nZeta,rFinal,la,lb,Slalbp,Slalbm,nComp)
!***********************************************************************
!                                                                      *
! Object: to assemble the pVp integrals                                *
!                                                                      *
!     Author: Bernd Hess, Institut fuer Physikalische und Theoretische *
!             Chemie, University of Bonn, Germany, August 1994         *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Constants, only: Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nZeta, la, lb, nComp
real(kind=wp) :: Beta(nZeta), rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),nComp), &
                 Slalbp(nZeta,nTri_Elem1(la),nTri_Elem1(lb+1),3,nComp), Slalbm(nZeta,nTri_Elem1(la),nTri_Elem1(lb-1),3,nComp)
#include "print.fh"
integer(kind=iwp) :: iComp, ipa, ipb, iPrint, iRout, ixa, ixb, iya, iyb, iza, izb, iZeta
character(len=80) :: Label
! Statement function for cartesian index
integer(kind=iwp) :: Ind, nElem, ixyz, ix, iz
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1
nElem(ix) = (ix+1)*(ix+2)/2

iRout = 211
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  write(u6,*)
  write(u6,*) ' In Ass_pXp la,lb,nComp,=',la,lb,nComp
  write(u6,*)
  call RecPrt('Beta','(10G15.8)',Beta,nZeta,1)
  do iComp=1,nComp
    write(u6,*) 'iComp=',iComp
    write(Label,'(A,I2,A)') ' Ass_pXp: Slalbp(1,iComp=',iComp,')'
    call RecPrt(Label,'(10G15.8)',Slalbp(1,1,1,1,iComp),nZeta,nElem(la)*nElem(lb+1))
    write(Label,'(A,I2,A)') ' Ass_pXp: Slalbp(2,iComp=',iComp,')'
    call RecPrt(Label,'(10G15.8)',Slalbp(1,1,1,2,iComp),nZeta,nElem(la)*nElem(lb+1))
    write(Label,'(A,I2,A)') ' Ass_pXp: Slalbp(3,iComp=',iComp,')'
    call RecPrt(Label,'(10G15.8)',Slalbp(1,1,1,3,iComp),nZeta,nElem(la)*nElem(lb+1))
    if (lb > 0) then
      write(Label,'(A,I2,A)') 'Ass_pXp: Slalbm(1,iComp=',iComp,')'
      call RecPrt(Label,'(10G15.8)',Slalbm(1,1,1,1,iComp),nZeta,nElem(la)*nElem(lb-1))
      write(Label,'(A,I2,A)') 'Ass_pXp: Slalbm(2,iComp=',iComp,')'
      call RecPrt(Label,'(10G15.8)',Slalbm(1,1,1,2,iComp),nZeta,nElem(la)*nElem(lb-1))
      write(Label,'(A,I2,A)') 'Ass_pXp: Slalbm(3,iComp=',iComp,')'
      call RecPrt(Label,'(10G15.8)',Slalbm(1,1,1,3,iComp),nZeta,nElem(la)*nElem(lb-1))
    end if
  end do
end if

do iComp=1,nComp

  do ixa=la,0,-1
    do iya=la-ixa,0,-1
      iza = la-ixa-iya
      ipa = Ind(la,ixa,iza)

      do ixb=lb,0,-1
        do iyb=lb-ixb,0,-1
          izb = lb-ixb-iyb
          ipb = Ind(lb,ixb,izb)

          do iZeta=1,nZeta
            rFinal(iZeta,ipa,ipb,iComp) = Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),1,iComp)+ &
                                          Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb),2,iComp)+ &
                                          Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),3,iComp)
          end do

          if (ixb > 0) then
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,iComp) = rFinal(iZeta,ipa,ipb,iComp)- &
                                            real(ixb,kind=wp)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),1,iComp)
            end do
          end if

          if (iyb > 0) then
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,iComp) = rFinal(iZeta,ipa,ipb,iComp)- &
                                            real(iyb,kind=wp)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb),2,iComp)
            end do
          end if

          if (izb > 0) then
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,iComp) = rFinal(iZeta,ipa,ipb,iComp)- &
                                            real(izb,kind=wp)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),3,iComp)
            end do
          end if

        end do
      end do

    end do
  end do

end do ! iComp

if (iPrint >= 49) then
  do iComp=1,nComp
    write(Label,'(A,I2,A,I2,A,I2,A)') ' Ass_pXp: pXp(iComp=',iComp,')'
    call RecPrt(Label,'(10G15.8)',rFinal(1,1,1,iComp),nZeta,nElem(la)*nElem(lb))
  end do
end if

return

end subroutine Ass_pXp
