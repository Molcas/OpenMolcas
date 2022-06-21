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

subroutine Util8(Beta,nZeta,final,la,lb,Slalbp,Slalbm)
!***********************************************************************
!                                                                      *
! Object: to assemble the Vp integrals from                            *
!         derivative integrals of the electric potential.              *
!                                                                      *
!     Author: Bernd Hess, Institut fuer Physikalische und Theoretische *
!             Chemie, University of Bonn, Germany, August 1994         *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,3), Slalbp(nZeta,(la+1)*(la+2)/2,(lb+2)*(lb+3)/2), &
       Slalbm(nZeta,(la+1)*(la+2)/2,lb*(lb+1)/2), Beta(nZeta)
character*80 Label
! Statement function for cartesian index
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1
nElem(ix) = (ix+1)*(ix+2)/2

iRout = 203
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  write(6,*) ' In util8 la,lb=',la,lb
  call RecPrt('Beta','(5f15.8)',Beta,nZeta,1)
  do ib=1,nElem(lb)
    write(Label,'(A,I2,A)') ' Slalbp(',la,ib,')'
    call RecPrt(Label,'(5f15.8)',Slalbp(1,1,ib),nZeta,nElem(la+1))
  end do
  if (lb > 0) then
    do ia=1,nElem(la)
      write(Label,'(A,I2,A)') ' Slalbm(',la,ib,')'
      call RecPrt(Label,'(5f15.8)',Slalbm(1,1,ib),nZeta,nElem(lb-1))
    end do
  end if
end if

do ixa=la,0,-1
  do iya=la-ixa,0,-1
    iza = la-ixa-iya
    ipa = Ind(la,ixa,iza)

    do ixb=lb,0,-1
      do iyb=lb-ixb,0,-1
        izb = lb-ixb-iyb
        ipb = Ind(lb,ixb,izb)

        if (ixb == 0) then
          do iZeta=1,nZeta
            final(iZeta,ipa,ipb,1) = Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb))
          end do
        else
          do iZeta=1,nZeta
            final(iZeta,ipa,ipb,1) = Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb))- &
                                     dble(ixb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb))
          end do
        end if

        if (iyb == 0) then

          do iZeta=1,nZeta
            final(iZeta,ipa,ipb,2) = Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb))
          end do
        else
          do iZeta=1,nZeta
            final(iZeta,ipa,ipb,2) = Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb))- &
                                     dble(iyb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb))
          end do
        end if

        if (izb == 0) then

          do iZeta=1,nZeta
            final(iZeta,ipa,ipb,3) = Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1))
          end do
        else
          do iZeta=1,nZeta
            final(iZeta,ipa,ipb,3) = Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1))- &
                                     dble(iza)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1))
          end do
        end if

      end do
    end do

  end do
end do

if (iPrint >= 49) then
  write(6,*) ' In UTIL8 la,lb=',la,lb
  do iiComp=1,3
    do jElem=1,nElem(lb)
      do iElem=1,nElem(la)
        do iiZeta=1,nZeta
          write(Label,'(A,I2,A,I2,A)') ' Final (',iElem,',',jElem,') '
          !call RecPrt(Label,'(5f15.8)',
          write(6,*) iiZeta,iElem,jElem,iiComp,final(iiZeta,iElem,jElem,iiComp)
        end do
      end do
    end do
  end do
end if

return

end subroutine Util8
