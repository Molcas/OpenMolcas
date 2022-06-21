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

subroutine Ass_pXp(Beta,nZeta,final,la,lb,Slalbp,Slalbm,nComp)
!***********************************************************************
!                                                                      *
! Object: to assemble the pVp integrals                                *
!                                                                      *
!     Author: Bernd Hess, Institut fuer Physikalische und Theoretische *
!             Chemie, University of Bonn, Germany, August 1994         *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
real*8 final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp), Slalbp(nZeta,(la+1)*(la+2)/2,(lb+2)*(lb+3)/2,3,nComp), &
       Slalbm(nZeta,(la+1)*(la+2)/2,lb*(lb+1)/2,3,nComp), Beta(nZeta)
character*80 Label
! Statement function for cartesian index
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1
nElem(ix) = (ix+1)*(ix+2)/2

iRout = 211
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  write(6,*)
  write(6,*) ' In Ass_pXp la,lb,nComp,=',la,lb,nComp
  write(6,*)
  call RecPrt('Beta','(10G15.8)',Beta,nZeta,1)
  do iComp=1,nComp
    write(6,*) 'iComp=',iComp
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
            final(iZeta,ipa,ipb,iComp) = Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb+1,izb),1,iComp)+ &
                                         Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb),2,iComp)+ &
                                         Two*Beta(iZeta)*Slalbp(iZeta,ipa,Ind(lb+1,ixb,izb+1),3,iComp)
          end do

          if (ixb > 0) then
            do iZeta=1,nZeta
              final(iZeta,ipa,ipb,iComp) = final(iZeta,ipa,ipb,iComp)-dble(ixb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb-1,izb),1,iComp)
            end do
          end if

          if (iyb > 0) then
            do iZeta=1,nZeta
              final(iZeta,ipa,ipb,iComp) = final(iZeta,ipa,ipb,iComp)-dble(iyb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb),2,iComp)
            end do
          end if

          if (izb > 0) then
            do iZeta=1,nZeta
              final(iZeta,ipa,ipb,iComp) = final(iZeta,ipa,ipb,iComp)-dble(izb)*Slalbm(iZeta,ipa,Ind(lb-1,ixb,izb-1),3,iComp)
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
    call RecPrt(Label,'(10G15.8)',final(1,1,1,iComp),nZeta,nElem(la)*nElem(lb))
  end do
end if

return

end subroutine Ass_pXp
