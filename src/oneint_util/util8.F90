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

subroutine Util8(Beta,nZeta,rFinal,la,lb,Slalbp,Slalbm)
!***********************************************************************
!                                                                      *
! Object: to assemble the Vp integrals from                            *
!         derivative integrals of the electric potential.              *
!                                                                      *
!     Author: Bernd Hess, Institut fuer Physikalische und Theoretische *
!             Chemie, University of Bonn, Germany, August 1994         *
!***********************************************************************

use Index_Functions, only: C_Ind, nTri_Elem1
use Constants, only: Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb
real(kind=wp), intent(in) :: Beta(nZeta), Slalbp(nZeta,nTri_Elem1(la),nTri_Elem1(lb+1)), &
                             Slalbm(nZeta,nTri_Elem1(la),nTri_Elem1(lb-1))
real(kind=wp), intent(out) :: rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),3)
#include "print.fh"
integer(kind=iwp) :: ia, ib, iElem, iiComp, iiZeta, ipa, ipb, iPrint, iRout, ixa, ixb, iya, iyb, iza, izb, jElem
character(len=80) :: Label

iRout = 203
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  write(u6,*) ' In util8 la,lb=',la,lb
  call RecPrt('Beta','(5f15.8)',Beta,nZeta,1)
  do ib=1,nTri_Elem1(lb)
    write(Label,'(A,I2,A)') ' Slalbp(',la,ib,')'
    call RecPrt(Label,'(5f15.8)',Slalbp(:,:,ib),nZeta,nTri_Elem1(la+1))
  end do
  if (lb > 0) then
    do ia=1,nTri_Elem1(la)
      write(Label,'(A,I2,A)') ' Slalbm(',la,ib,')'
      call RecPrt(Label,'(5f15.8)',Slalbm(:,:,ib),nZeta,nTri_Elem1(lb-1))
    end do
  end if
end if

do ixa=la,0,-1
  do iya=la-ixa,0,-1
    iza = la-ixa-iya
    ipa = C_Ind(la,ixa,iza)

    do ixb=lb,0,-1
      do iyb=lb-ixb,0,-1
        izb = lb-ixb-iyb
        ipb = C_Ind(lb,ixb,izb)

        if (ixb == 0) then
          rFinal(:,ipa,ipb,1) = Two*Beta*Slalbp(:,ipa,C_Ind(lb+1,ixb+1,izb))
        else
          rFinal(:,ipa,ipb,1) = Two*Beta*Slalbp(:,ipa,C_Ind(lb+1,ixb+1,izb))-real(ixb,kind=wp)*Slalbm(:,ipa,C_Ind(lb-1,ixb-1,izb))
        end if

        if (iyb == 0) then
          rFinal(:,ipa,ipb,2) = Two*Beta*Slalbp(:,ipa,C_Ind(lb+1,ixb,izb))
        else
          rFinal(:,ipa,ipb,2) = Two*Beta*Slalbp(:,ipa,C_Ind(lb+1,ixb,izb))-real(iyb,kind=wp)*Slalbm(:,ipa,C_Ind(lb-1,ixb,izb))
        end if

        if (izb == 0) then
          rFinal(:,ipa,ipb,3) = Two*Beta*Slalbp(:,ipa,C_Ind(lb+1,ixb,izb+1))
        else
          rFinal(:,ipa,ipb,3) = Two*Beta*Slalbp(:,ipa,C_Ind(lb+1,ixb,izb+1))-real(iza,kind=wp)*Slalbm(:,ipa,C_Ind(lb-1,ixb,izb-1))
        end if

      end do
    end do

  end do
end do

if (iPrint >= 49) then
  write(u6,*) ' In UTIL8 la,lb=',la,lb
  do iiComp=1,3
    do jElem=1,nTri_Elem1(lb)
      do iElem=1,nTri_Elem1(la)
        do iiZeta=1,nZeta
          write(Label,'(A,I2,A,I2,A)') ' rFinal (',iElem,',',jElem,') '
          !call RecPrt(Label,'(5f15.8)',
          write(u6,*) iiZeta,iElem,jElem,iiComp,rFinal(iiZeta,iElem,jElem,iiComp)
        end do
      end do
    end do
  end do
end if

return

end subroutine Util8
