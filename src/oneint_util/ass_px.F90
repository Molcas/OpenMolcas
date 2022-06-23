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

subroutine Ass_pX(Alpha,nZeta,rFinal,la,lb,Slaplb,Slamlb,nComp)
!***********************************************************************
!                                                                      *
! Object: to assemble the pV integrals from                            *
!         derivative integrals of the electric potential.              *
!                                                                      *
!     Author: Bernd Hess, Institut fuer Physikalische und Theoretische *
!             Chemie, University of Bonn, Germany, August 1994         *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Constants, only: Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nZeta, la, lb, nComp
real(kind=wp) :: Alpha(nZeta), rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),3,nComp),  &
                 Slaplb(nZeta,nTri_Elem1(la+1),nTri_Elem1(lb),nComp), Slamlb(nZeta,nTri_Elem1(la-1),nTri_Elem1(lb),nComp)
#include "print.fh"
integer(kind=iwp) :: iComp, ipa, ipb, iPrint, iRout, ixa, ixb, ixm, ixp, iya, iyb, iym, iyp, iza, izb, iZeta, izm, izp
character(len=80) :: Label
! Statement function for cartesian index
integer(kind=iwp) :: Ind, nElem, ixyz, ix, iz
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1
nElem(ix) = (ix+1)*(ix+2)/2

iRout = 203
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  write(u6,*)
  write(u6,*) ' In Ass_pX la,lb,nComp=',la,lb,nComp
  write(u6,*)
  call RecPrt('Alpha','(10G15.8)',Alpha,nZeta,1)
  do iComp=1,nComp
    write(u6,*)
    write(u6,*) 'iComp=',iComp
    write(u6,*)
    write(Label,'(A,I2,A)') 'Ass_pX:  Slaplb(iComp=',iComp,')'
    call RecPrt(Label,'(10f15.8)',Slaplb(1,1,1,iComp),nZeta,nElem(la+1)*nElem(lb))
    if (la > 0) then
      write(Label,'(A,I2,A)') 'Ass_pX: Slamlb(iComp=,',iComp,')'
      call RecPrt(Label,'(10G15.8)',Slamlb(1,1,1,iComp),nZeta,nElem(la-1)*nElem(lb))
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
          !                                                            *
          !*************************************************************
          !                                                            *
          if (ixa == 0) then
            ixp = Ind(la+1,ixa+1,iza)
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,1,iComp) = Two*Alpha(iZeta)*Slaplb(iZeta,ixp,ipb,iComp)
            end do
          else
            ixp = Ind(la+1,ixa+1,iza)
            ixm = Ind(la-1,ixa-1,iza)
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,1,iComp) = Two*Alpha(iZeta)*Slaplb(iZeta,ixp,ipb,iComp)- &
                                              real(ixa,kind=wp)*Slamlb(iZeta,ixm,ipb,iComp)
            end do
          end if

          if (iya == 0) then
            iyp = Ind(la+1,ixa,iza)
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,2,iComp) = Two*Alpha(iZeta)*Slaplb(iZeta,iyp,ipb,iComp)
            end do
          else
            iyp = Ind(la+1,ixa,iza)
            iym = Ind(la-1,ixa,iza)
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,2,iComp) = Two*Alpha(iZeta)*Slaplb(iZeta,iyp,ipb,iComp)- &
                                              real(iya,kind=wp)*Slamlb(iZeta,iym,ipb,iComp)
            end do
          end if

          if (iza == 0) then
            izp = Ind(la+1,ixa,iza+1)
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,3,iComp) = Two*Alpha(iZeta)*Slaplb(iZeta,izp,ipb,iComp)
            end do
          else
            izp = Ind(la+1,ixa,iza+1)
            izm = Ind(la-1,ixa,iza-1)
            do iZeta=1,nZeta
              rFinal(iZeta,ipa,ipb,3,iComp) = Two*Alpha(iZeta)*Slaplb(iZeta,izp,ipb,iComp)- &
                                              real(iza,kind=wp)*Slamlb(iZeta,izm,ipb,iComp)
            end do
          end if
          !                                                            *
          !*************************************************************
          !                                                            *
        end do
      end do

    end do
  end do

end do

if (iPrint >= 49) then
  write(u6,*) ' In Ass_pX la,lb,nComp=',la,lb,nComp
  do iComp=1,nComp
    write(u6,*)
    write(u6,*) 'iComp=',iComp
    write(u6,*)

    write(Label,'(A,I2,A)') ' Ass_pX: pX( 1,iComp=',iComp,')'
    call RecPrt(Label,'(10G15.8)',rFinal(1,1,1,1,iComp),nZeta,nElem(la)*nElem(lb))

    write(Label,'(A,I2,A)') ' Ass_pX: pX( 2,iComp=',iComp,')'
    call RecPrt(Label,'(10G15.8)',rFinal(1,1,1,2,iComp),nZeta,nElem(la)*nElem(lb))

    write(Label,'(A,I2,A)') ' Ass_pX: pX( 3,iComp=',iComp,')'
    call RecPrt(Label,'(10G15.8)',rFinal(1,1,1,3,iComp),nZeta,nElem(la)*nElem(lb))
  end do
end if

return

end subroutine Ass_pX
