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

use Index_Functions, only: C_Ind, nTri_Elem1
use Constants, only: Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb, nComp
real(kind=wp), intent(in) :: Alpha(nZeta), Slaplb(nZeta,nTri_Elem1(la+1),nTri_Elem1(lb),nComp), &
                             Slamlb(nZeta,nTri_Elem1(la-1),nTri_Elem1(lb),nComp)
real(kind=wp), intent(out) :: rFinal(nZeta,nTri_Elem1(la),nTri_Elem1(lb),3,nComp)
#include "print.fh"
integer(kind=iwp) :: iComp, ipa, ipb, iPrint, iRout, ixa, ixb, ixm, ixp, iya, iyb, iym, iyp, iza, izb, izm, izp
character(len=80) :: Label

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
    call RecPrt(Label,'(10f15.8)',Slaplb(:,:,:,iComp),nZeta,nTri_Elem1(la+1)*nTri_Elem1(lb))
    if (la > 0) then
      write(Label,'(A,I2,A)') 'Ass_pX: Slamlb(iComp=,',iComp,')'
      call RecPrt(Label,'(10G15.8)',Slamlb(:,:,:,iComp),nZeta,nTri_Elem1(la-1)*nTri_Elem1(lb))
    end if
  end do
end if

do iComp=1,nComp

  do ixa=la,0,-1
    do iya=la-ixa,0,-1
      iza = la-ixa-iya
      ipa = C_Ind(la,ixa,iza)

      do ixb=lb,0,-1
        do iyb=lb-ixb,0,-1
          izb = lb-ixb-iyb
          ipb = C_Ind(lb,ixb,izb)
          !                                                            *
          !*************************************************************
          !                                                            *
          if (ixa == 0) then
            ixp = C_Ind(la+1,ixa+1,iza)
            rFinal(:,ipa,ipb,1,iComp) = Two*Alpha*Slaplb(:,ixp,ipb,iComp)
          else
            ixp = C_Ind(la+1,ixa+1,iza)
            ixm = C_Ind(la-1,ixa-1,iza)
            rFinal(:,ipa,ipb,1,iComp) = Two*Alpha*Slaplb(:,ixp,ipb,iComp)-real(ixa,kind=wp)*Slamlb(:,ixm,ipb,iComp)
          end if

          if (iya == 0) then
            iyp = C_Ind(la+1,ixa,iza)
            rFinal(:,ipa,ipb,2,iComp) = Two*Alpha*Slaplb(:,iyp,ipb,iComp)
          else
            iyp = C_Ind(la+1,ixa,iza)
            iym = C_Ind(la-1,ixa,iza)
            rFinal(:,ipa,ipb,2,iComp) = Two*Alpha*Slaplb(:,iyp,ipb,iComp)-real(iya,kind=wp)*Slamlb(:,iym,ipb,iComp)
          end if

          if (iza == 0) then
            izp = C_Ind(la+1,ixa,iza+1)
            rFinal(:,ipa,ipb,3,iComp) = Two*Alpha*Slaplb(:,izp,ipb,iComp)
          else
            izp = C_Ind(la+1,ixa,iza+1)
            izm = C_Ind(la-1,ixa,iza-1)
            rFinal(:,ipa,ipb,3,iComp) = Two*Alpha*Slaplb(:,izp,ipb,iComp)-real(iza,kind=wp)*Slamlb(:,izm,ipb,iComp)
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
    call RecPrt(Label,'(10G15.8)',rFinal(:,:,:,1,iComp),nZeta,nTri_Elem1(la)*nTri_Elem1(lb))

    write(Label,'(A,I2,A)') ' Ass_pX: pX( 2,iComp=',iComp,')'
    call RecPrt(Label,'(10G15.8)',rFinal(:,:,:,2,iComp),nZeta,nTri_Elem1(la)*nTri_Elem1(lb))

    write(Label,'(A,I2,A)') ' Ass_pX: pX( 3,iComp=',iComp,')'
    call RecPrt(Label,'(10G15.8)',rFinal(:,:,:,3,iComp),nZeta,nTri_Elem1(la)*nTri_Elem1(lb))
  end do
end if

return

end subroutine Ass_pX
