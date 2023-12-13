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
! Copyright (C) 1996, Per Ake Malmqvist                                *
!***********************************************************************

subroutine AMPr(Beta,nZeta,Rslt,la,lb,Tabpp,Tabp,Tab0,Tabm,Tabmm)
!***********************************************************************
!                                                                      *
! Object: Compute matrix elements of hermitized products of angular    *
!         moment operators, using elementary overlaps, dipole, and     *
!         quadrupole integrals.                                        *
!                                                                      *
!     Author: Per-AAke Malmqvist, Dept. of Theoretical Chemistry,      *
!             University of Lund, SWEDEN                               *
!             November '96                                             *
!***********************************************************************

use Index_Functions, only: C_Ind3, nTri_Elem1
use Constants, only: Two, Four, Eight, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nZeta, la, lb
real(kind=wp), intent(in) :: Beta(nZeta), Tabpp(nZeta,nTri_Elem1(la),nTri_Elem1(lb+2),6), &
                             Tabp(nZeta,nTri_Elem1(la),nTri_Elem1(lb+1),3), Tab0(nZeta,nTri_Elem1(la),nTri_Elem1(lb),6), &
                             Tabm(nZeta,nTri_Elem1(la),nTri_Elem1(lb-1),3), Tabmm(nZeta,nTri_Elem1(la),nTri_Elem1(lb-2),6)
real(kind=wp), intent(out) :: Rslt(nZeta,nTri_Elem1(la),nTri_Elem1(lb),6)
#include "print.fh"
integer(kind=iwp) :: ia, ib, iElem, ipa, ipb, iPrint, iRout, ix, ixa, ixb, iy, iya, iyb, iz, iza, izb, iZeta, jElem
real(kind=wp) :: B, B2, Bx2, By2, Bz2, Term1, Term2, Term3, Term4, Term5, Term6
character(len=80) :: Label
integer(kind=iwp), parameter :: kx = 1, ky = 2, kz = 3, &
                                kxx = 1, kxy = 2, kxz = 3, &
                                kyx = 2, kyy = 4, kyz = 5, &
                                kzx = 3, kzy = 5, kzz = 6

iRout = 221
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  write(u6,*) ' In AMPr la,lb=',la,lb
  call RecPrt('Beta',' ',Beta,nZeta,1)
  do ia=1,nTri_Elem1(la)
    do ib=1,nTri_Elem1(lb+2)
      write(Label,'(A,I2,A,I2,A)') ' Tabpp(',ia,',',ib,'xx)'
      call RecPrt(Label,' ',Tabpp(:,ia,ib,1),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tabpp(',ia,',',ib,'xy)'
      call RecPrt(Label,' ',Tabpp(:,ia,ib,2),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tabpp(',ia,',',ib,'xz)'
      call RecPrt(Label,' ',Tabpp(:,ia,ib,3),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tabpp(',ia,',',ib,'yy)'
      call RecPrt(Label,' ',Tabpp(:,ia,ib,4),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tabpp(',ia,',',ib,'yz)'
      call RecPrt(Label,' ',Tabpp(:,ia,ib,5),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tabpp(',ia,',',ib,'zz)'
      call RecPrt(Label,' ',Tabpp(:,ia,ib,6),nZeta,1)
    end do
  end do
  do ia=1,nTri_Elem1(la)
    do ib=1,nTri_Elem1(lb+1)
      write(Label,'(A,I2,A,I2,A)') ' Tabp(',ia,',',ib,'x)'
      call RecPrt(Label,' ',Tabp(:,ia,ib,1),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tabp(',ia,',',ib,'y)'
      call RecPrt(Label,' ',Tabp(:,ia,ib,2),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tabp(',ia,',',ib,'z)'
      call RecPrt(Label,' ',Tabp(:,ia,ib,3),nZeta,1)
    end do
  end do
  do ia=1,nTri_Elem1(la)
    do ib=1,nTri_Elem1(lb)
      write(Label,'(A,I2,A,I2,A)') ' Tab0',ia,',',ib,'xx)'
      call RecPrt(Label,' ',Tab0(:,ia,ib,1),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tab0',ia,',',ib,'xy)'
      call RecPrt(Label,' ',Tab0(:,ia,ib,2),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tab0',ia,',',ib,'xz)'
      call RecPrt(Label,' ',Tab0(:,ia,ib,3),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tab0',ia,',',ib,'yy)'
      call RecPrt(Label,' ',Tab0(:,ia,ib,4),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tab0',ia,',',ib,'yz)'
      call RecPrt(Label,' ',Tab0(:,ia,ib,5),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tab0',ia,',',ib,'zz)'
      call RecPrt(Label,' ',Tab0(:,ia,ib,6),nZeta,1)
    end do
  end do
  if (lb > 0) then
    do ia=1,nTri_Elem1(la)
      do ib=1,nTri_Elem1(lb-1)
        write(Label,'(A,I2,A,I2,A)') ' Tabm(',ia,',',ib,'x)'
        call RecPrt(Label,' ',Tabm(:,ia,ib,1),nZeta,1)
        write(Label,'(A,I2,A,I2,A)') ' Tabm(',ia,',',ib,'y)'
        call RecPrt(Label,' ',Tabm(:,ia,ib,2),nZeta,1)
        write(Label,'(A,I2,A,I2,A)') ' Tabm(',ia,',',ib,'z)'
        call RecPrt(Label,' ',Tabm(:,ia,ib,3),nZeta,1)
      end do
    end do
  end if
  if (lb > 1) then
    do ia=1,nTri_Elem1(la)
      do ib=1,nTri_Elem1(lb-2)
        write(Label,'(A,I2,A,I2,A)') ' Tabmm(',ia,',',ib,'xx)'
        call RecPrt(Label,' ',Tabmm(:,ia,ib,1),nZeta,1)
        write(Label,'(A,I2,A,I2,A)') ' Tabmm(',ia,',',ib,'xy)'
        call RecPrt(Label,' ',Tabmm(:,ia,ib,2),nZeta,1)
        write(Label,'(A,I2,A,I2,A)') ' Tabmm(',ia,',',ib,'xz)'
        call RecPrt(Label,' ',Tabmm(:,ia,ib,3),nZeta,1)
        write(Label,'(A,I2,A,I2,A)') ' Tabmm(',ia,',',ib,'yy)'
        call RecPrt(Label,' ',Tabmm(:,ia,ib,4),nZeta,1)
        write(Label,'(A,I2,A,I2,A)') ' Tabmm(',ia,',',ib,'yz)'
        call RecPrt(Label,' ',Tabmm(:,ia,ib,5),nZeta,1)
        write(Label,'(A,I2,A,I2,A)') ' Tabmm(',ia,',',ib,'zz)'
        call RecPrt(Label,' ',Tabmm(:,ia,ib,6),nZeta,1)
      end do
    end do
  end if
end if

do ixa=la,0,-1
  do iya=la-ixa,0,-1
    iza = la-ixa-iya
    ipa = C_Ind3(ixa,iya,iza)

    do ixb=lb,0,-1
      do iyb=lb-ixb,0,-1
        izb = lb-ixb-iyb
        ipb = C_Ind3(ixb,iyb,izb)
        ix = ixb
        iy = iyb
        iz = izb

        do iZeta=1,nZeta
          B = Beta(iZeta)
          B2 = B**2
          Bx2 = real(2*ix,kind=wp)*B
          By2 = real(2*iy,kind=wp)*B
          Bz2 = real(2*iz,kind=wp)*B
          ! First compute Lx**2, Ly**2, and Lz**2:
          !------------------
          Term1 = Two*Bz2*Tab0(iZeta,ipa,C_Ind3(ix,iy,iz),kyy)-Four*B2*Tabpp(iZeta,ipa,C_Ind3(ix,iy,iz+2),kyy)+ &
                  Two*B*Tab0(iZeta,ipa,C_Ind3(ix,iy,iz),kyy)
          Term2 = -Two*By2*Tab0(iZeta,ipa,C_Ind3(ix,iy-1,iz+1),kyz)-Two*Bz2*Tab0(iZeta,ipa,C_Ind3(ix,iy+1,iz-1),kyz)+ &
                  Eight*B2*Tabpp(iZeta,ipa,C_Ind3(ix,iy+1,iz+1),kyz)
          Term3 = Two*By2*Tab0(iZeta,ipa,C_Ind3(ix,iy,iz),kzz)-Four*B2*Tabpp(iZeta,ipa,C_Ind3(ix,iy+2,iz),kzz)+ &
                  Two*B*Tab0(iZeta,ipa,C_Ind3(ix,iy,iz),kzz)
          Term4 = -Two*B*Tabp(iZeta,ipa,C_Ind3(ix,iy+1,iz),ky)-Two*B*Tabp(iZeta,ipa,C_Ind3(ix,iy,iz+1),kz)
          if (lb >= 1) then
            Term4 = Term4+real(iy,kind=wp)*Tabm(iZeta,ipa,C_Ind3(ix,iy-1,iz),ky)+ &
                    real(iz,kind=wp)*Tabm(iZeta,ipa,C_Ind3(ix,iy,iz-1),kz)
            if (lb >= 2) then
              Term1 = Term1-real(iz*(iz-1),kind=wp)*Tabmm(iZeta,ipa,C_Ind3(ix,iy,iz-2),kyy)
              Term2 = Term2+real(2*iy*iz,kind=wp)*Tabmm(iZeta,ipa,C_Ind3(ix,iy-1,iz-1),kyz)
              Term3 = Term3-real(iy*(iy-1),kind=wp)*Tabmm(iZeta,ipa,C_Ind3(ix,iy-2,iz),kzz)
            end if
          end if
          Rslt(iZeta,ipa,ipb,kxx) = Term1+Term2+Term3+Term4
          !------------------
          Term1 = Two*Bx2*Tab0(iZeta,ipa,C_Ind3(ix,iy,iz),kzz)-Four*B2*Tabpp(iZeta,ipa,C_Ind3(ix+2,iy,iz),kzz)+ &
                  Two*B*Tab0(iZeta,ipa,C_Ind3(ix,iy,iz),kzz)
          Term2 = -Two*Bz2*Tab0(iZeta,ipa,C_Ind3(ix+1,iy,iz-1),kzx)-Two*Bx2*Tab0(iZeta,ipa,C_Ind3(ix-1,iy,iz+1),kzx)+ &
                  Eight*B2*Tabpp(iZeta,ipa,C_Ind3(ix+1,iy,iz+1),kzx)
          Term3 = Two*Bz2*Tab0(iZeta,ipa,C_Ind3(ix,iy,iz),kxx)-Four*B2*Tabpp(iZeta,ipa,C_Ind3(ix,iy,iz+2),kxx)+ &
                  Two*B*Tab0(iZeta,ipa,C_Ind3(ix,iy,iz),kxx)
          Term4 = -Two*B*Tabp(iZeta,ipa,C_Ind3(ix,iy,iz+1),kz)-Two*B*Tabp(iZeta,ipa,C_Ind3(ix+1,iy,iz),kx)
          if (lb >= 1) then
            Term4 = Term4+real(iz,kind=wp)*Tabm(iZeta,ipa,C_Ind3(ix,iy,iz-1),kz)+ &
                    real(ix,kind=wp)*Tabm(iZeta,ipa,C_Ind3(ix-1,iy,iz),kx)
            if (lb >= 2) then
              Term1 = Term1-real(ix*(ix-1),kind=wp)*Tabmm(iZeta,ipa,C_Ind3(ix-2,iy,iz),kzz)
              Term2 = Term2+real(2*iz*ix,kind=wp)*Tabmm(iZeta,ipa,C_Ind3(ix-1,iy,iz-1),kzx)
              Term3 = Term3-real(iz*(iz-1),kind=wp)*Tabmm(iZeta,ipa,C_Ind3(ix,iy,iz-2),kxx)
            end if
          end if
          Rslt(iZeta,ipa,ipb,kyy) = Term1+Term2+Term3+Term4
          !------------------
          Term1 = Two*By2*Tab0(iZeta,ipa,C_Ind3(ix,iy,iz),kxx)-Four*B2*Tabpp(iZeta,ipa,C_Ind3(ix,iy+2,iz),kxx)+ &
                  Two*B*Tab0(iZeta,ipa,C_Ind3(ix,iy,iz),kxx)
          Term2 = -Two*Bx2*Tab0(iZeta,ipa,C_Ind3(ix-1,iy+1,iz),kxy)-Two*By2*Tab0(iZeta,ipa,C_Ind3(ix+1,iy-1,iz),kxy)+ &
                  Eight*B2*Tabpp(iZeta,ipa,C_Ind3(ix+1,iy+1,iz),kxy)
          Term3 = Two*Bx2*Tab0(iZeta,ipa,C_Ind3(ix,iy,iz),kyy)-Four*B2*Tabpp(iZeta,ipa,C_Ind3(ix+2,iy,iz),kyy)+ &
                  Two*B*Tab0(iZeta,ipa,C_Ind3(ix,iy,iz),kyy)
          Term4 = -Two*B*Tabp(iZeta,ipa,C_Ind3(ix+1,iy,iz),kx)-Two*B*Tabp(iZeta,ipa,C_Ind3(ix,iy+1,iz),ky)
          if (lb >= 1) then
            Term4 = Term4+real(ix,kind=wp)*Tabm(iZeta,ipa,C_Ind3(ix-1,iy,iz),kx)+ &
                    real(iy,kind=wp)*Tabm(iZeta,ipa,C_Ind3(ix,iy-1,iz),ky)
            if (lb >= 2) then
              Term1 = Term1-real(iy*(iy-1),kind=wp)*Tabmm(iZeta,ipa,C_Ind3(ix,iy-2,iz),kxx)
              Term2 = Term2+real(2*ix*iy,kind=wp)*Tabmm(iZeta,ipa,C_Ind3(ix-1,iy-1,iz),kxy)
              Term3 = Term3-real(ix*(ix-1),kind=wp)*Tabmm(iZeta,ipa,C_Ind3(ix-2,iy,iz),kyy)
            end if
          end if
          Rslt(iZeta,ipa,ipb,kzz) = Term1+Term2+Term3+Term4
          !------------------
          ! Compute (Lx*Ly+Ly*Lx)/2, etc cyclical.
          ! With Term5, (Lx*Ly) obtains. With Term6, (Ly*Lx).
          ! We want the hermitian average.
          !------------------ (Lx,Ly)
          Term1 = Bx2*Tab0(iZeta,ipa,C_Ind3(ix-1,iy,iz+1),kyz)+Bz2*Tab0(iZeta,ipa,C_Ind3(ix+1,iy,iz-1),kyz)- &
                  Four*B2*Tabpp(iZeta,ipa,C_Ind3(ix+1,iy,iz+1),kyz)
          Term2 = -Two*Bz2*Tab0(iZeta,ipa,C_Ind3(ix,iy,iz),kxy)+Four*B2*Tabpp(iZeta,ipa,C_Ind3(ix,iy,iz+2),kxy)- &
                  Two*B*Tab0(iZeta,ipa,C_Ind3(ix,iy,iz),kxy)
          Term3 = -Bx2*Tab0(iZeta,ipa,C_Ind3(ix-1,iy+1,iz),kzz)-By2*Tab0(iZeta,ipa,C_Ind3(ix+1,iy-1,iz),kzz)+ &
                  Four*B2*Tabpp(iZeta,ipa,C_Ind3(ix+1,iy+1,iz),kzz)
          Term4 = By2*Tab0(iZeta,ipa,C_Ind3(ix,iy-1,iz+1),kxz)+Bz2*Tab0(iZeta,ipa,C_Ind3(ix,iy+1,iz-1),kxz)- &
                  Four*B2*Tabpp(iZeta,ipa,C_Ind3(ix,iy+1,iz+1),kxz)
          Term5 = Two*B*Tabp(iZeta,ipa,C_Ind3(ix+1,iy,iz),ky)
          Term6 = Two*B*Tabp(iZeta,ipa,C_Ind3(ix,iy+1,iz),kx)
          if (lb >= 1) then
            Term5 = Term5-real(ix,kind=wp)*Tabm(iZeta,ipa,C_Ind3(ix-1,iy,iz),ky)
            Term6 = Term6-real(iy,kind=wp)*Tabm(iZeta,ipa,C_Ind3(ix,iy-1,iz),kx)
            if (lb >= 2) then
              Term1 = Term1-real(ix*iz,kind=wp)*Tabmm(iZeta,ipa,C_Ind3(ix-1,iy,iz-1),kyz)
              Term2 = Term2+real(iz*(iz-1),kind=wp)*Tabmm(iZeta,ipa,C_Ind3(ix,iy,iz-2),kxy)
              Term3 = Term3+real(ix*iy,kind=wp)*Tabmm(iZeta,ipa,C_Ind3(ix-1,iy-1,iz),kzz)
              Term4 = Term4-real(iy*iz,kind=wp)*Tabmm(iZeta,ipa,C_Ind3(ix,iy-1,iz-1),kxz)
            end if
          end if
          Rslt(iZeta,ipa,ipb,kxy) = Term1+Term2+Term3+Term4+Half*(Term5+Term6)
          ! (Ly,Lz)
          Term1 = By2*Tab0(iZeta,ipa,C_Ind3(ix+1,iy-1,iz),kzx)+Bx2*Tab0(iZeta,ipa,C_Ind3(ix-1,iy+1,iz),kzx)- &
                  Four*B2*Tabpp(iZeta,ipa,C_Ind3(ix+1,iy+1,iz),kzx)
          Term2 = -Two*Bx2*Tab0(iZeta,ipa,C_Ind3(ix,iy,iz),kyz)+Four*B2*Tabpp(iZeta,ipa,C_Ind3(ix+2,iy,iz),kyz)- &
                  Two*B*Tab0(iZeta,ipa,C_Ind3(ix,iy,iz),kyz)
          Term3 = -By2*Tab0(iZeta,ipa,C_Ind3(ix,iy-1,iz+1),kxx)-Bz2*Tab0(iZeta,ipa,C_Ind3(ix,iy+1,iz-1),kxx)+ &
                  Four*B2*Tabpp(iZeta,ipa,C_Ind3(ix,iy+1,iz+1),kxx)
          Term4 = Bz2*Tab0(iZeta,ipa,C_Ind3(ix+1,iy,iz-1),kyx)+Bx2*Tab0(iZeta,ipa,C_Ind3(ix-1,iy,iz+1),kyx)- &
                  Four*B2*Tabpp(iZeta,ipa,C_Ind3(ix+1,iy,iz+1),kyx)
          Term5 = Two*B*Tabp(iZeta,ipa,C_Ind3(ix,iy+1,iz),kz)
          Term6 = Two*B*Tabp(iZeta,ipa,C_Ind3(ix,iy,iz+1),ky)
          if (lb >= 1) then
            Term5 = Term5-real(iy,kind=wp)*Tabm(iZeta,ipa,C_Ind3(ix,iy-1,iz),kz)
            Term6 = Term6-real(iz,kind=wp)*Tabm(iZeta,ipa,C_Ind3(ix,iy,iz-1),ky)
            if (lb >= 2) then
              Term1 = Term1-real(iy*ix,kind=wp)*Tabmm(iZeta,ipa,C_Ind3(ix-1,iy-1,iz),kzx)
              Term2 = Term2+real(ix*(ix-1),kind=wp)*Tabmm(iZeta,ipa,C_Ind3(ix-2,iy,iz),kyz)
              Term3 = Term3+real(iy*iz,kind=wp)*Tabmm(iZeta,ipa,C_Ind3(ix,iy-1,iz-1),kxx)
              Term4 = Term4-real(iz*ix,kind=wp)*Tabmm(iZeta,ipa,C_Ind3(ix-1,iy,iz-1),kyx)
            end if
          end if
          Rslt(iZeta,ipa,ipb,kyz) = Term1+Term2+Term3+Term4+Half*(Term5+Term6)
          ! (Lx,Lz)
          Term1 = Bz2*Tab0(iZeta,ipa,C_Ind3(ix,iy+1,iz-1),kxy)+By2*Tab0(iZeta,ipa,C_Ind3(ix,iy-1,iz+1),kxy)- &
                  Four*B2*Tabpp(iZeta,ipa,C_Ind3(ix,iy+1,iz+1),kxy)
          Term2 = -Two*By2*Tab0(iZeta,ipa,C_Ind3(ix,iy,iz),kzx)+Four*B2*Tabpp(iZeta,ipa,C_Ind3(ix,iy+2,iz),kzx)- &
                  Two*B*Tab0(iZeta,ipa,C_Ind3(ix,iy,iz),kzx)
          Term3 = -Bz2*Tab0(iZeta,ipa,C_Ind3(ix+1,iy,iz-1),kyy)-Bx2*Tab0(iZeta,ipa,C_Ind3(ix-1,iy,iz+1),kyy)+ &
                  Four*B2*Tabpp(iZeta,ipa,C_Ind3(ix+1,iy,iz+1),kyy)
          Term4 = Bx2*Tab0(iZeta,ipa,C_Ind3(ix-1,iy+1,iz),kzy)+By2*Tab0(iZeta,ipa,C_Ind3(ix+1,iy-1,iz),kzy)- &
                  Four*B2*Tabpp(iZeta,ipa,C_Ind3(ix+1,iy+1,iz),kzy)
          Term5 = Two*B*Tabp(iZeta,ipa,C_Ind3(ix,iy,iz+1),kx)
          Term6 = Two*B*Tabp(iZeta,ipa,C_Ind3(ix+1,iy,iz),kz)
          if (lb >= 1) then
            Term5 = Term5-real(iz,kind=wp)*Tabm(iZeta,ipa,C_Ind3(ix,iy,iz-1),kx)
            Term6 = Term6-real(ix,kind=wp)*Tabm(iZeta,ipa,C_Ind3(ix-1,iy,iz),kz)
            if (lb >= 2) then
              Term1 = Term1-real(iz*iy,kind=wp)*Tabmm(iZeta,ipa,C_Ind3(ix,iy-1,iz-1),kxy)
              Term2 = Term2+real(iy*(iy-1),kind=wp)*Tabmm(iZeta,ipa,C_Ind3(ix,iy-2,iz),kzx)
              Term3 = Term3+real(iz*ix,kind=wp)*Tabmm(iZeta,ipa,C_Ind3(ix-1,iy,iz-1),kyy)
              Term4 = Term4-real(ix*iy,kind=wp)*Tabmm(iZeta,ipa,C_Ind3(ix-1,iy-1,iz),kzy)
            end if
          end if
          Rslt(iZeta,ipa,ipb,kxz) = Term1+Term2+Term3+Term4+Half*(Term5+Term6)

        end do

      end do
    end do

  end do
end do

if (iPrint >= 49) then
  write(u6,*) ' In AMPr la,lb=',la,lb
  do iElem=1,nTri_Elem1(la)
    do jElem=1,nTri_Elem1(lb)
      write(Label,'(A,I2,A,I2,A)') ' Rslt (',iElem,',',jElem,',xx) '
      call RecPrt(Label,' ',Rslt(:,iElem,jElem,kxx),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Rslt (',iElem,',',jElem,',xy) '
      call RecPrt(Label,' ',Rslt(:,iElem,jElem,kxy),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Rslt (',iElem,',',jElem,',xz) '
      call RecPrt(Label,' ',Rslt(:,iElem,jElem,kxz),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Rslt (',iElem,',',jElem,',yy) '
      call RecPrt(Label,' ',Rslt(:,iElem,jElem,kyy),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Rslt (',iElem,',',jElem,',yz) '
      call RecPrt(Label,' ',Rslt(:,iElem,jElem,kyz),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Rslt (',iElem,',',jElem,',zz) '
      call RecPrt(Label,' ',Rslt(:,iElem,jElem,kzz),nZeta,1)
    end do
  end do
  write(u6,*) ' Leaving AMPr.'
end if

return

end subroutine AMPr
