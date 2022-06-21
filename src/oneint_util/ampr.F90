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

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
real*8 Rslt(nZeta,((la+1)*(la+2))/2,((lb+1)*(lb+2))/2,6), Tabpp(nZeta,((la+1)*(la+2))/2,((lb+3)*(lb+4))/2,6), &
       Tabp(nZeta,((la+1)*(la+2))/2,((lb+2)*(lb+3))/2,3), Tab0(nZeta,((la+1)*(la+2))/2,((lb+1)*(lb+2))/2,6), &
       Tabm(nZeta,((la+1)*(la+2))/2,(lb*(lb+1))/2,3), Tabmm(nZeta,((la+1)*(la+2))/2,((lb-1)*lb)/2,6), Beta(nZeta)
character*80 Label
data kx,ky,kz/1,2,3/
data kxx,kxy,kxz,kyy,kyz,kzz/1,2,3,4,5,6/
data kyx,kzx,kzy/2,3,5/
! Statement function for cartesian index
Ind(j,k) = ((j+k)*(j+k+1))/2+k+1
nElem(ix) = (ix+1)*(ix+2)/2

iRout = 221
iPrint = nPrint(iRout)

if (iPrint >= 99) then
  write(6,*) ' In AMPr la,lb=',la,lb
  call RecPrt('Beta',' ',Beta,nZeta,1)
  do ia=1,nElem(la)
    do ib=1,nElem(lb+2)
      write(Label,'(A,I2,A,I2,A)') ' Tabpp(',ia,',',ib,'xx)'
      call RecPrt(Label,' ',Tabpp(1,ia,ib,1),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tabpp(',ia,',',ib,'xy)'
      call RecPrt(Label,' ',Tabpp(1,ia,ib,2),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tabpp(',ia,',',ib,'xz)'
      call RecPrt(Label,' ',Tabpp(1,ia,ib,3),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tabpp(',ia,',',ib,'yy)'
      call RecPrt(Label,' ',Tabpp(1,ia,ib,4),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tabpp(',ia,',',ib,'yz)'
      call RecPrt(Label,' ',Tabpp(1,ia,ib,5),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tabpp(',ia,',',ib,'zz)'
      call RecPrt(Label,' ',Tabpp(1,ia,ib,6),nZeta,1)
    end do
  end do
  do ia=1,nElem(la)
    do ib=1,nElem(lb+1)
      write(Label,'(A,I2,A,I2,A)') ' Tabp(',ia,',',ib,'x)'
      call RecPrt(Label,' ',Tabp(1,ia,ib,1),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tabp(',ia,',',ib,'y)'
      call RecPrt(Label,' ',Tabp(1,ia,ib,2),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tabp(',ia,',',ib,'z)'
      call RecPrt(Label,' ',Tabp(1,ia,ib,3),nZeta,1)
    end do
  end do
  do ia=1,nElem(la)
    do ib=1,nElem(lb)
      write(Label,'(A,I2,A,I2,A)') ' Tab0',ia,',',ib,'xx)'
      call RecPrt(Label,' ',Tab0(1,ia,ib,1),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tab0',ia,',',ib,'xy)'
      call RecPrt(Label,' ',Tab0(1,ia,ib,2),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tab0',ia,',',ib,'xz)'
      call RecPrt(Label,' ',Tab0(1,ia,ib,3),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tab0',ia,',',ib,'yy)'
      call RecPrt(Label,' ',Tab0(1,ia,ib,4),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tab0',ia,',',ib,'yz)'
      call RecPrt(Label,' ',Tab0(1,ia,ib,5),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Tab0',ia,',',ib,'zz)'
      call RecPrt(Label,' ',Tab0(1,ia,ib,6),nZeta,1)
    end do
  end do
  if (lb > 0) then
    do ia=1,nElem(la)
      do ib=1,nElem(lb-1)
        write(Label,'(A,I2,A,I2,A)') ' Tabm(',ia,',',ib,'x)'
        call RecPrt(Label,' ',Tabm(1,ia,ib,1),nZeta,1)
        write(Label,'(A,I2,A,I2,A)') ' Tabm(',ia,',',ib,'y)'
        call RecPrt(Label,' ',Tabm(1,ia,ib,2),nZeta,1)
        write(Label,'(A,I2,A,I2,A)') ' Tabm(',ia,',',ib,'z)'
        call RecPrt(Label,' ',Tabm(1,ia,ib,3),nZeta,1)
      end do
    end do
  end if
  if (lb > 1) then
    do ia=1,nElem(la)
      do ib=1,nElem(lb-2)
        write(Label,'(A,I2,A,I2,A)') ' Tabmm(',ia,',',ib,'xx)'
        call RecPrt(Label,' ',Tabmm(1,ia,ib,1),nZeta,1)
        write(Label,'(A,I2,A,I2,A)') ' Tabmm(',ia,',',ib,'xy)'
        call RecPrt(Label,' ',Tabmm(1,ia,ib,2),nZeta,1)
        write(Label,'(A,I2,A,I2,A)') ' Tabmm(',ia,',',ib,'xz)'
        call RecPrt(Label,' ',Tabmm(1,ia,ib,3),nZeta,1)
        write(Label,'(A,I2,A,I2,A)') ' Tabmm(',ia,',',ib,'yy)'
        call RecPrt(Label,' ',Tabmm(1,ia,ib,4),nZeta,1)
        write(Label,'(A,I2,A,I2,A)') ' Tabmm(',ia,',',ib,'yz)'
        call RecPrt(Label,' ',Tabmm(1,ia,ib,5),nZeta,1)
        write(Label,'(A,I2,A,I2,A)') ' Tabmm(',ia,',',ib,'zz)'
        call RecPrt(Label,' ',Tabmm(1,ia,ib,6),nZeta,1)
      end do
    end do
  end if
end if

do ixa=la,0,-1
  do iya=la-ixa,0,-1
    iza = la-ixa-iya
    ipa = Ind(iya,iza)

    do ixb=lb,0,-1
      do iyb=lb-ixb,0,-1
        izb = lb-ixb-iyb
        ipb = Ind(iyb,izb)
        ix = ixb
        iy = iyb
        iz = izb

        do iZeta=1,nZeta
          B = Beta(iZeta)
          B2 = B**2
          Bx2 = dble(2*ix)*B
          By2 = dble(2*iy)*B
          Bz2 = dble(2*iz)*B
          ! First compute Lx**2, Ly**2, and Lz**2:
          !------------------
          Term1 = 2d0*Bz2*Tab0(iZeta,ipa,Ind(iy,iz),kyy)-4d0*B2*Tabpp(iZeta,ipa,Ind(iy,iz+2),kyy)+ &
                  2d0*B*Tab0(iZeta,ipa,Ind(iy,iz),kyy)
          Term2 = -2d0*By2*Tab0(iZeta,ipa,Ind(iy-1,iz+1),kyz)-2d0*Bz2*Tab0(iZeta,ipa,Ind(iy+1,iz-1),kyz)+ &
                  8d0*B2*Tabpp(iZeta,ipa,Ind(iy+1,iz+1),kyz)
          Term3 = 2d0*By2*Tab0(iZeta,ipa,Ind(iy,iz),kzz)-4d0*B2*Tabpp(iZeta,ipa,Ind(iy+2,iz),kzz)+ &
                  2d0*B*Tab0(iZeta,ipa,Ind(iy,iz),kzz)
          Term4 = -2d0*B*Tabp(iZeta,ipa,Ind(iy+1,iz),ky)-2d0*B*Tabp(iZeta,ipa,Ind(iy,iz+1),kz)
          if (lb >= 1) then
            Term4 = Term4+dble(iy)*Tabm(iZeta,ipa,Ind(iy-1,iz),ky)+dble(iz)*Tabm(iZeta,ipa,Ind(iy,iz-1),kz)
            if (lb >= 2) then
              Term1 = Term1-dble(iz*(iz-1))*Tabmm(iZeta,ipa,Ind(iy,iz-2),kyy)
              Term2 = Term2+dble(2*iy*iz)*Tabmm(iZeta,ipa,Ind(iy-1,iz-1),kyz)
              Term3 = Term3-dble(iy*(iy-1))*Tabmm(iZeta,ipa,Ind(iy-2,iz),kzz)
            end if
          end if
          Rslt(iZeta,ipa,ipb,kxx) = Term1+Term2+Term3+Term4
          !------------------
          Term1 = 2d0*Bx2*Tab0(iZeta,ipa,Ind(iy,iz),kzz)-4d0*B2*Tabpp(iZeta,ipa,Ind(iy,iz),kzz)+ &
                  2d0*B*Tab0(iZeta,ipa,Ind(iy,iz),kzz)
          Term2 = -2d0*Bz2*Tab0(iZeta,ipa,Ind(iy,iz-1),kzx)-2d0*Bx2*Tab0(iZeta,ipa,Ind(iy,iz+1),kzx)+ &
                  8d0*B2*Tabpp(iZeta,ipa,Ind(iy,iz+1),kzx)
          Term3 = 2d0*Bz2*Tab0(iZeta,ipa,Ind(iy,iz),kxx)-4d0*B2*Tabpp(iZeta,ipa,Ind(iy,iz+2),kxx)+ &
                  2d0*B*Tab0(iZeta,ipa,Ind(iy,iz),kxx)
          Term4 = -2d0*B*Tabp(iZeta,ipa,Ind(iy,iz+1),kz)-2d0*B*Tabp(iZeta,ipa,Ind(iy,iz),kx)
          if (lb >= 1) then
            Term4 = Term4+dble(iz)*Tabm(iZeta,ipa,Ind(iy,iz-1),kz)+dble(ix)*Tabm(iZeta,ipa,Ind(iy,iz),kx)
            if (lb >= 2) then
              Term1 = Term1-dble(ix*(ix-1))*Tabmm(iZeta,ipa,Ind(iy,iz),kzz)
              Term2 = Term2+dble(2*iz*ix)*Tabmm(iZeta,ipa,Ind(iy,iz-1),kzx)
              Term3 = Term3-dble(iz*(iz-1))*Tabmm(iZeta,ipa,Ind(iy,iz-2),kxx)
            end if
          end if
          Rslt(iZeta,ipa,ipb,kyy) = Term1+Term2+Term3+Term4
          !------------------
          Term1 = 2d0*By2*Tab0(iZeta,ipa,Ind(iy,iz),kxx)-4d0*B2*Tabpp(iZeta,ipa,Ind(iy+2,iz),kxx)+ &
                  2d0*B*Tab0(iZeta,ipa,Ind(iy,iz),kxx)
          Term2 = -2d0*Bx2*Tab0(iZeta,ipa,Ind(iy+1,iz),kxy)-2d0*By2*Tab0(iZeta,ipa,Ind(iy-1,iz),kxy)+ &
                  8d0*B2*Tabpp(iZeta,ipa,Ind(iy+1,iz),kxy)
          Term3 = 2d0*Bx2*Tab0(iZeta,ipa,Ind(iy,iz),kyy)-4d0*B2*Tabpp(iZeta,ipa,Ind(iy,iz),kyy)+ &
                  2d0*B*Tab0(iZeta,ipa,Ind(iy,iz),kyy)
          Term4 = -2d0*B*Tabp(iZeta,ipa,Ind(iy,iz),kx)-2d0*B*Tabp(iZeta,ipa,Ind(iy+1,iz),ky)
          if (lb >= 1) then
            Term4 = Term4+dble(ix)*Tabm(iZeta,ipa,Ind(iy,iz),kx)+dble(iy)*Tabm(iZeta,ipa,Ind(iy-1,iz),ky)
            if (lb >= 2) then
              Term1 = Term1-dble(iy*(iy-1))*Tabmm(iZeta,ipa,Ind(iy-2,iz),kxx)
              Term2 = Term2+dble(2*ix*iy)*Tabmm(iZeta,ipa,Ind(iy-1,iz),kxy)
              Term3 = Term3-dble(ix*(ix-1))*Tabmm(iZeta,ipa,Ind(iy,iz),kyy)
            end if
          end if
          Rslt(iZeta,ipa,ipb,kzz) = Term1+Term2+Term3+Term4
          !------------------
          ! Compute (Lx*Ly+Ly*Lx)/2, etc cyclical.
          ! With Term5, (Lx*Ly) obtains. With Term6, (Ly*Lx).
          ! We want the hermitian average.
          !------------------ (Lx,Ly)
          Term1 = Bx2*Tab0(iZeta,ipa,Ind(iy,iz+1),kyz)+Bz2*Tab0(iZeta,ipa,Ind(iy,iz-1),kyz)- &
                  4d0*B2*Tabpp(iZeta,ipa,Ind(iy,iz+1),kyz)
          Term2 = -2d0*Bz2*Tab0(iZeta,ipa,Ind(iy,iz),kxy)+4d0*B2*Tabpp(iZeta,ipa,Ind(iy,iz+2),kxy)- &
                  2d0*B*Tab0(iZeta,ipa,Ind(iy,iz),kxy)
          Term3 = -Bx2*Tab0(iZeta,ipa,Ind(iy+1,iz),kzz)-By2*Tab0(iZeta,ipa,Ind(iy-1,iz),kzz)+ &
                  4d0*B2*Tabpp(iZeta,ipa,Ind(iy+1,iz),kzz)
          Term4 = By2*Tab0(iZeta,ipa,Ind(iy-1,iz+1),kxz)+Bz2*Tab0(iZeta,ipa,Ind(iy+1,iz-1),kxz)- &
                  4d0*B2*Tabpp(iZeta,ipa,Ind(iy+1,iz+1),kxz)
          Term5 = +2d0*B*Tabp(iZeta,ipa,Ind(iy,iz),ky)
          Term6 = +2d0*B*Tabp(iZeta,ipa,Ind(iy+1,iz),kx)
          if (lb >= 1) then
            Term5 = Term5-dble(ix)*Tabm(iZeta,ipa,Ind(iy,iz),ky)
            Term6 = Term6-dble(iy)*Tabm(iZeta,ipa,Ind(iy-1,iz),kx)
            if (lb >= 2) then
              Term1 = Term1-dble(ix*iz)*Tabmm(iZeta,ipa,Ind(iy,iz-1),kyz)
              Term2 = Term2+dble(iz*(iz-1))*Tabmm(iZeta,ipa,Ind(iy,iz-2),kxy)
              Term3 = Term3+dble(ix*iy)*Tabmm(iZeta,ipa,Ind(iy-1,iz),kzz)
              Term4 = Term4-dble(iy*iz)*Tabmm(iZeta,ipa,Ind(iy-1,iz-1),kxz)
            end if
          end if
          Rslt(iZeta,ipa,ipb,kxy) = Term1+Term2+Term3+Term4+0.5d0*(Term5+Term6)
          !------------------ (Ly,Lz)
          Term1 = By2*Tab0(iZeta,ipa,Ind(iy-1,iz),kzx)+Bx2*Tab0(iZeta,ipa,Ind(iy+1,iz),kzx)- &
                  4d0*B2*Tabpp(iZeta,ipa,Ind(iy+1,iz),kzx)
          Term2 = -2d0*Bx2*Tab0(iZeta,ipa,Ind(iy,iz),kyz)+4d0*B2*Tabpp(iZeta,ipa,Ind(iy,iz),kyz)- &
                  2d0*B*Tab0(iZeta,ipa,Ind(iy,iz),kyz)
          Term3 = -By2*Tab0(iZeta,ipa,Ind(iy-1,iz+1),kxx)-Bz2*Tab0(iZeta,ipa,Ind(iy+1,iz-1),kxx)+ &
                  4d0*B2*Tabpp(iZeta,ipa,Ind(iy+1,iz+1),kxx)
          Term4 = Bz2*Tab0(iZeta,ipa,Ind(iy,iz-1),kyx)+Bx2*Tab0(iZeta,ipa,Ind(iy,iz+1),kyx)- &
                  4d0*B2*Tabpp(iZeta,ipa,Ind(iy,iz+1),kyx)
          Term5 = +2d0*B*Tabp(iZeta,ipa,Ind(iy+1,iz),kz)
          Term6 = +2d0*B*Tabp(iZeta,ipa,Ind(iy,iz+1),ky)
          if (lb >= 1) then
            Term5 = Term5-dble(iy)*Tabm(iZeta,ipa,Ind(iy-1,iz),kz)
            Term6 = Term6-dble(iz)*Tabm(iZeta,ipa,Ind(iy,iz-1),ky)
            if (lb >= 2) then
              Term1 = Term1-dble(iy*ix)*Tabmm(iZeta,ipa,Ind(iy-1,iz),kzx)
              Term2 = Term2+dble(ix*(ix-1))*Tabmm(iZeta,ipa,Ind(iy,iz),kyz)
              Term3 = Term3+dble(iy*iz)*Tabmm(iZeta,ipa,Ind(iy-1,iz-1),kxx)
              Term4 = Term4-dble(iz*ix)*Tabmm(iZeta,ipa,Ind(iy,iz-1),kyx)
            end if
          end if
          Rslt(iZeta,ipa,ipb,kyz) = Term1+Term2+Term3+Term4+0.5d0*(Term5+Term6)
          !------------------ (Lx,Lz)
          Term1 = Bz2*Tab0(iZeta,ipa,Ind(iy+1,iz-1),kxy)+By2*Tab0(iZeta,ipa,Ind(iy-1,iz+1),kxy)- &
                  4d0*B2*Tabpp(iZeta,ipa,Ind(iy+1,iz+1),kxy)
          Term2 = -2d0*By2*Tab0(iZeta,ipa,Ind(iy,iz),kzx)+4d0*B2*Tabpp(iZeta,ipa,Ind(iy+2,iz),kzx)- &
                  2d0*B*Tab0(iZeta,ipa,Ind(iy,iz),kzx)
          Term3 = -Bz2*Tab0(iZeta,ipa,Ind(iy,iz-1),kyy)-Bx2*Tab0(iZeta,ipa,Ind(iy,iz+1),kyy)+ &
                  4d0*B2*Tabpp(iZeta,ipa,Ind(iy,iz+1),kyy)
          Term4 = Bx2*Tab0(iZeta,ipa,Ind(iy+1,iz),kzy)+By2*Tab0(iZeta,ipa,Ind(iy-1,iz),kzy)- &
                  4d0*B2*Tabpp(iZeta,ipa,Ind(iy+1,iz),kzy)
          Term5 = +2d0*B*Tabp(iZeta,ipa,Ind(iy,iz+1),kx)
          Term6 = +2d0*B*Tabp(iZeta,ipa,Ind(iy,iz),kz)
          if (lb >= 1) then
            Term5 = Term5-dble(iz)*Tabm(iZeta,ipa,Ind(iy,iz-1),kx)
            Term6 = Term6-dble(ix)*Tabm(iZeta,ipa,Ind(iy,iz),kz)
            if (lb >= 2) then
              Term1 = Term1-dble(iz*iy)*Tabmm(iZeta,ipa,Ind(iy-1,iz-1),kxy)
              Term2 = Term2+dble(iy*(iy-1))*Tabmm(iZeta,ipa,Ind(iy-2,iz),kzx)
              Term3 = Term3+dble(iz*ix)*Tabmm(iZeta,ipa,Ind(iy,iz-1),kyy)
              Term4 = Term4-dble(ix*iy)*Tabmm(iZeta,ipa,Ind(iy-1,iz),kzy)
            end if
          end if
          Rslt(iZeta,ipa,ipb,kxz) = Term1+Term2+Term3+Term4+0.5d0*(Term5+Term6)

        end do

      end do
    end do

  end do
end do

if (iPrint >= 49) then
  write(6,*) ' In AMPr la,lb=',la,lb
  do iElem=1,nElem(la)
    do jElem=1,nElem(lb)
      write(Label,'(A,I2,A,I2,A)') ' Rslt (',iElem,',',jElem,',xx) '
      call RecPrt(Label,' ',Rslt(1,iElem,jElem,kxx),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Rslt (',iElem,',',jElem,',xy) '
      call RecPrt(Label,' ',Rslt(1,iElem,jElem,kxy),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Rslt (',iElem,',',jElem,',xz) '
      call RecPrt(Label,' ',Rslt(1,iElem,jElem,kxz),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Rslt (',iElem,',',jElem,',yy) '
      call RecPrt(Label,' ',Rslt(1,iElem,jElem,kyy),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Rslt (',iElem,',',jElem,',yz) '
      call RecPrt(Label,' ',Rslt(1,iElem,jElem,kyz),nZeta,1)
      write(Label,'(A,I2,A,I2,A)') ' Rslt (',iElem,',',jElem,',zz) '
      call RecPrt(Label,' ',Rslt(1,iElem,jElem,kzz),nZeta,1)
    end do
  end do
  write(6,*) ' Leaving AMPr.'
end if

return

end subroutine AMPr
