************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1996, Per Ake Malmqvist                                *
************************************************************************
      SubRoutine AMPr(Beta,nZeta,Rslt,la,lb,Tabpp,Tabp,Tab0,
     &                Tabm,Tabmm)
************************************************************************
*                                                                      *
* Object: Compute matrix elements of hermitized products of angular    *
*         moment operators, using elementary overlaps, dipole, and     *
*         quadrupole integrals.                                        *
* Called from: AMPInt                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Per-Ake Malmqvist, Dept. of Theoretical Chemistry,       *
*             University of Lund, SWEDEN                               *
*             November '96                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8  Rslt  (nZeta,((la+1)*(la+2))/2,((lb+1)*(lb+2))/2,6),
     &        Tabpp (nZeta,((la+1)*(la+2))/2,((lb+3)*(lb+4))/2,6),
     &        Tabp  (nZeta,((la+1)*(la+2))/2,((lb+2)*(lb+3))/2,3),
     &        Tab0  (nZeta,((la+1)*(la+2))/2,((lb+1)*(lb+2))/2,6),
     &        Tabm  (nZeta,((la+1)*(la+2))/2,( lb   *(lb+1))/2,3),
     &        Tabmm (nZeta,((la+1)*(la+2))/2,((lb-1)* lb   )/2,6),
     &        Beta(nZeta)
      Character*80 Label
      data kx, ky, kz  / 1, 2, 3 /
      data kxx,kxy,kxz,kyy,kyz,kzz / 1,2,3,4,5,6 /
      data kyx,kzx,kzy/ 2,3,5 /
*
*     Statement function for cartesian index
      Ind(j,k)=((j+k)*(j+k+1))/2+k+1
      nElem(ix) = (ix+1)*(ix+2)/2
*
      iRout = 221
      iPrint = nPrint(iRout)
      Call qEnter('AMPr ')
*
      If (iPrint.ge.99) Then
          Write (6,*) ' In AMPr la,lb=',la,lb
          Call RecPrt('Beta',' ',Beta,nZeta,1)
          Do ia = 1, nElem(la)
             Do ib = 1, nElem(lb+2)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Tabpp(',ia,',',ib,'xx)'
                Call RecPrt(Label,' ',Tabpp(1,ia,ib,1),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Tabpp(',ia,',',ib,'xy)'
                Call RecPrt(Label,' ',Tabpp(1,ia,ib,2),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Tabpp(',ia,',',ib,'xz)'
                Call RecPrt(Label,' ',Tabpp(1,ia,ib,3),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Tabpp(',ia,',',ib,'yy)'
                Call RecPrt(Label,' ',Tabpp(1,ia,ib,4),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Tabpp(',ia,',',ib,'yz)'
                Call RecPrt(Label,' ',Tabpp(1,ia,ib,5),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Tabpp(',ia,',',ib,'zz)'
                Call RecPrt(Label,' ',Tabpp(1,ia,ib,6),nZeta,1)
             End Do
          End Do
          Do ia = 1, nElem(la)
             Do ib = 1, nElem(lb+1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Tabp(',ia,',',ib,'x)'
                Call RecPrt(Label,' ',Tabp(1,ia,ib,1),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Tabp(',ia,',',ib,'y)'
                Call RecPrt(Label,' ',Tabp(1,ia,ib,2),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Tabp(',ia,',',ib,'z)'
                Call RecPrt(Label,' ',Tabp(1,ia,ib,3),nZeta,1)
             End Do
          End Do
          Do ia = 1, nElem(la)
             Do ib = 1, nElem(lb)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Tab0',ia,',',ib,'xx)'
                Call RecPrt(Label,' ',Tab0(1,ia,ib,1),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Tab0',ia,',',ib,'xy)'
                Call RecPrt(Label,' ',Tab0(1,ia,ib,2),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Tab0',ia,',',ib,'xz)'
                Call RecPrt(Label,' ',Tab0(1,ia,ib,3),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Tab0',ia,',',ib,'yy)'
                Call RecPrt(Label,' ',Tab0(1,ia,ib,4),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Tab0',ia,',',ib,'yz)'
                Call RecPrt(Label,' ',Tab0(1,ia,ib,5),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                 ' Tab0',ia,',',ib,'zz)'
                Call RecPrt(Label,' ',Tab0(1,ia,ib,6),nZeta,1)
             End Do
          End Do
          If (lb.gt.0) Then
             Do ia = 1, nElem(la)
                Do ib = 1, nElem(lb-1)
                   Write (Label,'(A,I2,A,I2,A)')
     &                    ' Tabm(',ia,',',ib,'x)'
                   Call RecPrt(Label,' ',Tabm(1,ia,ib,1),nZeta,1)
                   Write (Label,'(A,I2,A,I2,A)')
     &                    ' Tabm(',ia,',',ib,'y)'
                   Call RecPrt(Label,' ',Tabm(1,ia,ib,2),nZeta,1)
                   Write (Label,'(A,I2,A,I2,A)')
     &                    ' Tabm(',ia,',',ib,'z)'
                   Call RecPrt(Label,' ',Tabm(1,ia,ib,3),nZeta,1)
                End Do
             End Do
          End If
          If (lb.gt.1) Then
             Do ia = 1, nElem(la)
                Do ib = 1, nElem(lb-2)
                   Write (Label,'(A,I2,A,I2,A)')
     &                    ' Tabmm(',ia,',',ib,'xx)'
                   Call RecPrt(Label,' ',Tabmm(1,ia,ib,1),nZeta,1)
                   Write (Label,'(A,I2,A,I2,A)')
     &                    ' Tabmm(',ia,',',ib,'xy)'
                   Call RecPrt(Label,' ',Tabmm(1,ia,ib,2),nZeta,1)
                   Write (Label,'(A,I2,A,I2,A)')
     &                    ' Tabmm(',ia,',',ib,'xz)'
                   Call RecPrt(Label,' ',Tabmm(1,ia,ib,3),nZeta,1)
                   Write (Label,'(A,I2,A,I2,A)')
     &                    ' Tabmm(',ia,',',ib,'yy)'
                   Call RecPrt(Label,' ',Tabmm(1,ia,ib,4),nZeta,1)
                   Write (Label,'(A,I2,A,I2,A)')
     &                    ' Tabmm(',ia,',',ib,'yz)'
                   Call RecPrt(Label,' ',Tabmm(1,ia,ib,5),nZeta,1)
                   Write (Label,'(A,I2,A,I2,A)')
     &                    ' Tabmm(',ia,',',ib,'zz)'
                   Call RecPrt(Label,' ',Tabmm(1,ia,ib,6),nZeta,1)
                End Do
             End Do
          End If
      End If
*
      Do 100 ixa = la, 0, -1
         Do 101 iya = la-ixa, 0, -1
            iza = la-ixa-iya
            ipa = Ind(iya,iza)
*
      Do 200 ixb = lb, 0, -1
         Do 201 iyb = lb-ixb, 0, -1
            izb = lb-ixb-iyb
            ipb = Ind(iyb,izb)
            ix=ixb
            iy=iyb
            iz=izb
*
            Do 300 iZeta = 1, nZeta
               B=Beta(iZeta)
               B2=B**2
               Bx2=Dble(2*ix)*B
               By2=Dble(2*iy)*B
               Bz2=Dble(2*iz)*B
C First compute Lx**2, Ly**2, and Lz**2:
C------------------
               Term1=    2d0*Bz2*Tab0 (iZeta,ipa,Ind(iy  ,iz  ),kyy)
     &                    -4d0*B2*Tabpp(iZeta,ipa,Ind(iy  ,iz+2),kyy)
     &                     +2d0*B*Tab0 (iZeta,ipa,Ind(iy  ,iz  ),kyy)
               Term2=   -2d0*By2*Tab0 (iZeta,ipa,Ind(iy-1,iz+1),kyz)
     &                  -2d0*Bz2*Tab0 (iZeta,ipa,Ind(iy+1,iz-1),kyz)
     &                    +8d0*B2*Tabpp(iZeta,ipa,Ind(iy+1,iz+1),kyz)
               Term3=    2d0*By2*Tab0 (iZeta,ipa,Ind(iy  ,iz  ),kzz)
     &                    -4d0*B2*Tabpp(iZeta,ipa,Ind(iy+2,iz  ),kzz)
     &                     +2d0*B*Tab0 (iZeta,ipa,Ind(iy  ,iz  ),kzz)
               Term4=      -2d0*B*Tabp (iZeta,ipa,Ind(iy+1,iz  ),ky )
     &                     -2d0*B*Tabp (iZeta,ipa,Ind(iy  ,iz+1),kz )
               if(lb.ge.1) then
               Term4=Term4 +Dble(iy)*Tabm (iZeta,ipa,Ind(iy-1,iz  ),ky )
     &                     +Dble(iz)*Tabm (iZeta,ipa,Ind(iy  ,iz-1),kz )
               if(lb.ge.2) then
               Term1=Term1
     &              -Dble(iz*(iz-1))*Tabmm(iZeta,ipa,Ind(iy  ,iz-2),kyy)
               Term2=Term2
     &              +Dble(2*iy*iz)*Tabmm(iZeta,ipa,Ind(iy-1,iz-1),kyz)
               Term3=Term3
     &              -Dble(iy*(iy-1))*Tabmm(iZeta,ipa,Ind(iy-2,iz  ),kzz)
               end if
               end if
               Rslt(iZeta,ipa,ipb,kxx) = Term1+Term2+Term3+Term4
C------------------
               Term1=    2d0*Bx2*Tab0 (iZeta,ipa,Ind(iy,iz  ),kzz)
     &                    -4d0*B2*Tabpp(iZeta,ipa,Ind(iy,iz  ),kzz)
     &                     +2d0*B*Tab0 (iZeta,ipa,Ind(iy,iz  ),kzz)
               Term2=   -2d0*Bz2*Tab0 (iZeta,ipa,Ind(iy,iz-1),kzx)
     &                  -2d0*Bx2*Tab0 (iZeta,ipa,Ind(iy,iz+1),kzx)
     &                    +8d0*B2*Tabpp(iZeta,ipa,Ind(iy,iz+1),kzx)
               Term3=    2d0*Bz2*Tab0 (iZeta,ipa,Ind(iy,iz  ),kxx)
     &                    -4d0*B2*Tabpp(iZeta,ipa,Ind(iy,iz+2),kxx)
     &                     +2d0*B*Tab0 (iZeta,ipa,Ind(iy,iz  ),kxx)
               Term4=      -2d0*B*Tabp (iZeta,ipa,Ind(iy,iz+1),kz )
     &                     -2d0*B*Tabp (iZeta,ipa,Ind(iy,iz  ),kx )
               if(lb.ge.1) then
               Term4=Term4  +Dble(iz)*Tabm (iZeta,ipa,Ind(iy,iz-1),kz )
     &                      +Dble(ix)*Tabm (iZeta,ipa,Ind(iy,iz  ),kx )
               if(lb.ge.2) then
               Term1=
     &           Term1-Dble(ix*(ix-1))*Tabmm(iZeta,ipa,Ind(iy,iz  ),kzz)
               Term2=
     &           Term2+  Dble(2*iz*ix)*Tabmm(iZeta,ipa,Ind(iy,iz-1),kzx)
               Term3=
     &           Term3-Dble(iz*(iz-1))*Tabmm(iZeta,ipa,Ind(iy,iz-2),kxx)
               end if
               end if
               Rslt(iZeta,ipa,ipb,kyy) = Term1+Term2+Term3+Term4
C------------------
               Term1=    2d0*By2*Tab0 (iZeta,ipa,Ind(iy  ,iz),kxx)
     &                    -4d0*B2*Tabpp(iZeta,ipa,Ind(iy+2,iz),kxx)
     &                     +2d0*B*Tab0 (iZeta,ipa,Ind(iy  ,iz),kxx)
               Term2=   -2d0*Bx2*Tab0 (iZeta,ipa,Ind(iy+1,iz),kxy)
     &                  -2d0*By2*Tab0 (iZeta,ipa,Ind(iy-1,iz),kxy)
     &                    +8d0*B2*Tabpp(iZeta,ipa,Ind(iy+1,iz),kxy)
               Term3=    2d0*Bx2*Tab0 (iZeta,ipa,Ind(iy  ,iz),kyy)
     &                    -4d0*B2*Tabpp(iZeta,ipa,Ind(iy  ,iz),kyy)
     &                     +2d0*B*Tab0 (iZeta,ipa,Ind(iy  ,iz),kyy)
               Term4=      -2d0*B*Tabp (iZeta,ipa,Ind(iy  ,iz),kx )
     &                     -2d0*B*Tabp (iZeta,ipa,Ind(iy+1,iz),ky )
               if(lb.ge.1) then
               Term4=Term4+  Dble(ix)*Tabm (iZeta,ipa,Ind(iy  ,iz),kx )
     &                      +Dble(iy)*Tabm (iZeta,ipa,Ind(iy-1,iz),ky )
               if(lb.ge.2) then
               Term1=
     &           Term1-Dble(iy*(iy-1))*Tabmm(iZeta,ipa,Ind(iy-2,iz),kxx)
               Term2=
     &           Term2+  Dble(2*ix*iy)*Tabmm(iZeta,ipa,Ind(iy-1,iz),kxy)
               Term3=
     &           Term3-Dble(ix*(ix-1))*Tabmm(iZeta,ipa,Ind(iy  ,iz),kyy)
               end if
               end if
               Rslt(iZeta,ipa,ipb,kzz) = Term1+Term2+Term3+Term4
C------------------
C Compute (Lx*Ly+Ly*Lx)/2, etc cyclical.
C With Term5, (Lx*Ly) obtains. With Term6, (Ly*Lx).
C We want the hermitian average.
C------------------ (Lx,Ly)
               Term1=       Bx2*Tab0 (iZeta,ipa,Ind(iy  ,iz+1),kyz)
     &                     +Bz2*Tab0 (iZeta,ipa,Ind(iy  ,iz-1),kyz)
     &                    -4d0*B2*Tabpp(iZeta,ipa,Ind(iy  ,iz+1),kyz)
               Term2=    -2d0*Bz2*Tab0 (iZeta,ipa,Ind(iy  ,iz  ),kxy)
     &                    +4d0*B2*Tabpp(iZeta,ipa,Ind(iy  ,iz+2),kxy)
     &                     -2d0*B*Tab0 (iZeta,ipa,Ind(iy  ,iz  ),kxy)
               Term3=      -Bx2*Tab0 (iZeta,ipa,Ind(iy+1,iz  ),kzz)
     &                     -By2*Tab0 (iZeta,ipa,Ind(iy-1,iz  ),kzz)
     &                    +4d0*B2*Tabpp(iZeta,ipa,Ind(iy+1,iz  ),kzz)
               Term4=       By2*Tab0 (iZeta,ipa,Ind(iy-1,iz+1),kxz)
     &                     +Bz2*Tab0 (iZeta,ipa,Ind(iy+1,iz-1),kxz)
     &                    -4d0*B2*Tabpp(iZeta,ipa,Ind(iy+1,iz+1),kxz)
               Term5=      +2d0*B*Tabp (iZeta,ipa,Ind(iy  ,iz  ),ky )
               Term6=      +2d0*B*Tabp (iZeta,ipa,Ind(iy+1,iz  ),kx )
               if(lb.ge.1) then
               Term5=Term5-Dble(ix)*Tabm (iZeta,ipa,Ind(iy  ,iz  ),ky )
               Term6=Term6-Dble(iy)*Tabm (iZeta,ipa,Ind(iy-1,iz  ),kx )
               if(lb.ge.2) then
               Term1=Term1
     &              -Dble(ix*iz)*Tabmm(iZeta,ipa,Ind(iy  ,iz-1),kyz)
               Term2=Term2
     &              +Dble(iz*(iz-1))*Tabmm(iZeta,ipa,Ind(iy  ,iz-2),kxy)
               Term3=Term3
     &              +Dble(ix*iy)*Tabmm(iZeta,ipa,Ind(iy-1,iz  ),kzz)
               Term4=Term4
     &              -Dble(iy*iz)*Tabmm(iZeta,ipa,Ind(iy-1,iz-1),kxz)
               end if
               end if
               Rslt(iZeta,ipa,ipb,kxy) = Term1+Term2+Term3+Term4
     &                     +0.5D0*(Term5+Term6)
C------------------ (Ly,Lz)
               Term1=       By2*Tab0 (iZeta,ipa,Ind(iy-1,iz  ),kzx)
     &                     +Bx2*Tab0 (iZeta,ipa,Ind(iy+1,iz  ),kzx)
     &                    -4d0*B2*Tabpp(iZeta,ipa,Ind(iy+1,iz  ),kzx)
               Term2=    -2d0*Bx2*Tab0 (iZeta,ipa,Ind(iy  ,iz  ),kyz)
     &                    +4d0*B2*Tabpp(iZeta,ipa,Ind(iy  ,iz  ),kyz)
     &                     -2d0*B*Tab0 (iZeta,ipa,Ind(iy  ,iz  ),kyz)
               Term3=      -By2*Tab0 (iZeta,ipa,Ind(iy-1,iz+1),kxx)
     &                     -Bz2*Tab0 (iZeta,ipa,Ind(iy+1,iz-1),kxx)
     &                    +4d0*B2*Tabpp(iZeta,ipa,Ind(iy+1,iz+1),kxx)
               Term4=       Bz2*Tab0 (iZeta,ipa,Ind(iy  ,iz-1),kyx)
     &                     +Bx2*Tab0 (iZeta,ipa,Ind(iy  ,iz+1),kyx)
     &                    -4d0*B2*Tabpp(iZeta,ipa,Ind(iy  ,iz+1),kyx)
               Term5=      +2d0*B*Tabp (iZeta,ipa,Ind(iy+1,iz  ),kz )
               Term6=      +2d0*B*Tabp (iZeta,ipa,Ind(iy  ,iz+1),ky )
               if(lb.ge.1) then
               Term5=Term5-Dble(iy)*Tabm (iZeta,ipa,Ind(iy-1,iz  ),kz )
               Term6=Term6-Dble(iz)*Tabm (iZeta,ipa,Ind(iy  ,iz-1),ky )
               if(lb.ge.2) then
               Term1=
     &          Term1 -Dble(iy*ix)*Tabmm(iZeta,ipa,Ind(iy-1,iz  ),kzx)
               Term2=Term2
     &              +Dble(ix*(ix-1))*Tabmm(iZeta,ipa,Ind(iy  ,iz  ),kyz)
               Term3=
     &          Term3 +Dble(iy*iz)*Tabmm(iZeta,ipa,Ind(iy-1,iz-1),kxx)
               Term4=
     &          Term4 -Dble(iz*ix)*Tabmm(iZeta,ipa,Ind(iy  ,iz-1),kyx)
               end if
               end if
               Rslt(iZeta,ipa,ipb,kyz) = Term1+Term2+Term3+Term4
     &                     +0.5D0*(Term5+Term6)
C------------------ (Lx,Lz)
               Term1=       Bz2*Tab0 (iZeta,ipa,Ind(iy+1,iz-1),kxy)
     &                     +By2*Tab0 (iZeta,ipa,Ind(iy-1,iz+1),kxy)
     &                    -4d0*B2*Tabpp(iZeta,ipa,Ind(iy+1,iz+1),kxy)
               Term2=    -2d0*By2*Tab0 (iZeta,ipa,Ind(iy  ,iz  ),kzx)
     &                    +4d0*B2*Tabpp(iZeta,ipa,Ind(iy+2,iz  ),kzx)
     &                     -2d0*B*Tab0 (iZeta,ipa,Ind(iy  ,iz  ),kzx)
               Term3=      -Bz2*Tab0 (iZeta,ipa,Ind(iy  ,iz-1),kyy)
     &                     -Bx2*Tab0 (iZeta,ipa,Ind(iy  ,iz+1),kyy)
     &                    +4d0*B2*Tabpp(iZeta,ipa,Ind(iy  ,iz+1),kyy)
               Term4=       Bx2*Tab0 (iZeta,ipa,Ind(iy+1,iz  ),kzy)
     &                     +By2*Tab0 (iZeta,ipa,Ind(iy-1,iz  ),kzy)
     &                    -4d0*B2*Tabpp(iZeta,ipa,Ind(iy+1,iz  ),kzy)
               Term5=      +2d0*B*Tabp (iZeta,ipa,Ind(iy  ,iz+1),kx )
               Term6=      +2d0*B*Tabp (iZeta,ipa,Ind(iy  ,iz  ),kz )
               if(lb.ge.1) then
               Term5=Term5 -Dble(iz)*Tabm (iZeta,ipa,Ind(iy  ,iz-1),kx )
               Term6=Term6 -Dble(ix)*Tabm (iZeta,ipa,Ind(iy  ,iz  ),kz )
               if(lb.ge.2) then
               Term1=Term1
     &              -Dble(iz*iy)*Tabmm(iZeta,ipa,Ind(iy-1,iz-1),kxy)
               Term2=Term2
     &              +Dble(iy*(iy-1))*Tabmm(iZeta,ipa,Ind(iy-2,iz  ),kzx)
               Term3=Term3
     &              +Dble(iz*ix)*Tabmm(iZeta,ipa,Ind(iy  ,iz-1),kyy)
               Term4=Term4
     &              -Dble(ix*iy)*Tabmm(iZeta,ipa,Ind(iy-1,iz  ),kzy)
               end if
               end if
               Rslt(iZeta,ipa,ipb,kxz) = Term1+Term2+Term3+Term4
     &                     +0.5D0*(Term5+Term6)

 300       Continue
*
 201     Continue
 200  Continue
*
 101     Continue
 100  Continue
*
      If (iPrint.ge.49) Then
          Write (6,*) ' In AMPr la,lb=',la,lb
          Do iElem = 1, nElem(la)
             Do jElem = 1, nElem(lb)
                Write (Label,'(A,I2,A,I2,A)')
     &                ' Rslt (',iElem,',',jElem,',xx) '
                Call RecPrt(Label,' ',Rslt(1,iElem,jElem,kxx),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                ' Rslt (',iElem,',',jElem,',xy) '
                Call RecPrt(Label,' ',Rslt(1,iElem,jElem,kxy),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                ' Rslt (',iElem,',',jElem,',xz) '
                Call RecPrt(Label,' ',Rslt(1,iElem,jElem,kxz),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                ' Rslt (',iElem,',',jElem,',yy) '
                Call RecPrt(Label,' ',Rslt(1,iElem,jElem,kyy),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                ' Rslt (',iElem,',',jElem,',yz) '
                Call RecPrt(Label,' ',Rslt(1,iElem,jElem,kyz),nZeta,1)
                Write (Label,'(A,I2,A,I2,A)')
     &                ' Rslt (',iElem,',',jElem,',zz) '
                Call RecPrt(Label,' ',Rslt(1,iElem,jElem,kzz),nZeta,1)
             End Do
          End Do
          Write (6,*) ' Leaving AMPr.'
      End If
*
      Call qExit('AMPr ')
      Return
      End
