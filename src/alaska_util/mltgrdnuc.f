************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine MltGrdNuc(Grad,nGrad,nOrdOp)
      use Basis_Info
      use Center_Info
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "print.fh"
#include "real.fh"
#include "disp.fh"
      Real*8 Grad(nGrad),C(3)
      parameter (lforce=20)
      real*8 Force(lforce)
      common /finfld/Force
*
      logical TstFnc,TF
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
      TF(mdc,iIrrep,iComp) = TstFnc(dc(mdc)%iCoSet,
     &                              iIrrep,iComp,dc(mdc)%nStab)
*
      iIrrep=0
      do 800 ixop=0,nOrdOp
      do 801 iyop=0,nOrdOp-ixop
      izop=nOrdOp-ixop-iyop
      icomp=Ind(nOrdOp,ixop,izop)
      ff=Force(icomp)
      if(ff.eq.0.d0) goto 801
        kdc = 0
        Do kCnttp = 1, nCnttp
           If (dbsc(kCnttp)%Charge.eq.0.d0) Go To 411
           Do kCnt = 1, dbsc(kCnttp)%nCntr
              C(1:3)=dbsc(kCnttp)%Coor(1:3,kCnt)
              ndc=kdc+kCnt
              Fact=-dbsc(kCnttp)%Charge*ff
              nDisp = IndDsp(ndc,iIrrep)
              Do iCar = 0, 2
                iComp = 2**iCar
                If ( TF(ndc,iIrrep,iComp) .and.
     &              .Not.dbsc(kCnttp)%pChrg ) Then
                  nDisp = nDisp + 1
                  If (Direct(nDisp)) Then
                   XGrad=0.d0
                   if(iCar.eq.0) then
                     if(ixop.gt.0) XGrad=Fact*Dble(ixop)*C(1)**(ixop-1)*
     &                              C(2)**iyop*C(3)**izop
                   else if(iCar.eq.1) then
                     if(iyop.gt.0) XGrad=Fact*Dble(iyop)*C(1)**ixop*
     &                              C(2)**(iyop-1)*C(3)**izop
                   else
                    if(izop.gt.0) XGrad=Fact*Dble(izop)*C(1)**ixop*
     &                              C(2)**iyop*C(3)**(izop-1)
                   endif
                   Grad(nDisp)=Grad(nDisp)+XGrad
                  End If
                End If
             Enddo
           Enddo
411     kdc = kdc + dbsc(kCnttp)%nCntr
        Enddo
801   continue
800   continue

      Return
      End
