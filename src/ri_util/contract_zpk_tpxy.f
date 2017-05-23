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
      subroutine contract_Zpk_Tpxy(Zpk,nZpk,Txy,nTxy,Scrt,nScrt,
     &                             Diag,nDiag,nnP,nBas_Aux,
     &                             nAdens,nAvec,nAct,nIrrep)
      Implicit none
      Integer nZpk,nTxy,nScrt,nDiag,nIrrep,nnP(0:nIrrep-1),
     &        nBas_Aux(0:nIrrep-1),nact(0:nIrrep-1),nAdens,nAvec,
     &        i,j,k,l,nCumnnP,nCumnnP2,nCumnnP3,
     &        ip,jp,kp,iSym,jSym,kSym,iDen
      Real*8 Zpk(nZpk,nAVec),Txy(nTxy,nAdens),Scrt(nScrt),
     &       Diag(nDiag,nADens)
************************************************************************
*
      Do l=1,nAVec
        iDen=1
        If (l.gt.1) iDen=2
*
        nCumnnP=0
        nCumnnP2=0
        nCumnnP3=0
        Do iSym=0,nIrrep-1
          Do i=1,nBas_Aux(iSym)
             ip=nCumnnP2+(i-1)*nnP(iSym)
             Do j=1,nnP(iSym)
                Scrt(j)=0.0d0
                Do k=1,nnP(iSym)
                  kp=nCumnnP3+(k-1)*nnP(iSym)
                  Scrt(j)=Scrt(j)+Sign(1.0d0,Diag(nCumnnP+k,iDen))
     &                    *Zpk(k+ip,l)*Txy(j+kp,iDen)
                End Do
             End Do
             Do j=1,nnP(iSym)
               Zpk(j+ip,l)=Scrt(j)
             End Do
          End Do
*Now correct for the 2 factor
          Do i=1,nBas_Aux(iSym)
             ip=nCumnnP2+(i-1)*nnP(iSym)
             Do jSym=0,nIrrep-1
               kSym=iEOR(jsym,isym)
               If (kSym.le.jSym) Then
                 Do j=1,nAct(jSym)
                   If (kSym.eq.jSym) Then
                     jp=ip+j*(j-1)/2
                     Do k=1,j-1
                       Zpk(jp+k,l)=Zpk(jp+k,l)/2.0d0
                     End Do
                   Else
                     jp=ip+(j-1)*nAct(kSym)
                     Do k=1,nAct(kSym)
                       Zpk(jp+k,l)=Zpk(jp+k,l)/2.0d0
                     End Do
                   EndIf
                 End Do
                 If (kSym.eq.jSym) Then
                   ip=ip+nAct(jSym)*(nAct(jSym)+1)/2
                 Else
                   ip=ip+nAct(jSym)*nAct(kSym)
                 End If
               End If
             End Do
          End Do
*
*       Call RecPrt('Zpk',' ',Zpk(nCumnnP2+1,l),nnP(iSym),nBas_Aux(iSym))
          nCumnnP=nCumnnP+nnP(iSym)
          nCumnnP2=nCumnnP2+nnP(iSym)*nBas_Aux(iSym)
          nCumnnP3=nCumnnP3+nnP(iSym)**2
        End Do
      End Do
      end
