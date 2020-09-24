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
* Copyright (C) 1993, Roland Lindh                                     *
*               1993, Per Boussard                                     *
************************************************************************
      SubRoutine PAM2Int(
#define _CALLING_
#include "int_interface.fh"
     &                  )
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of PAM integrals used in  *
*         PAM calculations. The operator is a gaussian type function   *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              DCopy   (ESSL)                                          *
*              DCR                                                     *
*              CrtCmp                                                  *
*              Assmbl                                                  *
*              CmbnMP                                                  *
*              DaXpY   (ESSL)                                          *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*      Alpha : exponents of bra gaussians                              *
*      nAlpha: number of primitives (exponents) of bra gaussians       *
*      Beta  : as Alpha but for ket gaussians                          *
*      nBeta : as nAlpha but for the ket gaussians                     *
*      Zeta  : sum of exponents (nAlpha x nBeta)                       *
*      ZInv  : inverse of Zeta                                         *
*      rKappa: gaussian prefactor for the products of bra and ket      *
*              gaussians.                                              *
*      P     : center of new gaussian from the products of bra and ket *
*              gaussians.                                              *
*      Final : array for computed integrals                            *
*      nZeta : nAlpha x nBeta                                          *
*      nComp : number of components in the operator (e.g. dipolmoment  *
*              operator has three components)                          *
*      la    : total angular momentum of bra gaussian                  *
*      lb    : total angular momentum of ket gaussian                  *
*      A     : center of bra gaussian                                  *
*      B     : center of ket gaussian                                  *
*      nRys  : order of Rys- or Hermite-Gauss polynomial               *
*      Array : Auxiliary memory as requested by ECPMem                 *
*      nArr  : length of Array                                         *
*      Ccoor : coordinates of the operator, zero for symmetric oper.   *
*      NOrdOp: Order of the operator                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden, and Per Boussard, Dept. of Theoretical  *
*             Physics, University of Stockholm, Sweden, October '93.   *
************************************************************************
      use Basis_Info
      use Center_Info
      use Her_RW
      use PAM2
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"

#include "int_interface.fh"

*     Local variables
      Real*8 TC(3), C(3)
      Character*80 Label
      Logical ABeq(3)
      Integer iDCRT(0:7)
*
*-----Statement function for Cartesian index
*
      nElem(k)=(k+1)*(k+2)/2
*
      iRout = 122
      iPrint = nPrint(iRout)
*     Call QEnter('PAM2Int')
*     Call GetMem(' Enter PAM2Int','LIST','REAL',iDum,iDum)
*
      nip = 1
      ipAxyz = nip
      nip = nip + nZeta*3*nHer*(la+1)
      ipBxyz = nip
      nip = nip + nZeta*3*nHer*(lb+1)
      ipRxyz = nip
      nip = nip + nZeta*3*nHer*(nOrdOp+1)
      ipQxyz = nip
      nip = nip + nZeta*3*(la+1)*(lb+1)*(nOrdOp+1)
      ipK = nip
      nip = nip + nZeta
      ipZ = nip
      nip = nip + nZeta
      ipPx= nip
      nip = nip + nZeta
      ipPy= nip
      nip = nip + nZeta
      ipPz= nip
      nip = nip + nZeta
      ipRes = nip
      nip = nip + nZeta*nComp*nElem(la)*nElem(lb)
      If (nip-1.gt.nArr*nZeta) Then
         Call WarningMessage(2,'PAM2Int: nip-1.gt.nArr*nZeta')
         Write (6,*) ' nArr is Wrong! ', nip-1,' > ',nArr*nZeta
         Write (6,*) ' Abend in PAM2Int'
         Call Abend()
      End If
*
      If (iPrint.ge.49) Then
         Call RecPrt(' In PAM2Int: A',' ',A,1,3)
         Call RecPrt(' In PAM2Int: RB',' ',RB,1,3)
         Call RecPrt(' In PAM2Int: Ccoor',' ',Ccoor,1,3)
         Call RecPrt(' In PAM2Int: Kappa',' ',rKappa,nAlpha,nBeta)
         Call RecPrt(' In PAM2Int: Zeta',' ',Zeta,nAlpha,nBeta)
         Call RecPrt(' In PAM2Int: P',' ',P,nZeta,3)
         Write (6,*) ' In PAM2Int: la,lb,nHer=',la,lb,nHer
      End If
*
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,Final,1)
*
*-----Loop over nuclear centers
*
      kdc=0
      if (kCnttpPAM.gt.1) Then
         do ikdc=1,kCnttpPAM-1
            kdc = kdc + dbsc(ikdc)%nCntr
         end do
      end if


      kCnttp = kCnttpPAM
         If (.Not.dbsc(kCnttp)%lPAM2) Go To 111
         If (dbsc(kCnttp)%nPAM2.eq.-1) Go To 111
*
         Do 101 kCnt = 1, dbsc(kCnttp)%nCntr
            C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)
*
            Call DCR(LmbdT,iStabM,nStabM,
     &               dc(kdc+kCnt)%iStab, dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
            Fact = DBLE(nStabM) / DBLE(LmbdT)
*
            Do 102 lDCRT = 0, nDCRT-1
               Call OA(iDCRT(lDCRT),C,TC)
*
                  Call GetMem(' Scr','ALLO','REAL',ipScr,
     &                       nZeta*nElem(la)*nElem(lb)*nComp)
                  call dcopy_(nZeta*nElem(la)*nElem(lb)*nComp,
     &                       [Zero],0,Work(ipScr),1)
               Do 1011 iM2xp = 1, iPAMPrim
                  Gamma = PAMexp(iM2xp,1)


                  If (iPrint.ge.99) Write (6,*) ' Gamma=',Gamma
*
*-----------------Modify the original basis.
*
                  Do 1012 iZeta = 1, nZeta
                     PTC2 = (P(iZeta,1)-TC(1))**2
     &                    + (P(iZeta,2)-TC(2))**2
     &                    + (P(iZeta,3)-TC(3))**2
                     Tmp0 = Zeta(iZeta)+Gamma
                     Tmp1 = Exp(-Zeta(iZeta)*Gamma*PTC2/Tmp0)
                     Array(ipK+iZeta-1)    = rKappa(iZeta) * Tmp1
                     Array(ipZ+iZeta-1)    = Tmp0
                     Array(ipPx+iZeta-1)   =
     &                  (Zeta(iZeta)*P(iZeta,1)+Gamma*TC(1))/Tmp0
                     Array(ipPy+iZeta-1)   =
     &                  (Zeta(iZeta)*P(iZeta,2)+Gamma*TC(2))/Tmp0
                     Array(ipPz+iZeta-1)   =
     &                  (Zeta(iZeta)*P(iZeta,3)+Gamma*TC(3))/Tmp0
 1012             Continue
                  If (iPrint.ge.99) Then
                     Write (6,*) ' The modified basis set'
                     Call RecPrt(' In PAM2Int: Kappa',' ',
     &                            Array(ipK),nAlpha,nBeta)
                     Call RecPrt(' In PAM2Int: Zeta',' ',
     &                            Array(ipZ),nAlpha,nBeta)
                 Call RecPrt(' In PAM2Int: P',' ',Array(ipPx),nZeta,3)
                  End If
*
*-----------------Compute the cartesian values of the basis functions
*                 angular part
*
                  ABeq(1) = A(1).eq.RB(1) .and. A(1).eq.TC(1)
                  ABeq(2) = A(2).eq.RB(2) .and. A(2).eq.TC(2)
                  ABeq(3) = A(3).eq.RB(3) .and. A(3).eq.TC(3)
                  Call CrtCmp(Array(ipZ),Array(ipPx),nZeta,A,
     &                        Array(ipAxyz),la,HerR(iHerR(nHer)),
     &                        nHer,ABeq)
                  Call CrtCmp(Array(ipZ),Array(ipPx),nZeta,RB,
     &                        Array(ipBxyz),lb,HerR(iHerR(nHer)),
     &                        nHer,ABeq)
*
*-----------------Compute the contribution from the multipole moment
*                 operator
*
                  ABeq(1) = .False.
                  ABeq(2) = .False.
                  ABeq(3) = .False.
                  Call CrtCmp(Array(ipZ),Array(ipPx),nZeta,TC,
     &                        Array(ipRxyz),nOrdOp,HerR(iHerR(nHer)),
     &                        nHer,ABeq)
*
*-----------------Compute the cartesian components for the multipole
*                 moment integrals. The integrals are factorized into
*                 components.
*
                  Call Assmbl(Array(ipQxyz),
     &                        Array(ipAxyz),la,
     &                        Array(ipRxyz),nOrdOp,
     &                        Array(ipBxyz),lb,
     &                        nZeta,HerW(iHerW(nHer)),nHer)
*
*-----------------Combine the cartesian components to the full one
*                 electron integral.
*
                  Call CmbnMP(Array(ipQxyz),nZeta,la,lb,nOrdOp,
     &                        Array(ipZ),Array(ipK),Array(ipRes),nComp)
                  If (iPrint.ge.99) Then
                     Write (6,*) ' Intermediate result in PAM2Int'
                     Do 9101 ia = 1, nElem(la)
                        Do 9201 ib = 1, nElem(lb)
                           iab = (ib-1)*nElem(la) + ia
                           ipab = (iab-1)*nZeta + ipRes
                           Write (Label,'(A,I2,A,I2,A)')
     &                          ' Array(',ia,',',ib,')'
                           If (nComp.ne.1) Then
                              Call RecPrt(Label,' ',
     &                                    Array(ipab),nZeta,nComp)
                           Else
                              Call RecPrt(Label,' ',
     &                                    Array(ipab),nAlpha,nBeta)
                           End If
 9201                   Continue
 9101                Continue
                  End If
*
*-----------------Multiply result by Zeff*Const
*
                  Factor = -dbsc(kCnttp)%Charge*PAMexp(iM2xp,2)
     &                   * Fact
*
*                 FOR DMFT calculation!!!
*
c                  write(6,*) ' Cff',PAMexp(iM2xp,2)
                  Factor = 1.00d0*Fact*PAMexp(iM2xp,2)
                  If (iPrint.ge.99) Write (6,*) ' Factor=',Factor
                  Call DaXpY_(nZeta*nElem(la)*nElem(lb)*nComp,Factor,
     &                       Array(ipRes),1,Work(ipScr),1)
*
 1011          Continue
*
*-----------------Accumulate contributions
*
            nOp = NrOpr(iDCRT(lDCRT))
            Call SymAdO(Work(ipScr),nZeta,la,lb,nComp,Final,
     &                 nIC,nOp,lOper,iChO,One)
            Call GetMem(' Scr','FREE','REAL',ipScr,
     &                       nZeta*nElem(la)*nElem(lb)*nComp)
*
 102        Continue
 101     Continue
 111     kdc = kdc + dbsc(kCnttp)%nCntr
*
c      If (nOrdOp.eq.1) Then
      If (iPrint.ge.99) Then
         Write (6,*) ' Result in PAM2Int'
         Do 9100 ia = 1, nElem(la)
            Do 9200 ib = 1, nElem(lb)
               Write (Label,'(A,I2,A,I2,A)')
     &              ' Final(ia=',ia,',ib=',ib,')'
               Call RecPrt(Label,' ',Final(1,ia,ib,1),nAlpha,nBeta)
 9200       Continue
 9100    Continue
      End If
*
*     Call GetMem(' Exit PAM2Int','LIST','REAL',iDum,iDum)
*     Call QExit('PAM2Int')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
         Call Unused_real_array(Alpha)
         Call Unused_real_array(Beta)
         Call Unused_real_array(ZInv)
      End If
      End
