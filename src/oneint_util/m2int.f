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
      SubRoutine M2Int(
#define _CALLING_
#include "int_interface.fh"
     &                )
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of M2 integrals used in   *
*         ECP calculations. The operator is a s-type gaussian          *
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
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"

#include "int_interface.fh"

*     Local variables.
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
*     Call GetMem(' Enter M2Int','LIST','REAL',iDum,iDum)
*
      nip = 1
      ipAxyz = nip
      nip = nip + nZeta*3*nHer*(la+1)
      ipBxyz = nip
      nip = nip + nZeta*3*nHer*(lb+1)
      ipRxyz = nip
      nip = nip + nZeta*3*nHer
      ipQxyz = nip
      nip = nip + nZeta*3*(la+1)*(lb+1)
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
         Call WarningMessage(2,'M2Int: nip-1.gt.nArr*nZeta')
         Write (6,*) ' nArr is Wrong! ', nip-1,' > ',nArr*nZeta
         Write (6,*) ' Abend in M2Int'
         Call Abend()
      End If
*
      If (iPrint.ge.49) Then
         Call RecPrt(' In M2Int: A',' ',A,1,3)
         Call RecPrt(' In M2Int: RB',' ',RB,1,3)
         Call RecPrt(' In M2Int: Ccoor',' ',Ccoor,1,3)
         Call RecPrt(' In M2Int: Kappa',' ',rKappa,nAlpha,nBeta)
         Call RecPrt(' In M2Int: Zeta',' ',Zeta,nAlpha,nBeta)
         Call RecPrt(' In M2Int: P',' ',P,nZeta,3)
         Write (6,*) ' In M2Int: la,lb,nHer=',la,lb,nHer
      End If
*
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,Final,1)
*
*-----Loop over nuclear centers
*
      kdc=0
      Do 100 kCnttp = 1, nCnttp
         If (.Not.dbsc(kCnttp)%ECP) Go To 111
         If (dbsc(kCnttp)%nM2.eq.0) Go To 111
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
               Do 1011 iM2xp = 1, dbsc(kCnttp)%nM2
                  Gamma = dbsc(kCnttp)%M2xp(iM2xp)
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
                     Call RecPrt(' In M2Int: Kappa',' ',
     &                            Array(ipK),nAlpha,nBeta)
                     Call RecPrt(' In M2Int: Zeta',' ',
     &                            Array(ipZ),nAlpha,nBeta)
                     Call RecPrt(' In M2Int: P',' ',Array(ipPx),nZeta,3)
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
                     Write (6,*) ' Intermediate result in M2Int'
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
                  Factor = -dbsc(kCnttp)%Charge*dbsc(kCnttp)%M2cf(iM2xp)
     &                   * Fact
                  If (iPrint.ge.99) Write (6,*) ' Factor=',Factor
                  Call DaXpY_(nZeta*nElem(la)*nElem(lb)*nIC,Factor,
     &                       Array(ipRes),1,Final,1)
*
 1011          Continue
*
 102        Continue
 101     Continue
 111     kdc = kdc + dbsc(kCnttp)%nCntr
*
 100  Continue
*
      If (iPrint.ge.99) Then
         Write (6,*) ' Result in M2Int'
         Do 9100 ia = 1, nElem(la)
            Do 9200 ib = 1, nElem(lb)
               Write (Label,'(A,I2,A,I2,A)')
     &              ' Final(ia=',ia,',ib=',ib,')'
               Call RecPrt(Label,' ',Final(1,ia,ib,1),nAlpha,nBeta)
 9200       Continue
 9100    Continue
      End If
*
*     Call GetMem(' Exit M2Int','LIST','REAL',iDum,iDum)
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Alpha)
         Call Unused_real_array(Beta)
         Call Unused_real_array(ZInv)
         Call Unused_integer_array(lOper)
         Call Unused_integer_array(iChO)
         Call Unused_real_array(PtChrg)
         Call Unused_integer(iAddPot)
      End If
      End
