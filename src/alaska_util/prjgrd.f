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
      SubRoutine PrjGrd(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,la,lb,A,RB,nRys,
     &                  Array,nArr,Ccoor,nOrdOp,Grad,nGrad,
     &                  IfGrad,IndGrd,DAO,mdc,ndc,kOp,lOper,nComp,
     &                  iStabM,nStabM)
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of ECP integrals.         *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              DCopy   (ESSL)                                          *
*              ZXia                                                    *
*              SetUp1                                                  *
*              Mlt1                                                    *
*              DGeTMO  (ESSL)                                          *
*              DGEMM_  (ESSL)                                          *
*              DScal   (ESSL)                                          *
*              DGEMM_  (ESSL)                                          *
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
      use Her_RW
      use Real_Spherical
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "disp.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3), Grad(nGrad),
     &       Array(nZeta*nArr), Ccoor(3), C(3), TC(3),
     &       DAO(nZeta,(la+1)*(la+2)/2*(lb+1)*(lb+2)/2)
      Integer iStabM(0:nStabM-1), iDCRT(0:7), lOper(nComp),
     &          iuvwx(4), kOp(2), lOp(4),
     &          IndGrd(3,2), JndGrd(3,4)
      Character*80 Label
      Logical IfGrad(3,2), JfGrad(3,4), TstFnc, TF, ABeq(3), EQ
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      TF(mdc,iIrrep,iComp) = TstFnc(iOper,nIrrep,iCoSet(0,0,mdc),
     &                       nIrrep/nStab(mdc),iChTbl,iIrrep,iComp,
     &                       nStab(mdc))

*
*     Call qEnter('PrjGrd')
      iRout = 192
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.49) Then
         Call RecPrt(' In PrjGrd: Grad',' ',Grad,1,nGrad)
         Call RecPrt(' In PrjGrd: A',' ',A,1,3)
         Call RecPrt(' In PrjGrd: RB',' ',RB,1,3)
         Call RecPrt(' In PrjGrd: Ccoor',' ',Ccoor,1,3)
         Call RecPrt(' In PrjGrd: P',' ',P,nZeta,3)
         Call RecPrt(' In PrjGrd: Alpha',' ',Alpha,nAlpha,1)
         Call RecPrt(' In PrjGrd: Beta',' ',Beta,nBeta,1)
         Write (6,*) ' In PrjGrd: la,lb=',' ',la,lb
      End If
*
      nDAO = nElem(la)*nElem(lb)
      iIrrep = 0
      iuvwx(1) = nStab(mdc)
      iuvwx(2) = nStab(ndc)
      lOp(1) = iOper(kOp(1))
      lOp(2) = iOper(kOp(2))
*
      iComp = 1
      kdc = 0
      Do 1960 kCnttp = 1, nCnttp
         If (.Not.ECP(kCnttp)) Go To 1961
         Do 1965 kCnt = 1,nCntr(kCnttp)
            ixyz = ipCntr(kCnttp) + (kCnt-1)*3
            call dcopy_(3,Work(ixyz),1,C,1)
            If (iPrint.ge.49) Call RecPrt(' In PrjGrd: C',' ',C,1,3)
*
            Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,
     &               jStab(0,kdc+kCnt),nStab(kdc+kCnt),iDCRT,nDCRT)
            Fact = DBLE(nStabM) / DBLE(LmbdT)
*
            iuvwx(3) = nStab(kdc+kCnt)
            iuvwx(4) = nStab(kdc+kCnt)
            Call ICopy(6,IndGrd,1,JndGrd,1)
            Do 10 i = 1, 3
               Do 11 j = 1, 2
                  JfGrad(i,j) = IfGrad(i,j)
 11            Continue
 10         Continue
*
            nDisp = IndDsp(kdc+kCnt,iIrrep)
            Do 220 iCar = 0, 2
               JfGrad(iCar+1,3) = .False.
               iCmp = 2**iCar
               If ( TF(kdc+kCnt,iIrrep,iCmp) .and.
     &              .Not.pChrg(kCnttp) ) Then
                  nDisp = nDisp + 1
                  If (Direct(nDisp)) Then
                     JndGrd(iCar+1,1) = Abs(JndGrd(iCar+1,1))
                     JndGrd(iCar+1,2) = Abs(JndGrd(iCar+1,2))
                     JndGrd(iCar+1,3) = -nDisp
                     JfGrad(iCar+1,1) = .True.
                     JfGrad(iCar+1,2) = .True.
                  Else
                     JndGrd(iCar+1,3) = 0
                  End If
               Else
                  JndGrd(iCar+1,3) = 0
               End If
 220        Continue
            Call ICopy(3,[0],0,JndGrd(1,4),1)
            JfGrad(1,4) = .False.
            JfGrad(2,4) = .False.
            JfGrad(3,4) = .False.
            mGrad = 0
            Do 231 iCar = 1, 3
               Do 232 i = 1, 2
                  If (JfGrad(iCar,i)) mGrad = mGrad + 1
 232           Continue
 231        Continue
            If (iPrint.ge.99) Write (6,*) ' mGrad=',mGrad
            If (mGrad.eq.0) Go To 1965
*
         Do 1967 lDCRT = 0, nDCRT-1
            lOp(3) = iDCRT(lDCRT)
            lOp(4) = lOp(3)
            TC(1) = DBLE(iPhase(1,iDCRT(lDCRT)))*C(1)
            TC(2) = DBLE(iPhase(2,iDCRT(lDCRT)))*C(2)
            TC(3) = DBLE(iPhase(3,iDCRT(lDCRT)))*C(3)
            If (EQ(A,RB).and.EQ(A,TC)) Go To 1967
            Do 1966 iAng = 0, nPrj_Shells(kCnttp)-1
               iShll = ipPrj(kCnttp) + iAng
               If (iPrint.ge.49) Then
                  Write (6,*) 'nExp(iShll)=',nExp(iShll)
                  Write (6,*) 'nBasis(iShll)=',nBasis(iShll)
                  Write (6,*) ' iAng=',iAng
                  Call RecPrt('TC',' ',TC,1,3)
               End If
               If (nExp(iShll).eq.0 .or. nBasis(iShll).eq.0) Go To 1966
*
               ip = 1
               ipF1 = ip
               nac = nElem(la)*nElem(iAng)*4
               ip = ip + nAlpha*nExp(iShll)*nac
               ipP1 = ip
               ip = ip + 3 * nAlpha*nExp(iShll)
               ipZ1 = ip
               ip = ip + nAlpha*nExp(iShll)
               ipK1 = ip
               ip = ip + nAlpha*nExp(iShll)
               ipZI1 = ip
               ip = ip + nAlpha*nExp(iShll)
               If (ip-1.gt.nArr*nZeta) Then
                  Write (6,*) '  ip-1.gt.nArr*nZeta(1) in PrjGrd'
                  Call Abend()
               End If
*
*--------------Calculate Effective center and exponent for <A|core>
*
               Call ZXia(Array(ipZ1),Array(ipZI1),nAlpha,nExp(iShll),
     &                   Alpha,Work(ipExp(iShll)))
               Call SetUp1(Alpha,nAlpha,Work(ipExp(iShll)),nExp(iShll),
     &                     A,TC,Array(ipK1),Array(ipP1),Array(ipZI1))
*
*--------------Calculate Overlap <A|core> and derivative <A'|core>
*
               nHer = ((la+1)+iAng+2)/2
               ipAxyz = ip
               ip = ip + nAlpha*nExp(iShll)*3*nHer*(la+2)
               ipCxyz = ip
               ip = ip + nAlpha*nExp(iShll)*3*nHer*(iAng+1)
               ipRxyz = ip
               ip = ip + nAlpha*nExp(iShll)*3*nHer*(nOrdOp+1)
               ipQ1 = ip
               ip = ip +
     &               nAlpha*nExp(iShll)*3*(la+2)*(iAng+1)*(nOrdOp+1)
               ipA = ip
               ip = ip + nAlpha*nExp(iShll)
               If (ip-1.gt.nArr*nZeta) Then
                  Write (6,*) '  ip-1.gt.nArr*nZeta(1b) in PrjGrd'
                  Call Abend()
               End If
               ABeq(1) = A(1).eq.TC(1)
               ABeq(2) = A(2).eq.TC(2)
               ABeq(3) = A(3).eq.TC(3)
               Call CrtCmp(Array(ipZ1),Array(ipP1),nAlpha*nExp(iShll),
     &                     A,Array(ipAxyz),la+1,HerR(iHerR(nHer)),
     &                     nHer,ABeq)
               Call CrtCmp(Array(ipZ1),Array(ipP1),nAlpha*nExp(iShll),
     &                     TC,Array(ipCxyz),iAng,HerR(iHerR(nHer)),
     &                     nHer,ABeq)
               ABeq(1) = .False.
               ABeq(2) = .False.
               ABeq(3) = .False.
               Call CrtCmp(Array(ipZ1),Array(ipP1),nAlpha*nExp(iShll),
     &                     Ccoor,Array(ipRxyz),nOrdOp,HerR(iHerR(nHer)),
     &                     nHer,ABeq)
               If (iPrint.ge.49) Then
                  Write (6,*) ' Array(ipAxyz)=',
     &             DNrm2_(nAlpha*nExp(iShll)*3*nHer*(la+2),
     &             Array(ipAxyz),1)
                  Write (6,*) ' Array(ipCxyz)=',
     &             DNrm2_(nAlpha*nExp(iShll)*3*nHer*(iAng+1),
     &             Array(ipCxyz),1)
                  Write (6,*) ' Array(ipRxyz)=',
     &             DNrm2_(nAlpha*nExp(iShll)*3*nHer*(nOrdOp+1),
     &             Array(ipRxyz),1)
               End If
               Call Assmbl(Array(ipQ1),
     &                     Array(ipAxyz),la+1,
     &                     Array(ipRxyz),nOrdOp,
     &                     Array(ipCxyz),iAng,
     &                     nAlpha*nExp(iShll),HerW(iHerW(nHer)),nHer)
               iStrt = ipA
               Do 20 iGamma = 1, nExp(iShll)
                  call dcopy_(nAlpha,Alpha,1,Array(iStrt),1)
                  iStrt = iStrt + nAlpha
 20            Continue
               If (iPrint.ge.49) Then
                  Write (6,*) ' Array(ipA)=',
     &            DNrm2_(nAlpha*nExp(iShll),Array(ipA),1)
               End If
               Call rKappa_Zeta(Array(ipK1),Array(ipZ1),
     &                          nExp(iShll)*nAlpha)
               ld=1
               Call CmbnAC(Array(ipQ1),nAlpha*nExp(iShll),la,iAng,
     &                     Array(ipK1),Array(ipF1),
     &                     Array(ipA),JfGrad(1,1),ld,nVecAC)
               If (iPrint.ge.49) Then
                Write (6,*) ' Array(ipQ1)=',
     &          DNrm2_(nAlpha*nExp(iShll)*3*(la+2)*(iAng+1)*(nOrdOp+1),
     &          Array(ipQ1),1)
                Write (6,*) ' Array(ipA)=',
     &          DNrm2_(nAlpha*nExp(iShll),Array(ipA),1)
               End If
               ip = ip - nAlpha*nExp(iShll)
     &            * ( 6 + 3*nHer*(la+2) + 3*nHer*(iAng+1)
     &            + 3*nHer*(nOrdOp+1) + 3*(la+2)*(iAng+1)*(nOrdOp+1) +1)
*
               ipF2 = ip
               ncb = nElem(iAng)*nElem(lb)*4
               ip = ip + nExp(iShll)*nBeta*ncb
               ipP2 = ip
               ip = ip + 3 * nExp(iShll)*nBeta
               ipZ2 = ip
               ip = ip + nExp(iShll)*nBeta
               ipK2 = ip
               ip = ip + nExp(iShll)*nBeta
               ipZI2 = ip
               ip = ip + nExp(iShll)*nBeta
               If (ip-1.gt.nArr*nZeta) Then
                  Write (6,*) '  ip-1.gt.nArr*nZeta(2) in PrjGrd'
                  Call Abend()
               End If
*
*--------------Calculate Effective center and exponent for <core|B>
*
               Call ZXia(Array(ipZ2),Array(ipZI2),nExp(iShll),nBeta,
     &                   Work(ipExp(iShll)),Beta)
               Call SetUp1(Work(ipExp(iShll)),nExp(iShll),Beta,nBeta,
     &                    TC,RB,Array(ipK2),Array(ipP2),Array(ipZI2))
*
*--------------Calculate Overlap <core|B> and <core|B'>
*
               nHer = (iAng+(lb+1)+2)/2
               ipCxyz = ip
               ip = ip + nBeta*nExp(iShll)*3*nHer*(iAng+1)
               ipBxyz = ip
               ip = ip + nBeta*nExp(iShll)*3*nHer*(lb+2)
               ipRxyz = ip
               ip = ip + nBeta*nExp(iShll)*3*nHer*(nOrdOp+1)
               ipQ1 = ip
               ip = ip +
     &               nBeta*nExp(iShll)*3*(iAng+1)*(lb+2)*(nOrdOp+1)
               ipB = ip
               ip = ip + nBeta*nExp(iShll)
               If (ip-1.gt.nArr*nZeta) Then
                  Write (6,*) '  ip-1.gt.nArr*nZeta(2b) in PrjGrd'
                  Call Abend()
               End If
               ABeq(1) = TC(1).eq.RB(1)
               ABeq(2) = TC(2).eq.RB(2)
               ABeq(3) = TC(3).eq.RB(3)
               Call CrtCmp(Array(ipZ2),Array(ipP2),nExp(iShll)*nBeta,
     &                     TC,Array(ipCxyz),iAng,HerR(iHerR(nHer)),
     &                     nHer,ABeq)
               Call CrtCmp(Array(ipZ2),Array(ipP2),nExp(iShll)*nBeta,
     &                     RB,Array(ipBxyz),lb+1,HerR(iHerR(nHer)),
     &                     nHer,ABeq)
               ABeq(1) = .False.
               ABeq(2) = .False.
               ABeq(3) = .False.
               Call CrtCmp(Array(ipZ2),Array(ipP2),nExp(iShll)*nBeta,
     &                     Ccoor,Array(ipRxyz),nOrdOp,HerR(iHerR(nHer)),
     &                     nHer,ABeq)
               If (iPrint.ge.49) Then
                  Write (6,*) ' Array(ipCxyz)=',
     &             DNrm2_(nBeta*nExp(iShll)*3*nHer*(iAng+1),
     &             Array(ipCxyz),1)
                  Write (6,*) ' Array(ipBxyz)=',
     &             DNrm2_(nBeta*nExp(iShll)*3*nHer*(lb+2),
     &             Array(ipBxyz),1)
                  Write (6,*) ' Array(ipRxyz)=',
     &             DNrm2_(nBeta*nExp(iShll)*3*nHer*(nOrdOp+1),
     &             Array(ipRxyz),1)
               End If
               Call Assmbl(Array(ipQ1),
     &                     Array(ipCxyz),iAng,
     &                     Array(ipRxyz),nOrdOp,
     &                     Array(ipBxyz),lb+1,
     &                     nExp(iShll)*nBeta,HerW(iHerW(nHer)),nHer)
               iStrt = ipB
               Do 21 iGamma = 1, nExp(iShll)
                  call dcopy_(nBeta,Beta,1,Array(iStrt),nExp(iShll))
                  iStrt = iStrt + 1
 21            Continue
               If (iPrint.ge.49) Then
                  Write (6,*) ' Array(ipB)=',
     &            DNrm2_(nExp(iShll)*nBeta,Array(ipB),1)
               End If
               Call rKappa_Zeta(Array(ipK2),Array(ipZ2),
     &                          nExp(iShll)*nBeta)
               ld=1
               Call CmbnCB(Array(ipQ1),nExp(iShll)*nBeta,iAng,lb,
     &                     Array(ipK2),Array(ipF2),
     &                     Array(ipB),JfGrad(1,2),ld,nVecCB)
               If (iPrint.ge.49) Then
                 Write (6,*) ' Array(ipQ1)=',
     &           DNrm2_(nExp(iShll)*nBeta*3*(la+2)*(iAng+1)*(nOrdOp+1),
     &           Array(ipQ1),1)
                 Write (6,*) ' Array(ipB)=',
     &           DNrm2_(nExp(iShll)*nBeta,Array(ipB),1)
               End If
               ip = ip - nBeta*nExp(iShll)
     &            * ( 6 + 3*nHer*(lb+2) + 3*nHer*(iAng+1)
     &            + 3*nHer*(nOrdOp+1) + 3*(lb+2)*(iAng+1)*(nOrdOp+1) +1)
               nac = nElem(la)*nElem(iAng)*nVecAC
               ncb = nElem(iAng)*nElem(lb)*nVecCB
               ipTmp = ip
               ip = ip +
     &              Max(nAlpha*nExp(iShll)*nac,nBeta*ncb*nBasis(iShll))
               If (ip-1.gt.nArr*nZeta) Then
                  Write (6,*) '  ip-1.gt.nArr*nZeta(3) in PrjGrd'
                  Call Abend()
               End If
               nac = nElem(la)*nElem(iAng)
               ncb = nElem(iAng)*nElem(lb)
*
*--------------Calculate Contraction over components of the core
*              orbitals of type <A|core>Bc<core|B> where we now have in
*              Array(ipF1) the cartesian components of <A|core>, and
*              similarily, in Array(ipF2), we have stored the cartesian
*              components of <core|B>. Observe that the core orbitals
*              orthonomal atomic orbitals. Hence, the transformation
*              to the spherical harmonics has to be for normilized
*              spherical harminics.
*
*--------------From the lefthandside overlap, form iKaC from ikac by
*              1) i,kac -> k,aci
*
               Call DgeTMo(Array(ipF1),nAlpha,nAlpha,
     &                     nExp(iShll)*nac*nVecAC,Array(ipTmp),
     &                     nExp(iShll)*nac*nVecAC)
*
*--------------2) aciK =  k,aci * k,K (Contract over core orbital)
*
               Call DGEMM_('T','N',
     &                     nac*nVecAC*nAlpha,nBasis(iShll),nExp(iShll),
     &                     1.0d0,Array(ipTmp),nExp(iShll),
     &                     Work(ipCff(iShll)),nExp(iShll),
     &                     0.0d0,Array(ipF1),nac*nVecAC*nAlpha)
*
*--------------3) Mult by shiftoperators aci,K -> Bk(K) * aci,K
*
               Do 1955 iBk = 0, nBasis(iShll)-1
                  Call DYaX(nac*nVecAC*nAlpha,Work(ipBk(iShll)+iBk),
     &                       Array(iBk*nac*nVecAC*nAlpha+ipF1),1,
     &                       Array(iBk*nac*nVecAC*nAlpha+ipTmp),1)
 1955          Continue
*
*--------------4) a,ciK -> ciKa
*
               Call DgeTMo(Array(ipTmp),nElem(la),nElem(la),
     &                     nElem(iAng)*nVecAC*nAlpha*nBasis(iShll),
     &                     Array(ipF1),
     &                     nElem(iAng)*nVecAC*nAlpha*nBasis(iShll))
*
*--------------5) iKa,C = c,iKa * c,C
*
               Call DGEMM_('T','N',
     &              nVecAC*nAlpha*nBasis(iShll)*nElem(la),
     &              (2*iAng+1),nElem(iAng),
     &              1.0d0,Array(ipF1),nElem(iAng),
     &              RSph(ipSph(iAng)),nElem(iAng),
     &              0.0d0,Array(ipTmp),
     &              nVecAC*nAlpha*nBasis(iShll)*nElem(la))
*
               Call DgeTMo(Array(ipTmp),nVecAC,nVecAC,
     &                     nAlpha*nBasis(iShll)*nElem(la)*(2*iAng+1),
     &                     Array(ipF1),
     &                     nAlpha*nBasis(iShll)*nElem(la)*(2*iAng+1))
*
*--------------And (almost) the same thing for the righthand side, form
*              KjCb from kjcb
*              1) jcb,K = k,jcb * k,K
*
               Call DGEMM_('T','N',
     &                     nBeta*ncb*nVecCB,nBasis(iShll),nExp(iShll),
     &                     1.0d0,Array(ipF2),nExp(iShll),
     &                     Work(ipCff(iShll)),nExp(iShll),
     &                     0.0d0,Array(ipTmp),nBeta*ncb*nVecCB)
*
*--------------2)  j,cbK -> cbK,j
*
               Call DgeTMo(Array(ipTmp),nBeta,nBeta,
     &                     ncb*nVecCB*nBasis(iShll),Array(ipF2),
     &                     ncb*nVecCB*nBasis(iShll))
*
*--------------3) bKj,C = c,bKj * c,C
*
               Call DGEMM_('T','N',
     &                     nElem(lb)*nVecCB*nBasis(iShll)*nBeta,
     &                     (2*iAng+1),nElem(iAng),
     &                     1.0d0,Array(ipF2),nElem(iAng),
     &                     RSph(ipSph(iAng)),nElem(iAng),
     &                     0.0d0,Array(ipTmp),
     &                     nElem(lb)*nVecCB*nBasis(iShll)*nBeta)
*
*--------------4) b,KjC -> KjC,b
*
               Call DgeTMo(Array(ipTmp),nElem(lb)*nVecCB,
     &                     nElem(lb)*nVecCB,
     &                     nBasis(iShll)*nBeta*(2*iAng+1),Array(ipF2),
     &                     nBasis(iShll)*nBeta*(2*iAng+1))
*
*--------------Next Contract (iKaC)*(KjCb) over K and C, producing ijab,
*              by the following procedure:
*              Loop over a and b
*                Loop over C
*                  Contract iK(aC)*Kj(Cb), over K producing ij(aCb),
*                    accumulate to ij(ab)
*                End loop C
*              End Loop b and a
*
               call dcopy_(nZeta*nElem(la)*nElem(lb)*6,[Zero],0,Final,1)
*
               mVec = 0
               mVecAC = 1
               mVecCB = 1
               Do 900 iCar = 1, 3
                  Do 901 iCent = 1, 2
                     If (JfGrad(iCar,iCent)) Then
                        mVec = mVec + 1
                        If (iCent.eq.1) Then
                           mVecAC = mVecAC+1
                           ipF1a = ipF1 + (mVecAC-1) *
     &                        nAlpha*nBasis(iShll)*nElem(la)*(2*iAng+1)
                           ipF2a = ipF2
                        Else
                           ipF1a = ipF1
                           mVecCB = mVecCB+1
                           ipF2a = ipF2 + (mVecCB-1) *
     &                       nBasis(iShll)*nBeta*(2*iAng+1)*nElem(lb)
                        End If
*
               Do 1030 ib = 1, nElem(lb)
                  Do 1031 ia = 1, nElem(la)
*
                     Do 1032 iC = 1, (2*iAng+1)
                        iaC = (iC-1)*nElem(la) + ia
                        ipaC = (iaC-1)*nAlpha*nBasis(iShll) + ipF1a
                        iCb = (ib-1)*(2*iAng+1) + iC
                        ipCb = (iCb-1)*nBasis(iShll)*nBeta  + ipF2a
*
                        Call DGEMM_('N','N',
     &                             nAlpha,nBeta,nBasis(iShll),
     &                             Fact,Array(ipaC),nAlpha,
     &                                 Array(ipCb),nBasis(iShll),
     &                             One,Final(1,ia,ib,mVec),nAlpha)
*
 1032                Continue
 1031             Continue
 1030          Continue
*
                     End If
 901              Continue
 900           Continue
*
               If (iPrint.ge.49) Then
                  Do 1111 iVec = 1, mVec
                     Write (6,*) iVec,
     &                  Sqrt(DNrm2_(nZeta*nElem(la)*nElem(lb),
     &                  Final(1,1,1,iVec),1))
 1111             Continue
               End If
               If (iPrint.ge.99) Then
                  Write (6,*) ' Result in PrjGrd'
                  Do 100 ia = 1, nElem(la)
                     Do 200 ib = 1, nElem(lb)
                        Do 300 iVec = 1, mVec
                           Write (Label,'(A,I2,A,I2,A)')
     &                           ' Final(',ia,',',ib,')'
                           Call RecPrt(Label,' ',Final(1,ia,ib,iVec),
     &                                 nAlpha,nBeta)
 300                    Continue
 200                 Continue
 100              Continue
               End If
*
*--------------Distribute contributions to the gradient
*
               Call Distg1X(Final,DAO,nZeta,nDAO,mVec,Grad,nGrad,
     &                      JfGrad,JndGrd,iuvwx,lOp,iChBas,MxFnc,nIrrep)
*
 1966       Continue
 1967    Continue
 1965    Continue
 1961    Continue
         kdc = kdc + nCntr(kCnttp)
 1960 Continue
*
*     Call QExit('PrjGrd')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Zeta)
         Call Unused_real_array(ZInv)
         Call Unused_real_array(rKappa)
         Call Unused_integer(nRys)
         Call Unused_integer_array(lOper)
      End If
      End
