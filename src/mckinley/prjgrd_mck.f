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
************************************************************************
      SubRoutine PrjGrd_mck(Alpha,nAlpha,Beta, nBeta,
     &                  Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,la,lb,A,RB,nHer,
     &                  Array,nArr,Ccoor,nOrdOp,
     &                  IfGrad,IndGrd,nop,
     &                  loper,iu,iv,nrop,idcar,idcnt,
     &                  iStabM,nStabM,ldum)
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
      use Real_Spherical
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "disp.fh"

      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nrop),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3),
     &       Array(nZeta*nArr), Ccoor(3), C(3), TC(3)
      Integer iStabM(0:nStabM-1), iDCRT(0:7), lOper,
     &          iuvwx(4), nOp(2), mOp(4),index(3,4),
     &          indgrd(0:7), JndGrd(3,4,0:7)
      Logical IfGrad(3,2), JfGrad(3,4), EQ,
     &        DiffCnt,tr(4),ifg(4),ifhess_dum(3,4,3,4)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*                                                                      *
************************************************************************
*                                                                      *
      iprint = 0
      iuvwx(1) = iu
      iuvwx(2) = iv
      mop(1) = nOp(1)
      mop(2) = nOp(2)
      DiffCnt=(IfGrad(iDCar,1).or.IfGrad(iDCar,2))
      iComp = 1

      If (iPrint.ge.49) Then
            Call RecPrt(' In PrjGrd: A',' ',A,1,3)
            Call RecPrt(' In PrjGrd: RB',' ',RB,1,3)
            Call RecPrt(' In PrjGrd: P',' ',P,nZeta,3)
            Call RecPrt(' In PrjGrd: Alpha',' ',Alpha,nAlpha,1)
            Call RecPrt(' In PrjGrd: Beta',' ',Beta,nBeta,1)
            Write (6,*) ' In PrjGrd: la,lb=',' ',la,lb
            Write (6,*) ' In PrjGrd: Diffs=',' ',
     &                    IfGrad(iDCar,1),IfGrad(iDCar,2)
            Write (6,*) ' In PrjGrd: Center=',' ',iDCNT
      End If

      kdc = 0
      Do 1960 kCnttp = 1, nCnttp
         If (.Not.ECP(kCnttp)) Go To 1961
         If (nSRO_Shells(kCnttp).le.0) Go To 1961

         Do 1965 kCnt = 1,nCntr(kCnttp)
            If ((.not.DiffCnt).and.((kdc+kCnt).ne.iDCnt)) Goto 1965

            ixyz = ipCntr(kCnttp) + (kCnt-1)*3
            call dcopy_(3,Work(ixyz),1,C,1)
*
            Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,
     &               jStab(0,kdc+kCnt),nStab(kdc+kCnt),iDCRT,nDCRT)
            Fact = DBLE(nStabM) / DBLE(LmbdT)
            iuvwx(3) = nStab(kdc+kCnt)
            iuvwx(4) = nStab(kdc+kCnt)

            Call LCopy(12,.false.,0,JFgrad,1)
            Call LCopy(4 ,.false.,0,iFg,1)
            Call LCopy(4 ,.false.,0,tr,1)
            Call ICopy(12*nIrrep,0,0,jndGrd,1)
            Do iCnt = 1, 2
                  JfGrad(iDCar,iCnt) = IfGrad(iDCar,iCnt)
            End Do
            Do ICnt=1,2
               If (ifgrad(idcar,iCnt)) Then
               ifg(icnt)=.true.
                 Do iIrrep=0,nIrrep-1
                   jndGrd(iDCar,iCnt,iIrrep)=IndGrd(iIrrep)
                 End Do
               End IF
            End Do

*
            If ((kdc+kCnt).eq.iDCnt) Then
                ifg(1)=.true.
                ifg(2)=.true.
                 Tr(3)=.true.
                 JfGrad(iDCar,1) = .true.
                 JfGrad(iDCar,2) = .true.
                 Do iIrrep=0,nIrrep-1
                 jndGrd(iDCar,3,iIrrep) = - IndGrd(iIrrep)
                 End Do
            End If

*
         Do 1967 lDCRT = 0, nDCRT-1

            mop(3)=nropr(iDCRT(lDCRT),ioper,nirrep)
            mop(4)=mop(3)
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
               ip=1

               ipFin= ip
               ip=ip+nZeta*(la+1)*(la+2)/2*(lb+1)*(lb+2)/2*6

               ipFA1 = ip
               ip = ip + nAlpha*nExp(iShll)*nElem(la)*nElem(iAng)*4

               ipFB1 = ip
               ip = ip + nExp(iShll)*nBeta*nElem(iAng)*nElem(lb)*4

               ipFB2 = ip
               ipFA2 = ip
               if (ip.ge.narr) then
                 write(6,*) 'No mem in prjgrd',ip,narr
                 call abend()
                 endif

               call dcopy_(nArr,0.0d0,0,Array,1)

*
               Call Acore(iang,la,ishll,nordop,TC,A,Array(ip),
     &                     narr-ip+1,Alpha,nalpha,Array(ipFA1),
     &                     array(ipFA2),jfgrad(1,1),ifhess_dum,
     &                     1,iprint.ge.49)
               call LToCore(Array(ipFA1),nalpha,ishll,la,iAng, 4)



               Call coreB(iang,lb,ishll,nordop,TC,RB,Array(ip),
     &                    narr-ip+1,Beta,nbeta,Array(ipFB1),
     &                    array(ipFB2),jfgrad(1,2),ifhess_dum,1,
     &                    iprint.ge.49)
               call RToCore(Array(ipFB1),nBeta,ishll,lb,iAng,4)



               call  CmbnACB1(Array(ipFA1),Array(ipFB1),Array(ipFin),
     &                        Fact,nAlpha,nBeta,Dum,nBasis(ishll),
     &                        la,lb,iang,jfgrad,Dum,.false.,
     &                        index,mvec,idcar)


               nt=nAlpha*nBeta*nElem(lb)*nElem(la)
               Call SmAdNa(Array(ipFin),nt,Final,
     &                     mop,loper,JndGrd,iuvwx,JfGrad,index,
     &                     idcar,1.0d0,iFG,tr)


 1966       Continue
 1967    Continue
 1965    Continue
 1961    Continue
         kdc = kdc + nCntr(kCnttp)
 1960 Continue
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Zeta)
         Call Unused_real_array(ZInv)
         Call Unused_real_array(rKappa)
         Call Unused_integer(nHer)
         Call Unused_real_array(Ccoor)
         Call Unused_integer(ldum)
      End If
      End
