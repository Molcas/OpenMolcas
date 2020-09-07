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
      SubRoutine PrjHss(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                 Final,nZeta,la,lb,A,RB,nRys,
     &                 Array,nArr,Ccoor,nOrdOp,Hess,nHess,
     &                 IfHss,IndHss,ifgrd,IndGrd,DAO,mdc,ndc,nOp,
     &                 lOper,nComp,iStabM,nStabM)
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
*             Physics, University of Stockholm, Sweden, October 1993.  *
************************************************************************
      use Basis_Info
      use Real_Spherical
      use Phase_Info
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "disp.fh"
#include "disp2.fh"
      Real*8 Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3),
     &       Array(nArr), Ccoor(3), C(3), TC(3),Coor(3,4),
     &       DAO(nZeta,(la+1)*(la+2)/2*(lb+1)*(lb+2)/2),Hess(nHess),
     &       g2(78)
      Integer iStabM(0:nStabM-1), iDCRT(0:7), lOper(nComp),
     &          iuvwx(4), nOp(2), kOp(4),mop(4),
     &          IndGrd(3,2,0:7), JndGrd(3,4,0:7),jndhss(4,3,4,3,0:7),
     &          indhss(2,3,2,3,0:7)
      Logical  ifgrd(3,2),JfGrd(3,4),  EQ,
     &         jfhss(4,3,4,3),ifhss(2,3,2,3) ,ifg(4),tr(4)
      Dimension Dum(1)

      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

*
      iuvwx(1) = nStab(mdc)
      iuvwx(2) = nStab(ndc)
      call icopy(2,nop,1,mop,1)
      kOp(1) = iOper(nOp(1))
      kOp(2) = iOper(nOp(2))

*
      call dcopy_(3,A,1,Coor(1,1),1)
      call dcopy_(3,RB,1,Coor(1,2),1)

      kdc = 0
      Do 1960 kCnttp = 1, nCnttp
         If (.Not.dbsc(kCnttp)%ECP) Go To 1961
         Do 1965 kCnt = 1,dbsc(kCnttp)%nCntr
            C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)
*
            Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,
     &               jStab(0,kdc+kCnt),nStab(kdc+kCnt),iDCRT,nDCRT)
            Fact = DBLE(nStabM) / DBLE(LmbdT)
*
            iuvwx(3) = nStab(kdc+kCnt)
            iuvwx(4) = nStab(kdc+kCnt)


*
         Do 1967 lDCRT = 0, nDCRT-1
            kOp(3) = iDCRT(lDCRT)
            kOp(4) = kOp(3)
            mop(3) = nropr(kop(3),ioper,nirrep)
            mop(4) = mop(3)
            TC(1) = DBLE(iPhase(1,iDCRT(lDCRT)))*C(1)
            TC(2) = DBLE(iPhase(2,iDCRT(lDCRT)))*C(2)
            TC(3) = DBLE(iPhase(3,iDCRT(lDCRT)))*C(3)
            call dcopy_(3,TC,1,Coor(1,3),1)

            If (EQ(A,RB).and.EQ(A,TC)) Go To 1967
            Call NucInd(coor,kdc+kCnt,ifgrd,ifhss,indgrd,indhss,
     &                  jfgrd,jfhss,jndgrd,jndhss,tr,ifg)
            Do 1966 iAng = 0, dbsc(kCnttp)%nPrj-1
               iShll = dbsc(kCnttp)%iPrj + iAng
               nExpi=Shells(iShll)%nExp
               nBasisi=Shells(iShll)%nBasis
               If (nExpi.eq.0 .or. nBasisi.eq.0) Go To 1966
*
               ip = 1

               ipFin = ip
               ip = ip + nZeta*nElem(la)*nElem(lb)*21

               ipFA1 = ip
               ip = ip + nAlpha*nExpi*nElem(la)*nElem(iAng)*4

               ipFA2 = ip
               ip = ip + nAlpha*nExpi*nElem(la)*nElem(iAng)*6

               ipFB1 = ip
               ip = ip + nExpi*nBeta*nElem(iAng)*nElem(lb)*4

               ipFB2 = ip
               ip = ip + nExpi*nBeta*nElem(iAng)*nElem(lb)*6

               call dcopy_(nArr,[0.0d0],0,Array,1)
*              <a|c>,<a'|c>,<a"|c>
               Call Acore(iang,la,ishll,nordop,TC,A,Array(ip),
     &                     narr-ip+1,Alpha,nalpha,Array(ipFA1),
     &                     array(ipfa2),jfgrd(1,1),jfhss,
     &                     2,.false.)
*              Transform to core orbital
               call LToCore(Array(ipFA1),nalpha,ishll,la,iAng, 4)
               call LToCore(Array(ipFA2),nalpha,ishll,la,iAng, 6)
*              <c|b>,<c|b'>,<c|b">
               Call coreB(iang,lb,ishll,nordop,TC,RB,Array(ip),
     &                    narr-ip+1,Beta,nbeta,Array(ipFB1),
     &                    array(ipfb2),jfgrd(1,2),jfhss,
     &                    2,.false.)
*              Transform to core orbital
               call RToCore(Array(ipFB1),nbeta,ishll,lb,iAng, 4)
               call RToCore(Array(ipFB2),nbeta,ishll,lb,iAng, 6)
*              Construct complete derivative (contract core)
               Call CmbnACB2(Array(ipFa1),Array(ipFa2),Array(ipFb1),
     &                        Array(ipFb2),Array(ipFin),Fact,
     &                        nalpha,nbeta,
     &                        Dum,nBasisi,
     &                        la,lb,iang,jfhss,dum,.false.)

*              contract density
               nt=nZeta*(la+1)*(la+2)/2*(lb+1)*(lb+2)/2
               call dcopy_(78,[0.0d0],0,g2,1)
               Call dGeMV_('T',nT,21,
     &                    One,Array(ipFin),nT,
     &                    DAO,1,
     &                    Zero,g2,1)

*              distribute in hessian
               Call Distg2(g2,Hess,nHess,JndGrd,
     &                  JfHss,JndHss,iuvwx,kOp,mop,
     &                  Tr,IfG)

*
 1966       Continue !iang
 1967    Continue !DCR
 1965    Continue !cnt
 1961    Continue !cont
         kdc = kdc + dbsc(kCnttp)%nCntr
 1960 Continue !cnttp
         Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Zeta)
         Call Unused_real_array(ZInv)
         Call Unused_real_array(rKappa)
         Call Unused_real_array(P)
         Call Unused_real(Final)
         Call Unused_integer(nRys)
         Call Unused_real_array(Ccoor)
         Call Unused_integer_array(lOper)
      End If
         End
