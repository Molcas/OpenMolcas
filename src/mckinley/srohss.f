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
      SubRoutine SroHss(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
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
*      p     : center of new gaussian from the products of bra and ket *
*              gaussians.                                              *
*      final : array for computed integrals                            *
*      nzeta : nalpha x nbeta                                          *
*      ncomp : number of components in the operator (e.g. dipolmoment  *
*              operator has three components)                          *
*      la    : total angular momentum of bra gaussian                  *
*      lb    : total angular momentum of ket gaussian                  *
*      a     : center of bra gaussian                                  *
*      b     : center of ket gaussian                                  *
*      nrys  : order of rys- or Hermite-Gauss polynomial               *
*      array : auxiliary memory as requested by ECPMem                 *
*      narr  : length of array                                         *
*      ccoor : coordinates of the operator, zero for symmetric oper.   *
*      nordop: order of the operator                                   *
************************************************************************
      use Basis_Info
      use Real_Spherical
      implicit real*8 (a-h,o-z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "disp.fh"
#include "disp2.fh"
      real*8 zeta(nzeta), zinv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rkappa(nzeta), p(nZeta,3), A(3), RB(3),
     &       array(narr), Ccoor(3), C(3), TC(3),Coor(3,4),
     &       dao(nzeta,(la+1)*(la+2)/2*(lb+1)*(lb+2)/2),Hess(nHess),
     &       g2(78)
      integer istabm(0:nstabm-1), iDCRT(0:7), lOper(nComp),
     &          iuvwx(4), nop(2), kOp(4),mop(4),
     &          indgrd(3,2,0:7), JndGrd(3,4,0:7),jndhss(4,3,4,3,0:7),
     &          indhss(2,3,2,3,0:7)
      logical  ifgrd(3,2),jfgrd(3,4),  EQ,
     & jfhss(4,3,4,3),ifhss(2,3,2,3) ,ifg(4),tr(4)
       nelem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
*
      iuvwx(1) = nstab(mdc)
      iuvwx(2) = nstab(ndc)
      call icopy(2,nop,1,mop,1)
      kop(1) = ioper(nop(1))
      kop(2) = ioper(nop(2))
      call dcopy_(3,a,1,coor(1,1),1)
      call dcopy_(3,rb,1,coor(1,2),1)

*
      kdc = 0
      do 1960 kcnttp = 1, ncnttp
         if (.not.ecp(kcnttp)) Go To 1961
         if (nsro_shells(kcnttp).le.0) Go To 1961
         do 1965 kcnt = 1,ncntr(kCnttp)
            ixyz = dbsc(kCnttp)%ipcntr + (kCnt-1)*3
            call dcopy_(3,work(ixyz),1,C,1)
*
            call dcr(lmbdt,ioper,nIrrep,iStabM,nStabM,
     &               jstab(0,kdc+kCnt),nStab(kdc+kCnt),iDCRT,nDCRT)
            fact = dble(nstabm) / DBLE(LmbdT)
*
            iuvwx(3) = nstab(kdc+kCnt)
            iuvwx(4) = nstab(kdc+kCnt)

*
         do 1967 ldcrt = 0, ndcRT-1


            kop(3) = idcrt(ldcrT)
            kop(4) = kop(3)
            mop(3) = nropr(kop(3),ioper,nirrep)
            mop(4) = mop(3)

            tc(1) = DBLE(iphase(1,idCRT(lDCRT)))*C(1)
            tc(2) = DBLE(iphase(2,idCRT(lDCRT)))*C(2)
            tc(3) = DBLE(iphase(3,idCRT(lDCRT)))*C(3)
            call dcopy_(3,tc,1,coor(1,3),1)

            if (eq(a,rb).and.eq(A,TC)) Go To 1967
            call nucind(coor,kdc+kCnt,ifgrd,ifhss,indgrd,indhss,
     &                  jfgrd,jfhss,jndgrd,jndhss,tr,ifg)
            do 1966 iang = 0, nSRO_Shells(kCnttp)-1
               ishll = ipsro(kcnttp) + iAng
               if (nexp(ishll).eq.0) Go To 1966
*
               ip = 1
               ipfin = ip
               ip = ip + nzeta*nElem(la)*nElem(lb)*21
               ipfa1 = ip
               ip = ip + nalpha*nExp(iShll)*nElem(la)*nElem(iAng)*4
               iptmp = ip
               ip = ip + nalpha*nExp(iShll)
               ipfa2 = ip
               ip = ip + nalpha*nExp(iShll)*nElem(la)*nElem(iAng)*6
               ipfb1 = ip
               ip = ip + nexp(iShll)*nBeta*nElem(iAng)*nElem(lb)*4
               ipfb2 = ip
               ip = ip + nexp(iShll)*nBeta*nElem(iAng)*nElem(lb)*6

               call dcopy_(narr,[Zero],0,Array,1)
*              <a|c>, <a'|c>, <a",c>
               Call Acore(iang,la,ishll,nordop,TC,A,Array(ip),
     &                     narr-ip+1,Alpha,nalpha,Array(ipFA1),
     &                     array(ipfa2),jfgrd(1,1),jfhss,
     &                     2,.false.)
*              Transform core orbital to spherical harmonics
               call LToSph(Array(ipFA1),nAlpha,ishll,la,iAng,4)
               call LToSph(Array(ipFA2),nAlpha,ishll,la,iAng,6)

*              <c|b>,<c,b'>,<c|b">
               Call coreB(iang,lb,ishll,nordop,TC,RB,Array(ip),
     &                    narr-ip+1,Beta,nbeta,Array(ipFB1),
     &                    array(ipfb2),jfgrd(1,2),jfhss,
     &                    2,.false.)
*              Transform core orbital to spherical harmonics
               call RToSph(Array(ipFB1),nBeta,ishll,lb,iAng,4)
               call RToSph(Array(ipFB2),nBeta,ishll,lb,iAng,6)

*              Construct complete derivatives (contracting core)
               Call CmbnACB2(Array(ipFa1),Array(ipFa2),Array(ipFb1),
     &                        Array(ipFb2),Array(ipFin),Fact,
     &                        nalpha,nbeta,
     &                        Work(ipAkl(iShll)),nexp(ishll),
     &                        la,lb,iang,jfhss,Array(ipTmp),.true.)


*              contract density
               nt=nZeta*(la+1)*(la+2)/2*(lb+1)*(lb+2)/2
               mvec=21
               call dcopy_(78,[Zero],0,g2,1)
               Call dGeMV_('T',nT,21,
     &                    One,Array(ipFin),nT,
     &                    DAO,1,
     &                    Zero,g2,1)

*              distribute in hessian
               Call Distg2(g2,Hess,nHess,JndGrd,
     &                  JfHss,JndHss,iuvwx,kOp,mop,Tr,IfG)

*
 1966       Continue
 1967    Continue
 1965    Continue
 1961    Continue
         kdc = kdc + dbsc(kCnttp)%nCntr
 1960    Continue
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
