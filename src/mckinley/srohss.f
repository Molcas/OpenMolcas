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
      SubRoutine SroHss(
#define _CALLING_
#include "hss_interface.fh"
     &                 )
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
      use Center_Info
      use Real_Spherical
      use Symmetry_Info, only: iOper
      implicit real*8 (a-h,o-z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "disp.fh"
#include "disp2.fh"

#include "hss_interface.fh"

*     Local variables
      Real*8 C(3), TC(3), Coor(3,4),  g2(78)
      Integer iDCRT(0:7), iuvwx(4), kOp(4),mop(4),
     &        JndGrd(3,4,0:7),jndhss(4,3,4,3,0:7)
      logical jfgrd(3,4),  EQ, jfhss(4,3,4,3), ifg(4),tr(4)
*
       nelem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      nRys=nHer
*
      iuvwx(1) = dc(mdc)%nStab
      iuvwx(2) = dc(ndc)%nStab
      call icopy(2,nop,1,mop,1)
      kop(1) = ioper(nop(1))
      kop(2) = ioper(nop(2))
      call dcopy_(3,A,1,coor(1,1),1)
      call dcopy_(3,RB,1,coor(1,2),1)

*
      kdc = 0
      do 1960 kcnttp = 1, ncnttp
         if (.not.dbsc(kcnttp)%ECP) Go To 1961
         if (dbsc(kcnttp)%nSRO.le.0) Go To 1961
         do 1965 kcnt = 1,dbsc(kCnttp)%nCntr
            C(1:3)=dbsc(kCnttp)%Coor(1:3,kCnt)
*
            call dcr(lmbdt,iStabM,nStabM,
     &               dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
            fact = dble(nstabm) / DBLE(LmbdT)
*
            iuvwx(3) = dc(kdc+kCnt)%nStab
            iuvwx(4) = dc(kdc+kCnt)%nStab

*
         do 1967 ldcrt = 0, ndcRT-1


            kop(3) = idcrt(ldcrT)
            kop(4) = kop(3)
            mop(3) = nropr(kop(3))
            mop(4) = mop(3)

            Call OA(iDCRT(lDCRT),C,TC)
            call dcopy_(3,TC,1,Coor(1,3),1)

            if (eq(a,rb).and.eq(A,TC)) Go To 1967
            call nucind(coor,kdc+kCnt,ifgrd,ifhss,indgrd,indhss,
     &                  jfgrd,jfhss,jndgrd,jndhss,tr,ifg)
            do 1966 iang = 0, dbsc(kCnttp)%nSRO-1
               ishll = dbsc(kcnttp)%iSRO + iAng
               nExpi=Shells(iShll)%nExp
               if (nExpi.eq.0) Go To 1966
*
               ip = 1
               ipfin = ip
               ip = ip + nzeta*nElem(la)*nElem(lb)*21
               ipfa1 = ip
               ip = ip + nalpha*nExpi*nElem(la)*nElem(iAng)*4
               iptmp = ip
               ip = ip + nalpha*nExpi
               ipfa2 = ip
               ip = ip + nalpha*nExpi*nElem(la)*nElem(iAng)*6
               ipfb1 = ip
               ip = ip + nExpi*nBeta*nElem(iAng)*nElem(lb)*4
               ipfb2 = ip
               ip = ip + nExpi*nBeta*nElem(iAng)*nElem(lb)*6

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
     &                        Shells(iShll)%Akl,nExpi,
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
         Call Unused_real_array(Final)
         Call Unused_integer(nRys)
         Call Unused_real_array(Ccoor)
         Call Unused_integer_array(lOper)
      End If
      End
