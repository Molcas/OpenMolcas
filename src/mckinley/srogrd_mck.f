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
      SubRoutine SroGrd_mck(
#define _CALLING_
#include "grd_mck_interface.fh"
     &                     )
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
      use Basis_Info
      use Center_Info
      use Real_Spherical
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "disp.fh"

#include "grd_mck_interface.fh"

*     Local variables
      Real*8 C(3), TC(3)
      Integer iDCRT(0:7), iuvwx(4), mOp(4),index(3,4), JndGrd(3,4,0:7)
      Logical JfGrad(3,4), EQ, DiffCnt,tr(4),ifg(4),ifhess_dum(3,4,3,4)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*                                                                      *
************************************************************************
*                                                                      *
      iuvwx(1) = iu
      iuvwx(2) = iv
      mOp(1) = nOp(1)
      mOp(2) = nOp(2)
*
      iprint = 0
      DiffCnt=(IfGrad(iDCar,1).or.IfGrad(iDCar,2))



      If (iPrint.ge.49) Then
            Call RecPrt(' In SROGrd: A',' ',A,1,3)
            Call RecPrt(' In SROGrd: RB',' ',RB,1,3)
            Call RecPrt(' In SROGrd: P',' ',P,nZeta,3)
            Call RecPrt(' In SROGrd: Alpha',' ',Alpha,nAlpha,1)
            Call RecPrt(' In SROGrd: Beta',' ',Beta,nBeta,1)
            Write (6,*) ' In SROGrd: la,lb=',' ',la,lb
            Write (6,*) ' In SROGrd: Diffs=',' ',
     &                    IfGrad(iDCar,1),IfGrad(iDCar,2)
            Write (6,*) ' In SROGrd: Center=',' ',iDCNT
      End If

      kdc = 0
      Do 1960 kCnttp = 1, nCnttp

         If (.Not.dbsc(kCnttp)%ECP) Go To 1961
         If (dbsc(kCnttp)%nSRO.le.0) Go To 1961
         Do 1965 kCnt = 1,dbsc(kCnttp)%nCntr

            If ((.not.DiffCnt).and.((kdc+kCnt).ne.iDCnt)) Goto 1965

            C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)
*
            Call DCR(LmbdT,iStabM,nStabM,
     &               dc(kdc+kCnt)%iStab,dc(kdc+kCnt)%nStab,iDCRT,nDCRT)
            Fact = DBLE(nStabM) / DBLE(LmbdT)
            iuvwx(3) = dc(kdc+kCnt)%nStab
            iuvwx(4) = dc(kdc+kCnt)%nStab

            Call LCopy(12,[.false.],0,JFgrad,1)
            Call LCopy(4,[.false.],0,tr,1)
            Call LCopy(4,[.false.],0,ifg,1)
            Call ICopy(12*nIrrep,[0],0,jndGrd,1)

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
                 Tr(3)=.true.
                 ifg(1)=.true.
                 ifg(2)=.true.
                 JfGrad(iDCar,1) = .true.
                 JfGrad(iDCar,2) = .true.
                 Do iIrrep=0,nIrrep-1
                 jndGrd(iDCar,3,iIrrep) = - IndGrd(iIrrep)
                 End Do
            End If

*
         Do 1967 lDCRT = 0, nDCRT-1

            mOp(3) = nropr(iDCRT(lDCRT))
            mOp(4) = mOp(3)

            Call OA(iDCRT(lDCRT),C,TC)

            If (EQ(A,RB).and.EQ(A,TC)) Go To 1967

            Do 1966 iAng = 0, dbsc(kCnttp)%nSRO-1
               iShll = dbsc(kCnttp)%iSRO + iAng
               nExpi=Shells(iShll)%nExp
               nBasisi=Shells(iShll)%nBasis
               If (iPrint.ge.49) Then
                  Write (6,*) 'nExpi=',nExpi
                  Write (6,*) 'nBasis(iShll)=',nBasisi
                  Write (6,*) ' iAng=',iAng
                  Call RecPrt('TC',' ',TC,1,3)
               End If

               If (nExpi.eq.0) Go To 1966
*
               ip = 1

               ipFin= ip
               ip=ip+nZeta*(la+1)*(la+2)/2*(lb+1)*(lb+2)/2*6

               ipTmp = ip
               ip = ip + MAX(nBeta,nAlpha)*nExpi

               ipFA1 = ip
               ip = ip + nAlpha*nExpi*nElem(la)*nElem(iAng)*2
               ipFA2 = ip ! Not in use for 1st derivative

               ipFB1 = ip
               ip = ip + nExpi*nBeta*nElem(iAng)*nElem(lb)*2

               ipFB2 = ip ! Not in use for 1st derivatives

               call dcopy_(nArr,[Zero],0,Array,1)

               Call Acore(iang,la,ishll,nordop,TC,A,Array(ip),
     &                     narr-ip+1,Alpha,nalpha,Array(ipFA1),
     &                     array(ipFA2),jfgrad(1,1),ifhess_dum,
     &                     1,iprint.ge.49)
               call LToSph(Array(ipFA1),nalpha,ishll,la,iAng,2)



               call dcopy_(nBeta*nExpi*nElem(lb)*nElem(iAng)*2,
     &                    [Zero],0,Array(ipFB1),1)
               Call coreB(iang,lb,ishll,nordop,TC,RB,Array(ip),
     &                    narr-ip+1,Beta,nbeta,Array(ipFB1),
     &                    array(ipFB2),jfgrad(1,2),ifhess_dum,1,
     &                    iprint.ge.49)
               call RToSph(Array(ipFB1),nBeta,ishll,lb,iAng,2)


*
               call CmbnACB1(Array(ipFA1),Array(ipFB1),Array(ipFin),
     &                  Fact,nAlpha,nBeta,Shells(iShll)%Akl,
     &                  nExpi,la,lb,iang,jfgrad,Array(ipTmp),
     &                  .true.,index,mvec,idcar)

*
               nt=nAlpha*nBeta*nElem(lb)*nElem(la)
               Call SmAdNa(Array(ipFin),nt,Final,
     &                     mop,loper,JndGrd,iuvwx,JfGrad,index,
     &                     idcar,1.0d0,iFG,tr)

 1966       Continue
 1967    Continue
 1965    Continue
 1961    Continue
         kdc = kdc + dbsc(kCnttp)%nCntr
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
         Call Unused_logical_array(Trans)
      End If
      End
