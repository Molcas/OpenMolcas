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
* Copyright (C) 1993, Bernd Artur Hess                                 *
*               1999, Roland Lindh                                     *
************************************************************************
      SubRoutine pXpInt(
#define _CALLING_
#include "int_interface.fh"
     &                 )
************************************************************************
*                                                                      *
* Object: kernel routine for the comutation of pXp integrals           *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : qEnter                                                  *
*              RecPrt                                                  *
*              pvint                                                   *
*              Util2                                                   *
*              DCopy   (ESSL)                                          *
*              GetMem                                                  *
*              qExit                                                   *
*                                                                      *
*     Author: Bernd Hess, Institut fuer Physikalische und Theoretische *
*             Chemie, University of Bonn, Germany, April 1993          *
*             R. Lindh, modified to molcas 4.1 form, Oct 1999          *
************************************************************************
      use Symmetry_Info, only: nIrrep, iChBas
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"

#include "int_interface.fh"

*     Local variables
      Parameter (mComp=200)
      Integer kOper(mComp), kChO(mComp)
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = ((ixyz+1)*(ixyz+2))/2
*
      iRout = 220
      iPrint = nPrint(iRout)
*
      iSize=nZeta*nElem(la)*nElem(lb)*nComp
      call dcopy_(iSize,[Zero],0,Final,1)
      call dcopy_(nZeta*nArr,[Zero],0,Array,1)
      nip = 1
      ipB = nip
      nip = nip + nZeta
      ipS1 = nip
      nip = nip + nZeta*nElem(la)*nElem(lb+1)*3*nIC
      If (lb.gt.0) Then
         ipS2 = nip
         nip = nip + nZeta*nElem(la)*nElem(lb-1)*3*nIC
      Else
         ipS2=ipS1
      End If
      ipArr = nip
      mArr=nArr-(nip-1)/nZeta
      If (mArr.lt.0) Then
          Call WarningMessage(2,'pXpInt: mArr<0!')
          Call Abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     nIC: the number of blocks of the symmetry adapted operator pXp
*     nComp: number of components of the operator pXp
*
*     pXp = d/dx X d/dx + d/dy X d/dy + d/dz X d/dz
*
      kIC=nIC*3
      kComp=nComp*3
      kOrdOp = nOrdOp-1
      If (kComp.gt.mComp) Then
         Write (6,*) 'pxpint: kComp.gt.mComp'
         Call Abend()
      End If
C     Write (6,*)
C     Write (6,*) 'pXpInt:**********'
*
      iSym_p1 = IrrFnc(1)
      iSym_p2 = IrrFnc(2)
      iSym_p3 = IrrFnc(4)
C     Write (6,*) 'iSym_p=',iSym_p1,iSym_p2,iSym_p3
      ipar_p1 = iChBas(2)
      ipar_p2 = iChBas(3)
      ipar_p3 = iChBas(4)
C     Write (6,*) 'ipar_p=',ipar_p1,ipar_p2,ipar_p3
      Do iComp = 1, nComp
         jComp1 = (iComp-1)*3 + 1
         jComp2 = (iComp-1)*3 + 2
         jComp3 = (iComp-1)*3 + 3
         iTemp = lOper(iComp)
         ipar  = iChO(iComp)
*
         jTemp1= 0
         jTemp2= 0
         jTemp3= 0
         Do iSym_pXp = 0, nIrrep-1
            If (iAnd(2**iSym_pXp,iTemp).ne.0) Then
               iSym_pX = iEor(iSym_pXp,iSym_p1)
C              Write (6,*) 'iSym_pXp,iSym_pX=',iSym_pXp,iSym_pX
               jTemp1=iOr(jTemp1,2**iSym_pX)
            End If
            If (iAnd(2**iSym_pXp,iTemp).ne.0) Then
               iSym_pX = iEor(iSym_pXp,iSym_p2)
C              Write (6,*) 'iSym_pXp,iSym_pX=',iSym_pXp,iSym_pX
               jTemp2=iOr(jTemp2,2**iSym_pX)
            End If
            If (iAnd(2**iSym_pXp,iTemp).ne.0) Then
               iSym_pX = iEor(iSym_pXp,iSym_p3)
C              Write (6,*) 'iSym_pXp,iSym_pX=',iSym_pXp,iSym_pX
               jTemp3=iOr(jTemp3,2**iSym_pX)
            End If
         End Do
         kOper(jComp1)=jTemp1
         kOper(jComp2)=jTemp2
         kOper(jComp3)=jTemp3
*
         kChO(jComp1)=iEOr(ipar,ipar_p1)
         kChO(jComp2)=iEOr(ipar,ipar_p2)
         kChO(jComp3)=iEOr(ipar,ipar_p3)
*
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Call pXint(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &           Array(ipS1),nZeta,kIC,kComp,la,lb+1,A,RB,iDum,
     &           Array(ipArr),mArr,CCoor,kOrdOp,kOper,kChO,
     &           iStabM,nStabM,
     &           PtChrg,nGrid,iAddPot)
*                                                                      *
************************************************************************
*                                                                      *
      If (lb.gt.0) Then
         Call pXint(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &              Array(ipS2),nZeta,kIC,kComp,la,lb-1,A,RB,iDum,
     &              Array(ipArr),mArr,CCoor,kOrdOp,kOper,kChO,
     &              iStabM,nStabM,
     &              PtChrg,nGrid,iAddPot)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      ipOff = ipB
      Do iAlpha = 1, nAlpha
         call dcopy_(nBeta,Beta,1,Array(ipOff),nAlpha)
         ipOff = ipOff + 1
      End Do
*
      If (iPrint.ge.99) Then
         Call RecPrt(' In pXpint: Beta (expanded)','(5D20.13)',
     *         Array(ipB),nZeta,1)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Combine pX integrals to generate the pXp integrals.
*
*     Note that the pX integrals have 3*nComp components.
*
      Call Ass_pXp(Array(ipB),nZeta,Final,la,lb,Array(ipS1),Array(ipS2),
     &             nComp)
*                                                                      *
************************************************************************
*                                                                      *
      If (iPrint.ge.49) Call RecPrt('pXpInt: Final',' ',Final,
     &                               nZeta,nElem(la)*nElem(lb))
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nHer)
      End
