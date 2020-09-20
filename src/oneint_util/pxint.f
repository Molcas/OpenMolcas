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
* Copyright (C) 2006, Roland Lindh                                     *
************************************************************************
      SubRoutine PXInt(
#define _CALLING_
#include "int_interface.fh"
     &                )
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of  pX integrals          *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Author: Roland Lindh, Dept. Chem. Phys., Lund University,            *
*         June 2006                                                    *
************************************************************************
      use Symmetry_Info, only: iChBas
      Implicit Real*8 (A-H,O-Z)
      External NAInt, MltInt, EFInt, CntInt
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "property_label.fh"

#include "int_interface.fh"

*     Local variables
      Parameter (mComp=200)
      Integer kOper(mComp), kChO(mComp)
*                                                                      *
************************************************************************
*                                                                      *
*      Interface
*      Subroutine PVINT(
*#define _CALLING_
*#include "int_interface.fh"
*     &                , Kernel)
*#include "int_interface.fh"
*      External Kernel
*      End Subroutine PVINT
*      End Interface
*                                                                      *
************************************************************************
*                                                                      *
*
      Call QEnter('PXInt')
*                                                                      *
************************************************************************
*                                                                      *
*     nIC: number of symmetry adapted blocks in total for the nComp
*          elements of the compund operator, pX.
*     nComp: is the number of elements of the compund operator
*
*     kIC: number of symmetry adapted blocks in total for the kComp
*          elements of the operator X
*     kComp: is the number of elements of the operator X.
*                                                                      *
************************************************************************
*                                                                      *
*     Note that the p operator's each element is only a basis function
*     of a single irreducible representation. Hence, 3*kIC=nIC.
*
*     In addition if the operator X has kComp elements pX has 3*kComp
*     elements.
*                                                                      *
************************************************************************
*                                                                      *
      nRys = nHer
      kIC=nIC/3
      kComp=nComp/3
      kOrdOp = nOrdOp-1
*                                                                      *
************************************************************************
*                                                                      *
*     Now produce the kOper array with kComp elements from the lOper
*     array. Dito kChO/iChO.
*
*     lOper is an integer which bit pattern indicate to which irreps
*     the component of the operator is a basis function. Note that the
*     operator is not symmetry adapted, i.e. it can be a basis function
*     in more than one irrep.
*
*     iChO is an integer which describe the parity character of the
*     operator with respect to X, Y, and Z coordinates. For example,
*     if the first bit is set this means that the operator change sign
*     under a reflextion  in the yz-plane, etc.
*
      If (kComp.gt.mComp) Then
         Call WarningMessage(2,'PXInt: kComp.gt.mComp')
         Write (6,*) 'kComp=',kComp
         Write (6,*) 'mComp=',mComp
         Call Abend()
      End If
*
*     As we remove the p operator (three of them) X should be the same
*     regardless of if we remove d/dx, d/dy, or d/dz.
*
      iSym_p1 = IrrFnc(1)   ! d/dx
      iSym_p2 = IrrFnc(2)   ! d/dy
      iSym_p3 = IrrFnc(4)   ! d/dz
C     Write (6,*)
C     Write (6,*) 'pXInt******'
C     Write (*,*) 'iSym_p=',iSym_p1,iSym_p2,iSym_p3
      ipar_p1 = iChBas(2)
      ipar_p2 = iChBas(3)
      ipar_p3 = iChBas(4)
C     Write (6,*) 'ipar_p=',ipar_p1,ipar_p2,ipar_p3
      Do iComp = 1, kComp
         jComp1 = (iComp-1)*3 + 1
         jComp2 = (iComp-1)*3 + 2
         jComp3 = (iComp-1)*3 + 3
         iTemp1= lOper(jComp1)
         iTemp2= lOper(jComp2)
         iTemp3= lOper(jComp3)
         jpar_p1 = iChO(jComp1)
         jpar_p2 = iChO(jComp2)
         jpar_p3 = iChO(jComp3)
*
*        Look thru all irreps and check if pX is a basis function in
*        irrep iSym_pX. If so find the symmetry to which X is a basis
*        function
*
         jTemp1= 0
         jTemp2= 0
         jTemp3= 0
C        Write (6,*) 'lOper=',lOper(jComp1),lOper(jComp2),lOper(jComp3)
         Do iSym_pX = 0, nIrrep-1
            If (iAnd(2**iSym_pX,lOper(jComp1)).ne.0) Then
                iSym_X = iEOr(iSym_pX,iSym_p1)
C               Write (6,*) 'iSym_pX,iSym_X=',iSym_pX,iSym_X
                jTemp1=iOr(jTemp1,2**iSym_X)
C               Write (6,*) 'jTemp1=',jTemp1
            End If
            If (iAnd(2**iSym_pX,lOper(jComp2)).ne.0) Then
                iSym_X = iEOr(iSym_pX,iSym_p2)
C               Write (6,*) 'iSym_pX,iSym_X=',iSym_pX,iSym_X
                jTemp2=iOr(jTemp2,2**iSym_X)
C               Write (6,*) 'jTemp2=',jTemp2
            End If
            If (iAnd(2**iSym_pX,lOper(jComp3)).ne.0) Then
                iSym_X = iEOr(iSym_pX,iSym_p3)
C               Write (6,*) 'iSym_pX,iSym_X=',iSym_pX,iSym_X
                jTemp3=iOr(jTemp3,2**iSym_X)
C               Write (6,*) 'jTemp3=',jTemp3
            End If
         End Do
*
*        Check for consistency!
*
         If (jTemp1.ne.jTemp2.or.jTemp1.ne.jTemp3) Then
            Call WarningMessage(2,'PXInt: corrupted jTemps!')
            Write (6,*) 'jTemp1,jTemp2,jTemp3=',
     &                   jTemp1,jTemp2,jTemp3
            Call Abend()
         End If
*
*        Compute the parity of X
*
         jpar_p1=iEOr(jpar_p1,ipar_p1)
         jpar_p2=iEOr(jpar_p2,ipar_p2)
         jpar_p3=iEOr(jpar_p3,ipar_p3)
*
         If (jpar_p1.ne.jpar_p2.or.jpar_p1.ne.jpar_p3) Then
            Call WarningMessage(2,'PXInt: corrupted jpars!')
            Call Abend()
         End If
*
*        Store the data
*
         kOper(iComp)=jTemp1
         kChO(iComp)=jpar_p1
      End Do
*
C     Write (6,*) 'pXpInt'
C     Do iComp = 1, nComp
C        Write (6,*) lOper(iComp), iChO(iComp)
C     End Do
C     Write (6,*)
C     Do iComp = 1, kComp
C        Write (6,*) kOper(iComp), kChO(iComp)
C     End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Compute now the integrals
*
      If (PLabel.eq.'NAInt ') Then
         Call PVInt(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &              Final,nZeta,kIC,kComp,la,lb,A,RB,nRys,
     &              Array,nArr,CCoor,kOrdOp,kOper,kChO,
     &              iStabM,nStabM,
     &              PtChrg,nGrid,iAddPot,
     &              NAInt)
      Else If (PLabel.eq.'MltInt') Then
         Call PVInt(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &              Final,nZeta,kIC,kComp,la,lb,A,RB,nRys,
     &              Array,nArr,CCoor,kOrdOp,kOper,kChO,
     &              iStabM,nStabM,
     &              PtChrg,nGrid,iAddPot,
     &              MltInt)
      Else If (PLabel.eq.'EFInt ') Then
         Call PVInt(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &              Final,nZeta,kIC,kComp,la,lb,A,RB,nRys,
     &              Array,nArr,CCoor,kOrdOp,kOper,kChO,
     &              iStabM,nStabM,
     &              PtChrg,nGrid,iAddPot,
     &              EFInt)
      Else If (PLabel.eq.'CntInt') Then
         Call PVInt(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &              Final,nZeta,kIC,kComp,la,lb,A,RB,nRys,
     &              Array,nArr,CCoor,kOrdOp,kOper,kChO,
     &              iStabM,nStabM,
     &              PtChrg,nGrid,iAddPot,
     &              CntInt)
      Else
         Call WarningMessage(2,'PXInt: Illegal type!')
         Write(6,*) '       PLabel=',PLabel
         Call Abend()
      End If
*
      Call QExit('PXInt')
      Return
      End
