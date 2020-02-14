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
* Copyright (C) 1991,1997, Roland Lindh                                *
************************************************************************
      SubRoutine Cllct2(Strng,Vector,dVector,Value,Names,nAtom,Coor,
     &                  nCntr,mCntr,xyz,Grad,Ind,Type,rMss,qMss,Lbl,
     &                  nSym,lWrite,iOper,jStab,nStab,mxdc,Deg,Smmtrc,
     &                  Hess,lIter)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from: DefInt2                                                 *
*                                                                      *
* Calling    : QEnter                                                  *
*              DCopy  (ESSL)                                           *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
*             University of Lund, SWEDEN                               *
*             May '91                                                  *
*                                                                      *
*             Modified to be used in optimizations with constraints,   *
*             June '97 (R. Lindh)                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "Molcas.fh"
      Character*(*) Strng
      Character*(LENIN) Names(nAtom)
      Character*(LENIN5)Label
      Character*(LENIN) Name
      Character Oper*3, Type*6, Lbl*8
      Real*8 Coor(3,nAtom), Vector(3,nAtom), xyz(3,nCntr+mCntr),
     &       Grad(3,nCntr+mCntr), dVector(3,nAtom,3,nAtom),
     &       Axis(3),
     &       Perp_Axis(3,2),rMss(nAtom), qMss(nCntr+mCntr),
     &       Hess(3,nCntr+mCntr,3,nCntr+mCntr)
      Integer   Ind(nCntr+mCntr,2), iOper(0:nSym-1), nStab(mxdc),
     &          jStab(0:7,mxdc), iDCR(MxAtom)

      Logical lWrite, ldB, lWarn, PSPrint, Smmtrc(3,nAtom)
*
      iRout = 50
      iPrint = nPrint(iRout)
      Call qEnter('Cllct2')
      ldB=.True.
      lWarn  = lWrite
      If (iPrint.gt.20) lWrite = .True.
      PSPrint=.False.
      If (iPrint.gt.20) PSPrint = .True.
      If (iPrint.ge.99) Then
         Call RecPrt(' In Cllct2: Coor',' ',
     &                               Coor,3,nAtom)
         Call RecPrt('rMss',' ',rMss,1,nAtom)
      End If
*
      iFrst = 1
      iEnd  = 1
      lStrng=LEN(Strng)
*
*     Pick up cartesian coordinates associated with the
*     internal coordinate
*
      nCent=nCntr+mCntr
      Do ixyz = 1, nCent
         Call NxtWrd(Strng,iFrst,iEnd)
         If (iEnd.ge.iFrst) Then
            Label = Strng(iFrst:iEnd)
            nPar1 = Index(Label,'(')
            nPar2 = Index(Label,')')
         Else
            Label = ' '
            nPar1 = 0
            nPar2 = 0
         End If
         iPhase = 0
         If (nPar1.ne.0 .and. nPar2.ne.0) Then
            Name = '    '
            Name = Label(1:nPar1-1)
            Oper = Label(nPar1+1:nPar2-1)
            Call UpCase(Oper)
            If (Index(Oper,'X').ne.0) iPhase=iEor(iPhase,1)
            If (Index(Oper,'Y').ne.0) iPhase=iEor(iPhase,2)
            If (Index(Oper,'Z').ne.0) iPhase=iEor(iPhase,4)

*
*---------- Check if operator belongs to the current point group
*
            i = 0
            Do j = 1, nSym-1
               If (iPhase.eq.iOper(j)) i = j
            End Do
            iDCR(ixyz)=iOper(i)
            If (i.eq.0) Then
               Call WarningMessage(2,' Undefined symmetry operator')
               Write (6,'(A)') Oper
               Call Quit_OnUserError()
            End If
            iFrst = iEnd + 1
         Else If (nPar1.eq.0 .and. nPar2.eq.0) Then
            Name = '    '
            Oper = ' '
            If (iEnd.ge.iFrst) Then
               Name = Strng(iFrst:iEnd)
               If (iEnd.ge.1.and.iEnd.lt.lStrng) iFrst = iEnd + 1
            End If
            iDCR(ixyz)=iOper(0)
         Else
            Call WarningMessage(2,' Syntax error in:'//Label)
            Call Quit_OnUserError()
         End If
*
*------- Find corresponding coordinate
*
         If (Type(1:5).ne.'EDIFF'   .and.
     &       Type(1:3).ne.'NAC'     .and.
     &       Type(1:6).ne.'SPHERE'  .and.
     &       Type(1:6).ne.'TRANSV'        ) Then
            jsAtom = 0
            Do isAtom = 1, nAtom
               If (Name.eq.Names(isAtom)) jsAtom = isAtom
            End Do
            If (jsAtom.eq.0) Then
               Call WarningMessage(2,
     &                 ' Unrecognizable atom label '//Name)
               Call Quit_OnUserError()
            End If
         Else
            jsAtom=ixyz
         End If
*
*--------Store away the unique center index and the operator
*
         Ind(ixyz,1) = jsAtom
         Ind(ixyz,2) = iPhase
         call dcopy_(3,Coor(1,jsAtom),1,xyz(1,ixyz),1)
*--------Generate actual coordinate
         If (iAnd(iPhase,1).ne.0) xyz(1,ixyz) = - xyz(1,ixyz)
         If (iAnd(iPhase,2).ne.0) xyz(2,ixyz) = - xyz(2,ixyz)
         If (iAnd(iPhase,4).ne.0) xyz(3,ixyz) = - xyz(3,ixyz)
         If (Type.eq.'DISSOC') qMss(ixyz) = rMss(jsAtom)
*
      End Do  ! Do ixyz = 1, nCntr+mCntr
*
      If (iPrint.ge.99) Then
         Call RecPrt(' Coordinates',' ',
     &                              xyz,3,nCntr+mCntr)
         Call RecPrt('qMss',' ',qMss,1,nCntr+mCntr)
       End If
*                                                                      *
************************************************************************
*                                                                      *
*
*---- Process the internal coordinate
*
      If (Type.eq.'X     ') Then
         Value = xyz(1,1)
         call dcopy_(3,[Zero],0,Grad,1)
         call dcopy_(9,[Zero],0,Hess,1)
         Grad(1,1) = One
         If (lWrite) Write (6,'(1X,A,A,2X,F10.4,A)') Lbl,
     &          ' : x-component=',Value,'/ bohr'
         Deg=D_Cart(Ind(1,1),nStab,jStab,mxdc,nSym)
      Else If (Type.eq.'Y     ') Then
         Value = xyz(2,1)
         call dcopy_(3,[Zero],0,Grad,1)
         call dcopy_(9,[Zero],0,Hess,1)
         Grad(2,1) = One
         If (lWrite) Write (6,'(1X,A,A,2X,F10.4,A)') Lbl,
     &          ' : y-component=',Value,'/ bohr'
         Deg=D_Cart(Ind(1,1),nStab,jStab,mxdc,nSym)
      Else If (Type.eq.'Z     ') Then
         Value = xyz(3,1)
         call dcopy_(3,[Zero],0,Grad,1)
         call dcopy_(9,[Zero],0,Hess,1)
         Grad(3,1) = One
         If (lWrite) Write (6,'(1X,A,A,2X,F10.4,A)') Lbl,
     &          ' : z-component=',Value,'/ bohr'
         Deg=D_Cart(Ind(1,1),nStab,jStab,mxdc,nSym)
      Else If (Type.eq.'STRTCH') Then
         Call Strtch(xyz,nCntr,Value,Grad,lWrite,Lbl,Hess,ldB)
         Deg=D_Bond(Ind,Ind(1,2),nStab,jStab,mxdc,nSym)
      Else If (Type.eq.'LBEND1')Then
         Call CoSys(xyz,Axis,Perp_Axis)
         Call LBend(xyz,nCntr,Value,Grad,lWrite,lWarn,Lbl,Hess,ldB,
     &              Axis,Perp_Axis(1,1),.False.)
         Deg=D_Bend(Ind,Ind(1,2),nStab,jStab,mxdc,nSym)
      Else If (Type.eq.'LBEND2')Then
         Call CoSys(xyz,Axis,Perp_Axis)
         Call LBend(xyz,nCntr,Value,Grad,lWrite,lWarn,Lbl,Hess,ldB,
     &              Axis,Perp_Axis(1,2),.True.)
         Deg=D_Bend(Ind,Ind(1,2),nStab,jStab,mxdc,nSym)
      Else If (Type.eq.'BEND  ')Then
         Call Bend(xyz,nCntr,Value,Grad,lWrite,lWarn,Lbl,Hess,ldB)
         Deg=D_Bend(Ind,Ind(1,2),nStab,jStab,mxdc,nSym)
      Else If (Type.eq.'TRSN  ')Then
         Call Trsn(xyz,nCntr,Value,Grad,lWrite,lWarn,Lbl,Hess,ldB)
         Deg=D_Trsn(Ind,Ind(1,2),nStab,jStab,mxdc,nSym)
      Else If (Type.eq.'OUTOFP')Then
         Call OutOfP(xyz,nCntr,Value,Grad,lWrite,lWarn,Lbl,Hess,ldB)
         Deg=D_Trsn(Ind,Ind(1,2),nStab,jStab,mxdc,nSym)
      Else If (Type(1:3).eq.'NAC')Then
         Call NACInt(xyz,nCntr,Value,Grad,lWrite,lWarn,Lbl,Hess,ldB,
     &               lIter)
         Deg=One
      Else If (Type(1:5).eq.'EDIFF')Then
         Call ConInt(xyz,nCntr,Value,Grad,lWrite,lWarn,Lbl,Hess,ldB,
     &               lIter)
         Deg=One
      Else If (Type(1:6).eq.'SPHERE')Then
         Call SphInt(xyz,nCntr,Value,Grad,lWrite,lWarn,Lbl,Hess,ldB)
         Deg=One
      Else If (Type(1:6).eq.'TRANSV')Then
         Call Transverse(xyz,nCntr,Value,Grad,lWrite,lWarn,Lbl,Hess,ldB)
         Deg=One
      Else If (Type.eq.'DISSOC')Then
         Call Dissoc(xyz,nCntr,mCntr,qMss,nAtom,Value,Grad,lWrite,
     &               Lbl,Hess,ldB)
         Deg=One
      Else
         Call WarningMessage(2,
     &                     ' Type declaration is not supported:'//Type)
         Call Quit_OnUserError()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Deg=Sqrt(Deg)
*
      Call ProjSym2(nAtom,nCent,Ind,nStab,jStab,xyz,
     &             iDCR,Grad,Vector,Smmtrc,PSPrint,Hess,dVector)
      If (iPrint.ge.99) Then
         Call RecPrt(' symmetry adapted vector',
     &                              ' ',Vector,3,nAtom)
         Call RecPrt(' symmetry adapted dvector',
     &                              ' ',dVector,3*nAtom,3*nAtom)
      End If
*
      Call qExit('Cllct2')
      Return
      End
