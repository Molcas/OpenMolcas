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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      SubRoutine Cllct(Strng,Vector,Value,nAtom,Coor,nCntr,mCntr,
     &                 xyz,Temp,Ind,Type,qMss,TMtrx,First,Lbl,
     &                 lWrite,Deg,lAtom)
************************************************************************
*     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
*             University of Lund, SWEDEN                               *
*             May '91                                                  *
************************************************************************
      use Symmetry_Info, only: nIrrep, iOper
      use Slapaf_Info, only: dMass, AtomLbl
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "Molcas.fh"
      Character*(*) Strng
      Character Label*(LENIN5), Name*(LENIN)
      Character Oper*3, Type*6, Lbl*8
      Real*8 Coor(3,nAtom), Vector(3,nAtom), xyz(3,nCntr+mCntr),
     &       Temp(3,nCntr+mCntr), qMss(nCntr+mCntr),
     &       TMtrx(3,nAtom,3,(nCntr+mCntr)), Axis(3), Perp_Axis(3,2)
      Integer   Ind(nCntr+mCntr,2)
      Logical First, lWrite, ldB, lWarn, lAtom(nAtom)
      Dimension Dummy(1)
*
      iRout = 50
      iPrint = nPrint(iRout)
      ldB=.False.
      lWrite = First
      lWarn  = First
      If (iPrint.gt.20) lWrite = .True.
      If (iPrint.ge.99) Call RecPrt(' In Cllct: Coor',' ',
     &                               Coor,3,nAtom)
*
      iFrst = 1
      nCent = nCntr+mCntr
*
*     Pick up cartesian coordinates associated with the
*     internal coordinate
*
      Do 10 ixyz = 1, nCent
         Call NxtWrd(Strng,iFrst,iEnd)
         Label = Strng(iFrst:iEnd)
         nPar1 = Index(Label,'(')
         nPar2 = Index(Label,')')
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
            Do 11 j = 1, nIrrep-1
               If (iPhase.eq.iOper(j)) i = j
 11         Continue
            If (i.eq.0) Then
               Call WarningMessage(2,'Error in Cllclt')
               Write (6,*) '*********** ERROR ***********'
               Write (6,*) ' Undefined symmetry operator '
               Write (6,*) '*****************************'
               Write (6,*) ' ',Oper
               Write (6,*)
               Call Quit_OnUserError()
            End If
         Else If (nPar1.eq.0 .and. nPar2.eq.0) Then
            Name = '    '
            Name = Strng(iFrst:iEnd)
            Oper = ' '
         Else
            Call WarningMessage(2,'Error in Cllclt')
            Write (6,*) '********** ERROR **********'
            Write (6,*) ' Syntax error in:'
            Write (6,*) ' ',Label
            Write (6,*) '***************************'
            Write (6,*)
            Call Quit_OnUserError()
         End If
*
*------- Find corresponding coordinate
*
         jsAtom = 0
         Do 20 isAtom = 1, nAtom
            If (Name.eq.AtomLbl(isAtom)) jsAtom = isAtom
 20      Continue
         If (jsAtom.eq.0) Then
            Call WarningMessage(2,'Error in Cllclt')
            Write (6,*) '********** ERROR **********'
            Write (6,*) ' Unrecognizable atom label '
            Write (6,*) ' ',Name
            Write (6,*) '***************************'
            Write (6,*)
            Call Quit_OnUserError()
         End If
         lAtom(jsAtom) = .False.
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
         If (Type.eq.'DISSOC') qMss(ixyz) = dMass(jsAtom)
         iFrst = iEnd + 1
 10   Continue
*
      If (iPrint.ge.99) Call RecPrt(' Coordinates',' ',
     &                              xyz,3,nCent)
*                                                                      *
************************************************************************
*                                                                      *
*
*---- Process the internal coordinate
*
      If (Type.eq.'X     ') Then
         Value = xyz(1,1)
         call dcopy_(3,[Zero],0,Temp,1)
         Temp(1,1) = One
         If (lWrite) Write (6,'(1X,A,A,2X,F10.4,A)') Lbl,
     &          ' : x-component=',Value,'/ bohr'
         Deg=D_Cart(Ind,nIrrep)
      Else If (Type.eq.'Y     ') Then
         Value = xyz(2,1)
         call dcopy_(3,[Zero],0,Temp,1)
         Temp(2,1) = One
         If (lWrite) Write (6,'(1X,A,A,2X,F10.4,A)') Lbl,
     &          ' : y-component=',Value,'/ bohr'
         Deg=D_Cart(Ind,nIrrep)
      Else If (Type.eq.'Z     ') Then
         Value = xyz(3,1)
         call dcopy_(3,[Zero],0,Temp,1)
         Temp(3,1) = One
         If (lWrite) Write (6,'(1X,A,A,2X,F10.4,A)') Lbl,
     &          ' : z-component=',Value,'/ bohr'
         Deg=D_Cart(Ind,nIrrep)
      Else If (Type.eq.'STRTCH') Then
         Call Strtch(xyz,nCent,Value,Temp,lWrite,Lbl,Dummy,ldB)
         Deg=D_Bond(Ind,Ind(1,2),nIrrep)
      Else If (Type.eq.'LBEND1')Then
         Call CoSys(xyz,Axis,Perp_Axis)
         Call LBend(xyz,nCent,Value,Temp,lWrite,Lbl,Dummy,ldB,
     &              Axis,Perp_Axis(1,1),.False.)
         Deg=D_Bend(Ind,Ind(1,2),nIrrep)
      Else If (Type.eq.'LBEND2')Then
         Call CoSys(xyz,Axis,Perp_Axis)
         Call LBend(xyz,nCent,Value,Temp,lWrite,Lbl,Dummy,ldB,
     &              Axis,Perp_Axis(1,2),.True.)
         Deg=D_Bend(Ind,Ind(1,2),nIrrep)
      Else If (Type.eq.'BEND  ')Then
         Call Bend(xyz,nCent,Value,Temp,lWrite,lWarn,Lbl,Dummy,ldB)
         Deg=D_Bend(Ind,Ind(1,2),nIrrep)
      Else If (Type.eq.'TRSN  ')Then
         Call Trsn(xyz,nCent,Value,Temp,lWrite,lWarn,Lbl,Dummy,ldB)
         Deg=D_Trsn(Ind,Ind(1,2),nIrrep)
      Else If (Type.eq.'OUTOFP')Then
         Call OutOfP(xyz,nCent,Value,Temp,lWrite,lWarn,Lbl,Dummy,ldB)
         Deg=D_Trsn(Ind,Ind(1,2),nIrrep)
      Else If (Type.eq.'DISSOC')Then
         Call Dissoc(xyz,nCntr,mCntr,qMss,Value,Temp,lWrite,
     &               Lbl,Dummy,ldB)
         Deg=One
      Else
         Call WarningMessage(2,'Error in Cllclt')
         Write (6,'(A,A)') ' Type declaration is corrupted:',Type
         Call Quit_OnUserError()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Deg=Sqrt(Deg)
*
*-----Symmetrize and normalize
*     Observe that the B matrix is defined in both symmetry adapted
*     internal coordinates and symmetry adapted cartesian coordinates.
*
      If (iPrint.ge.99) Call RecPrt(' Transformation vector',' ',
     &   Temp,3,nCent)
*
*
*     Symmetrization of P: 1/g Sum over G of X(G)*GPG
*
*     i) project away symmetry breaking displacements
*     ii) divide the contribution for each center with
*         the multiplicity, g/u
*     iii) sum over unique centers.
*
      call dcopy_(3*nAtom,[Zero],0,Vector,1)
      call dcopy_(3*nAtom*3*nCent,[Zero],0,TMtrx,1)
      Do 30 ixyz = 1, nCent
         jsAtom = Ind(ixyz,1)
         iPhase = Ind(ixyz,2)
         tx=One
         ty=One
         tz=One
*
*--------Step i
*--------Project away nonsymmetric displacements
*
*--------Restrict loop to the stabilizers of the center.
         Do 350 iIrrep= 0, nIrrep-1
            If (Coor(1,jsAtom).ne.Zero.and.iAnd(iOper(iIrrep),
     &          1).ne.0) Go To 350
            If (Coor(2,jsAtom).ne.Zero.and.iAnd(iOper(iIrrep),
     &          2).ne.0) Go To 350
            If (Coor(3,jsAtom).ne.Zero.and.iAnd(iOper(iIrrep),
     &          4).ne.0) Go To 350
*
            If (iAnd(iOper(iIrrep),1).ne.0) tx=Zero
            If (iAnd(iOper(iIrrep),2).ne.0) ty=Zero
            If (iAnd(iOper(iIrrep),4).ne.0) tz=Zero
 350     Continue
*
*--------Step ii
*
*--------Rotate vector back to the unique center
*
         If (iAnd(iPhase,1).ne.0) tx = - tx
         If (iAnd(iPhase,2).ne.0) ty = - ty
         If (iAnd(iPhase,4).ne.0) tz = - tz
         TMtrx(1,jsAtom,1,ixyz)=tx
         TMtrx(2,jsAtom,2,ixyz)=ty
         TMtrx(3,jsAtom,3,ixyz)=tz
 30   Continue
      Call dGeMV_('N',3*nAtom,3*nCent,
     &           One,TMtrx,3*nAtom,
     &               Temp,1,
     &           Zero,Vector,1)
      If (iPrint.ge.99) Then
         Call RecPrt('TMtrx',' ',TMtrx,3*nAtom,3*nCent)
         Call RecPrt(' symmetry adapted vector',
     &                              ' ',Vector,3,nAtom)
      End If
*
      Return
      End
