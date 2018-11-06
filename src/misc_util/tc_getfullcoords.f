************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine GetFullCoord(Coor,FMass,FAtLbl,nFAtoms,mxAtm,lSlapaf)
      Implicit Real*8 ( a-h,o-z )
      Implicit Integer (i-n)
#include "Molcas.fh"
#include "real.fh"
#include "constants2.fh"
#include "WrkSpc.fh"
      Dimension iOper(8)
      Dimension RotVec(3)
      Character*(LENIN) AtomLbl(mxAtom)
      Character*(LENIN) FAtLbl(mxAtom), Byte4
      Real*8  Coor(3,mxAtom), Mass(MxAtom), FMass(mxAtom), AMass
      Logical lSlapaf
*
      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('Symmetry operations',iOper,nSym)
      Call Get_iScalar('Unique atoms',nAtoms)
      Call Get_cArray('Unique Atom Names',AtomLbl,LENIN*nAtoms)
      Call GetMem('Coor','ALLO','REAL',lw1,3*8*nAtoms)
      If (lSlapaf) then
        Call Get_dArray('Initial Coordinates',Work(lw1),3*nAtoms)
      else
        Call Get_dArray('Unique Coordinates',Work(lw1),3*nAtoms)
      EndIf
*
      Call Get_Mass(Mass,nAtoms)
      Call dScal_(nAtoms,One/uToau,Mass,1)
*
      If (nSym.EQ.1) then
        nFAtoms = nAtoms
        Do i = 0, nFAtoms-1
          FMass(i+1) = Mass(i+1)
          FAtLbl(i+1) = AtomLbl(i+1)
          Coor(1,i+1) = Work(lw1+3*i  )
          Coor(2,i+1) = Work(lw1+3*i+1)
          Coor(3,i+1) = Work(lw1+3*i+2)
        EndDo
        GoTo 9999
      EndIf
*
      lw2=1
      nOper=0
      If ( nSym.eq.2 ) nOper=1
      If ( nSym.eq.4 ) nOper=2
      If ( nSym.eq.8 ) nOper=3
      nCenter=nAtoms
      Do i=1,nOper
        jOper=i+1
        If ( i.eq.3 ) jOper=5
        RotVec(1)=One
        If ( IAND(iOper(jOper),1).eq.1 ) RotVec(1)=-One
        RotVec(2)=One
        If ( IAND(iOper(jOper),2).eq.2 ) RotVec(2)=-One
        RotVec(3)=One
        If ( IAND(iOper(jOper),4).eq.4 ) RotVec(3)=-One
        mCenter=nCenter
        Do iAt=0,mCenter-1
          Xold=Work(lw1+iAt*3+0)
          Yold=Work(lw1+iAt*3+1)
          Zold=Work(lw1+iAt*3+2)
          Byte4=AtomLbl(lw2+iAt)
          AMass=Mass(lw2+iAt)
          FMass(lw2+iAt)=AMass
          Xnew=RotVec(1)*Xold
          Ynew=RotVec(2)*Yold
          Znew=RotVec(3)*Zold
          Do jAt=0,nCenter-1
             If (Byte4.eq.AtomLbl(Lw2+jAt)) Then
                Xold2=Work(lw1+jAt*3+0)
                Yold2=Work(lw1+jAt*3+1)
                Zold2=Work(lw1+jAt*3+2)
                If ( Xnew.eq.Xold2 .and.
     &               Ynew.eq.Yold2 .and.
     &               Znew.eq.Zold2      ) Go To 999
             End If
          End Do
          nCenter=nCenter+1
          Work(lw1+(nCenter-1)*3+0)=Xnew
          Work(lw1+(nCenter-1)*3+1)=Ynew
          Work(lw1+(nCenter-1)*3+2)=Znew
          AtomLbl(lw2+nCenter-1)=Byte4
          Mass(lw2+nCenter-1)=AMass
 999      Continue
        End Do
      End Do
      nFAtoms = nCenter

      Do iAt=0,nCenter-1
        FAtLbl(iAt+1) =  AtomLbl(lw2+iAt)
        FMass(iAt+1)  =  Mass(lw2+iAt)
        Coor(1,iAt+1) =  Work(lw1+3*iAt)
        Coor(2,iAt+1) =  Work(lw1+3*iAt+1)
        Coor(3,iAt+1) =  Work(lw1+3*iAt+2)
      End Do
*
 9999 Call GetMem('Coor','FREE','REAL',lw1,3*8*nAtoms)
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(mxAtm)
      End
