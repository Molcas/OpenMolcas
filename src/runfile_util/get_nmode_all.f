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
      Subroutine Get_NMode_All(Vectors,nVectors,nFreq,nUnique_Atoms,
     &                          Vectors_All,nAll_Atoms,mDisp)
      use Symmetry_Info, only: iChTbl, nIrrep, iOper, Symmetry_Info_Get
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8  Vectors(nVectors), Vectors_All(3*nAll_Atoms*nFreq)
      Integer iGen(3), iCoSet(0:7,0:7), mDisp(0:7),
     &        iChCar(3), nDisp(0:7), iStab(0:7)
#ifdef _DEBUGPRINT_
      Logical Temp
#endif
      Integer, Save:: Active=0
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Call RecPrt('Vectors',' ',Vectors,1,nVectors)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      If (Active.eq.0) Then
         Call Symmetry_Info_Get()
         Active=1
      End If
*                                                                      *
************************************************************************
*                                                                      *
      nGen=0
      If (nIrrep.eq.2) nGen=1
      If (nIrrep.eq.4) nGen=2
      If (nIrrep.eq.8) nGen=3
      If (nGen.ge.1) iGen(1)=iOper(1)
      If (nGen.ge.2) iGen(2)=iOper(2)
      If (nGen.eq.3) iGen(3)=iOper(4)
#ifdef _DEBUGPRINT_
      Write (6,*) 'nGen=',nGen
      Write (6,*) 'iGen=',(iGen(i),i=1,nGen)
#endif
      Call ChCar(iChCar,iGen,nGen)
*                                                                      *
************************************************************************
*                                                                      *
*     Generate list of symmetry adapted displacements
*
      Call Get_iScalar('Unique atoms',mUnique_Atoms)
      If (mUnique_Atoms.ne.nUnique_Atoms) Then
         Write (6,*) 'Get_NMode_All: mUnique_Atoms.ne.nUnique_Atoms'
         Call Abend()
      End If
      Call Allocate_Work(ipCoor,3*mUnique_Atoms)
      Call Get_dArray('Unique Coordinates',Work(ipCoor),3*mUnique_Atoms)
#ifdef _DEBUGPRINT_
      Write(6,*) 'nVectors,nAll_Atoms,nFreq=',
     &            nVectors,nAll_Atoms,nFreq
      Write (6,*)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*---- Find the number of displacements in each irrep
*
      MaxDCR=0
      Do iIrrep = 0, nIrrep-1
#ifdef _DEBUGPRINT_
         Write (6,*)
         Write (6,*) 'iIrrep=',iIrrep
         Write (6,*)
#endif
         nDisp(iIrrep)=0
         ipTmp=ipCoor
         Do iUnique_Atom = 1, nUnique_Atoms
#ifdef _DEBUGPRINT_
            Write (6,*) 'iUnique_Atom=',iUnique_Atom
#endif
            iChAtom=iChxyz(Work(ipTmp),iGen,nGen)
            ipTmp=ipTmp+3
            Call Stblz(iChAtom,nStab,iStab,MaxDCR,iCoSet)
            nCoSet=nIrrep/nStab
#ifdef _DEBUGPRINT_
            Write (6,*) 'nCoSet=',nCoSet
            Write (6,*) 'iCoSet=',(iCoSet(i,0),i=0,nCoset-1)
            Write (6,*) 'iChAtom=',iChAtom
#endif
            Do  iCar = 0, 2
                iComp=2**iCar
#ifdef _DEBUGPRINT_
                Write (6,*) 'iComp=',iComp
                Temp=TF(iIrrep,iComp)
                Write (6,*) 'TF(iIrrep,iComp)=',Temp
#endif
                If (TF(iIrrep,iComp)) nDisp(iIrrep)=nDisp(iIrrep)+1
            End Do
         End Do
      End Do
#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*) 'Grand Total'
      Write (6,*) 'nDisp=',(nDisp(i),i=0,nIrrep-1)
      Write (6,*)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Loop Irreps
*
      iVector=0
      iVector_all=0
      iFreq=0
      Do iIrrep = 0, nIrrep-1
#ifdef _DEBUGPRINT_
         Write (6,*) 'iIrrep,nDisp(iIrrep)=',iIrrep,nDisp(iIrrep)
#endif
*
      Do iMode = 1, mDisp(iIrrep)
         iFreq=iFreq+1
#ifdef _DEBUGPRINT_
         Write (6,*) 'iMode=',iMode
#endif
*
*        Loop over symmetry unique centers
*
         ipTmp=ipCoor
         Do iUnique_Atom = 1, nUnique_Atoms
*
#ifdef _DEBUGPRINT_
            Write (6,*) 'iUnique_Atom=',iUnique_Atom
#endif
*           Get permutational character of the center
*
#ifdef _DEBUGPRINT_
            Write (6,*) Work(ipTmp),Work(ipTmp+1),Work(ipTmp+2)
#endif
            iChAtom=iChxyz(Work(ipTmp),iGen,nGen)
            ipTmp=ipTmp+3
            Call Stblz(iChAtom,nStab,iStab,MaxDCR,iCoSet)
            nCoSet=nIrrep/nStab
#ifdef _DEBUGPRINT_
            Write (6,*) 'nCoSet=',nCoSet
            Write (6,*) 'iCoSet=',(iCoSet(i,0),i=0,nCoset-1)
            Write (6,*) 'iChAtom=',iChAtom
#endif
*
            iVec=0 ! dummy initialize
            Do iCo = 0, nCoSet-1
               kOp=iCoSet(iCo,0)
#ifdef _DEBUGPRINT_
               Write (6,*) 'iVector_All=',iVector_All
               Write (6,*) 'iVector=',iVector
               Write (6,*) 'iCo,kOp=',iCo,kOp
#endif
               iVec=0
               Do  iCar = 0, 2
                   iComp=2**iCar
                   iVector_All=iVector_All+1
#ifdef _DEBUGPRINT_
                   Write (6,*) 'iCar=',iCar
                   Write (6,*) 'iVector_All=',iVector_All
                   Write (6,*) 'iComp=',iComp
                   Temp=TF(iIrrep,iComp)
C                  Write (6,*) 'TF(iIrrep,iComp)=',Temp
#endif
                   If (TF(iIrrep,iComp)) Then
C                     Write (*,*) 'Belong!'
                      iVec=iVec+1
*
*     In some cases we only want the normal modes of the first irrep! In
*     that case branch out.
*
                      If (iVector+iVec.gt.nVectors) Then
                         Go To 999
                      End If
                      Vec=Vectors(iVector+iVec)
                      XR=DBLE(iPrmt(NrOpr(kOp),iComp))
                      XY=DBLE(iChTbl(iIrrep,NrOpr(kOp)))
                      Vectors_All(iVector_All)=Vec*XR*XY
                   Else
C                     Write (*,*) 'Doesn''t belong!'
                      Vectors_All(iVector_All)=Zero
                   End If
#ifdef _DEBUGPRINT_
                   Write (6,*) 'iVec=',iVec
#endif
               End Do   ! iCar
            End Do      ! iCo
            iVector=iVector+iVec
*
         End Do ! iUnique_Atom
*
#ifdef _DEBUGPRINT_
         Call RecPrt('Normal mode',' ',
     &               Vectors_All((iFreq-1)*3*nAll_Atoms+1),
     &               3,nAll_Atoms)
#endif
      End Do  ! iMode
      End Do  ! iIrrep
 999  Continue
      Call Free_Work(ipCoor)
#ifdef _DEBUGPRINT_
      Call RecPrt('Normal mode',' ',
     &             Vectors_All,
     &             3*nAll_Atoms,nFreq)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
*
      Contains
      Logical Function TF(iIrrep,iComp)
      Implicit Real*8 (a-h,o-z)
      Logical, External :: TstFnc
      TF = TstFnc(iCoSet,iIrrep,iComp,nIrrep/nCoSet)
      End Function TF
      End
