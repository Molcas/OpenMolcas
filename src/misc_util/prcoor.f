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
      Subroutine PrCoor
************************************************************************
*                                                                      *
*     purpose: Write coordinates.                                      *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "angstr.fh"
      Integer iGen(3), iPhase(3,0:7), iCoSet(0:7,0:7), iStab(0:7),
     &        iOper(0:7)
c      Real*8 RotVec(3)
      Character*(LENIN) AtomLbl(MxAtom), Byte4
      Data iPhase/ 1, 1, 1,   -1, 1, 1,   1,-1, 1,  -1,-1, 1,
     &             1, 1,-1,   -1, 1,-1,   1,-1,-1,  -1,-1,-1/
*----------------------------------------------------------------------*
*     Prologue                                                         *
*----------------------------------------------------------------------*
      Call qEnter('PrCoor')
*----------------------------------------------------------------------*
*     Read no.of symm. species                                         *
*----------------------------------------------------------------------*
      Call Get_iScalar('nSym',nSym)
*----------------------------------------------------------------------*
*     Read symm. oper per symm. species                                *
*----------------------------------------------------------------------*
      Call Get_iArray('Symmetry operations',iOper,nSym)
*----------------------------------------------------------------------*
*     Read no. of unique atoms in the system                           *
*----------------------------------------------------------------------*
      Call Get_iScalar('Unique atoms',nAtoms)
*----------------------------------------------------------------------*
*     Read atom labels                                                 *
*----------------------------------------------------------------------*
      Call Get_cArray('Unique Atom Names',AtomLbl,LENIN*nAtoms)
*----------------------------------------------------------------------*
*     Read coordinates of atoms                                        *
*----------------------------------------------------------------------*
      Call GetMem('Coor','ALLO','REAL',lw1,3*8*nAtoms)
      Call Get_dArray('Unique Coordinates',Work(lw1),3*nAtoms)
*----------------------------------------------------------------------*
*     Read nuclear repulsion energy                                    *
*----------------------------------------------------------------------*
      Call Get_dScalar('PotNuc',PotNuc)
*----------------------------------------------------------------------*
*     Apply the symmetry operations                                    *
*----------------------------------------------------------------------*
      nGen=0
      If (nSym.eq.2) nGen=1
      If (nSym.eq.4) nGen=2
      If (nSym.eq.8) nGen=3
      If (nGen.ge.1) iGen(1)=iOper(1)
      If (nGen.ge.2) iGen(2)=iOper(2)
      If (nGen.ge.3) iGen(3)=iOper(4)
*
      iAll_Atom=0
      MaxDCR=0
      iAll_Atom=nAtoms
      Do iAtom = 1, nAtoms
         iChAtom=iChxyz(Work(lw1+(iAtom-1)*3),iGen,nGen)
         Call Stblz(iChAtom,iOper,nSym,nStab,iStab,MaxDCR,iCoSet)
         nCoSet=nSym/nStab
         XOld=Work(lw1+(iAtom-1)*3 )
         YOld=Work(lw1+(iAtom-1)*3+1)
         ZOld=Work(lw1+(iAtom-1)*3+2)
         Byte4=AtomLbl(iAtom)

*
         Do iCo = 1, nCoSet-1
*
            iAll_Atom = iAll_Atom + 1
            Work(lw1+(iAll_Atom-1)*3  )=
     &           XOld*DBLE(iPhase(1,iCoSet(iCo,0)))
            Work(lw1+(iAll_Atom-1)*3+1)=
     &           YOld*DBLE(iPhase(2,iCoSet(iCo,0)))
            Work(lw1+(iAll_Atom-1)*3+2)=
     &           ZOld*DBLE(iPhase(3,iCoSet(iCo,0)))
            AtomLbl(iAll_Atom)=Byte4
*
         End Do
*
      End Do
*----------------------------------------------------------------------*
*     Print coordinates of the system                                  *
*----------------------------------------------------------------------*
      Write(6,*)
      Write(6,'(6X,A)')'Cartesian coordinates in Angstrom:'
      Write(6,'(6X,A)')'--------------------------------------------'//
     &                 '---------'
      Write(6,'(6X,A)')'No.  Label        X            Y            '//
     &                 'Z        '
      Write(6,'(6X,A)')'--------------------------------------------'//
     &                 '---------'
      Do iAt=0,iAll_Atom-1
        Write(6,'(4X,I4,3X,A,2X,3F13.8)')
     &  iAt+1,AtomLbl(iAt+1),
     &  Angstr*Work(lw1+3*iAt),
     &  Angstr*Work(lw1+3*iAt+1),
     &  Angstr*Work(lw1+3*iAt+2)
      End Do
      Write(6,'(6X,A)')'--------------------------------------------'//
     &                 '---------'
      Write(6,'(6X,A,F14.8)')'Nuclear repulsion energy =',PotNuc
      Call GetMem('Coor','FREE','REAL',lw1,3*8*nAtoms)
*----------------------------------------------------------------------*
*     Normal exit                                                      *
*----------------------------------------------------------------------*
      Call qExit('PrCoor')
      Return
      End
