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
* Copyright (C) 2005, Giovanni Ghigo                                   *
************************************************************************
      Subroutine Cho_GenC(iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV,
     &                      iAddCou,LenCou, LenEx)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden & Torino University, Italy         *
*           July 2005                                                  *
*----------------------------------------------------------------------*
* This is the routine that really generates the A,B block of coulomb   *
* integrals (in symmetries iSymA,iSymB) for occupied MO iI,iJ in       *
* symmetries iSymI,iSymJ. The 3 x 3 A,B block is built gatering 9      *
* sub-blocks. These are combination of inactive, active, and secondary *
* A,B MO                                                               *
* OBS !!!!!  By now, it works only for iSymA .EQ. iSymB  !!!           *
************************************************************************
*
*   <DOC>
*     <Name>Cho\_GenC</Name>
*     <Syntax>Call Cho\_GenC(iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV,
*     iAddCou,LenCou, LenEx)
*     </Syntax>
*     <Arguments>
*      \Argument{iSymI,iSymJ,iSymA,iSymB}{Symmetry block of the
*      Two-electrons integrals}{Integers}{in}
*      \Argument{NumV}{Number of Cholesky vectors to transform in the
*      current batch}{Integer}{in}
*      \Argument{iAddCou}{Memory pointer of the A,B integrals block}
*      {Integer}{in}
*      \Argument{LenCou}{Length of the A,B integrals block}{Integer}{in}
*     </Arguments>
*     <Purpose>
*      Routine for the generation of the A,B block of coulomb integrals
*      (in symmetries iSymA,iSymB) for occupied MO iI,iJ in symmetries
*      iSymI,iSymJ.\\
*      Called by Cho\_TwoEl
*     </Purpose>
*     <Dependencies>
*     </Dependencies>
*     <Author>
*      G. Ghigo
*     </Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*      The routine generates the A,B block of integrals gatering 9
*      sub-blocks. These are combination of inactive, active, and
*      secondary A,B MO.
*     </Description>
*    </DOC>
*
******************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "rasdim.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "cho_tra.fh"
      Dimension iAddSB(3,3), LenSB(3,3)

CGG   ------------------------------------------------------------------
c      IfTest=.True.
CGG   ------------------------------------------------------------------

      Do iSB_A = 1, 3
        Do iSB_B = 1, 3
          iAddSB(iSB_A,iSB_B)= 0  ! Mem Address of the SubBlocks
          LenSB (iSB_A,iSB_B)= 0  ! Length of the SubBlocks
        EndDo
      EndDo

* --- GENERATION of SubBlocks
CGG   ------------------------------------------------------------------
      If(IfTest .and.SubBlocks(1,1)) then
      Write(6,*)'       SB_11 :',nIsh(iSymA),' x',nIsh(iSymB)
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------
      If (SubBlocks(1,1)) Call MkCouSB11(iAddSB(1,1),LenSB(1,1),
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
CGG   ------------------------------------------------------------------
      If(IfTest .and.SubBlocks(1,1)) then
      Write(6,'(8F10.6)')(Work(iAddSB(1,1)+k),k=0,LenSB(1,1)-1)
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------

CGG   ------------------------------------------------------------------
      If(IfTest .and.SubBlocks(1,2)) then
      Write(6,*)'       SB_12 :',nIsh(iSymA),' x',nAsh(iSymB)
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------
      If (SubBlocks(1,2)) then
        If (iSymA.NE.iSymB) then
          Continue ! CGG
c excluded
c          Call MkCouSB12(iAddSB(1,2),LenSB(1,2),
c     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
        else
          LenSB(1,2)=nIsh(iSymA)*nAsh(iSymB)
          Call GetMem('SB','Allo','Real',iAddSB(1,2),LenSB(1,2))
        Endif
      Endif
CGG   ------------------------------------------------------------------
      If(IfTest .and.SubBlocks(1,2)) then
      Write(6,'(8F10.6)')(Work(iAddSB(1,2)+k),k=0,LenSB(1,2)-1)
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------

CGG   ------------------------------------------------------------------
      If(IfTest .and.SubBlocks(1,3)) then
      Write(6,*)'       SB_13 :',nIsh(iSymA),' x',nSsh(iSymB)
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------
      If (SubBlocks(1,3)) then
        If (iSymA.NE.iSymB) then
          Continue ! CGG
c excluded
c          Call MkCouSB13(iAddSB(1,3),LenSB(1,3),
c     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
        else
          LenSB(1,3)=nIsh(iSymA)*nSsh(iSymB)
          Call GetMem('SB','Allo','Real',iAddSB(1,3),LenSB(1,3))
        Endif
      Endif
CGG   ------------------------------------------------------------------
      If(IfTest .and.SubBlocks(1,3)) then
      Write(6,'(8F10.6)')(Work(iAddSB(1,3)+k),k=0,LenSB(1,3)-1)
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------

CGG   ------------------------------------------------------------------
      If(IfTest .and.SubBlocks(2,1)) then
      Write(6,*)'       SB_21 :',nAsh(iSymA),' x',nIsh(iSymB)
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------
      If (SubBlocks(2,1)) Call MkCouSB21(iAddSB(2,1),LenSB(2,1),
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
CGG   ------------------------------------------------------------------
      If(IfTest .and.SubBlocks(2,1)) then
      Write(6,'(8F10.6)')(Work(iAddSB(2,1)+k),k=0,LenSB(2,1)-1)
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------

CGG   ------------------------------------------------------------------
      If(IfTest .and.SubBlocks(2,2)) then
      Write(6,*)'       SB_22 :',nAsh(iSymA),' x',nAsh(iSymB)
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------
      If (SubBlocks(2,2)) Call MkCouSB22(iAddSB(2,2),LenSB(2,2),
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
CGG   ------------------------------------------------------------------
      If(IfTest .and.SubBlocks(2,2)) then
      Write(6,'(8F10.6)')(Work(iAddSB(2,2)+k),k=0,LenSB(2,2)-1)
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------

CGG   ------------------------------------------------------------------
      If(IfTest .and.SubBlocks(2,3)) then
      Write(6,*)'       SB_23 :',nAsh(iSymA),' x',nSsh(iSymB)
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------
      If (SubBlocks(2,3)) then
        If (iSymA.NE.iSymB) then
          Continue ! CGG
c excluded
c          Call MkCouSB23(iAddSB(2,3),LenSB(2,3),
c     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
        else
          LenSB(2,3)=nAsh(iSymA)*nSsh(iSymB)
          Call GetMem('SB','Allo','Real',iAddSB(2,3),LenSB(2,3))
        Endif
      Endif
CGG   ------------------------------------------------------------------
      If(IfTest .and.SubBlocks(2,3)) then
      Write(6,'(8F10.6)')(Work(iAddSB(2,3)+k),k=0,LenSB(2,3)-1)
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------

CGG   ------------------------------------------------------------------
      If(IfTest .and.SubBlocks(3,1)) then
      Write(6,*)'       SB_31 :',nSsh(iSymA),' x',nIsh(iSymB)
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------
      If (SubBlocks(3,1)) Call MkCouSB31(iAddSB(3,1),LenSB(3,1),
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
CGG   ------------------------------------------------------------------
      If(IfTest .and.SubBlocks(3,1)) then
      Write(6,'(8F10.6)')(Work(iAddSB(3,1)+k),k=0,LenSB(3,1)-1)
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------

CGG   ------------------------------------------------------------------
      If(IfTest .and.SubBlocks(3,2)) then
      Write(6,*)'       SB_32 :',nSsh(iSymA),' x',nAsh(iSymB)
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------
      If (SubBlocks(3,2)) Call MkCouSB32(iAddSB(3,2),LenSB(3,2),
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
CGG   ------------------------------------------------------------------
      If(IfTest .and.SubBlocks(3,2)) then
      Write(6,'(8F10.6)')(Work(iAddSB(3,2)+k),k=0,LenSB(3,2)-1)
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------

CGG   ------------------------------------------------------------------
      If(IfTest .and.SubBlocks(3,3)) then
      Write(6,*)'       SB_33 :',nSsh(iSymA),' x',nSsh(iSymB)
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------
      If(SubBlocks(3,3)) Call MkCouSB33(iAddSB(3,3),LenSB(3,3),
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
CGG   ------------------------------------------------------------------
      If(IfTest .and.SubBlocks(3,3)) then
      Write(6,'(8F10.6)')(Work(iAddSB(3,3)+k),k=0,LenSB(3,3)-1)
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------

CGG   ------------------------------------------------------------------
      If(IfTest) then
      Write(6,*)'     END GENERATION of SubBlocks'
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------
* --- END GENERATION of SubBlocks

* --- GATERING of SubBlocks
      Call GetMem('CSq','Allo','Real',iAddSq,LenEx)
      iAddCouSB = iAddSq
      iLenAi = nIsh(iSymA)
      iLenAt = nAsh(iSymA)
      iLenAa = nSsh(iSymA)
      Do iSB_B = 1, 3
        iLenB = 0
        If (iSB_B.EQ.1) iLenB = nIsh(iSymB)
        If (iSB_B.EQ.2) iLenB = nAsh(iSymB)
        If (iSB_B.EQ.3) iLenB = nSsh(iSymB)
        Do iB = 1,iLenB
          If (iLenAi.GT.0) then                          ! SB(1,iSB_B)
            iAddSBi = iAddSB(iSB_B,1) + iLenAi * (iB-1)
            Call dCopy_(iLenAi,Work(iAddSBi),1,Work(iAddCouSB),1)
            iAddCouSB = iAddCouSB + iLenAi
          EndIf
          If (iLenAt.GT.0) then                          ! SB(2,iSB_B)
            iAddSBi = iAddSB(iSB_B,2) + iLenAt * (iB-1)
            Call dCopy_(iLenAt,Work(iAddSBi),1,Work(iAddCouSB),1)
            iAddCouSB = iAddCouSB + iLenAt
          EndIf
          If (iLenAa.GT.0) then                          ! SB(3,iSB_B)
            iAddSBi = iAddSB(iSB_B,3) + iLenAa * (iB-1)
            Call dCopy_(iLenAa,Work(iAddSBi),1,Work(iAddCouSB),1)
            iAddCouSB = iAddCouSB + iLenAa
          EndIf
        EndDo  ! iB
      EndDo ! iSB_B
      nOrbA = nOrb(iSymA)
CGG   ------------------------------------------------------------------
      If(IfTest) then
      Write(6,*)
      Write(6,*)'        The Square Gatered matrix'
      Call PrintSquareMat(nOrbA,Work(iAddSq))
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------

      Call Local_Triang(nOrbA, Work(iAddSq))
      call daxpy_(LenCou,1.0d0,Work(iAddSq),1,Work(iAddCou),1)
      Call GetMem('CSq','Free','Real',iAddSq,LenEx)
CGG   ------------------------------------------------------------------
      If(IfTest) then
      Write(6,*)
      Call XFlush(6)
      Write(6,*)'        The Triangular Integrals matrix'
      Call PrintDiagMat(nOrbA,Work(iAddCou))
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------

* --- END GATERING of SubBlocks

      Do iSB_A = 1, 3
        Do iSB_B = 1, 3
          If (iAddSB(iSB_A,iSB_B).GT.0) Call GetMem('SB','Free','Real',
     &                         iAddSB(iSB_A,iSB_B), LenSB(iSB_A,iSB_B))
        EndDo
      EndDo

CGG   ------------------------------------------------------------------
c      IfTest=.False.
CGG   ------------------------------------------------------------------
      Return
      End
