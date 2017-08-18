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
* Copyright (C) 2004, Giovanni Ghigo                                   *
************************************************************************
*  Cho_TwoEl
*
*> @brief
*>   Driver for the generation of the two-electrons integrals file (``MOLINT``) from the Transformed Cholesky vectors TCVx
*> @author Giovanni Ghigo
*>
*> @details
*> Within the symmetry \p iSymL, the routine generates the symmetry
*> block \p iSymI, \p iSymJ, \p iSymA, \p iSymB of two-electrons integrals. The
*> number of integrals to generate for each  \p iSymI, \p iSymJ couple is
*> defined by \c LenInt. The exch-1 and exch-2 integrals are then
*> generated for the occupied MO calling ::Cho_GenE.
*>
*> @param[in]     iBatch                  Main batch current number
*> @param[in]     nBatch                  Main batch total number
*> @param[in]     NumV                    Number of Cholesky vectors to transform in the current batch
*> @param[in]     LUINTM                  Unit number of two-electrons integrals file (``MOLINT``)
*> @param[in,out] iAddrIAD2M              Current disk address in ``MOLINT``
*> @param[in]     iSymI,iSymJ,iSymA,iSymB Symmetry block of the two-electrons integrals
*> @param[in]     iSymL                   Symmetry of the Cholesky vector
************************************************************************
      Subroutine Cho_TwoEl(iBatch,nBatch,numV, LUINTM,iAddrIAD2M,
     &                                   iSymI,iSymJ,iSymA,iSymB, iSymL)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden                                    *
*           October-November 2004                                      *
*----------------------------------------------------------------------*
* Routine for the generation of the new Two-electrons integrals files  *
* (MOLINT) from the Transformed Cholesky Full Vectors (TCVX).          *
*                                                                      *
* A,B are MO indices, counting only non-frozen and non-deleted.        *
* I,J are occupied MO indices, only non-frozen and non-deleted.        *
* (AB/IJ) ARE ALWAYS GENERATED                                         *
* (AI/BJ) IF ISP.GE.ISR                                                *
* (AI/JB) IF ISP.GT.ISS AND ISP.NE.ISQ                                 *
* (IA/BJ) IF ISQ.GT.ISR AND ISP.NE.ISQ                                 *
* (IA/JB) IF ISQ.GE.ISS AND ISP.NE.ISQ                                 *
*                                                                      *
*   IAD2M CONTAINS START ADRESS FOR EACH TYPE OF INTEGRALS:            *
*    IAD2M(1,ISPQRS)   COULOMB INTEGRALS <AB|IJ>                       *
*    IAD2M(2,ISPQRS)   EXCHANGE INTEGRALS <AB|IJ> FOR iSymI > iSymJ    *
*    IAD2M(3,ISPQRS)   EXCHANGE INTEGRALS <AB|IJ> FOR iSymI < iSymJ    *
*   THE LAST ADRESS IS ZERO IF iSymI = iSymJ                           *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "rasdim.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "cho_tra.fh"

      nSymP=(nSym**2+nSym)/2
      Call LenInt(iSymI,iSymJ,iSymA,iSymB,nN_IJ,nN_AB,nN_Ex1,nN_Ex2)
CGG   ------------------------------------------------------------------
      If(IfTest) then
      Write(6,*)
      Write(6,'(A,4I3,A,I8,A,I9,A,I9)')
     &'    * [CGG:Cho_TwoEl]: SYMMETRY BLOCK < A B | I J >',iSymA,iSymB,
     &iSymI,iSymJ,': nN_AB=',nN_AB,', nN_Ex1=',nN_Ex1,', nN_Ex2=',nN_Ex2
      If(nN_IJ*(nN_AB+nN_Ex1+nN_Ex2).EQ.0)
     &Write(6,*)'                      Nothing to do!'
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------
      If (nN_IJ*(nN_AB+nN_Ex1+nN_Ex2).GT.0) then                  ! If !
        iIJAB = ( (iSymI**2-iSymI)/2 + iSymJ-1 ) * nSymP +
     &                 (iSymA**2-iSymA)/2 + iSymB

*  ***  START GENERATION of COULOMB INTEGRALS  *************************
        IF (DoCoul .and. nN_AB.GT.0) THEN
C          Call Def_SubBlockC(iSymA,iSymB)
          Call Def_SubBlockE(iSymA,iSymB)
CGG   ------------------------------------------------------------------
      If(IfTest) then
      Write(6,*)
      Write(6,*) '    Generation of Coulomb Integrals'
      Call XFlush(6)
      Write(6,*)'       Available TCVx for Cou: '
      Do iType=1,6
      If(TCVXist(iType,iSymA,iSymI))
     &Write(6,*)'       -TCV x=',iType,' Sym=',iSymA,iSymI
      If(TCVXist(iType,iSymB,iSymJ) .and. iSymA.NE.iSymB)
     &Write(6,*)'       -TCV x=',iType,' Sym=',iSymB,iSymJ
      EndDo
      Write(6,*)
      Write(6,*)'       SubBlocks to create for Cou: '
      Do i = 1, 3
      Do j = 1, 3
      If(SubBlocks(i,j)) Write(6,*) '       -SB(',i,',',j,')'
      EndDo
      EndDo
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------
          If (iBatch.EQ.1) then
            IAD2M(1,iIJAB)=iAddrIAD2M
          else
            iAddrIAD2M=IAD2M(1,iIJAB)
          EndIf
*   ---   Start Loop on i, j
          iAddrIAD2Mij=iAddrIAD2M
          Do iI=1,nOsh(iSymI)
            If(iSymI.EQ.iSymJ) then
              iEndJ=iI
            else
              iEndJ=nOsh(iSymJ)
            EndIf
            Do iJ=1,iEndJ
CGG   ------------------------------------------------------------------
      If(IfTest) then
      Write(6,*)
      Write(6,*) '   Coulomb Integrals for |ij> pair',
     &  iI,iJ,'  iAddrIAD2Mij=',iAddrIAD2Mij
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------
              Call GetMem('Coul','Allo','Real',iAddCou,nN_AB)
              If (iBatch.GT.1) then
                Call dDaFile(LUINTM,2,Work(iAddCou),nN_AB,
     &                                     iAddrIAD2Mij)   ! Reload Int
                iAddrIAD2Mij=iAddrIAD2Mij-nN_AB
              else
                Call dCopy_(nN_AB,0.0d0,0,Work(iAddCou),1)
              EndIf
              Call Cho_GenC(iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV,
     &                                          iAddCou,nN_AB, nN_Ex1 )
              Call GAdSum(Work(iAddCou),nN_AB)
              Call dDaFile(LUINTM,1,Work(iAddCou),nN_AB,iAddrIAD2Mij)
              Call GetMem('Coul','Free','Real',iAddCou,nN_AB)
            EndDo
          EndDo
*   ---   End Loop on i, j
          iAddrIAD2M=iAddrIAD2Mij ! Last written+1 in iAddrIAD2M
        ENDIF
*  ***  END GENERATION of COULOMB INTEGRALS  ***************************


*  ***  START GENERATION of EXCHANGE-1 INTEGRALS  **********************
        IF (nN_Ex1.GT.0) THEN
          Call Def_SubBlockE(iSymA,iSymB)
CGG   ------------------------------------------------------------------
c      If(IfTest) then
c      Write(6,*)
c      Write(6,*) '    Generation of Exchange-1 Integrals'
c      Write(6,*)'       Available TCVx for Ex-1: '
c      Do iType=1,MxTCVx
c      If(TCVXist(iType,iSymA,iSymI))
c     &Write(6,*)'       -TCV x=',iType,' Sym=',iSymA,iSymI
c      If(TCVXist(iType,iSymB,iSymJ) .and. iSymA.NE.iSymB)
c     &Write(6,*)'       -TCV x=',iType,' Sym=',iSymB,iSymJ
c      EndDo
c      Write(6,*)
c      Write(6,*)'       SubBlocks to create for Ex-1: '
c      Do i = 1, 3
c      Do j = 1, 3
c      If(SubBlocks(i,j)) Write(6,*) '       -SB(',i,',',j,')'
c      EndDo
c      EndDo
c      Call XFlush(6)
c      EndIf
CGG   ------------------------------------------------------------------
          If (iBatch.EQ.1) then
            IAD2M(2,iIJAB)=iAddrIAD2M
          else
            iAddrIAD2M=IAD2M(2,iIJAB)
          EndIf
*   ---   Start Loop on i, j
          iAddrIAD2Mij=iAddrIAD2M
          Do iI=1,nOsh(iSymI)
            If(iSymI.EQ.iSymJ) then
              iEndJ=iI
            else
              iEndJ=nOsh(iSymJ)
            EndIf
            Do iJ=1,iEndJ
CGG   ------------------------------------------------------------------
c      If(IfTest) then
c      Write(6,*)
c      Write(6,*) '   Excha-1 Integrals for |ij> pair',iI,iJ,
c     &  '  iAddrIAD2Mij=',iAddrIAD2Mij,'   #'
c      Call XFlush(6)
c      EndIf
CGG   ------------------------------------------------------------------
              Call GetMem('Ex1','Allo','Real',iAddEx1,nN_Ex1)
              If (iBatch.GT.1) then
                Call dDaFile(LUINTM,2,Work(iAddEx1),nN_Ex1,
     &                                     iAddrIAD2Mij)   ! Reload Int
                iAddrIAD2Mij=iAddrIAD2Mij-nN_Ex1
              else
                call dcopy_(nN_Ex1,0.0d0,0,Work(iAddEx1),1)
              EndIf
              Call Cho_GenE(iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV,
     &                                             iAddEx1,nN_Ex1 )
              Call GAdSum(Work(iAddEx1),nN_Ex1)
              Call dDaFile(LUINTM,1,Work(iAddEx1),nN_Ex1,iAddrIAD2Mij)
              Call GetMem('Ex1','Free','Real',iAddEx1,nN_Ex1)
            EndDo
          EndDo
*   ---   End Loop on i, j
          iAddrIAD2M=iAddrIAD2Mij
        ENDIF
*  ***  END GENERATION of EXCHANGE-1 INTEGRALS  ************************


*  ***  START GENERATION of EXCHANGE-2 INTEGRALS  **********************
        IF (nN_Ex2.GT.0 .and. DoExc2) THEN
          iIJAB = ( (iSymI**2-iSymI)/2 + iSymJ-1 ) * nSymP +
     &                 (iSymB**2-iSymB)/2 + iSymA
          Call Def_SubBlockE(iSymA,iSymB)
CGG   ------------------------------------------------------------------
c      If(IfTest) then
c      Write(6,*)
c      Write(6,*) '    Generation of Exchange-2 Integrals'
c      Write(6,*)'       Available TCVx for Ex-2: '
c      Do iType=1,MxTCVx
c      If(TCVXist(iType,iSymA,iSymI))
c     &Write(6,*)'       -TCV x=',iType,' Sym=',iSymA,iSymI
c      If(TCVXist(iType,iSymB,iSymJ) .and. iSymA.NE.iSymB)
c     &Write(6,*)'       -TCV x=',iType,' Sym=',iSymB,iSymJ
c      EndDo
c      Write(6,*)
c      Write(6,*)'       SubBlocks to create for Ex-2: '
c      Do i = 1, 3
c      Do j = 1, 3
c      If(SubBlocks(i,j)) Write(6,*) '       -SB(',i,',',j,')'
c      EndDo
c      EndDo
c      Call XFlush(6)
c      EndIf
CGG   ------------------------------------------------------------------
          If (iBatch.EQ.1) then
            IAD2M(3,iIJAB)=iAddrIAD2M
          else
            iAddrIAD2M=IAD2M(3,iIJAB)
          EndIf
*   ---   Start Loop on i, j
          iAddrIAD2Mij=iAddrIAD2M
          Do iI=1,nOsh(iSymI)
            If(iSymI.EQ.iSymJ) then
              iEndJ=iI
            else
              iEndJ=nOsh(iSymJ)
            EndIf
            Do iJ=1,iEndJ
CGG   ------------------------------------------------------------------
c      If(IfTest) then
c      Write(6,*)
c      Write(6,*) '   Excha-2 Integrals for |ij> pair',iI,iJ,
c     &  '  iAddrIAD2Mij=',iAddrIAD2Mij,'   #'
c      Call XFlush(6)
c      EndIf
CGG   ------------------------------------------------------------------
              If (DoTCVA) then
                nA=nOrb(iSymA)
                nB=nOrb(iSymB)
              else
                nA=nSsh(iSymA)
                nB=nSsh(iSymB)
              EndIf
              Call GetMem('Ex2','Allo','Real',iAddEx2,nN_Ex2)
              Call GetMem('Ex2t','Allo','Real',iAddEx2t,nN_Ex2)
              If (iBatch.GT.1) then
                Call dDaFile(LUINTM,2,Work(iAddEx2),nN_Ex2,
     &                                     iAddrIAD2Mij)   ! Reload Int
                iAddrIAD2Mij=iAddrIAD2Mij-nN_Ex2
                Call Trnsps(nA,nB,Work(iAddEx2),Work(iAddEx2t))
              else
                call dcopy_(nN_Ex2,0.0d0,0,Work(iAddEx2t),1)
              EndIf
              Call Cho_GenE(iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV,

     &                                            iAddEx2t,nN_Ex2 )
              Call Trnsps(nB,nA,Work(iAddEx2t),Work(iAddEx2))
CGG   ------------------------------------------------------------------
c      If(IfTest) then
c      Write(6,*) '    Integrals from Production code:'
c      Write(6,'(8F10.6)') (Work(iAddEx2+i),i=0,nN_Ex2-1)
c      Call XFlush(6)
c      EndIf
CGG   ------------------------------------------------------------------
              Call GAdSum(Work(iAddEx2),nN_Ex2)
              Call dDaFile(LUINTM,1,Work(iAddEx2),nN_Ex2,iAddrIAD2Mij)
              Call GetMem('Ex2t','Free','Real',iAddEx2t,nN_Ex2)
              Call GetMem('Ex2','Free','Real',iAddEx2,nN_Ex2)
            EndDo
          EndDo
*   ---   End Loop on i, j
          iAddrIAD2M=iAddrIAD2Mij
        ENDIF
*  ***  END GENERATION of EXCHANGE-2 INTEGRALS  ************************
      EndIf                                                       ! If !

      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(nBatch)
         Call Unused_integer(iSymL)
      End If
      End
