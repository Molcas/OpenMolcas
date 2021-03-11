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
* Copyright (C) 2004,2005, Giovanni Ghigo                              *
************************************************************************
      Subroutine ChoMP2_TwoEl(iBatch,nBatch,numV, LUINTM,iAddrIAD2M,
     &                                   iSymI,iSymJ,iSymA,iSymB, iSymL)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden                                    *
* Written:  October-November 2004                                      *
* Modified for Cholesky-MP2 May 2005                                   *
************************************************************************
      use Cho_Tra
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "rasdim.fh"
#include "stdalloc.fh"
#include "SysDef.fh"

      Real*8, Allocatable:: AddEx1(:), AddEx2(:), AddEx2t(:)

      nSymP=(nSym**2+nSym)/2
      Call LenInt(iSymI,iSymJ,iSymA,iSymB,nN_IJ,nN_AB,nN_Ex1,nN_Ex2)

* *** START GENERATION of EXCHANGE-1 INTEGRALS  **********************
      IF (nN_IJ*nN_Ex1.GT.0) THEN
        iIJAB = ( (iSymI**2-iSymI)/2 + iSymJ-1 ) * nSymP +
     &                 (iSymA**2-iSymA)/2 + iSymB
        SubBlocks(3,3)=.True.
        If (iBatch.EQ.1) then
          IAD2M(2,iIJAB)=iAddrIAD2M
        else
          iAddrIAD2M=IAD2M(2,iIJAB)
        EndIf
*  ---  Start Loop on i, j
        iAddrIAD2Mij=iAddrIAD2M
        Do iI=1,nOsh(iSymI)
          If(iSymI.EQ.iSymJ) then
            iEndJ=iI
          else
            iEndJ=nOsh(iSymJ)
          EndIf
          Do iJ=1,iEndJ
            Call mma_allocate(AddEx1,nN_Ex1,Label='AddEx1')
            If (iBatch.GT.1) then
              ! Reload Int
              Call dDaFile(LUINTM,2,AddEx1,nN_Ex1,iAddrIAD2Mij)
              iAddrIAD2Mij=iAddrIAD2Mij-nN_Ex1
            else
              AddEx1(:)=0.0D0
            EndIf
            Call ChoMP2_GenE(iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV,
     &                                            AddEx1,nN_Ex1 )
            Call dDaFile(LUINTM,1,AddEx1,nN_Ex1,iAddrIAD2Mij)
            Call mma_deallocate(AddEx1)
          EndDo
        EndDo
*  ---  End Loop on i, j
        iAddrIAD2M=iAddrIAD2Mij
      ENDIF
* *** END GENERATION of EXCHANGE-1 INTEGRALS  ************************


* *** START GENERATION of EXCHANGE-2 INTEGRALS  **********************
      IF (nN_IJ*nN_Ex2.GT.0) THEN
        iIJAB = ( (iSymI**2-iSymI)/2 + iSymJ-1 ) * nSymP +
     &                 (iSymB**2-iSymB)/2 + iSymA
        SubBlocks(3,3)=.True.
        If (iBatch.EQ.1) then
          IAD2M(3,iIJAB)=iAddrIAD2M
        else
          iAddrIAD2M=IAD2M(3,iIJAB)
        EndIf
*  ---  Start Loop on i, j
        iAddrIAD2Mij=iAddrIAD2M
        Do iI=1,nOsh(iSymI)
          If(iSymI.EQ.iSymJ) then
            iEndJ=iI
          else
            iEndJ=nOsh(iSymJ)
          EndIf
          Do iJ=1,iEndJ
            nA=nSsh(iSymA)
            nB=nSsh(iSymB)
            Call mma_allocate(AddEx2,nN_Ex2,Label='AddEx2')
            Call mma_allocate(AddEx2t,nN_Ex2,Label='AddEx2t')
            If (iBatch.GT.1) then
              ! Reload Int
              Call dDaFile(LUINTM,2,AddEx2,nN_Ex2,iAddrIAD2Mij)
              iAddrIAD2Mij=iAddrIAD2Mij-nN_Ex2
              Call Trnsps(nA,nB,AddEx2,AddEx2t)
            else
              AddEx2t(:)=0.0D0
            EndIf
            Call ChoMP2_GenE(iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV,
     &                                           AddEx2t,nN_Ex2 )
            Call Trnsps(nB,nA,AddEx2t,AddEx2)
            Call dDaFile(LUINTM,1,AddEx2,nN_Ex2,iAddrIAD2Mij)
            Call mma_deallocate(AddEx2t)
            Call mma_deallocate(AddEx2)
          EndDo
        EndDo
*  ---  End Loop on i, j
        iAddrIAD2M=iAddrIAD2Mij
      ENDIF
* *** END GENERATION of EXCHANGE-2 INTEGRALS  ************************

      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(nBatch)
         Call Unused_integer(iSymL)
      End If
      End
