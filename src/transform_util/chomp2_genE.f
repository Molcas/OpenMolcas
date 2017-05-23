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
      Subroutine ChoMP2_GenE(iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV,
     &                      iAddEx,LenE)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden & Torino University, Italy         *
*           Jannuary-February 2005                                     *
* Modified for Cholesky-MP2 May 2005                                   *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "rasdim.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "cho_tra.fh"
      Dimension iAddSB(3,3), LenSB(3,3)

      iAddSB(3,3)= 0  ! Mem Address of the SubBlocks
      LenSB (3,3)= 0  ! Length of the SubBlocks

      Call MkExMP2(iAddSB(3,3),LenSB(3,3),
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV, iAddEx,LenE)

      Call GetMem('SB','Free','Real',iAddSB(3,3), LenSB(3,3))

      Return
      End

      Subroutine MkExMP2(iAddSB,LenSB,
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV, iAddEx,LenE)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden & Torino University, Italy         *
*           Jannuary February 2005                                     *
* Modified for Cholesky-MP2 May 2005                                   *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "rasdim.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "cho_tra.fh"
      Logical SameLx

*   - SubBlock 3 3
      LenSB = nSsh(iSymA) * nSsh(iSymB)
      Call GetMem('SB','Allo','Real',iAddSB,LenSB)

*     Build Lx
      Call GetMem('Lx','Allo','Real',iAddLx0,nSsh(iSymA)*numV)
      LxType=0
      iIx=0
      SameLx=.False.
      Call MkL3(iSymA,iSymI,iI,numV, LxType,iIx, iAddLx0,SameLx)

*     Build Ly
      Call GetMem('Ly','Allo','Real',iAddLy0,nSsh(iSymB)*numV)
      If(iSymA.EQ.iSymB) SameLx=.True.
      Call MkL3(iSymB,iSymJ,iJ,numV, LxType,iIx, iAddLy0,SameLx)

*     Generate the SubBlock
      If (.NOT.SameLx) then
        Call DGEMM_('N','T',nSsh(iSymB),nSsh(iSymA),numV,1.0d0,
     &   Work(iAddLy0),nSsh(iSymB), Work(iAddLx0),nSsh(iSymA),
     &                        1.0d0,Work(iAddEx),nSsh(iSymB) )
      else
        Call DGEMM_('N','T',nSsh(iSymA),nSsh(iSymA),numV,1.0d0,
     &   Work(iAddLx0),nSsh(iSymA), Work(iAddLx0),nSsh(iSymA),
     &                        1.0d0,Work(iAddEx),nSsh(iSymA) )
      EndIf

      Call GetMem('Ly','Free','Real',iAddLy0,nSsh(iSymB)*numV)
      Call GetMem('Lx','Free','Real',iAddLx0,nSsh(iSymA)*numV)

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(LenE)
      End
