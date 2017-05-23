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
      Subroutine MkExSB13(iAddSB,LenSB,
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden & Torino University, Italy         *
*           February 2005                                              *
*----------------------------------------------------------------------*
* Purpuse:  Generation of the SubBlock(1,3) (p,q inactive,secondary)   *
*           two-electron integral matrix for each i,j occupied couple. *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "rasdim.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "cho_tra.fh"
      Logical SameLx

*   - SubBlock 1 3
      LenSB = nIsh(iSymA) * nSsh(iSymB)
      Call GetMem('SB','Allo','Real',iAddSB,LenSB)

*     Build Lx
      Call GetMem('Lx','Allo','Real',iAddLx0,nIsh(iSymA)*numV)
      LxType=0
      iIx=0
      SameLx=.False.
      Call MkL1(iSymA,iSymI,iI,numV, LxType,iIx, iAddLx0,SameLx)

*     Build Ly
      Call GetMem('Ly','Allo','Real',iAddLy0,nSsh(iSymB)*numV)
      Call MkL3(iSymB,iSymJ,iJ,numV, LxType,iIx, iAddLy0,SameLx)

*     Generate the SubBlock
      Call DGEMM_('N','T',nSsh(iSymB),nIsh(iSymA),numV,1.0d0,
     &    Work(iAddLy0),nSsh(iSymB), Work(iAddLx0),nIsh(iSymA),
     &                     0.0d0,Work(iAddSB),nSsh(iSymB) )

      Call GetMem('Ly','Free','Real',iAddLy0,nSsh(iSymB)*numV)
      Call GetMem('Lx','Free','Real',iAddLx0,nIsh(iSymA)*numV)

      Return
      End
