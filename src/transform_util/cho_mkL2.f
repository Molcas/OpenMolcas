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
      Subroutine MkL2(iSymA,iSymI,iI, numV, LyType,iJy, iAddLx0, SameLx)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden & Torino University, Italy         *
*           February 2005                                              *
*----------------------------------------------------------------------*
* Purpuse:  Generation of the Cholesky matrix of Active(iSymA) for     *
*           occupied iI(iSymI) for numV vectors.                       *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "rasdim.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "cho_tra.fh"
      Logical SameLx

*     Build Lx
      If (iI.LE.nIsh(iSymI)) then
        LxType = 2
        iIx = iI
        nIx = nIsh(iSymI)
      else
        LxType = 4
        iIx = iI - nIsh(iSymI)
        nIx = nAsh(iSymI)
      EndIf

      If (.NOT.SameLx) then
        LyType=LxType
        iJy   =iIx
      else
        If (LyType.EQ.LxType .and. iIx.EQ.iJy) then
          Return
        else
          SameLx=.False.
        EndIf
      EndIf

      iAddLx  = iAddLx0
      iAddTCVX= iMemTCVX(LxType,iSymA,iSymI,1)+nAsh(iSymA)*(iIx-1)
      Do iV=1,numV
        Call dCopy_(nAsh(iSymA),Work(iAddTCVX),1,Work(iAddLx),1)
        iAddTCVX= iAddTCVX +  nAsh(iSymA) * nIx
        iAddLx  = iAddLx + nAsh(iSymA)
      EndDo

      Return
      End
