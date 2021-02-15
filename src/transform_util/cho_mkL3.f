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
      Subroutine MkL3(iSymA,iSymI,iI, numV, LyType,iJy, AddLx0, SameLx)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden & Torino University, Italy         *
*           February 2005                                              *
*----------------------------------------------------------------------*
* Purpuse:  Generation of the Cholesky matrix of Secondary(iSymA) for  *
*           occupied iI(iSymI) for numV vectors.                       *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Integer iSymA,iSymI,iI, numV, LyType,iJy
      Real*8 AddLx0(*)
      Logical SameLx
#include "rasdim.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "cho_tra.fh"

*     Build Lx
      If (iI.LE.nIsh(iSymI)) then
        LxType = 3
        iIx = iI
        nIx = nIsh(iSymI)
      else
        LxType = 5
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

      iAddLx  = 1
      iAddTCVX= iMemTCVX(LxType,iSymA,iSymI,1)+nSsh(iSymA)*(iIx-1)
      Do iV=1,numV
        Call dCopy_(nSsh(iSymA),Work(iAddTCVX),1,AddLx0(iAddLx),1)
        iAddTCVX= iAddTCVX +  nSsh(iSymA) * nIx
        iAddLx  = iAddLx + nSsh(iSymA)
      EndDo

      Return
      End
