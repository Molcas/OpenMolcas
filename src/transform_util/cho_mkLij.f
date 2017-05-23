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
      Subroutine MkLij(iSymI,iSymJ,iI,iJ,numV, iAddLij)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Torino University, Italy                                   *
*           July 2005                                                  *
*----------------------------------------------------------------------*
* Purpuse:  Generation of the Cholesky vector for occupied iI(iSymI),  *
*           iJ(iSymJ) for numV vectors.                                *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "rasdim.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "cho_tra.fh"

      If (iI.LE.nIsh(iSymI)) then
        iIx = iI
        nIx = nIsh(iSymI)
        If (iJ.LE.nIsh(iSymJ)) then
          LijType = 1
          iJy = iJ
          nJy = nIsh(iSymJ)
        else
          LijType = 7
          iJy = iJ - nIsh(iSymJ)
          nJy = nAsh(iSymJ)
        EndIf
      else
        iIx = iI - nIsh(iSymI)
        nIx = nAsh(iSymI)
        If (iJ.LE.nIsh(iSymJ)) then
          LijType = 2
          iJy = iJ
          nJy = nIsh(iSymJ)
        else
          LijType = 4
          iJy = iJ - nIsh(iSymJ)
          nJy = nAsh(iSymJ)
        EndIf
      EndIf

CGG   ------------------------------------------------------------------
      If(IfTest) then
      Write(6,*)'     Cho_MkLij: TCVx(',LijType,': ',iSymI,',',iSymJ,')'
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------
      iAddTCVX= iMemTCVX(LijType,iSymI,iSymJ,1) + nIx*(iJy-1) + iIx - 1
      Call dCopy_(numV,Work(iAddTCVX),nIx*nJy,Work(iAddLij),1)

      Return
      End
