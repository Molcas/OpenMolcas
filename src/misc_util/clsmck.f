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
* Copyright (C) 1993, Markus P. Fuelscher                              *
*               1993, Per-Olof Widmark                                 *
************************************************************************
      Subroutine ClsMCK(rc,option)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Close the one-electron integral file.                            *
*                                                                      *
*     input:                                                           *
*     option : Switch to set options                                   *
*              (not used at present)                                   *
*                                                                      *
*     output:                                                          *
*     rc     : Return code.                                            *
*              A value of 0 (zero) is returned upon successful         *
*              completion of the request. A nonzero value indi-        *
*              cates an error.                                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. P. Fuelscher and P.O. Widmark                                 *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Integer (A-Z)
*
#include "MckRc.fh"
#include "MckFlags.fh"
#include "SysDef.fh"
#include "MckDat.fh"
*----------------------------------------------------------------------*
*     Check file status                                                *
*----------------------------------------------------------------------*
*     Call qEnter('ClsMCK')
      If( AuxMCK(pOpen).ne.1 ) Then
         rc=rcCL01
      Call SysAbendMsg('ClsMCK',
     *  'The MCK file has not been opened',' ')
      End If
      If ( iAnd(Option,1024).ne.0 ) Then
         Write (6,'(i6,z8)') pFID,TocOne(pFID)
         Write (6,'(i6,z8)') pVersN,TocOne(pVersN)
         Write (6,'(i6,z8)') pTitle,TocOne(pTitle)
         Write (6,'(i6,z8)') pOp,TocOne(pOp)
         Write (6,'(i6,z8)') pSym,TocOne(pSym)
         Write (6,'(i6,z8)') pSymOp,TocOne(pSymOp)
         Write (6,'(i6,z8)') pBas,TocOne(pBas)
         Write (6,'(i6,z8)') pNext,TocOne(pNext)
         Write (6,'(i6,z8)') pEnd,TocOne(pEnd)
      End If
*----------------------------------------------------------------------*
*     Reset error code,open flag and unit number. Close file.          *
*----------------------------------------------------------------------*
      LuMCK=AuxMCK(pLu)
      Call DaClos(LuMCK)
      AuxMCK(pLu   ) = 0
      AuxMCK(pOpen ) = 0
      rc=rc0000
*----------------------------------------------------------------------*
*     Terminate procedure                                              *
*----------------------------------------------------------------------*
*     Call qExit('ClsMCK')
      Return
      End
