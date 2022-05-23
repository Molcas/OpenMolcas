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
      SubRoutine Cho_RI_PutInfo(iPass,iRed)
#if defined (_MOLCAS_MPP_)
      Use Para_Info, Only: nProcs, Is_Real_Par
#endif
      Implicit None
      Integer iPass, iRed
#include "cholesky.fh"
#include "choglob.fh"
      Logical doSwap
      Integer iTmp
#if defined (_MOLCAS_MPP_)
      doSwap = nProcs.gt.1 .and. Is_Real_Par()
#else
      doSwap = .False.
#endif

      If (doSwap) Then
         iTmp = LuRed
         LuRed = LuRed_G
         Call Cho_PutRed(iPass,iRed) ! save reduced set indices on disk
         LuRed = iTmp
         iTmp = LuRst
         LuRst = LuRst_G
         Call Cho_WrRstC(iPass) ! save disk addresses etc. on disk
         LuRst = iTmp
      Else
         Call Cho_PutRed(iPass,iRed) ! save reduced set indices on disk
         Call Cho_WrRstC(iPass) ! save disk addresses etc. on disk
      End If
*
      Return
      End
