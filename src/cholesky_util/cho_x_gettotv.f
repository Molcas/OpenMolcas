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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine Cho_X_GetTotV(nV,n)
C
C     Thomas Bondo Pedersen, April 2010.
C
C     Purpose: return total number of vectors (over all nodes).
C
      Implicit None
      Integer n
      Integer nV(n)
#include "cholesky.fh"
#include "chopar.fh"
#include "cho_para_info.fh"

      Integer iSym

#if defined (_DEBUGPRINT_)
      If (n .lt. nSym) Then
         Call Cho_Quit('Illegal input variable in Cho_X_GetTotV',104)
      End If
#endif

      If (Cho_Real_Par) Then
         Do iSym = 1,nSym
            nV(iSym)=NumCho_Bak(iSym)
         End Do
      Else
         Do iSym = 1,nSym
            nV(iSym)=NumCho(iSym)
         End Do
      End If

      End
