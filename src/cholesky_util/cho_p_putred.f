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
      SubRoutine Cho_P_PutRed(iPass,iLoc)
C
C     Purpose: write global and local reduced set index arrays
C              to disk and set
C              address for next write. iLoc specifies the location in
C              the index arrays to write, and iPass identifies the
C              reduced set (i.e. the integral pass).
C
      Implicit None
      Integer iPass, iLoc
#include "cho_para_info.fh"
#include "cholesky.fh"
#include "choglob.fh"

      Integer iTmp

      Real*8 c1, c2, w1, w2

      Call Cho_Timer(c1,w1)

      If (Cho_Real_Par) Then

C        Swap local and global reduced set indices and use original serial
C        routine to write global index arrays.
C        -----------------------------------------------------------------

         Call Cho_P_IndxSwp()
         iTmp = LuRed
         LuRed = LuRed_G
         Call Cho_PutRed(iPass,iLoc)
         LuRed = iTmp
         Call Cho_P_IndxSwp()

      End If

C     Write local index arrays.
C     -------------------------

      Call Cho_PutRed(iPass,iLoc)

      Call Cho_Timer(c2,w2)
      tMisc(1,2)=tMisc(1,2)+c2-c1
      tMisc(2,2)=tMisc(2,2)+w2-w1

      End
