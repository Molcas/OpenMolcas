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
* Copyright (C) 2004, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine ChoMP2_Tra(COcc,CVir,Diag,DoDiag)
C
C     Thomas Bondo Pedersen, Dec. 2004.
C
C     Purpose: transform Cholesky vectors to (ai) MO basis.
C
      Implicit None
      Real*8  COcc(*), CVir(*), Diag(*)
      Logical DoDiag
#include "cholesky.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"

      Character*10 SecNam
      Parameter (SecNam = 'ChoMP2_Tra')

      Integer kOffD, iSym, ipW, lW

C     Allocate remaining memory.
C     --------------------------

      Call GetMem('TraGetMax','Max ','Real',ipW,lW)
      Call GetMem('TraMax','Allo','Real',ipW,lW)

      kOffD = 1
      Do iSym = 1,nSym

C        Open files for MO vectors.
C        --------------------------

         Call ChoMP2_OpenF(1,1,iSym)

C        Transform vectors.
C        ------------------

         Call ChoMP2_Tra_1(COcc,CVir,Diag(kOffD),DoDiag,Work(ipW),lW,
     &                     iSym)
         If (DoDiag) kOffD = kOffD + nT1am(iSym)

C        Close files for MO vectors.
C        ---------------------------

         Call ChoMP2_OpenF(2,1,iSym)

      End Do

C     Free memory.
C     ------------

      Call GetMem('TraMax','Free','Real',ipW,lW)

      End
