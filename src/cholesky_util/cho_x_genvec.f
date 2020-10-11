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
      SubRoutine Cho_X_GenVec(irc,Diag)

      Implicit None
      Integer irc
      Real*8  Diag(*)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Character*12 SecNam
      Parameter (SecNam = 'Cho_X_GenVec')

      Integer MaxQual_SAVE, ip_iQuAB_SAVE, l_iQuAB_SAVE
      Integer iSym


C     Set return code.
C     ----------------

      irc = 0

C     Re-allocate the iQuAB index array, save old allocation.
C     This is used to trick the integral extraction from Seward so as to
C     reduce the number of re-calculations of shell pairs.
C     ------------------------------------------------------------------

      ip_iQuAB_SAVE = ip_iQuAB
      l_iQuAB_SAVE  = l_iQuAB
      MaxQual_SAVE  = MaxQual

      MaxQual = NumCho(1)
      Do iSym = 2,nSym
         MaxQual = Max(MaxQual,NumCho(iSym))
      End Do

      l_iQuAB = MaxQual*nSym
      Call Cho_Mem('iQuAB_2','Allo','Inte',ip_iQuAB,l_iQuAB)

C     Read initial diagonal.
C     ----------------------

      Call Cho_IODiag(Diag,2)

C     Reinitialize the number of zeroed negative diagonals.
C     Turn on damped screening for second step.
C     -----------------------------------------------------

      nNZTot = 0
      MODE_SCREEN = 1

C     Generate vectors.
C     -----------------

      Call Cho_GnVc_Drv(irc,Diag)
      If (irc .ne. 0) Then
         Write(Lupri,*) SecNam,': Cho_GnVc_Drv returned ',irc
         Go To 100 ! exit
      End If

C     De-allocations.
C     Restore original iQuAB array.
C     -----------------------------

  100 Continue
      Call Cho_Mem('iQuAB_2','Free','Inte',ip_iQuAB,l_iQuAB)
      ip_iQuAB = ip_iQuAB_SAVE
      l_iQuAB  = l_iQuAB_SAVE
      MaxQual  = MaxQual_SAVE


      End
