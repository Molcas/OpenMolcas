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
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
*  Cho_X_nVecRS
*
*> @brief
*>   Find first vector and number of vectors in reduced set \p iRed, sym. block \p iSym
*> @author Thomas Bondo Pedersen
*>
*> @details
*> This routine finds the first vector and number of
*> vectors in reduced set \p iRed, sym. block \p iSym.
*> Note that \p iVec=\p nVec = ``0`` may be returned---this is
*> perfectly acceptable: a given reduced set may be
*> empty. However, if negative numbers are returned
*> (\p iVec < ``0`` and \p nVec < ``0``), an error has ocurred. This
*> should be tested by the caller!!
*>
*> @note
*> The Cholesky procedure must have been successfully initialized (by ::Cho_X_Init).
*>
*> @param[in]  iRed Reduced set
*> @param[in]  iSym Symmetry block (1-8)
*> @param[out] iVec First vector in red. set \p iRed, sym. \p iSym
*> @param[out] nVec Number of vectors in red. set \p iRed, sym. \p iSym
************************************************************************
      SubRoutine Cho_X_nVecRS(iRed,iSym,iVec,nVec)
      use ChoSwp, only: InfVec
      Implicit None
      Integer iRed, iSym, iVec, nVec
#include "cholesky.fh"
#include "choptr.fh"

      Character*12 SecNam
      Parameter (SecNam = 'Cho_X_nVecRS')

      Logical Found

      Integer irc, LastRed, jVec, jRed

C     Check input.
C     ------------

      irc = 0
      If (iSym.lt.1 .or. iSym.gt.nSym) Then
         irc = -1
      End If
      If (NumCho(iSym).lt.0 .or. NumCho(iSym).gt.MaxVec) Then
         irc = -2
      End If
      LastRed = InfVec(NumCho(iSym),2,iSym)
      If (LastRed .lt. 1) Then
         irc = -3
      End If
      If (iRed .lt. 1) Then
         irc = -4
      End If
      If (irc .ne. 0) Then
         iVec = irc
         nVec = irc
         Return
      End If
      If (iRed .gt. LastRed) Then
         iVec = 0
         nVec = 0
         Return
      End If

C     Find first vector in reduced set iRed.
C     --------------------------------------

      Found = .false.
      jVec  = 0
      Do While (jVec.lt.NumCho(iSym) .and. .not.Found)
         jVec = jVec + 1
         jRed = InfVec(jVec,2,iSym)
         If (jRed .eq. iRed) Then
            iVec  = jVec
            Found = .true.
         Else If (jRed .gt. iRed) Then
            jVec = NumCho(iSym)  ! break loop
         End If
      End Do

C     No first vector <=> 0 vectors in reduced set iRed.
C     --------------------------------------------------

      If (.not. Found) Then
         iVec = 0
         nVec = 0
         Return
      End If

C     Count number of vectors in reduced set iRed.
C     --------------------------------------------

      nVec = 1
      jVec = iVec
      Do While (jVec .lt. NumCho(iSym))
         jVec = jVec + 1
         jRed = InfVec(jVec,2,iSym)
         If (jRed .eq. iRed) Then
            nVec = nVec + 1
         Else
            jVec = NumCho(iSym) ! break loop
         End If
      End Do

#if defined (_DEBUGPRINT_)
C     Debug: print result.
C     --------------------

      Write(6,*) SecNam,': there are ',nVec,' vectors in reduced set ',
     &           iRed,' (sym. block ',iSym,')'
      Write(6,*) SecNam,': first vector is: ',iVec
#endif

      End
