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
*  Cho_X_NumRd
*
*> @brief
*>   Return the number of Cholesky vectors that may be read into \p Mem words of memory.
*> @author Thomas Bondo Pedersen
*>
*> @details
*> The count starts at vector \p iVec1 of symmetry \p iSym
*> (this is needed since the vectors are stored in
*> different reduces sets).
*> On exit, \p Cho_X_NumRd is negative if some error has
*> occurred (``-1``, ``-2``, and ``-3`` signify errors in input
*> variables, ``-4`` indicates an error in \p Cho_X_SetRed).
*>
*> @note
*> The Cholesky procedure must have been successfully initialized (by ::Cho_X_Init).
*>
*> @param[in] iVec1 First vector
*> @param[in] iSym  Symmetry
*> @param[in] iRedC Reduced set in core (location ``3``); ``0`` (or ``-1``) if unknown or undefined
*> @param[in] Mem   Memory available for read
************************************************************************
      Integer Function Cho_X_NumRd(iVec1,iSym,iRedC,Mem)
      use ChoArr, only: nDimRS
      use ChoSwp, only: InfVec
#include "implicit.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Integer iRed

      If (iSym.lt.1 .or. iSym.gt.nSym) Then
         Cho_X_NumRd = -1
      Else If (NumCho(iSym).lt.1) Then
         Cho_X_NumRd = 0
      Else If (NumCho(iSym).gt.MaxVec) Then
         Cho_X_NumRd = -2
      Else If (iVec1.lt.1 .or. iVec1.gt.NumCho(iSym)) Then
         Cho_X_NumRd = -3
      Else If (Mem .lt. 1) Then
         Cho_X_NumRd = 0
      Else
         NumRd = 0
         Need  = 0
         iVec  = iVec1 - 1
         If (.NOT.Allocated(nDimRS)) Then
            iLoc  = 3
            Do While (iVec.lt.NumCho(iSym) .and. Need.lt.Mem)
               iVec = iVec + 1
               iRed = InfVec(iVec,2,iSym)
               If (iRed .ne. iRedC) Then
                  irc = 0
                  Call Cho_X_SetRed(irc,iLoc,iRed)
                  If (irc .ne. 0) Then
                     Cho_X_NumRd = -4
                     iRedC = -1
                     Return
                  End If
                  iRedC = iRed
               End If
               Need = Need + nnBstR(iSym,iLoc)
               If (Need .le. Mem) Then
                  NumRd = NumRd + 1
               End If
            End Do
         Else
            Do While (iVec.lt.NumCho(iSym) .and. Need.lt.Mem)
               iVec = iVec + 1
               iRed = InfVec(iVec,2,iSym)
               Need = Need + nDimRS(iSym,iRed)
               If (Need .le. Mem) Then
                  NumRd = NumRd + 1
               End If
            End Do
         End If
         Cho_X_NumRd = NumRd
      End If

      End
