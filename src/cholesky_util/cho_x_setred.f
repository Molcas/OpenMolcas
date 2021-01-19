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
*  Cho_X_SetRed
*
*> @brief
*>   Read and set index arrays for reduced set \p iRed at location \p iLoc
*> @author Thomas Bondo Pedersen
*>
*> @details
*> Reads information for reduced set \p iRed (= ``1``, ``2``, ..., \c MaxRed)
*> and sets up the index arrays
*>
*> - \p nnBstRT(iLoc)      &rarr; stored in cholesky.fh
*> - \p nnBstR(:,iLoc)     &rarr; stored in cholesky.fh
*> - \p iiBstR(:,iLoc)     &rarr; stored in cholesky.fh
*> - \p nnBstRSh(:,:,iLoc) &rarr; accesible via choswp.f90
*> - \p iiBstRSh(:,:,iLoc) &rarr; accesible via choswp.f90
*> - \p IndRed(:,iLoc)     &rarr; accesible via choswp.f90
*>
*> On succesful completion, \p irc = ``0`` is returned.
*> Note that the only allowed \p iLoc values are ``2`` and ``3``; any other
*> value gives rise to error code \p irc = ``1`` and nothing is done!
*> If \p iRed is out of bounds, \p irc = ``2`` is returned and nothing is done!
*>
*> @note
*> The Cholesky procedure must have been successfully initialized (by ::Cho_X_Init).
*>
*> @param[out] irc  return code
*> @param[in]  iLoc location in index arrays
*> @param[in]  iRed reduced set on disk
************************************************************************
      Subroutine Cho_X_SetRed(irc,iLoc,iRed)
      use ChoArr, only: iSP2F
      use ChoSwp, only: nnBstRSh, iiBstRSh, IndRSh, InfRed, IndRed
#include "implicit.fh"
#include "cholesky.fh"

      If (iLoc.eq.2 .or. iLoc.eq.3) Then
         If (iRed.lt.1 .or. iRed.gt.MaxRed) Then
            irc = 2
         Else
            Call Cho_GetRed(InfRed,nnBstRSh(:,:,iLoc),
     &                      IndRed(:,iLoc),IndRSh,iSP2F,
     &                      MaxRed,nSym,nnShl,nnBstRT(1),iRed,.false.)
            Call Cho_SetRedInd(iiBstRSh,nnBstRSh,nSym,nnShl,iLoc)
            irc = 0
            If (iRed .eq. 1) Then ! set correct IndRed array
               Do iab = 1,nnBstRT(1)
                  IndRed(iab,iLoc) = iab
               End Do
            End If
         End If
      Else
         irc = 1
      End If

      End
