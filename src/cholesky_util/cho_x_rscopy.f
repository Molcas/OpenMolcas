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
*  Cho_X_RSCopy
*
*> @brief
*>   Copy reduced set index arrays
*> @author Thomas Bondo Pedersen
*>
*> @details
*> Copy index arrays from location \p iRS (``1``,``2``,``3``) to location
*> \p jRS (``1``,``2``,``3``). Special handling of the case \p iRS = ``1``:
*> ``IndRed(i,jRS) = i`` rather than ``IndRed(i,jRS)=IndRed(i,iRS)``.
*> If \p iRS and/or \p jRS are out of bounds, \p irc = ``1`` on exit and nothing
*> has been done!
*>
*> @note
*> cholesky.fh must have been initialized.
*>
*> @param[out] irc return code
*> @param[in]  iRS location of reference reduced set
*> @param[in]  jRS location of target reduced set
************************************************************************
      Subroutine Cho_X_RSCopy(irc,iRS,jRS)
      use ChoSwp, only: nnBstRsh, iiBstRSh
#include "implicit.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      If (iRS.lt.1 .or. iRS.gt.3 .or. jRS.lt.1 .or. jRS.gt.3) Then
         irc = 1
      Else
         Call Cho_RSCopy(iiBstRSh,nnBStRSh,
     &                   iWork(ip_IndRed),iRS,jRS,nSym,nnShl,nnBstRT(1),
     &                   3)
         irc = 0
      End If

      End
