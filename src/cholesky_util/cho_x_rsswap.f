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
*  Cho_X_RSSwap
*
*> @brief
*>   Swap reduced set index arrays
*> @author Thomas Bondo Pedersen
*>
*> @details
*> Swap index arrays at locations \p iRS (``1``,``2``,``3``) and \p jRS (``1``,``2``,``3``).
*> If \p iRS and/or \p jRS are out of bounds, \p irc = ``1`` on exit and nothing
*> has been done!
*>
*> @warning
*> No special action is taken to redefine the \c IndRed array for first reduced set.
*>
*> @note
*> cholesky.fh must have been initialized.
*>
*> @param[out] irc return code
*> @param[in]  iRS location of reduced set
*> @param[in]  jRS location of reduced set
************************************************************************
      Subroutine Cho_X_RSSwap(irc,iRS,jRS)
      use ChoSwp, only: nnBstRSh
      Implicit None
      Integer irc, iRS, jRS
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Integer N, iOff, jOff, iTemp

      If (iRS.lt.1 .or. iRS.gt.3 .or. jRS.lt.1 .or. jRS.gt.3) Then
         irc = 1
      Else
         If (iRS .ne. jRS) Then
            N = nSym*nnShl
            iOff = ip_iiBstRSh + N*(iRs-1)
            jOff = ip_iiBstRSh + N*(jRs-1)
            Call iSwap(N,iWork(iOff),1,iWork(jOff),1)
            Call iSwap(N,nnBstRsh(:,:,iRS),1,nnBstRSh(:,:,jRS),1)
            Call iSwap(N,iWork(iOff),1,iWork(jOff),1)
            Call iSwap(nSym,iiBstR(1,iRS),1,iiBstR(1,jRS),1)
            Call iSwap(nSym,nnBstR(1,iRS),1,nnBstR(1,jRS),1)
            iOff = ip_IndRed + nnBstRT(1)*(iRs-1)
            jOff = ip_IndRed + nnBstRT(1)*(jRs-1)
            Call iSwap(nnBstRT(1),iWork(iOff),1,iWork(jOff),1)
            iTemp = nnBstRT(iRS)
            nnBstRT(iRS) = nnBstRT(jRS)
            nnBstRT(jRS) = iTemp
         End If
         irc = 0
      End If

      End
