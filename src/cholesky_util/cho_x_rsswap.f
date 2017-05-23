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
      Subroutine Cho_X_RSSwap(irc,iRS,jRS)
************************************************************************
*
*   <DOC>
*     <Name>Cho\_X\_RSSwap</Name>
*     <Syntax>Call Cho\_X\_RSSwap(irc,iRS,jRS)</Syntax>
*     <Arguments>
*       \Argument{irc}{return code}{Integer}{out}
*       \Argument{iRS}{location of reduced set}{Integer}{in}
*       \Argument{jRS}{location of reduced set}{Integer}{in}
*     </Arguments>
*     <Purpose>Swap reduced set index arrays</Purpose>
*     <Dependencies>cholesky.fh initialized</Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*        Swap index arrays at locations iRS (1,2,3) and jRS (1,2,3).
*        If iRS and/or jRS are out of bounds, irc=1 on exit and nothing
*        has been done!
*        Warning: no special action is taken to redefine the IndRed
*        array for first reduced set.
*     </Description>
*    </DOC>
*
************************************************************************

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
            iOff = ip_nnBstRSh + N*(iRs-1)
            jOff = ip_nnBstRSh + N*(jRs-1)
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
