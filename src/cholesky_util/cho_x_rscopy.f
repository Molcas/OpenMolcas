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
      Subroutine Cho_X_RSCopy(irc,iRS,jRS)
************************************************************************
*
*   <DOC>
*     <Name>Cho\_X\_RSCopy</Name>
*     <Syntax>Call Cho\_X\_RSCopy(irc,iRS,jRS)</Syntax>
*     <Arguments>
*       \Argument{irc}{return code}{Integer}{out}
*       \Argument{iRS}{location of reference reduced set}{Integer}{in}
*       \Argument{jRS}{location of target reduced set}{Integer}{in}
*     </Arguments>
*     <Purpose>Copy reduced set index arrays</Purpose>
*     <Dependencies>cholesky.fh initialized</Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*        Copy index arrays from location iRS (1,2,3) to location
*        jRS (1,2,3). Special handling of the case iRS=1:
*        IndRed(i,jRS) = i rather than IndRed(i,jRS)=IndRed(i,iRS).
*        If iRS and/or jRS are out of bounds, irc=1 on exit and nothing
*        has been done!
*     </Description>
*    </DOC>
*
************************************************************************

#include "implicit.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      If (iRS.lt.1 .or. iRS.gt.3 .or. jRS.lt.1 .or. jRS.gt.3) Then
         irc = 1
      Else
         Call Cho_RSCopy(iWork(ip_iiBstRSh),iWork(ip_nnBStRSh),
     &                   iWork(ip_IndRed),iRS,jRS,nSym,nnShl,nnBstRT(1),
     &                   3)
         irc = 0
      End If

      End
