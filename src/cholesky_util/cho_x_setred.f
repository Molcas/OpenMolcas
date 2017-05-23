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
      Subroutine Cho_X_SetRed(irc,iLoc,iRed)
************************************************************************
*
*   <DOC>
*     <Name>Cho\_X\_SetRed</Name>
*     <Syntax>Call Cho\_X\_SetRed(irc,iLoc,iRed)</Syntax>
*     <Arguments>
*        \Argument{irc}{return code}{Integer}{out}
*        \Argument{iLoc}{location in index arrays}{Integer}{in}
*        \Argument{iRed}{reduced set on disk}{Integer}{in}
*     </Arguments>
*     <Purpose>
*        Read and set index arrays for reduced set iRed at location
*        iLoc
*     </Purpose>
*     <Dependencies>Call Cho\_X\_Init</Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*        Reads information for reduced set iRed (=1,2,...,MaxRed)
*        and sets up the index arrays
*        \begin{itemize}
*          \item nnBstRT(iLoc)      $\rightarrow$ stored in cholesky.fh
*          \item nnBstR(:,iLoc)     $\rightarrow$ stored in cholesky.fh
*          \item iiBstR(:,iLoc)     $\rightarrow$ stored in cholesky.fh
*          \item nnBstRSh(:,:,iLoc) $\rightarrow$ accesible via
*                                   ip\_nnBstRSh in choptr.fh
*          \item iiBstRSh(:,:,iLoc) $\rightarrow$ accesible via
*                                   ip\_iiBstRSh in choptr.fh
*          \item IndRed(:,iLoc)     $\rightarrow$ accesible via
*                                   ip\_IndRed in choptr.fh
*        \end{itemize}
*        On succesful completion, irc=0 is returned.
*        Note that the only allowed iLoc values are 2 and 3; any other
*        value gives rise to error code irc=1 and nothing is done!
*        If iRed is out of bounds, irc=2 is returned and nothing is
*        done!
*     </Description>
*    </DOC>
*
************************************************************************
#include "implicit.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      If (iLoc.eq.2 .or. iLoc.eq.3) Then
         If (iRed.lt.1 .or. iRed.gt.MaxRed) Then
            irc = 2
         Else
            kOff1 = ip_nnBstRSh + nSym*nnShl*(iLoc - 1)
            kOff2 = ip_IndRed   + nnBstRT(1)*(iLoc - 1)
            Call Cho_GetRed(iWork(ip_InfRed),iWork(kOff1),
     &                      iWork(kOff2),iWork(ip_IndRSh),
     &                      iWork(ip_iSP2F),
     &                      MaxRed,nSym,nnShl,nnBstRT(1),iRed,.false.)
            Call Cho_SetRedInd(iWork(ip_iiBstRSh),iWork(ip_nnBstRSh),
     &                         nSym,nnShl,iLoc)
            irc = 0
            If (iRed .eq. 1) Then ! set correct IndRed array
               kOff2 = kOff2 - 1
               Do iab = 1,nnBstRT(1)
                  kOff = kOff2 + iab
                  iWork(kOff) = iab
               End Do
            End If
         End If
      Else
         irc = 1
      End If

      End
