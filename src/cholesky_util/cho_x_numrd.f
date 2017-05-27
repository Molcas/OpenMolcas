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
      Integer Function Cho_X_NumRd(iVec1,iSym,iRedC,Mem)
************************************************************
*
*   <DOC>
*     <Name>Cho\_X\_NumRd</Name>
*     <Syntax>Cho\_X\_NumRd(iVec1,iSym,iRedC,Mem)</Syntax>
*     <Arguments>
*       \Argument{iVec1}{First vector}{Integer}{in}
*       \Argument{iSym}{Symmetry}{Integer}{in}
*       \Argument{iRedC}{Reduced set in core (location 3);
*                        0 (or -1) if unknown or undefined}
*                {Integer}{in}
*       \Argument{Mem}{Memory available for read}{Integer}{in}
*     </Arguments>
*     <Purpose>
*        Return the number of Cholesky vectors that may be read
*        into Mem words of memory.
*     </Purpose>
*     <Dependencies>
*        The Cholesky procedure must have been successfully
*        initialized (by Cho\_X\_Init).
*     </Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*        The count starts at vector iVec1 of symmetry iSym
*        (this is needed since the vectors are stored in
*        different reduces sets).
*        On exit, Cho\_X\_NumRd is negative if some error has
*        occurred (-1, -2, and -3 signify errors in input
*        variables, -4 indicates an error in Cho\_X\_SetRed).
*     </Description>
*    </DOC>
*
************************************************************

#include "implicit.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Parameter (N2 = InfVec_N2)
      Integer iRed

      InfVec(i,j,k)=iWork(ip_InfVec-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)
      nDimRS(i,j)=iWork(ip_nDimRS-1+nSym*(j-1)+i)

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
         If (l_nDimRS .lt. 1) Then
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
