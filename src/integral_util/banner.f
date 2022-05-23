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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      SubRoutine Banner(Lines,nLines,nWidth)
************************************************************************
*     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
*             University of Lund, SWEDEN                               *
*             May '91                                                  *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Parameter(MxWdth=132)
      Character*(*)   Lines(nLines)
      Character*(MxWdth-2) Line, Format*72
*
      mWidth = nWidth
      nChar = Len(Lines(1))
      If (nChar+2.gt.mWidth) mWidth = nChar + 2
      mWidth = Min(MxWdth-2,mWidth)
      Write (Format,'(A,i3,A)') '(1X,A',mWidth,')'
      Do 100 i = 1, mWidth
         Line(i:i) = '*'
 100  Continue
      Write (6,Format) Line
      Do 110 i = 2, mWidth-1
         Line(i:i) = ' '
 110  Continue
      Write (6,Format) Line
      Do 10 i = 1, nLines
         Do 20 j = 1, nChar
            If (Lines(i)(j:j).ne.' ') Go To 21
 20      Continue
 21      Continue
         iFrst = j
         Do 30 j = nChar, iFrst, -1
            If (Lines(i)(j:j).ne.' ') Go To 31
 30      Continue
 31      Continue
         iEnd = j
         Do 120 k = 2, mWidth-1
            Line(k:k) = ' '
 120     Continue
         Length = iEnd-iFrst+1
         nSplit = (mWidth-2-Length)/2
         jFrst = 1+nSplit+1
         jEnd = jFrst+Length-1
         Line(jFrst:jEnd) = Lines(i)(iFrst:iEnd)
         Write (6,Format) Line
 10   Continue
*
      Do 130 k = 2, mWidth-1
         Line(k:k) = ' '
 130  Continue
      Write (6,Format) Line
      Do 140 k = 2, mWidth-1
         Line(k:k) = '*'
 140  Continue
      Write (6,Format) Line
*
      Return
      End
