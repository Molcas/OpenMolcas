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
* Copyright (C) 1993, Markus P. Fuelscher                              *
************************************************************************
      Subroutine Center(String)
************************************************************
*
*   <DOC>
*     <Name>Center</Name>
*     <Syntax>Call Center(String)</Syntax>
*     <Arguments>
*       \Argument{String}{a string}{Character *(*)}{inout}
*     </Arguments>
*     <Purpose>Center a string for printing purpose</Purpose>
*     <Dependencies></Dependencies>
*     <Author>M.P. Fuelscher and P.O. Widmark</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*        Add spaces to the beginning and the end of a string
*     </Description>
*    </DOC>
*
************************************************************
      Character*(*) String
*----------------------------------------------------------------------*
*     get the length of the line                                       *
*----------------------------------------------------------------------*
      lString=0
      lString=Len(String)
*----------------------------------------------------------------------*
*     get the number of leading blanks                                 *
*----------------------------------------------------------------------*
      lLeft=0
      Do 100 i=lString,1,-1
        If( String(i:i).ne.' ' ) lLeft=i-1
100   Continue
*----------------------------------------------------------------------*
*     get the number of trailing blanks                                *
*----------------------------------------------------------------------*
      lRight=0
      Do 200 i=1,lString
        If( String(i:i).ne.' ' ) lRight=lString-i
200   Continue
*----------------------------------------------------------------------*
*     shift the line                                                   *
*----------------------------------------------------------------------*
      If ( lLeft+lRight.ne.0 ) Then
         lShift=(lRight-lLeft)/2
         If ( lShift.gt.0 ) Then
           Do 310 i=lString,lShift+1,-1
             String(i:i)=String(i-lShift:i-lShift)
310        Continue
           Do 320 i=1,lLeft+lShift
             String(i:i)=' '
320        Continue
         Else If ( lShift.lt.0 ) Then
           Do 410 i=1,lString-lShift
             String(i:i)=String(i-lShift:i-lShift)
410        Continue
           Do 420 i=lString,lString-lRight-lShift,-1
             String(i:i)=' '
420        Continue
         End If
      End If
*----------------------------------------------------------------------*
*     normal termination                                               *
*----------------------------------------------------------------------*
      Return
      End
