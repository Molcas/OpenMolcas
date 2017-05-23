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
*               1993, Per-Olof Widmark                                 *
************************************************************************
      Subroutine RightAd(String)
************************************************************************
*                                                                      *
*     Left adjust a character string.                                  *
*                                                                      *
*     calling argument:                                                *
*     String  : Character string                                       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. P. Fuelscher and P.O. Widmark                                 *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Character*(*) String
*----------------------------------------------------------------------*
*     Get the length of the line                                       *
*----------------------------------------------------------------------*
      lString=LEN(String)
*----------------------------------------------------------------------*
*     get the number of trailing blanks                                *
*----------------------------------------------------------------------*
      lRight=0
      Do 100 i=1,lString
         If( String(i:i).ne.' ' ) lRight=lString-i
100   Continue
*----------------------------------------------------------------------*
*     Shift the line as long as there are leading blanks               *
*----------------------------------------------------------------------*
      If ( lRight.eq.0 ) Return
      nShift=lRight
      Do 210 i=lString,nShift+1,-1
         String(i:i)=String(i-nShift:i-nShift)
210   Continue
      Do 220 i=nShift,1,-1
         String(i:i)=' '
220   Continue
*----------------------------------------------------------------------*
*     Normal termination                                               *
*----------------------------------------------------------------------*
      Return
      End
