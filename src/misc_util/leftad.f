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
      Subroutine LeftAd(String)
************************************************************************
*                                                                      *
*     Left adjust a character string.                                  *
*                                                                      *
*     calling argument                                                 *
*     String  : Type character string, input/output                    *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher and P.O. Widmark                                  *
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
*     get the number of leading blanks                                 *
*----------------------------------------------------------------------*
      lLeft=0
      Do 100 i=lString,1,-1
         If( String(i:i).ne.' ' ) lLeft=i-1
100   Continue
*----------------------------------------------------------------------*
*     Shift the line and insert trailing blanks                        *
*----------------------------------------------------------------------*
      If ( lLeft.eq.0 ) Return
      nShift=lLeft
      Do 210 i=1,lString-nShift
          String(i:i)=String(i+nShift:i+nShift)
210   Continue
      Do 220 i=lString-nShift+1,lString
          String(i:i)=' '
220   Continue
*----------------------------------------------------------------------*
*     Normal termination                                               *
*----------------------------------------------------------------------*
      Return
      End
