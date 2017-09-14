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
      Subroutine StdFmt(NameIn,NameUt)
************************************************************************
*                                                                      *
*     Convert the first word of the incoming string, NameIn, to        *
*     standard format (no leading blanks, all upper case letters).     *
*                                                                      *
*     calling arguments                                                *
*     NameIn  : Type character string, input                           *
*     NameUt  : Type character string, output                          *
*               First token of NameIn in upper case letters            *
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
      Character*(*) NameIn
      Character*(*) NameUt
*----------------------------------------------------------------------*
*     start                                                            *
*----------------------------------------------------------------------*
      lIn=Len(NameIn)
      lUt=Len(NameUt)
*----------------------------------------------------------------------*
*     set NameUt to a blank string                                     *
*----------------------------------------------------------------------*
      NameUt=' '
c      Do 100 i=1,lUt
c         NameUt(i:i)=' '
c100   Continue
*----------------------------------------------------------------------*
*     get the number of leading blanks                                 *
*----------------------------------------------------------------------*
c      lLeft=0
      do 200 i=1,lIn
       if(NameIn(i:i).ne.' ') goto 201
200   continue
201   lLeft=i

c      Do 200 i=lIn,1,-1
c         If( NameIn(i:i).ne.' ' ) lLeft=i-1
c200   Continue
*----------------------------------------------------------------------*
*     copy only the first token of NameIn to NameUt                    *
*----------------------------------------------------------------------*
      j=0
      Do 300 i=lLeft,lIn
         If ( NameIn(i:i).eq.' ' .or. j.eq.lUt ) Goto 400
         j=j+1
         NameUt(j:j)=NameIn(i:i)
300   Continue
*----------------------------------------------------------------------*
*     change lower case character to upper case                        *
*----------------------------------------------------------------------*
400   Call UpCase(NameUt)
*----------------------------------------------------------------------*
*     exit                                                             *
*----------------------------------------------------------------------*
      Return
      End
