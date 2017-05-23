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
* Copyright (C) Per Ake Malmqvist                                      *
************************************************************************
************************************************************************
*                                                                      *
* This routine uppercases a text string.                               *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-AAke Malmqvist                                          *
*          Lund University                                             *
*                                                                      *
************************************************************************
      Subroutine LoCase(string)
      Character*(*)  string
      Character*26   up,lw
      Dimension      itab(0:255)
      Save           up,lw,ifset,itab

      Data up    /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      Data lw    /'abcdefghijklmnopqrstuvwxyz'/
      Data ifset / 0 /

      If (ifset.eq.0) Then
        ifset=1
        Do i=0,255
          itab(i)=i
        End Do
        Do ii=1,26
          i=ichar(up(ii:ii))
          j=ichar(lw(ii:ii))
          itab(i)=j
        End Do
      End If

      Do ii=1,len(string)
        i=ichar(string(ii:ii))
        j=itab(i)
        string(ii:ii)=char(j)
      End Do

      Return
      End
