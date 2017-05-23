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
      Subroutine RdNLst(iUnit,NameIn)
      Character*(*) NameIn
      Logical No_Input_OK
*
      No_Input_OK=.False.
      Call RdNLst_(iUnit,NameIn,No_Input_OK)
*
      Return
      End
      Subroutine RdNLst_(iUnit,NameIn,No_Input_OK)
************************************************************************
*                                                                      *
*     Locate the beginning of an input stream                          *
*     (similar to FORTRAN NAMELIST read known to some systems)         *
*                                                                      *
*     calling arguments:                                               *
*     iUnit  : Type integer, input                                     *
*              FORTRAN unit number                                     *
*     NameIn : Type character string, input                            *
*              Character string marking the beginning of the input     *
*     No_Input_OK: Logical                                             *
*                  On input determines if an input has to be found     *
*                  On exit determines if an input was found.           *
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
      Character*8 StdNam
      Character*80 Line
      Integer StrnLn
      Logical No_Input_OK
      common /igetline/igetline,myunit
        igetline=0
*----------------------------------------------------------------------*
*     push the entry name on the calling stack                         *
*----------------------------------------------------------------------*
*     Call qEnter('RdNlst')
*----------------------------------------------------------------------*
*     convert the Name to internal standard format.                    *
*----------------------------------------------------------------------*
      Call StdFmt(NameIn,StdNam)
      lStdNam=StrnLn(StdNam)
*----------------------------------------------------------------------*
*     read until an input Line is located which starts with            *
*     the string, Name, not before the second column                   *
*----------------------------------------------------------------------*
      lLine=LEN(Line)
100   Read(iUnit,'(A)',End=900) Line
      Call LeftAd(Line)
      Call UpCase(Line)
      If ( Line(1:1).eq.'&' .and.
     &     Line(2:lStdNam+1).eq.StdNam(1:lStdNam) ) then
         Return
      End If
      Goto 100
*----------------------------------------------------------------------*
*     error exit                                                       *
*----------------------------------------------------------------------*
900   If (No_Input_OK) then
         No_Input_OK=.False.
         Return
      EndIf
      Write (6,*) 'RdNLst: Input section not found in input file'
      Write (6,*) '        Looking for:',StdNam(1:lStdNam)
      Call QTrace()
      Call Quit_OnUserError()
      End
