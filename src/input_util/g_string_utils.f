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
* Copyright (C) 2004,2006, Giovanni Ghigo                              *
************************************************************************
************************************************************************
* This routine check if -string- is an integer. If it is a number      *
* NaN=.False. and -iNumber- contains the number, otherwise, NaN=.True. *
*----------------------------------------------------------------------*
* Author:   Giovanni Ghigo - November 2004 - Lund(SE)                  *
* Author:   Giovanni Ghigo                                             *
************************************************************************
      Subroutine Get_iNumber(string,iNumber,iErr)
      Implicit Integer (i-n)
      Character*(*)  string
      Character*11  Chars
      Data Chars /' 1234567890' /
      Logical NaN
      iErr=0
      iNumber=0
      NaN=.True.
      Do i=1,len(string)
      NaN=.True.
        Do j=1,11
          If (string(i:i).EQ.Chars(j:j)) NaN=.False.
        EndDo
      If (NaN) GoTo 10
      End Do
10    If (.NOT.NaN) Read(string,*) iNumber
      If (NaN) iErr=1
      Return
      End

************************************************************************
* This routine check if -string- is a Real*8 number. If it is a number *
* NaN=.False. and -Number- contains the number, otherwise, NaN=.True.  *
*----------------------------------------------------------------------*
* Author:   Giovanni Ghigo - November 2004 - Lund(SE)                  *
* Author:   Giovanni Ghigo                                             *
************************************************************************
      Subroutine Get_dNumber(string,dNumber,iErr)
      Implicit Integer (i-n)
      Implicit Real*8 (a-h,o-z)
      Character*(*)  string
      Character*14  Chars
      Data Chars /' +-1234567890.' /
      Logical NaN
      iErr=0
      iNumber=0
      NaN=.True.
      Do i=1,len(string)
      NaN=.True.
        Do j=1,14
          If (string(i:i).EQ.Chars(j:j)) NaN=.False.
        EndDo
      If (NaN) GoTo 10
      End Do
10    If (.NOT.NaN) Read(string,*) dNumber
      If (NaN) iErr=1
      Return
      End

************************************************************************
* This routine check if -string- is empty (spaces are empty).          *
* If empty, Logical Empty=.True.                                       *
*----------------------------------------------------------------------*
* Author:   Giovanni Ghigo - November 2004 - Lund(SE)                  *
************************************************************************
      Subroutine Empty_(string,Empty)
      Implicit Integer (i-n)
      Character*(*)  string
      Logical Empty
      Empty=.True.
      Do i=1,len(string)
        If (string(i:i).NE.' ') Empty=.False.
      End Do
      Return
      End


************************************************************************
* This routine find in a -string- up to -MaxNWords- and put the words  *
* in the array -Words(MaxNwords)- and the number in -Nwords-           *
*----------------------------------------------------------------------*
* Author:   Giovanni Ghigo - October 2006 - Torino (IT)                *
************************************************************************
      Subroutine Pick_Words(string,MaxNwords,Nwords,Words)
      Implicit Integer (i-n)
      Character*(*) string
      Character*(*) Words(MaxNwords)
      Character*64  Word
      Logical Empty

      Call Empty_(string,Empty)
      If (Empty) then
        Nwords = 0
        GoTo 1000
      EndIf
      Do i = 1, MaxNwords
        Words(i) = ' '
      EndDo
      NWords = 0
      i = 0
10    i = i + 1
      Word = ''
      If (i.GT.Len(string)) GoTo 1000
      If (string(i:i).EQ.' ') GoTo 10
        NWords = NWords + 1
        ilett  = 1
        Word(ilett:ilett) = string(i:i) ! First letter of word
        If (i.EQ.Len(string)) GoTo 100
20        i = i + 1
          If (string(i:i).EQ.' ' .OR. i.GT.Len(string)) GoTo 100
          ilett  = ilett + 1
          Word(ilett:ilett) = string(i:i)
        Goto 20
100     Words(NWords) = Word
      If (i.EQ.Len(string) .OR. NWords.EQ.MaxNwords) GoTo 1000
      GoTo 10

1000  Return
      End


      Subroutine FoundAtomicNumber(LuWr,Symb,NAT,iErr)
      Implicit Integer (i-n)
      Character*6 Symb
#include "periodic_table.fh"
      If ((Symb(1:1).le.'z').and.(Symb(1:1).ge.'a'))
     &     Symb(1:1)=char(ichar(Symb(1:1))-32)
      If ((Symb(2:2).le.'Z').and.(Symb(2:2).ge.'A'))
     &     Symb(2:2)=char(ichar(Symb(2:2))+32)
      iErr = 1
* --- Z atoms are ghost atoms that will disappear in seward
      If (Symb(1:1).EQ.'Z' .AND.
     &    Symb(2:2).NE.'n' .AND. Symb(2:2).NE.'r') then
        iErr = 0
        NAT = -1
        Return
      EndIf
* --- X atoms are dummy atoms that will appear in seward
      If (Symb(1:1).EQ.'X' .AND. Symb(2:2).NE.'e') then
        iErr = 0
        NAT =  0
        Return
      EndIf

      Do i = 1, Num_Elem
         If (Symb(1:2).EQ.AdjustL(PTab(i))) then
          iErr = 0
          NAT = i
          Return
        EndIf
      EndDo
      Do i = 1, Num_Elem
        If (Symb(1:1)//' '.EQ.AdjustL(PTab(i))) then
          iErr = 0
          NAT = i
          Return
        EndIf
      EndDo
      Write(LuWr,*) '   [FoundAtomicNumber]: Wrong atomic symbol !'
      Return
      End
