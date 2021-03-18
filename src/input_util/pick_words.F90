!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2006, Giovanni Ghigo                                   *
!***********************************************************************
!***********************************************************************
! This routine finds in a -string- up to -MaxNWords- and put the words *
! in the array -Words(MaxNwords)- and the number in -Nwords-           *
!----------------------------------------------------------------------*
! Author:   Giovanni Ghigo - October 2006 - Torino (IT)                *
!***********************************************************************
      Subroutine Pick_Words(string,MaxNwords,Nwords,Words)
      Implicit Integer (i-n)
      Character*(*) string
      Character*(*) Words(MaxNwords)
      Character*64  Word

      If (string.eq.' ') then
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
