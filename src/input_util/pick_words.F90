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

subroutine Pick_Words(string,MaxNwords,Nwords,Words)

implicit integer(i-n)
character*(*) string
character*(*) Words(MaxNwords)
character*64 Word

if (string == ' ') then
  Nwords = 0
  goto 1000
end if
do i=1,MaxNwords
  Words(i) = ' '
end do
NWords = 0
i = 0
10 i = i+1
Word = ''
if (i > len(string)) goto 1000
if (string(i:i) == ' ') goto 10
NWords = NWords+1
ilett = 1
Word(ilett:ilett) = string(i:i) ! First letter of word
if (i == len(string)) goto 100
20 i = i+1
if (string(i:i) == ' ' .or. i > len(string)) goto 100
ilett = ilett+1
Word(ilett:ilett) = string(i:i)
goto 20
100 Words(NWords) = Word
if (i == len(string) .or. NWords == MaxNwords) goto 1000
goto 10

1000 return

end subroutine Pick_Words
