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

use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: string
integer(kind=iwp), intent(in) :: MaxNwords
integer(kind=iwp), intent(out) :: Nwords
character(len=*), intent(out) :: Words(MaxNwords)
integer(kind=iwp) :: i, ilett
character(len=64) :: Word

if (string == ' ') then
  Nwords = 0
  return
end if
Words(:) = ' '
NWords = 0
i = 0
do
  i = i+1
  Word = ''
  if (i > len(string)) exit
  if (string(i:i) == ' ') cycle
  NWords = NWords+1
  ilett = 1
  Word(ilett:ilett) = string(i:i) ! First letter of word
  if (i /= len(string)) then
    do
      i = i+1
      if ((string(i:i) == ' ') .or. (i > len(string))) exit
      ilett = ilett+1
      Word(ilett:ilett) = string(i:i)
    end do
  end if
  Words(NWords) = Word
  if ((i == len(string)) .or. (NWords == MaxNwords)) exit
end do

return

end subroutine Pick_Words
