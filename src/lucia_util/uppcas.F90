!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine UPPCAS(LINE,LENGTH)
! Convert letters in character string LINE to upper case
!
! very stupid and not vectorized!

character*(*) LINE
parameter(NCHAR=41)
character*1 LOWER(NCHAR)
character*1 UPPER(NCHAR)
data LOWER/'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','+','-','<', &
           '>','=','0','1','2','3','4','5','6','7','8','9'/
data UPPER/'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','+','-','<', &
           '>','=','0','1','2','3','4','5','6','7','8','9'/

do ICHA=1,LENGTH
  do I=1,NCHAR
    if (LINE(ICHA:ICHA) == LOWER(I)) LINE(ICHA:ICHA) = UPPER(I)
  end do
end do

end subroutine UPPCAS
