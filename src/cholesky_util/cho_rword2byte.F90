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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_RWord2Byte(Word,Byte,Unt)

implicit none
real*8 Word
real*8 Byte
character*2 Unt

Byte = Word*8.0d0
Unt = 'b '
if (abs(Byte) > 1.0d3) then
  Byte = Byte/1.024d3
  Unt = 'kb'
  if (abs(Byte) > 1.0d3) then
    Byte = Byte/1.024d3
    Unt = 'Mb'
    if (abs(Byte) > 1.0d3) then
      Byte = Byte/1.024d3
      Unt = 'Gb'
      if (abs(Byte) > 1.0d3) then
        Byte = Byte/1.024d3
        Unt = 'Tb'
      end if
    end if
  end if
end if

end subroutine Cho_RWord2Byte
