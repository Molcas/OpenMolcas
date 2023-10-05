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

use Constants, only: Eight
use Definitions, only: wp

implicit none
real(kind=wp), intent(in) :: Word
real(kind=wp), intent(out) :: Byte
character(len=2), intent(out) :: Unt

Byte = Word*Eight
Unt = 'b '
if (abs(Byte) > 1.0e3_wp) then
  Byte = Byte/1.024e3_wp
  Unt = 'kb'
  if (abs(Byte) > 1.0e3_wp) then
    Byte = Byte/1.024e3_wp
    Unt = 'Mb'
    if (abs(Byte) > 1.0e3_wp) then
      Byte = Byte/1.024e3_wp
      Unt = 'Gb'
      if (abs(Byte) > 1.0e3_wp) then
        Byte = Byte/1.024e3_wp
        Unt = 'Tb'
      end if
    end if
  end if
end if

end subroutine Cho_RWord2Byte
