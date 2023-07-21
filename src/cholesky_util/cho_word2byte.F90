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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************
!  Cho_Word2Byte
!
!> @brief
!>   Convert integer number of \p n -byte words to byte with a reasonable prefix (``-``, ``k``, ``M``, ``G``, or ``T``)
!> @author Thomas Bondo Pedersen
!>
!> @details
!> Convert integer number of \p n -byte words to byte, kilobyte,
!> megabyte, gigabyte, or terabyte.
!>
!> @param[in]  iWord Number of \p n -byte words
!> @param[in]  n     Number of byte per word
!> @param[out] Byte  \p iWord in bytes/kb/Mb/Gb/Tb
!> @param[out] Unt   Unit of Byte ['``b ``', '``kb``', '``Mb``', '``Gb``', or '``Tb``']
!***********************************************************************

subroutine Cho_Word2Byte(iWord,n,Byte,Unt)

implicit none
integer iWord, n
real*8 Byte
character*2 Unt

Byte = dble(iWord)*dble(n)
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

end subroutine Cho_Word2Byte
