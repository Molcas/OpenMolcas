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

subroutine setposition(lunit, keyin, line, iRC)
  ! Read until, and include, a line beginning with a particular string in an ASCII file, assumed already opened, with unit number LUnit. That line is returned.
  ! key lengths up to 16 bytes can be used. It is determined by the size of the input variable
  use mcpdft_output, only: terse, lf, iPrLoc
  implicit none

  integer, intent(in) :: lunit
  character*(*), intent(in) :: keyin
  character*(*), intent(out) :: line
  integer, intent(out) :: iRC

#include "warnings.h"
  integer :: iPrLev, key_len
  character*16 :: command, key

  intrinsic :: len, min

  command = ' '
  key = ' '

  iPrLev = iPrLoc(1)
  iRC = _RC_ALL_IS_WELL_
  key_len = min(16, len(keyin))

  rewind(lunit)

  key(1:key_len) = keyin(1:key_len)
  call upcase(key)
10 continue
  read(lunit,'(A)',END=9910,ERR=9910) line
  command(1:key_len)=line(1:key_len)
  call upcase(command)
  if(command /= key) goto 10
  return

9910 continue
  if(iPrLev >= terse) then
    write(lf,*) 'SETPOSITION: Attempt to find an input line beginning'
    write(lf,*) 'with the keyword ''',KeyIn,''' failed.'
  end if
  iRC=_RC_INPUT_ERROR_
  return
end subroutine setposition
