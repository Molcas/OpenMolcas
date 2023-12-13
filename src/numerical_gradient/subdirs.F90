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
! Copyright (C) 2015, Ignacio Fdez. Galvan                             *
!***********************************************************************
! Interfaces to C functions needed to create and remove
! a subdirectory in WorkDir

module subdirs

implicit none
private

character(len=1024) :: Sub, OldWorkDir, NewWorkDir

public :: f_setsubdir, Sub, OldWorkDir, NewWorkDir

contains

subroutine f_setsubdir(sub)
# ifdef _HAVE_EXTRA_
  use, intrinsic :: iso_c_binding, only: c_null_char
  implicit none
  character(len=*), intent(in) :: sub
  interface
    subroutine c_setsubdir(sub) bind(C,name="setsubdir")
      use, intrinsic :: iso_c_binding, only: c_char
      character(kind=c_char) :: sub(*)
    end subroutine c_setsubdir
  end interface
  if (trim(sub) == '') then
    call c_setsubdir(''//c_null_char)
  else
    call c_setsubdir('/'//trim(sub)//c_null_char)
  end if
# else
  use Prgm, only: SetSubDir
  implicit none
  character(len=*), intent(in) :: sub
  if (trim(sub) == '') then
    call SetSubDir('')
  else
    call SetSubDir('/'//trim(sub))
  end if
# endif
end subroutine f_setsubdir

end module subdirs
