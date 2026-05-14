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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine Copy_JobIph(InFile,OutFile)
! This is very nasty routine, actually a temporary hack!!!
! should be replaced to a proper copy of JobIph files!

use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: InFile, OutFile
integer(kind=iwp) :: ierr

call fcopy(InFile,OutFile,ierr)
if (iErr /= 0) call abend()
!Junk = 'cp -p '//File1(1:lFile1)//' '//File2(1:lFile2)

return

end subroutine Copy_JobIph
