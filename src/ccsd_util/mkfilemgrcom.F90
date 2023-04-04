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

subroutine mkfilemgrcom()
! this routine makes names of temp files

use ccsd_global, only: filename, filerst, maxfiles, minfiles
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nhelp

!3 def filenames

do nhelp=minfiles,maxfiles
  write(filename(nhelp),'("Temp",I2.2)') nhelp
end do
filename(10) = 'INTAB'
filename(11) = 'INTA1'
filename(12) = 'INTA2'
filename(13) = 'INTA3'
filename(14) = 'INTA4'
filename(15) = 'INTSTA'
filename(16) = filerst

return

end subroutine mkfilemgrcom
