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

subroutine SubWorkDir

use subdirs, only: f_setsubdir, Sub, OldWorkDir, NewWorkDir
use filesystem, only: getcwd_, chdir_, mkdir_
implicit none
integer :: i, Length, iErr
integer, parameter :: nFiles = 22
character(len=1024) :: Names(nFiles), OldFile(nFiles), NewFile(nFiles)
logical :: Found

! Define name of subdirectory and files that must be copied over
Sub = 'NG'
Names(1) = 'RUNFILE'
Names(2) = 'SEWARINP'
Names(3) = 'SCFINP'
Names(4) = 'RASSCINP'
Names(5) = 'CASPTINP'
Names(6) = 'MBPT2INP'
Names(7) = 'RASSIINP'
Names(8) = 'MOTRAINP'
Names(9) = 'CCSDTINP'
Names(10) = 'CHCCINP'
Names(11) = 'CHT3INP'
Names(12) = 'ESPFINP'
Names(13) = 'FALSEINP'
Names(14) = 'JOBIPH'
Names(15) = 'ESPF.SAV'
Names(16) = 'TINKER.XYZ'
Names(17) = 'TINKER.KEY'
Names(18) = 'MCPDFINP'
Names(19) = 'CHEMNATFIE'
Names(20) = 'CHEMCANFIE'
Names(21) = 'CHEMNATMPS0'
Names(22) = 'CHEMCANMPS0'

! Get real filenames to copy
do i=1,nFiles
  call prgmtranslate(Names(i),OldFile(i),Length)
end do

! Create the new directory and switch to it
call getcwd_(OldWorkDir)
NewWorkDir = trim(OldWorkDir)//'/'//trim(Sub)
call mkdir_(NewWorkDir)
call chdir_(NewWorkDir)
call f_setsubdir(Sub)

! Get real target filenames
do i=1,nFiles
  call prgmtranslate(Names(i),NewFile(i),Length)
  ! ESPF.SAV is copied to ESPF.DATA
  if (Names(i) == 'ESPF.SAV') then
    call prgmtranslate('ESPF.DATA',NewFile(i),Length)
  end if
end do

! Copy the files from the old directory to the new
do i=1,nFiles
  call f_inquire(OldFile(i),Found)
  if (Found) call fcopy(OldFile(i),NewFile(i),iErr)
end do

! The INPORB file is special...
call f_inquire('../INPORB',Found)
if (Found) call fcopy('../INPORB','./INPORB',iErr)

return

end subroutine SubWorkDir
