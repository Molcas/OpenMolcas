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
! Copyright (C) 2003, Per-Olof Widmark                                 *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine sets the name for the runfile.                          *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund university, Sweden                                     *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************

subroutine NameRun(fName)

use RunFile_data, only: RnNmStk, RunName

implicit none
character(len=*), intent(in) :: fName

if (fName == '#Pop') then
  RunName = RnNmStk(1)
  RnNmStk(1) = RnNmStk(2)
  RnNmStk(2) = RnNmStk(3)
  RnNmStk(3) = RnNmStk(4)
  RnNmStk(4) = ''
else
  RnNmStk(4) = RnNmStk(3)
  RnNmStk(3) = RnNmStk(2)
  RnNmStk(2) = RnNmStk(1)
  RnNmStk(1) = RunName
  RunName = fName
end if

call ClrRunCache()

! Do not create the run file when naming it, patch 6.7.263

!iRc = 0
!iOpt = 1
!call MkRun(iRc,iOpt)

return

end subroutine NameRun
