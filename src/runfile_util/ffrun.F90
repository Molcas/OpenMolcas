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
! This routine locates a field in the runfile.                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund university, Sweden                                     *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************

subroutine ffRun(Label,nData,RecTyp)

use RunFile_data, only: rcNotFound, rcOK, TypUnk
use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: Label
integer(kind=iwp), intent(out) :: nData, RecTyp
integer(kind=iwp) :: iOpt, iRc
character(len=64) :: ErrMsg

!----------------------------------------------------------------------*
! Call extended version routine.                                       *
!----------------------------------------------------------------------*
iRc = rcOK
iOpt = 0
call ffxRun(iRc,Label,nData,RecTyp,iOpt)
if (iRc == rcNotFound) then
  nData = 0
  RecTyp = TypUnk
else if (iRc /= rcOK) then
  write(ErrMsg,'(3a)') 'Error locating field "',Label,'" in runfile'
  call SysAbendMsg('ffRun',ErrMsg,' ')
end if
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine ffRun
