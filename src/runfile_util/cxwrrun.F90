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
! This routine writes a record into the runfile.                       *
! Data type is Character.                                              *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund university, Sweden                                     *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************

subroutine cxWrRun(iRc,Label,cData,nData,iOpt)

use RunFile_data, only: TypStr
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: iRc
character(len=*), intent(in) :: Label
character, intent(in) :: cData(*)
integer(kind=iwp), intent(in) :: nData, iOpt
character(len=64) :: ErrMsg

!----------------------------------------------------------------------*
! Check that arguments are ok.                                         *
!----------------------------------------------------------------------*
if (iOpt /= 0) then
  write(ErrMsg,*) 'Illegal option flag:',iOpt
  call SysAbendMsg('cxWrRun',ErrMsg,' ')
end if
iRc = 0
!----------------------------------------------------------------------*
! Call generic writing routine.                                        *
!----------------------------------------------------------------------*
call gxWrRun(iRc,Label,cData,nData,iOpt,TypStr)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine cxWrRun
