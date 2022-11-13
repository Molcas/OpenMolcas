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
! Data type is Real*8.                                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund university, Sweden                                     *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************

subroutine dWrRun(Label,rData,nData)

use Definitions, only: wp, iwp

implicit none
character(len=*), intent(in) :: Label
real(kind=wp), intent(in) :: rData(*)
integer(kind=iwp), intent(in) :: nData
integer(kind=iwp) :: iOpt, iRc
character(len=64) :: ErrMsg

!----------------------------------------------------------------------*
! Call extended writing routine.                                       *
!----------------------------------------------------------------------*
iRc = 0
iOpt = 0
call dxWrRun(iRc,Label,rData,nData,iOpt)
if (iRc /= 0) then
  write(ErrMsg,'(3a)') 'Error writing field "',Label,'" into runfile'
  call SysAbendMsg('dWrRun',ErrMsg,' ')
end if
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine dWrRun
