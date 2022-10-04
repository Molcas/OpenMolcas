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
! Data type is Integer.                                                *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund university, Sweden                                     *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************

subroutine iWrRun(Label,data,nData)

!----------------------------------------------------------------------*
! Declare arguments                                                    *
!----------------------------------------------------------------------*
character*(*) Label
integer data(*)
integer nData
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
character*64 ErrMsg
integer iRc
integer iOpt

!----------------------------------------------------------------------*
! Call extended writing routine.                                       *
!----------------------------------------------------------------------*
iRc = 0
iOpt = 0
call ixWrRun(iRc,Label,data,nData,iOpt)
if (iRc /= 0) then
  write(ErrMsg,'(3a)') 'Error writing field "',Label,'" into runfile'
  call SysAbendMsg('iWrRun',ErrMsg,' ')
end if
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine iWrRun
