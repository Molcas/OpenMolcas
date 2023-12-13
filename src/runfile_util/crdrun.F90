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
! This routine reads a record from the runfile.                        *
! Data type is Real*8.                                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund university, Sweden                                     *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************

subroutine cRdRun(Label,cData,nData)

use Definitions, only: iwp

#include "intent.fh"

implicit none
character(len=*), intent(in) :: Label
character, intent(_OUT_) :: cData(*)
integer(kind=iwp), intent(in) :: nData
integer(kind=iwp) :: iOpt, iRc
character(len=64) :: ErrMsg

!----------------------------------------------------------------------*
! Call extended reading routine.                                       *
!----------------------------------------------------------------------*
iRc = 0
iOpt = 0
call cxRdRun(iRc,Label,cData,nData,iOpt)
if (iRc /= 0) then
  write(ErrMsg,'(3a)') 'Error reading field "',Label,'" from runfile'
  call SysAbendMsg('cRdRun',ErrMsg,' ')
end if
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine cRdRun
