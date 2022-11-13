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

subroutine ixWrRun(iRc,Label,iData,nData,iOpt)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use RunFile_data, only: TypInt
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: iRc
character(len=*), intent(in) :: Label
integer(kind=iwp), intent(in) :: iData(*), nData, iOpt
character(len=64) :: Errmsg

call ixWrRun_Internal(iData)

! This is to allow type punning without an explicit interface
contains

subroutine ixWrRun_Internal(iData)

  integer(kind=iwp), target, intent(in) :: iData(*)
  character, pointer :: cData(:)

  !--------------------------------------------------------------------*
  ! Check that arguments are ok.                                       *
  !--------------------------------------------------------------------*
  if (iOpt /= 0) then
    write(ErrMsg,*) 'Illegal option flag:',iOpt
    call SysAbendMsg('ixWrRun',ErrMsg,' ')
  end if
  iRc = 0
  !--------------------------------------------------------------------*
  ! Call generic writing routine.                                      *
  !--------------------------------------------------------------------*
  call c_f_pointer(c_loc(iData(1)),cData,[1])
  call gxWrRun(iRc,Label,cData,nData,iOpt,TypInt)
  nullify(cData)
  !--------------------------------------------------------------------*
  !                                                                    *
  !--------------------------------------------------------------------*
  return

end subroutine ixWrRun_Internal

end subroutine ixWrRun
