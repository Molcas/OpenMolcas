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
! Data type is Integer.                                                *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund university, Sweden                                     *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************

subroutine ixRdRun(iRc,Label,iData,nData,iOpt)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use RunFile_data, only: TypInt
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: iRc
character(len=*), intent(in) :: Label
integer(kind=iwp), intent(_OUT_) :: iData(*)
integer(kind=iwp), intent(in) :: nData, iOpt
character(len=64) :: ErrMsg

call ixRdRun_Internal(iData)

! This is to allow type punning without an explicit interface
contains

subroutine ixRdRun_Internal(iData)

  integer(kind=iwp), target, intent(_OUT_) :: iData(*)
  character, pointer :: cData(:)

  !--------------------------------------------------------------------*
  ! Check that arguments are ok.                                       *
  !--------------------------------------------------------------------*
  if (iOpt /= 0) then
    write(ErrMsg,*) 'Illegal option flag:',iOpt
    call SysAbendMsg('ixRdRun',ErrMsg,' ')
  end if
  iRc = 0
  !--------------------------------------------------------------------*
  ! Call generic reading routine.                                      *
  !--------------------------------------------------------------------*
  call c_f_pointer(c_loc(iData(1)),cData,[nData])
  call gxRdRun(iRc,Label,cData,nData,iOpt,TypInt)
  nullify(cData)
  !--------------------------------------------------------------------*
  !                                                                    *
  !--------------------------------------------------------------------*
  return

end subroutine ixRdRun_internal

end subroutine ixRdRun
