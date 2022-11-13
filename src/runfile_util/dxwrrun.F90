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

subroutine dxWrRun(iRc,Label,rData,nData,iOpt)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use RunFile_data, only: TypDbl
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: iRc
character(len=*), intent(in) :: Label
real(kind=wp), intent(in) :: rData(*)
integer(kind=iwp), intent(in) :: nData, iOpt
character(len=64) :: ErrMsg

call dxWrRun_Internal(rData)

! This is to allow type punning without an explicit interface
contains

subroutine dxWrRun_Internal(rData)

  real(kind=wp), target, intent(in) :: rData(*)
  character, pointer :: cData(:)

  !--------------------------------------------------------------------*
  ! Check that arguments are ok.                                       *
  !--------------------------------------------------------------------*
  if (iOpt /= 0) then
    write(ErrMsg,*) 'Illegal option flag:',iOpt
    call SysAbendMsg('dxWrRun',ErrMsg,' ')
  end if
  iRc = 0
  !--------------------------------------------------------------------*
  ! Call generic writing routine.                                      *
  !--------------------------------------------------------------------*
  call c_f_pointer(c_loc(rData(1)),cData,[nData])
  call gxWrRun(iRc,Label,cData,nData,iOpt,TypDbl)
  nullify(cData)
  !--------------------------------------------------------------------*
  !                                                                    *
  !--------------------------------------------------------------------*
  return

end subroutine dxWrRun_Internal

end subroutine dxWrRun
