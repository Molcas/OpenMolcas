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
! Copyright (C) 2009, Per-Olof Widmark                                 *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine reads or writes a record from/to the runfile.           *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund university, Sweden                                     *
! Written: December 2009                                               *
!                                                                      *
!***********************************************************************

subroutine gzRWRun(Lu,icXX,cData,nData,iDisk,RecTyp)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use RunFile_data, only: TypDbl, TypInt, TypLgl, TypStr
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: Lu, icXX, nData, RecTyp
integer(kind=iwp), intent(inout) :: iDisk
!TODO: No "intent" possible, since it depends on icXX
character :: cData(*)

!----------------------------------------------------------------------*
! Read/write data from/to runfile.                                     *
!----------------------------------------------------------------------*
select case (RecTyp)
  case (TypInt)
    call c_iDaFile(Lu,icXX,cData,nData,iDisk)
  case (TypDbl)
    call c_dDaFile(Lu,icXX,cData,nData,iDisk)
  case (TypStr)
    call cDaFile(Lu,icXX,cData,nData,iDisk)
  case (TypLgl)
    call SysAbendMsg('gzRWRun','Records of logical type not implemented','Aborting')
  case default
    call SysAbendMsg('gzRWRun','Argument RecTyp is of wrong type','Aborting')
end select
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

! This is to allow type punning without an explicit interface
contains

subroutine c_iDaFile(Lu,iOpt,Buf,lBuf_,iDisk_)

  integer(kind=iwp), intent(in) :: Lu, iOpt, lBuf_
  character, target :: Buf(*)
  integer(kind=iwp), intent(inout) :: iDisk_
  integer(kind=iwp), pointer :: pBuf(:)

  call c_f_pointer(c_loc(Buf(1)),pBuf,[lBuf_])
  call iDaFile(Lu,iOpt,pBuf,lBuf_,iDisk_)
  nullify(pBuf)

end subroutine c_iDaFile

subroutine c_dDaFile(Lu,iOpt,Buf,lBuf_,iDisk_)

  integer(kind=iwp), intent(in) :: Lu, iOpt, lBuf_
  character, target :: Buf(*)
  integer(kind=iwp), intent(inout) :: iDisk_
  real(kind=wp), pointer :: pBuf(:)

  call c_f_pointer(c_loc(Buf(1)),pBuf,[lBuf_])
  call dDaFile(Lu,iOpt,pBuf,lBuf_,iDisk_)
  nullify(pBuf)

end subroutine c_dDaFile

end subroutine gzRWRun
