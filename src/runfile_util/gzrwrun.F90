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

subroutine gzRWRun(Lu,icXX,data,nData,iDisk,RecTyp)

#include "runinfo.fh"
#include "runtypes.fh"
!----------------------------------------------------------------------*
! Declare arguments                                                    *
!----------------------------------------------------------------------*
integer Lu
integer icXX
character data(*)
integer nData
integer iDisk
integer RecTyp

!----------------------------------------------------------------------*
! Read/write data from/to runfile.                                     *
!----------------------------------------------------------------------*
select case (RecTyp)
  case (TypInt)
    call c_iDaFile(Lu,icXX,data,nData,iDisk)
  case (TypDbl)
    call c_dDaFile(Lu,icXX,data,nData,iDisk)
  case (TypStr)
    call cDaFile(Lu,icXX,data,nData,iDisk)
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

  use iso_c_binding

  integer Lu, iOpt, lBuf_, iDisk_
  character, target :: Buf(*)
  integer, pointer :: pBuf(:)

  call c_f_pointer(c_loc(Buf(1)),pBuf,[lBuf_])
  call iDaFile(Lu,iOpt,pBuf,lBuf_,iDisk_)
  nullify(pBuf)

end subroutine c_iDaFile

subroutine c_dDaFile(Lu,iOpt,Buf,lBuf_,iDisk_)

  use iso_c_binding

  integer Lu, iOpt, lBuf_, iDisk_
  character, target :: Buf(*)
  real*8, pointer :: pBuf(:)

  call c_f_pointer(c_loc(Buf(1)),pBuf,[lBuf_])
  call dDaFile(Lu,iOpt,pBuf,lBuf_,iDisk_)
  nullify(pBuf)

end subroutine c_dDaFile

end subroutine gzRWRun
