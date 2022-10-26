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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!               1993, Per-Olof Widmark                                 *
!***********************************************************************

subroutine OpnMCK(rc,Option,Name,Lu)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Open the one-electron integral file and initialize/load          *
!     the table of contents.                                           *
!                                                                      *
!     input:                                                           *
!     option : Switch to set options                                   *
!              = 0 old file                                            *
!              = 1 new file                                            *
!     FnCom  : Logical file name                                       *
!     LuCom  : FORTRAN unit number                                     *
!                                                                      *
!     output:                                                          *
!     rc     : Return code.                                            *
!              A value of 0 (zero) is returned upon successful         *
!              completion of the request. A nonzero value indi-        *
!              cates an error.                                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M. P. Fuelscher and P.O. Widmark                                 *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

implicit integer(A-Z)
#include "FileIDs.fh"
#include "MckDat.fh"
character*(*) Name
character*8 FnMCK
logical exist, NewToc
character*16 TheName
data TheName/'OpnMck'/

!---------------------------------------------------------------------*
! Start procedure:                                                    *
!---------------------------------------------------------------------*
NewToc = iand(option,sNew) /= 0
rc = rc0000
!---------------------------------------------------------------------*
! Start procedure:                                                    *
!---------------------------------------------------------------------*
AuxMCK(pLu) = 0
AuxMCK(pOpen) = 0
call StdFmt(Name,FnMCK)
LuMCK = Lu
call f_Inquire(FnMCK,Exist)
!----------------------------------------------------------------------*
! Check the options                                                    *
!----------------------------------------------------------------------*
if (Option /= 0) then
  SumOpt = 0
  if (iand(Option,sNew) /= 0) SumOpt = SumOpt+sNew
  if (iand(Option,1024) /= 0) SumOpt = SumOpt+1024
  if (SumOpt /= Option) then
    call SysWarnMsg(TheName,'MSG: invalid option',' ')
    call SysCondMsg('SumOpt /= Option',SumOpt,'/=',Option)
  end if
end if
!----------------------------------------------------------------------*
! Compare file status with options                                     *
!----------------------------------------------------------------------*
if ((.not. NewToc) .and. (.not. exist)) then
  !--------------------------------------------------------------------*
  ! Old file did not exist                                             *
  !--------------------------------------------------------------------*
  call SysAbendMsg(TheName,'MCK file does not exist',' ')
else if (NewToc) then
  !--------------------------------------------------------------------*
  ! New toc                                                            *
  !--------------------------------------------------------------------*
  call DaName(LuMCK,FnMCK)
  call iCopy(lToc,[NaN],0,TocOne,1)
  TocOne(pFID) = IDone
  TocOne(pVersN) = VNone
  iDisk = 0
  call iDaFile(LuMCK,1,TocOne,lToc,iDisk)
  TocOne(pNext) = iDisk
  iDisk = 0
  call iDaFile(LuMCK,1,TocOne,lToc,iDisk)
  AuxMCK(pLu) = LuMCK
  AuxMCK(pOpen) = 1
else
  !--------------------------------------------------------------------*
  ! Old toc                                                            *
  !--------------------------------------------------------------------*
  call DaName(LuMCK,FnMCK)
  iDisk = 0
  call iDaFile(LuMCK,2,TocOne,lToc,iDisk)
  if ((TocOne(pFID) /= IDone) .or. (TocOne(pVersN) /= VNone)) then
    call SysFileMsg(TheName,'file version number is outdated',LuMCK,' ')
  end if
  AuxMCK(pLu) = LuMCK
  AuxMCK(pOpen) = 1
end if
! Report back the actual unit number
Lu = LuMCK

!----------------------------------------------------------------------*
! normal end                                                           *
!----------------------------------------------------------------------*
return

end subroutine OpnMCK
