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

subroutine OpnMCK(rc,Option,FName,Lu)
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

use MckDat, only: AuxMck, FInfoMck, lTocMck, NaN, pFID, pNext, pVersN, rcMck, sDmp, sNew, TocMck
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp), intent(in) :: Option
character(len=*), intent(in) :: FName
integer(kind=iwp), intent(inout) :: Lu
integer(kind=iwp) :: iDisk, LuMCK, SumOpt
logical(kind=iwp) :: Exists, NewToc
character(len=8) :: FnMCK
character(len=*), parameter :: TheName = 'OpnMck'

!---------------------------------------------------------------------*
! Start procedure:                                                    *
!---------------------------------------------------------------------*
NewToc = btest(option,sNew)
rc = rcMck%good
!---------------------------------------------------------------------*
! Start procedure:                                                    *
!---------------------------------------------------------------------*
AuxMck%Lu = 0
AuxMck%Opn = .false.
call StdFmt(FName,FnMCK)
LuMCK = Lu
call f_Inquire(FnMCK,Exists)
!----------------------------------------------------------------------*
! Check the options                                                    *
!----------------------------------------------------------------------*
if (Option /= 0) then
  SumOpt = 0
  if (btest(Option,sNew)) SumOpt = ibset(SumOpt,sNew)
  if (btest(Option,sDmp)) SumOpt = ibset(SumOpt,sDmp)
  if (SumOpt /= Option) then
    call SysWarnMsg(TheName,'MSG: invalid option',' ')
    call SysCondMsg('SumOpt /= Option',SumOpt,'/=',Option)
  end if
end if
!----------------------------------------------------------------------*
! Compare file status with options                                     *
!----------------------------------------------------------------------*
if ((.not. NewToc) .and. (.not. Exists)) then
  !--------------------------------------------------------------------*
  ! Old file did not exist                                             *
  !--------------------------------------------------------------------*
  call SysAbendMsg(TheName,'MCK file does not exist',' ')
else if (NewToc) then
  !--------------------------------------------------------------------*
  ! New toc                                                            *
  !--------------------------------------------------------------------*
  call DaName(LuMCK,FnMCK)
  TocMck(:) = NaN
  TocMck(pFID) = FInfoMck%ID
  TocMck(pVersN) = FInfoMck%VN
  iDisk = 0
  call iDaFile(LuMCK,1,TocMck,lTocMck,iDisk)
  TocMck(pNext) = iDisk
  iDisk = 0
  call iDaFile(LuMCK,1,TocMck,lTocMck,iDisk)
  AuxMck%Lu = LuMCK
  AuxMck%Opn = .true.
else
  !--------------------------------------------------------------------*
  ! Old toc                                                            *
  !--------------------------------------------------------------------*
  call DaName(LuMCK,FnMCK)
  iDisk = 0
  call iDaFile(LuMCK,2,TocMck,lTocMck,iDisk)
  if ((TocMck(pFID) /= FinfoMck%ID) .or. (TocMck(pVersN) /= FinfoMck%VN)) &
    call SysFileMsg(TheName,'file version number is outdated',LuMCK,' ')
  AuxMck%Lu = LuMCK
  AuxMck%Opn = .true.
end if
! Report back the actual unit number
Lu = LuMCK

!----------------------------------------------------------------------*
! normal end                                                           *
!----------------------------------------------------------------------*
return

end subroutine OpnMCK
