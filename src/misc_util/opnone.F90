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

subroutine OpnOne(rc,Option,FName,Lu)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Open the one electron integral file                              *
!                                                                      *
!     input:                                                           *
!     option : Switch to set options                                   *
!     FnOne  : Logical file name                                       *
!     LuOne  : FORTRAN unit number                                     *
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

use OneDat, only: AuxOne, FInfoOne, lTocOne, NaN, nBas, nSym, pFID, pNext, pVersN, rcOne, sDmp, sNew, TocOne
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp), intent(in) :: Option, Lu
character(len=*), intent(in) :: FName
integer(kind=iwp) :: iDisk, LuOne, SumOpt
logical(kind=iwp) :: Exists, NewToc
character(len=8) :: FnOne
character(len=*), parameter :: TheName = 'OpnOne'

!----------------------------------------------------------------------*
! Start procedure:                                                     *
!----------------------------------------------------------------------*
rc = rcOne%good
!----------------------------------------------------------------------*
! Get basis sets dimensions                                            *
!----------------------------------------------------------------------*
call Get_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)
!----------------------------------------------------------------------*
! Truncate the name to 8 characters and convert it to upper case       *
!----------------------------------------------------------------------*
LuOne = Lu
FnOne = FName
call UpCase(FnOne)
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
call f_Inquire(FnOne,Exists)
NewToc = btest(Option,sNew)
!----------------------------------------------------------------------*
! Compare file status with options                                     *
!----------------------------------------------------------------------*
if ((.not. Exists) .and. (.not. NewToc)) then
  !--------------------------------------------------------------------*
  ! Old file did not exist                                             *
  !--------------------------------------------------------------------*
  call SysAbendMsg(TheName,'The ONEINT file does not exist',' ')
else if (NewToc) then
  !--------------------------------------------------------------------*
  ! New toc                                                            *
  !--------------------------------------------------------------------*
  AuxOne%Lu = NaN
  AuxOne%Opn = .false.
  TocOne(:) = NaN
  call DaName_MF(LuOne,FnOne)
  TocOne(pFID) = FInfoOne%ID
  TocOne(pVersN) = FInfoOne%VN
  iDisk = 0
  call iDaFile(LuOne,1,TocOne,lTocOne,iDisk)
  TocOne(pNext) = iDisk
  iDisk = 0
  call iDaFile(LuOne,1,TocOne,lTocOne,iDisk)
  AuxOne%Lu = LuOne
  AuxOne%Opn = .true.
else
  !--------------------------------------------------------------------*
  ! Keep toc                                                           *
  !--------------------------------------------------------------------*
  call DaName_MF(LuOne,FnOne)
  iDisk = 0
  call iDaFile(LuOne,2,TocOne,lTocOne,iDisk)
  if ((TocOne(pFID) /= FInfoOne%ID) .or. (TocOne(pVersN) /= FInfoOne%VN)) &
    call SysFileMsg(TheName,'file version number is outdated',LuOne,' ')
  AuxOne%Lu = LuOne
  AuxOne%Opn = .true.
end if
!----------------------------------------------------------------------*
! Dump the TOC upon request                                            *
!----------------------------------------------------------------------*
if (btest(Option,sDmp)) call DmpOne()

!----------------------------------------------------------------------*
! exit                                                                 *
!----------------------------------------------------------------------*
return

end subroutine OpnOne
