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

subroutine OpnOne(rc,Option,Name,Lu)
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

implicit integer(A-Z)
#include "FileIDs.fh"
#include "OneDat.fh"
character*(*) Name
character*8 FnOne
logical Exist, NewToc
character*16 TheName
data TheName/'OpnOne'/

!----------------------------------------------------------------------*
! Start procedure:                                                     *
!----------------------------------------------------------------------*
rc = rc0000
!----------------------------------------------------------------------*
! Get basis sets dimensions                                            *
!----------------------------------------------------------------------*
call Get_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)
!----------------------------------------------------------------------*
! Truncate the name to 8 characters and convert it to upper case       *
!----------------------------------------------------------------------*
LuOne = Lu
FnOne = Name
call UpCase(FnOne)
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
call f_Inquire(FnOne,Exist)
NewToc = iand(Option,sNew) /= 0
!----------------------------------------------------------------------*
! Compare file status with options                                     *
!----------------------------------------------------------------------*
if ((.not. Exist) .and. (.not. NewToc)) then
  !--------------------------------------------------------------------*
  ! Old file did not exist                                             *
  !--------------------------------------------------------------------*
  call SysAbendMsg(TheName,'The ONEINT file does not exist',' ')
else if (NewToc) then
  !--------------------------------------------------------------------*
  ! New toc                                                            *
  !--------------------------------------------------------------------*
  call iCopy(lAux,[NaN],0,AuxOne,1)
  call iCopy(lToc,[NaN],0,TocOne,1)
  call DaName_MF(LuOne,FnOne)
  TocOne(pFID) = IDrlx
  TocOne(pVersN) = VNrlx
  iDisk = 0
  call iDaFile(LuOne,sWrite,TocOne,lToc,iDisk)
  TocOne(pNext) = iDisk
  iDisk = 0
  call iDaFile(LuOne,sWrite,TocOne,lToc,iDisk)
  AuxOne(pLu) = LuOne
  AuxOne(pOpen) = 1
else
  !--------------------------------------------------------------------*
  ! Keep toc                                                           *
  !--------------------------------------------------------------------*
  call DaName_MF(LuOne,FnOne)
  iDisk = 0
  call iDaFile(LuOne,sRead,TocOne,lToc,iDisk)
  if ((TocOne(pFID) /= IDrlx) .or. (TocOne(pVersN) /= VNrlx)) then
    call SysFileMsg(TheName,'file version number is outdated',LuOne,' ')
  end if
  AuxOne(pLu) = LuOne
  AuxOne(pOpen) = 1
end if
!----------------------------------------------------------------------*
! Dump the TOC upon request                                            *
!----------------------------------------------------------------------*
if (iand(Option,1024) /= 0) call DmpOne()

!----------------------------------------------------------------------*
! exit                                                                 *
!----------------------------------------------------------------------*
return

end subroutine OpnOne
