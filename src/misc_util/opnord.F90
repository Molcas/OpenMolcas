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
!***********************************************************************

subroutine OpnOrd(rc,Option,Name,Lu)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Open the two-electron integral file and initialize/load          *
!     the table of contents.                                           *
!                                                                      *
!     input:                                                           *
!     option : Switch to set options                                   *
!              = 0 old file                                            *
!              = 1 new file                                            *
!     Name   : Logical file name                                       *
!     Lu     : FORTRAN unit number                                     *
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
!     M. P. Fuelscher                                                  *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

implicit integer(A-Z)
#include "FileIDs.fh"
#include "Molcas.fh"
#include "TwoDat.fh"
parameter(NaN=-1)
character*(*) Name
character*8 FnTwo
logical Exist, NewToc
logical lDummy
integer nDummy1(8), nDummy2(8)
character*16 TheName
data TheName/'OpnOrd'/

!----------------------------------------------------------------------*
! Start procedure:                                                     *
!----------------------------------------------------------------------*
rc = rc0000
!----------------------------------------------------------------------*
! Check the file status                                                *
!----------------------------------------------------------------------*
AuxTwo(isUnit) = iNoNum
AuxTwo(isStat) = iNoNum
AuxTwo(isDaDa) = iNoNum
TocTwo(isPkPa) = iNoNum
TocTwo(isPkAs) = iNoNum
call StdFmt(Name,FnTwo)
LuTwo = Lu
call f_Inquire(FnTwo,Exist)
!----------------------------------------------------------------------*
! Check the options                                                    *
!----------------------------------------------------------------------*
if (Option /= 0) then
  SumOpt = 0
  if (iand(Option,sNew) /= 0) SumOpt = SumOpt+sNew
  if (SumOpt /= Option) then
    call SysWarnMsg(TheName,'MSG: invalid option',' ')
    call SysCondMsg('SumOpt /= Option',SumOpt,'/=',Option)
  end if
end if
NewToc = iand(sNew,Option) /= 0
!----------------------------------------------------------------------*
! Compare file status with options                                     *
!----------------------------------------------------------------------*
if ((.not. Exist) .and. (.not. NewToc)) then
  !--------------------------------------------------------------------*
  ! Old file did not exist                                             *
  !--------------------------------------------------------------------*
  call SysAbendMsg(TheName,'ORDINT file does not exist',' ')
else if (NewToc) then
  !--------------------------------------------------------------------*
  ! New Toc                                                            *
  !--------------------------------------------------------------------*
  call DaName_MF(LuTwo,FnTwo)
  call iCopy(lTocTwo,[NaN],0,TocTwo,1)
  TocTwo(isId) = IDtwo
  TocTwo(isVer) = VNtwo
  TocTwo(isForm) = 0
  iDisk = 0
  call iDaFile(LuTwo,1,TocTwo,lTocTwo,iDisk)
  AuxTwo(isUnit) = LuTwo
  AuxTwo(isStat) = 1
  AuxTwo(isDaDa) = 0
else
  !--------------------------------------------------------------------*
  ! Keep Toc                                                           *
  !--------------------------------------------------------------------*
  call DaName_MF(LuTwo,FnTwo)
  iDisk = 0
  call iDaFile(LuTwo,2,TocTwo,lTocTwo,iDisk)
  if ((TocTwo(isId) /= IDtwo) .or. (TocTwo(isVer) /= VNtwo)) then
    call SysFileMsg(TheName,'file version number is outdated',LuTwo,' ')
  end if
  AuxTwo(isUnit) = LuTwo
  AuxTwo(isStat) = 1
  AuxTwo(isDaDa) = iDisk
end if

! Call to GetOrd to fill nBatch etc.

if (Option == 0) then
  call GetOrd(rd_Dummy,lDummy,iDummy,nDummy1,nDummy2)
end if

!----------------------------------------------------------------------*
! normal end                                                           *
!----------------------------------------------------------------------*
return

end subroutine OpnOrd
