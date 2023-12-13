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

subroutine OpnOrd(rc,Option,FName,Lu)
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
!     FName  : Logical file name                                       *
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

use TwoDat, only: AuxTwo, FInfoTwo, iNoNum, lTocTwo, isId, isForm, isPkPa, isVer, rcTwo, sNew, TocTwo
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp), intent(in) :: Option, Lu
character(len=*), intent(in) :: FName
integer(kind=iwp) :: iDisk, iDummy, LuTwo, nDummy1(8), nDummy2(8), rd_Dummy, SumOpt
logical(kind=iwp) :: Exists, lDummy, NewToc
character(len=8) :: FnTwo
integer(kind=iwp), parameter :: NaN = -1
character(len=*), parameter :: TheName = 'OpnOrd'

!----------------------------------------------------------------------*
! Start procedure:                                                     *
!----------------------------------------------------------------------*
rc = rcTwo%good
!----------------------------------------------------------------------*
! Check the file status                                                *
!----------------------------------------------------------------------*
AuxTwo%Opn = .false.
AuxTwo%Unt = iNoNum
AuxTwo%DaDa = iNoNum
TocTwo(isPkPa) = iNoNum
call StdFmt(FName,FnTwo)
LuTwo = Lu
call f_Inquire(FnTwo,Exists)
!----------------------------------------------------------------------*
! Check the options                                                    *
!----------------------------------------------------------------------*
if (Option /= 0) then
  SumOpt = 0
  if (btest(Option,sNew)) SumOpt = ibset(SumOpt,sNew)
  if (SumOpt /= Option) then
    call SysWarnMsg(TheName,'MSG: invalid option',' ')
    call SysCondMsg('SumOpt /= Option',SumOpt,'/=',Option)
  end if
end if
NewToc = btest(Option,sNew)
!----------------------------------------------------------------------*
! Compare file status with options                                     *
!----------------------------------------------------------------------*
if ((.not. Exists) .and. (.not. NewToc)) then
  !--------------------------------------------------------------------*
  ! Old file did not exist                                             *
  !--------------------------------------------------------------------*
  call SysAbendMsg(TheName,'ORDINT file does not exist',' ')
else if (NewToc) then
  !--------------------------------------------------------------------*
  ! New Toc                                                            *
  !--------------------------------------------------------------------*
  call DaName_MF(LuTwo,FnTwo)
  TocTwo(:) = NaN
  TocTwo(isId) = FInfoTwo%ID
  TocTwo(isVer) = FInfoTwo%VN
  TocTwo(isForm) = 0
  iDisk = 0
  call iDaFile(LuTwo,1,TocTwo,lTocTwo,iDisk)
  AuxTwo%Opn = .true.
  AuxTwo%Unt = LuTwo
  AuxTwo%DaDa = 0
else
  !--------------------------------------------------------------------*
  ! Keep Toc                                                           *
  !--------------------------------------------------------------------*
  call DaName_MF(LuTwo,FnTwo)
  iDisk = 0
  call iDaFile(LuTwo,2,TocTwo,lTocTwo,iDisk)
  if ((TocTwo(isId) /= FInfoTwo%ID) .or. (TocTwo(isVer) /= FInfoTwo%VN)) &
    call SysFileMsg(TheName,'file version number is outdated',LuTwo,' ')
  AuxTwo%Opn = .true.
  AuxTwo%Unt = LuTwo
  AuxTwo%DaDa = iDisk
end if

! Call to GetOrd to fill nBatch etc.

if (Option == 0) call GetOrd(rd_Dummy,lDummy,iDummy,nDummy1,nDummy2)

!----------------------------------------------------------------------*
! normal end                                                           *
!----------------------------------------------------------------------*
return

end subroutine OpnOrd
