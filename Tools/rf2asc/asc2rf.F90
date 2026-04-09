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
! This program converts a runfile into an ascii file.                  *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
! Written: July 2003                                                   *
!          Lund University, Sweden                                     *
!                                                                      *
!***********************************************************************

program Asc2RF

use RunFile_data, only: lw, TypDbl, TypInt, TypLgl, TypStr
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, iOpt, iRc, istatus, Length, RUNASCII, typ
logical(kind=iwp) :: Found
character(len=lw+2) :: Record
character(len=lw) :: Label
integer(kind=iwp), allocatable :: iBuf(:)
real(kind=wp), allocatable :: dBuf(:)
character, allocatable :: cBuf(:)
integer(kind=iwp), external :: isFreeUnit

call IniMem()
call Init_LinAlg()
call PrgmInit('Asc2RF')
!----------------------------------------------------------------------*
! Open input file.                                                     *
!----------------------------------------------------------------------*
call f_inquire('RUNASCII',Found)
if (.not. Found) then
  call WarningMessage(2,'RUNASCII not found')
  stop
end if
RUNASCII = isFreeUnit(15)
call MOLCAS_OPEN(RUNASCII,'RUNASCII')
!----------------------------------------------------------------------*
! Create RunFile.                                                      *
!----------------------------------------------------------------------*
call NameRun('RUNFILE')
iOpt = 0
iRc = 0
call MkRun(iRc,iOpt)
!----------------------------------------------------------------------*
! Populate runfile.                                                    *
!----------------------------------------------------------------------*
do
  read(RUNASCII,'(a)',iostat=istatus) Record
  if (istatus /= 0) then
    if (istatus > 0) call WarningMessage(2,'An error occurred')
    exit
  end if
  if (Record(1:1) == '#') cycle
  if ((Record(1:1) /= '<') .or. (Record(lw+2:lw+2) /= '>')) then
    call WarningMessage(2,'An error occurred')
    exit
  end if
  Label = Record(2:lw+1)
  read(RUNASCII,*) Length,typ
  write(u6,*) 'Processing:',Record,Length,typ

  select case (typ)
    case (TypDbl)
      call mma_allocate(dBuf,Length,Label='dBuf')
      read(RUNASCII,*) dBuf(:)
      call dxWrRun(iRc,Label,dBuf,Length,iOpt)
      call mma_deallocate(dBuf)
    case (TypInt)
      call mma_allocate(iBuf,Length,Label='iBuf')
      read(RUNASCII,*) iBuf(:)
      call ixWrRun(iRc,Label,iBuf,Length,iOpt)
      call mma_deallocate(iBuf)
    case (TypStr)
      call mma_allocate(cBuf,Length,Label='cBuf')
      read(RUNASCII,'(64a1)') cBuf(:)
      call cxWrRun(iRc,Label,cBuf,Length,iOpt)
      call mma_deallocate(cBuf)
    case (TypLgl)
      write(u6,*) 'Cannot handle type logical'
      exit
    case default
      write(u6,*) 'Unknown data type:',typ
      exit
  end select
end do
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
close(RUNASCII)

call GetMem('Finish','LIST','REAL',I,0)
call GetMem('Finish','TERM','REAL',I,0)

end program Asc2RF
