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
! Copyright (C) 1991, Per-Olof Widmark                                 *
!               1993,1996,1997, Markus P. Fuelscher                    *
!               1996, Luis Serrano-Andres                              *
!               2012, Victor P. Vysotskiy                              *
!***********************************************************************

#if defined(_I8_) || defined(_OPENMP)
#define NO_SPLITTING
#endif

subroutine DaName_Main(Lu,String,mf,wa)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Open unit Lu for direct access I/O and link the data stream to   *
!     the logical file name Name.                                      *
!                                                                      *
!     calling arguments:                                               *
!     Lu      : integer, input                                         *
!               logical unit number                                    *
!     LuName  : character string, input                                *
!               logical file name                                      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark, IBM Sweden, 1991                                   *
!     M.P. Fuelscher, University of Lund, Sweden, 1993, 1996, 1997     *
!     L. Serrano-Andres, University of Lund, Sweden, 1996              *
!----------------------------------------------------------------------*
!                                                                      *
! History:                                                             *
!          V.P. Vysotskiy, University of Lund, Sweden, 2012            *
!                                                                      *
!***********************************************************************

#ifndef _HAVE_EXTRA_
use Prgm
#endif
implicit integer(A-Z)
#include "fio.fh"
#include "switch.fh"
#ifdef _OLD_IO_STAT_
#include "ofio.fh"
#else
#include "pfio.fh"
#endif
#include "blksize.fh"
#include "filesize.fh"
character*(*) String
logical mf, wa
character*8 StdNam
character*80 Text
character*16 TheName
data TheName/'DaName_Main'/

if (Trace) then
  write(6,*) ' >>> Enter DaName_Main <<<'
  write(6,*) ' unit :',Lu
  write(6,*) ' name :',String,mf,wa
end if

tmp = Lu
Lu = isfreeunit(tmp)
! Check calling arguments
if ((Lu <= 0) .or. (Lu > MxFile)) call SysFileMsg(TheName,'MSG: unit',Lu,String)

! Check for consistency
if (isOpen(Lu) /= 0) call SysFileMsg(TheName,'MSG: used',Lu,String)

! Reformat file name to standard notation
! (capital letters, no leading blank)
call StdFmt(String,StdNam)

! If no file name has been given link it to the default name
if (StdNam == '        ') write(StdNam,'(A,I2.2,A)') 'FT',Lu,'F001'

! Check the storage scheme
#ifndef _GA_
isFiM(Lu) = 0
isFiM(Lu) = isinmem(StdNam)
#ifdef _DEBUGPRINT_IO_
if (isFiM(Lu) > 0) write(6,*) 'The file ',StdNam,' will be kept in memory'
#endif
! Open file
temp = isFiM(Lu)
#else
temp = 0
#endif
iRc = AixOpn(temp,StdNam,.true.)
#ifndef _GA_
if (iRc == eFiMFo) then
# ifdef _DEBUGPRINT_IO_
  write(6,*) 'Failed to open file in memory'
# endif
  isFiM(Lu) = 0
  iRc = 0
end if
#endif
if (iRc /= 0) then
  iRc = AixErr(Text)
  call SysFileMsg(TheName,'MSG: open',Lu,Text)
end if
isOpen(Lu) = 1
FSCB(Lu) = temp
LuName(Lu) = StdNam
#ifndef _OLD_IO_STAT_
inUse = 0
do i=1,NProfFiles
  if (LuNameProf(i) == StdNam) then
    inUse = 1
  end if
end do

if (inUse == 0) then
  if (NProfFiles+1 <= MxFile) then
    NProfFiles = NProfFiles+1
    LuNameProf(NProfFiles) = StdNam
  else
    write(6,*) 'IO error: NProfFiles+1.gt.MxFile'
    write(6,*) 'Increase MxFile in src/Include/MxFile.fh!'
    call Abend()
  end if
end if
#endif
call SetLuMark(Lu)
Addr(Lu) = 0
MPUnit(0,Lu) = Lu
#ifdef _OLD_IO_STAT_
MxAddr(Lu) = 0
count(1,Lu) = 0
count(2,Lu) = 0
count(3,Lu) = 0
count(4,Lu) = 0
#endif
Multi_File(Lu) = .false.
MBL(Lu) = MBl_nwa
if (wa) MBL(Lu) = MBl_wa
#ifndef NO_SPLITTING
if (mf) then
  Multi_File(Lu) = .true.
  MaxFileSize = AllocDisk()
  MFMB = int(dble(Max_File_Length)/(1024**2))
  if (MaxFileSize > MFMB) then
    write(6,*)
    write(6,*) 'DANAME_MF: Requested MaxFileSize is too large!'
    write(6,*) ' Requested value of ',MaxFileSize
    MaxFileSize = MFMB
    write(6,*) ' has been reset to  ',MaxFileSize
  else if (MaxFileSize /= 0) then
    if (Trace) write(6,*) ' This is a partitioned data set'
    lName = StrnLn(String)
    if ((lName == 0) .or. (lName == 8)) &
      call SysFileMsg(TheName,'Invalid file name.\n File names used in multiple unit files must be less than 8 characters.', &
                      Lu,String)
    do i=1,MaxSplitFile-1
      MPUnit(i,Lu) = -99
    end do
  end if
end if
#endif
if (Trace) then
  write(6,*) ' >>> Exit DaName_Main <<<'
end if

return

end subroutine DaName_Main
