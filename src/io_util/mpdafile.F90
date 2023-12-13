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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!               1996, Luis Serrano-Andres                              *
!               2012, Victor P. Vysotskiy                              *
!               2016, Steven Vancoillie                                *
!***********************************************************************

subroutine MpDaFile(Lu,MaxFileSizel,iOpt,Buf,lBuf,iDisk)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     File systems or disks may be to small to keep a given large      *
!     dataset. Therefore, it may be split into subsets which then      *
!     may be linked to different physical units. This process is       *
!     controlled by the environment variable MOLCASDISK                *
!     (c.f. allocdisk). This routine computes the proper splitting     *
!     and calls the IO routines.                                       *
!                                                                      *
!     calling arguments:                                               *
!     Lu      : integer, input                                         *
!               logical unit number (Lu={1,2,...40,50,60,70,80,90}     *
!     MaxFileSize : integer, input                                     *
!               max length of a file.                                  *
!     iOpt    : integer, input                                         *
!               option code                                            *
!               iOpt= 0 dummy write                                    *
!               iOpt= 1 synchronous write                              *
!               iOpt= 2 synchronous read                               *
!               iOpt= 3 synchronous gather and write                   *
!               iOpt= 4 synchronous read and scatter                   *
!               iOpt= 5 synchronous rewind                             *
!               iOpt= 6 asynchronous write                             *
!               iOpt= 7 asynchronous read                              *
!               iOpt= 8 asynchronous gather and write                  *
!               iOpt= 9 asynchronous read and scatter                  *
!               iOpt=10 asynchronous rewind                            *
!               Note: At present the asynchronous modes are not        *
!                     supported and work identically the synchronous   *
!                     modes                                            *
!     Buf     : array of characters, input/output                      *
!               Buffer carrying/receiving the data to write/read       *
!     lBuf    : integer, input                                         *
!               length of the buffer Buf in bytes!!                    *
!     iDisk   : integer, input                                         *
!               disk address as byte offset!!!                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher, University of Lund, Sweden, 1996                 *
!     L. Serrano-Andres, University of Lund, Sweden, 1996              *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! History:                                                             *
!     V.P. Vysotskiy, University of Lund, Sweden, 2012                 *
!     Steven Vancoillie: use of byte lengths/offests (2016)            *
!                                                                      *
!***********************************************************************

use Fast_IO, only: Addr, FSCB, isOpen, LuName, Max_File_Length, MaxSplitFile, MBL, MPUnit, Multi_File
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: Lu, MaxFileSizel, iOpt, lBuf, iDisk
character, intent(inout) :: Buf(*)
integer(kind=iwp) :: iDisk_mod, iext, iRc, lStdNam, ltmp, Lu_mod, max_Bytes, max_File_Size, n_Bytes, n_Copy, n_Save, offset, &
                     p_Buf, start_char, temp
character(len=256) :: tmp
character(len=80) :: Text
character(len=8) :: Stdnam, ext
character(len=*), parameter :: TheName = 'MpDaFile'
integer(kind=iwp), external :: AixOpn, isFreeUnit, StrnLn
interface
  function AixErr(FileName) bind(C,name='aixerr_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: AixErr
    character(kind=c_char) :: FileName(*)
  end function AixErr
end interface

max_File_Size = MaxFileSizel*10**6
max_Bytes = min(max_File_Length,max_File_Size)
n_Bytes = lBuf
offset = iDisk/max_Bytes
! set correct char to append: 0-9, A-Z
start_char = 48

if (offset <= 9) then
  start_char = 48
else
  start_char = 55
end if

if ((offset < 0) .or. (offset > MaxSplitFile-1)) then
  StdNam = LuName(Lu)
  call PrgmTranslate(StdNam,tmp,ltmp)
  write(u6,*) '          Current I/O Status as follows'
  write(u6,*)
  call FASTIO('STATUS')
  call SysWarnFileMsg(TheName,StdNam,'Extensions out of range!','increase MOLCAS_DISK value or MaxSplitFile in Fast_IO')
  call Abend()
end if
Lu_mod = MPUnit(offset,LU)
iDisk_mod = iDisk-offset*max_Bytes

! If extension not opened before: open and initiate here!

StdNam = LuName(Lu)
call PrgmTranslate(StdNam,tmp,ltmp)
lStdNam = ltmp
if (Lu_Mod < 0) then
  Lu_Mod = isfreeunit(Lu)
  MPUnit(offset,LU) = Lu_Mod
  temp = 0
  tmp(1+lStdNam:1+lStdNam) = char(start_char+offset)
  iext = strnln(StdNam)
  ext = StdNam
  if (offset <= 9) then
    ext(1+iext:1+iext) = char(start_char+offset)
  else
    ext(1+iext:1+iext) = char(start_char+offset/10)
    ext(2+iext:2+iext) = char(start_char+offset-offset/10*10)
  end if

  iRc = AixOpn(temp,tmp,.false.)
  if (iRc /= 0) then
    iRc = AixErr(Text)
    call SysFileMsg(TheName,'MSG: open',Lu_Mod,Text)
  end if
  isOpen(Lu_Mod) = 1
  FSCB(Lu_Mod) = temp

  !vv LuName(Lu_Mod) = tmp
  LuName(Lu_Mod) = ext
  Addr(Lu_Mod) = 0
  Multi_File(Lu_Mod) = .true.
  MPUnit(0,LU_Mod) = Lu
  MBL(Lu_Mod) = MBL(Lu)

end if

if ((iDisk_mod+n_Bytes) <= max_Bytes) then

  ! The whole record in contained on a single file.

  p_Buf = 1
  n_Copy = lBuf
  call CHDAFILE(LU_mod,iOpt,Buf(p_Buf),n_Copy,iDisk_mod)
else

  ! The record must be split over several files.

  p_Buf = 1
  n_Save = lBuf
  n_Copy = (max_Bytes-iDisk_mod)
  do while (n_Save > 0)
    if (LU_MOD < 0) then
      Lu_Mod = isfreeunit(Lu)
      MPUnit(offset,LU) = Lu_Mod
      temp = 0
      start_char = 48
      if (offset <= 9) then
        start_char = 48
      else
        start_char = 55
      end if

      tmp(1+lStdNam:1+lStdNam) = char(start_char+offset)
      iext = strnln(StdNam)
      ext = StdNam
      if (offset <= 9) then
        ext(1+iext:1+iext) = char(start_char+offset)
      else
        ext(1+iext:1+iext) = char(start_char+offset/10)
        ext(2+iext:2+iext) = char(start_char+offset-offset/10*10)
      end if

      iRc = AixOpn(temp,tmp,.false.)
      if (iRc /= 0) then
        iRc = AixErr(Text)
        call SysFileMsg(TheName,'MSG: open',Lu_Mod,Text)
      end if
      isOpen(Lu_mod) = 1
      FSCB(Lu_mod) = temp

      !vv LuName(Lu_Mod) = tmp
      LuName(Lu_Mod) = ext

      Addr(Lu_Mod) = 0
      Multi_File(Lu_Mod) = .true.
      MPUnit(0,LU_Mod) = Lu
      MBL(Lu_Mod) = MBL(Lu)

    end if
    call CHDAFILE(LU_mod,iOpt,Buf(p_Buf),n_Copy,iDisk_mod)
    p_Buf = p_Buf+n_Copy
    n_Save = n_Save-n_Copy
    n_Copy = min(max_Bytes,n_Save)
    offset = offset+1
    if (offset > MaxSplitFile-1) then
      write(u6,*) '          Current I/O Status as follows'
      write(u6,*)
      call FASTIO('STATUS')
      call SysWarnFileMsg(TheName,StdNam,'Extensions out of range!','increase MOLCAS_DISK value or MaxSplitFile in Fast_IO')
      call Abend()
    end if

    Lu_Mod = MPUnit(offset,LU)
    iDisk_mod = 0
  end do
end if

return

end subroutine MpDaFile
