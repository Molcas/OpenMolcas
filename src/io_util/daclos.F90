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

subroutine DaClos(Lu)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Close unit Lu                                                    *
!                                                                      *
!     calling arguments:                                               *
!     Lu      : integer, input                                         *
!               logical unit number                                    *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark, IBM Sweden, 1991                                   *
!     M.P. Fuelscher, University of Lund, Sweden, 1993, 1996, 1997     *
!     L. Serrano-Andres, University of Lund, Sweden, 1996              *
!     V.P. Vysotskiy, University of Lund, Sweden, 2012                 *
!                                                                      *
!***********************************************************************

use Fast_IO, only: FlsSize, FSCB, isOpen, LuName, LuNameProf, MaxFileSize, MaxSplitFile, MBL, MPUnit, Multi_File, MxFile, &
                   NProfFiles, Trace
#if defined (_HAVE_EXTRA_) && ! defined (_GA_)
use Fast_IO, only: isFiM
#endif
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: Lu
integer(kind=iwp) :: i, iRc, Lu_, LuP
character(len=80) :: Text
character(len=*), parameter :: TheName = 'DaClos'
integer(kind=iwp), external :: AixCls, AixFsz
interface
  function AixErr(FileName) bind(C,name='aixerr_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: AixErr
    character(kind=c_char) :: FileName(*)
  end function AixErr
end interface

if (Trace) then
  write(u6,*) ' >>> Enter DaClos <<<'
  write(u6,*) ' unit :',Lu
  write(u6,*) ' name :',LuName(Lu)
end if
LuP = 0
do i=1,NProfFiles
  if (LuNameProf(i) == LuName(Lu)) then
    LuP = i
  end if
end do
#if defined (_HAVE_EXTRA_) && ! defined (_GA_)
if (isFiM(Lu) == 0) then
#endif
  FlsSize(LuP) = AixFsz(FSCB(Lu))
#if defined (_HAVE_EXTRA_) && ! defined (_GA_)
else
  FlsSize(LuP) = FimFsz(FSCB(Lu))
end if
#endif

if ((Lu <= 0) .or. (Lu > MxFile)) call SysFileMsg(TheName,'MSG: unit',Lu,' ')
if (isOpen(Lu) == 0) call SysFileMsg(TheName,'MSG: notopened',Lu,' ')
#if defined (_HAVE_EXTRA_) && ! defined (_GA_)
if (isFiM(Lu) == 0) then
#endif
  iRc = AixCls(FSCB(Lu))
#if defined (_HAVE_EXTRA_) && ! defined (_GA_)
else
  iRc = FimCls(FSCB(Lu))
  isFiM(Lu) = 0
end if
#endif
if (iRc /= 0) then
  iRc = AixErr(Text)
  call SysFileMsg(TheName,'MSG: close',Lu,Text)
end if
isOpen(Lu) = 0
MBL(Lu) = 0
if (Multi_File(Lu)) then
  if (MaxFileSize /= 0) then
    if (Trace) write(u6,*) ' This is a partitioned data set'
    do i=1,MaxSplitFile-1
      Lu_ = MPUnit(i,Lu)
      if (Lu_ >= 1) then
        irc = 0
        if (isOpen(Lu_) /= 0) iRc = AixCls(FSCB(Lu_))
        if (iRc /= 0) then
          iRc = AixErr(Text)
          call SysFileMsg(TheName,'MSG: close',Lu_,Text)
        end if
        isOpen(Lu_) = 0
        MPUnit(i,Lu) = -99
        Multi_File(Lu_) = .false.
        MBL(Lu_) = 0
      end if
    end do
  end if
  MPUnit(0,Lu) = 0
  Multi_File(Lu) = .false.
end if

if (Trace) then
  write(u6,*) ' >>> Exit DaClos <<<'
end if

return

end subroutine DaClos
