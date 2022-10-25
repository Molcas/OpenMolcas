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

subroutine DaEras(Lu)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Close unit Lu and remove the file                                *
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

use Fast_IO, only: FSCB, isOpen, LuName, MaxFileSize, MaxSplitFile, MPUnit, Multi_File, MxFile, Trace
#if defined (_HAVE_EXTRA_) && ! defined (_GA_)
use Fast_IO, only: isFiM
#endif
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: Lu
integer(kind=iwp) :: i, iRc, Lu_
character(len=80) :: Text
character(len=*), parameter :: TheName = 'DaEras'
integer(kind=iwp), external :: AixCls, AixRm
interface
  function AixErr(FileName) bind(C,name='aixerr_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: AixErr
    character(kind=c_char) :: FileName(*)
  end function AixErr
end interface

if (Trace) then
  write(u6,*) ' >>> Enter DaEras <<<'
  write(u6,*) ' unit :',Lu
end if

! Save calling arguments

if ((Lu <= 0) .or. (Lu > MxFile)) call SysFileMsg(TheName,'MSG: unit',Lu,' ')

if (isOpen(Lu) == 0) call SysFileMsg(TheName,'MSG: used',Lu,' ')
#if defined (_HAVE_EXTRA_) && ! defined (_GA_)
if (isFiM(Lu) == 0) then
#endif
  iRc = AixCls(FSCB(Lu))
#if defined (_HAVE_EXTRA_) && ! defined (_GA_)
else
  iRc = FimCls(FSCB(Lu))
end if
#endif
if (iRc /= 0) then
  iRc = AixErr(Text)
  call SysFileMsg(TheName,'MSG: close',Lu,Text)
end if
#if defined (_HAVE_EXTRA_) && ! defined (_GA_)
if (isFiM(Lu) == 0) then
#endif
  iRc = AixRm(LuName(Lu))
  if (iRc /= 0) then
    iRc = AixErr(Text)
    call SysFileMsg(TheName,'MSG: delete',Lu,Text)
  end if
#if defined (_HAVE_EXTRA_) && ! defined (_GA_)
end if
#endif
isOpen(Lu) = 0

if (Multi_File(Lu)) then
  if (MaxFileSize /= 0) then
    if (Trace) then
      write(u6,*) ' This is a partitioned data set'
    end if
    do i=1,MaxSplitFile-1
      Lu_ = MPUnit(i,Lu)
      if (Lu_ >= 1) then
        if (isOpen(Lu_) /= 0) then
          iRc = AixCls(FSCB(Lu_))
          if (iRc /= 0) then
            iRc = AixErr(Text)
            call SysFileMsg(TheName,'MSG: close',Lu_,Text)
          end if
          iRc = AixRm(LuName(Lu_))
          if (iRc /= 0) then
            iRc = AixErr(Text)
            call SysFileMsg(TheName,'MSG: delete',Lu_,Text)
          end if
          isOpen(Lu_) = 0
        end if
      end if
    end do
  end if
end if

if (Trace) then
  write(u6,*) ' >>> Exit DaEras <<<'
end if

return

end subroutine DaEras
