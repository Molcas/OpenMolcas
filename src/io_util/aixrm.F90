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
! Copyright (C) 1990, Per-Olof Widmark                                 *
!***********************************************************************
!***********************************************************************
!                                                                      *
!                             A I X - I / O                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! rc=AixRm(FileName)                                                   *
!                                                                      *
! Erase file with specified file name.                                 *
!                                                                      *
! Input:  FileName - This file will be erased. Given as a character    *
!                    string.                                           *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          S&TC, ACIS, IBM Sweden                                      *
! Written: November 1990                                               *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! History:                                                             *
!                                                                      *
!***********************************************************************

function AixRm(filename)

use Fast_IO, only: eBlNme, eTlFn
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: AixRm
character(len=*), intent(in) :: filename
integer(kind=iwp) :: n, ltmp, rc
character(len=256) :: tmp, outname
character(len=80) :: ErrTxt
interface
  function AixErr(FileName) bind(C,name='aixerr_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: AixErr
    character(kind=c_char) :: FileName(*)
  end function AixErr
  function c_remove(FileName) bind(C,name='c_remove_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: c_remove
    character(kind=c_char) :: FileName(*)
  end function c_remove
end interface

!----------------------------------------------------------------------*
! Entry to AixRm                                                       *
!----------------------------------------------------------------------*
AixRm = 0
!----------------------------------------------------------------------*
! Strip file name and append string terminator                         *
!----------------------------------------------------------------------*
n = len(filename)
do
  if (filename(n:n) /= ' ') exit
  n = n-1
  if (n <= 0) then
    AixRm = eBlNme
    return
  end if
end do
n = n+1
if (n >= len(tmp)) then
  AixRm = eTlFn
  return
end if
tmp = filename
tmp(n:n) = char(0)
!----------------------------------------------------------------------*
! erase file                                                           *
!----------------------------------------------------------------------*
outname = ' '

call PrgmTranslate(filename,outname,ltmp)
outname(ltmp+1:ltmp+1) = char(0)
rc = c_remove(outname)
if (rc /= 0) then
  AixRm = AixErr(ErrTxt)
  call SysAbendMsg('AixRm','MSG: delete',ErrTxt)
  return
end if
!----------------------------------------------------------------------*
! Finished so return to caller                                         *
!----------------------------------------------------------------------*
return

end function AixRm
