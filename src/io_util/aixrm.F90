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

integer function AixRm(name)
implicit integer(a-z)
#include "switch.fh"
#include "ctl.fh"
character*(*) name
character*256 tmp, out
character*80 ErrTxt

!----------------------------------------------------------------------*
! Entry to AixRm                                                       *
!----------------------------------------------------------------------*
AixRm = 0
!----------------------------------------------------------------------*
! Strip file name and append string terminator                         *
!----------------------------------------------------------------------*
n = len(name)
100 continue
if (name(n:n) == ' ') then
  n = n-1
  if (n <= 0) then
    AixRm = eBlNme
    return
  end if
  Go To 100
end if
n = n+1
if (n >= len(tmp)) then
  AixRm = eTlFn
  return
end if
tmp = name
tmp(n:n) = char(0)
!----------------------------------------------------------------------*
! erase file                                                           *
!----------------------------------------------------------------------*
out = ' '

call PrgmTranslate(Name,out,ltmp)
out(ltmp+1:ltmp+1) = char(0)
rc = c_remove(out)
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
