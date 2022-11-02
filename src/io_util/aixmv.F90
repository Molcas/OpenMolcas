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
!               2020, Ignacio Fdez. Galvan                             *
!***********************************************************************
!***********************************************************************
!                                                                      *
!                             A I X - I / O                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! rc=AixMv(FileName,NewName)                                           *
!                                                                      *
! Rename (move) file.                                                  *
!                                                                      *
! Input:  FileName - This file will be renamed. Given as a character   *
!                    string.                                           *
! Input:  NewName  - New name for the file. Given as a character       *
!                    string.                                           *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Ignacio Fdez. Galvan                                        *
! Written: April 2020                                                  *
!          (based on AixRm)                                            *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! History:                                                             *
!                                                                      *
!***********************************************************************

function AixMv(FileName,NewName)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: AixMv
character(len=*), intent(in) :: FileName, NewName
integer(kind=iwp) :: ltmp, rc
character(len=80) :: ErrTxt
character(len=256) :: out1, out2
interface
  function AixErr(FileName) bind(C,name='aixerr_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: AixErr
    character(kind=c_char) :: FileName(*)
  end function AixErr
  function c_rename(FileName,NewName) bind(C,name='c_rename_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    integer(kind=MOLCAS_C_INT) :: c_rename
    character(kind=c_char) :: FileName(*), NewName(*)
  end function c_rename
end interface

!----------------------------------------------------------------------*
! Entry to AixMv                                                       *
!----------------------------------------------------------------------*
AixMv = 0
!----------------------------------------------------------------------*
! rename file                                                          *
!----------------------------------------------------------------------*
out1 = ' '
out2 = ' '

call PrgmTranslate(FileName,out1,ltmp)
out1(ltmp+1:ltmp+1) = char(0)
call PrgmTranslate(NewName,out2,ltmp)
out2(ltmp+1:ltmp+1) = char(0)
rc = c_rename(out1,out2)
if (rc /= 0) then
  AixMv = AixErr(ErrTxt)
  call SysAbendMsg('AixMv','MSG: rename',ErrTxt)
  return
end if
!----------------------------------------------------------------------*
! Finished so return to caller                                         *
!----------------------------------------------------------------------*
return

end function AixMv
