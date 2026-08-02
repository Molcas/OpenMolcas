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
! Copyright (C) 2000-2016, Valera Veryazov                             *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Author:   Valera Veryazov 2000-2016                                  *
!           Theoretical Chemistry                                      *
!           Lund University                                            *
!           Sweden                                                     *
!                                                                      *
!***********************************************************************
!
!  This is a simple wrapper for getenv
!

subroutine getenvf(VarName,Val)

use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: VarName
character(len=*), intent(out) :: Val
integer(kind=iwp) :: ilen, irl, maxlen
interface
  subroutine getenvf2c(name_,ilen,value_,maxlen,irl) bind(C,name='getenvf2c_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    character(kind=c_char) :: name_(*), value_(*)
    integer(kind=MOLCAS_C_INT) :: ilen, maxlen, irl
  end subroutine getenvf2c
end interface

Val = ' '
ilen = len(VarName)
maxlen = len(Val)
call getenvf2c(VarName,ilen,Val,maxlen,irl)
if (irl == 0) then
  Val = ' '
else
  Val = Val(1:irl)
end if

return

end subroutine getenvf
