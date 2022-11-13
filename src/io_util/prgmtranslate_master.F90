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
! Copyright (C) 2001-2016, Valera Veryazov                             *
!***********************************************************************

subroutine prgmtranslate_master(namein,nameout,lout)

use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: namein
character(len=*), intent(out) :: nameout
integer(kind=iwp), intent(out) :: lout
integer(kind=iwp), external :: Strnln
integer(kind=iwp) :: lin
interface
  subroutine PrgmTranslateC(InStr,l1,OutStr,l2,Par) bind(C,name='prgmtranslatec_')
  use, intrinsic :: iso_c_binding, only: c_char
  use Definitions, only: MOLCAS_C_INT
# include "intent.fh"
  character(kind=c_char), intent(in) :: InStr(*)
  integer(kind=MOLCAS_C_INT), intent(in) :: l1, Par
  character(kind=c_char), intent(_OUT_) :: OutStr(*)
  integer(kind=MOLCAS_C_INT), intent(out) :: l2
  end subroutine PrgmTranslateC
end interface

lin = Strnln(namein)
nameout = ' '
if (index(namein,'/') /= 0) then
  ! just in case we process a translated name!
  nameout = namein
  lout = lin
else
  call prgmtranslatec(namein,lin,nameout,lout,0)
end if
nameout = nameout(1:lout)

return

end subroutine prgmtranslate_master
