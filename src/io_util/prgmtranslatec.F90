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
! Copyright (C) 2017, Ignacio Fdez. Galvan                             *
!***********************************************************************

! Wrapper for the OpenMolcas PrgmTranslate_Mod routine, so it works
! without using the module (i.e., from C). Some C-to-Fortran conversion
! needs to be done, which feels quite hackish.

#ifndef _HAVE_EXTRA_
#define MAXSTR 1024

subroutine PrgmTranslateC(InStr,l1,OutStr,l2,Par)

use iso_c_binding, only: c_char, c_null_char
use Prgm
implicit none
character(Kind=c_char), intent(In) :: InStr(*)
character(Kind=c_char), intent(Out) :: OutStr(*)
integer, intent(In) :: Par, l1
integer, intent(Out) :: l2
character(Len=MAXSTR) :: TmpStr, TmpStr2
integer :: i

TmpStr = ''
do i=1,l1
  TmpStr(i:i) = InStr(i)
end do
call PrgmTranslate_Mod(TmpStr,l1,TmpStr2,l2,Par)
do i=1,l2
  OutStr(i) = TmpStr2(i:i)
end do
OutStr(l2+1:l2+1) = c_null_char

end subroutine PrgmTranslateC

#elif defined (NAGFOR)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(PrgmTranslateC)

#endif
