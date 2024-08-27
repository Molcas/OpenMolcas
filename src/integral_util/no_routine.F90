!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

#include "compiler_features.h"
#ifdef _IN_MODULE_

subroutine No_Routine( &
#                     define _CALLING_
#                     include "int_wrout_interface.fh"
                     )
implicit none
#include "int_wrout_interface.fh"

! Avoid unused argument warnings
if (.false.) then
  call Unused_integer_array(iCmp)
  call Unused_integer_array(iShell)
  call Unused_integer(iBas)
  call Unused_integer(jBas)
  call Unused_integer(kBas)
  call Unused_integer(lBas)
  call Unused_integer_array(kOp)
  call Unused_logical(Shijij)
  call Unused_integer_array(iAO)
  call Unused_integer_array(iAOst)
  call Unused_integer(ijkl)
  call Unused_real_array(AOInt)
  call Unused_real_array(SOInt)
  call Unused_integer(nSOint)
  call Unused_integer_array(iSOSym)
  call Unused_real_array(TInt)
  call Unused_integer(mSym)
end if
end subroutine No_Routine

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(No_Routine)

#endif
