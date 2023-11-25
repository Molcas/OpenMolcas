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

      Subroutine No_Routine(                                            &
#define _FIXED_FORMAT_
#define _CALLING_
#include "int_wrout_interface.fh"
     &                     )
      Implicit None
#include "int_wrout_interface.fh"
! Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(iCmp)
         Call Unused_integer_array(iShell)
         Call Unused_integer(iBas)
         Call Unused_integer(jBas)
         Call Unused_integer(kBas)
         Call Unused_integer(lBas)
         Call Unused_integer_array(kOp)
         Call Unused_logical(Shijij)
         Call Unused_integer_array(iAO)
         Call Unused_integer_array(iAOst)
         Call Unused_integer(ijkl)
         Call Unused_real_array(AOInt)
         Call Unused_real_array(SOInt)
         Call Unused_integer(nSOint)
         Call Unused_integer_array(iSOSym)
         Call Unused_real_array(TInt)
         Call Unused_integer(mSym)
      End If
      End Subroutine No_Routine

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
      dummy_empty_procedure(No_Routine)

#endif

