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

      SubRoutine Integral_WrOut2(                                       &
#define _FIXED_FORMAT_
#define _CALLING_
#include "int_wrout_interface.fh"
     &                          )
!     calls the proper routines IndSft/PLF
!     if IntOrd_jikl==.TRUE. integral order within symblk: jikl
!                      else  integral order within symblk: ijkl
      Implicit None
!
#include "int_wrout_interface.fh"
!
      If (mSym.eq.1) Then
        Call PLF2(AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),           &
     &           iShell,iAO,iAOst,iBas,jBas,kBas,lBas,kOp)
      Else
        Call IndSft2(iCmp,iShell,                                       &
     &               iBas,jBas,kBas,lBas,Shijij,                        &
     &               iAO,iAOst,ijkl,SOInt,nSOint,iSOSym,nSOs)
      End If
!
! Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(TInt)
         Call Unused_integer(mSym)
      End If
      End SubRoutine Integral_WrOut2

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
      dummy_empty_procedure(Integral_WrOut2)

#endif

