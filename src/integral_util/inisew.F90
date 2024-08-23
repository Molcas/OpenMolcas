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
! Copyright (C) 1999, Roland Lindh                                     *
!***********************************************************************
      SubRoutine IniSew(DSCF,nDiff)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Chemical Physics,                 *
!             University of Lund, SWEDEN                               *
!***********************************************************************
      use Basis_Info, only: Seward_Activated
      Implicit None
#include "SysDef.fh"
      Logical DSCF
      Integer nDiff
!
      If (Seward_Activated) Then
         Call ClsSew()
         Call xRlsMem_Ints()
      End If
!
      Call Seward_Init()
!
      Call GetInf(DSCF,nDiff)
!
      Return
      End SubRoutine IniSew
