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
! Copyright (C) 2017,2022, Roland Lindh                                *
!***********************************************************************
      Subroutine TraClc_x_qNR(QNR1st)
      Implicit None
#include "real.fh"
#include "stdalloc.fh"
      Logical QNR1st

!     Extrapolation case.
!
!     displacements: del = -H^(-1)g, where g=dE/dX_m

      If (QNR1st) Then

!        compute actual gradient
         Call GrdClc('All',.True.)

         QNR1st=.FALSE.
      Else

!        Note that the required displacements are actually computed
!        in linser!

         Call GrdClc('Lst',.True.)

      End If

      Return
      End
