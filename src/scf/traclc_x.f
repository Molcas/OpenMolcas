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
      Subroutine TraClc_x(FrstDs,iOpt)
      Implicit None
#include "real.fh"
#include "stdalloc.fh"
      Integer iOpt
      Logical FrstDs

!     Extrapolation case.
!
!     gradients: energy derivatives w.r.t the antisymmetric matrix, X,
!                which defines the orbital rotations, exp[P(X)]

!---  only DIIS, compute gradients

      If (FrstDs) Then
!        On first iteration compute all gradients and put them on file.
         Call GrdClc('All',iOpt)
         FrstDs=.FALSE.
      Else
!        Compute just the last one.
         Call GrdClc('Lst',iOpt)
      End If

      Return
      End
