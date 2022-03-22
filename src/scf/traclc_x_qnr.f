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
      Subroutine TraClc_x_qNR(kOptim,QNR1st,CInter,nCI,nD,nOV,iter,LLx)
      Implicit None
#include "real.fh"
#include "stdalloc.fh"
      Integer kOptim,nCI,nD,nOV,iter,LLx
      Logical QNR1st
      Real*8 CInter(nCI,nD)
      Real*8, Dimension(:,:), Allocatable:: Xn

!     Extrapolation case.
!
!     displacements: del = -H^(-1)g, where g=dE/dX_m

      If (QNR1st) Then

!------  1st QNR step, reset kOptim to 1

         kOptim = 1
         CInter(1,1) = One
         CInter(1,nD) = One

!        init 1st orb rot parameter X1 (set it to zero)
         Call mma_allocate(Xn,nOV,nD,Label='Xn')
         Call FZero(Xn,nOV*nD)
!        and store it on appropriate LList
         Call PutVec(Xn,nOV*nD,iter,'NOOP',LLx)
         Call mma_deallocate(Xn)

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
