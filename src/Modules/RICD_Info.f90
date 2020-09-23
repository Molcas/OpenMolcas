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
!
Module RICD_Info
Private
Public :: iRI_Type, LDF, Do_RI, Cholesky, Do_acCD_Basis, &
          RICD_Info_Dmp, RICD_Info_Get
Integer :: iRI_Type=-1
Logical :: LDF=.False.
Logical :: Do_RI=.False.
Logical :: Cholesky=.False.
Logical :: Do_acCD_Basis=.True.
#include "stdalloc.fh"

Contains

Subroutine RICD_Info_Dmp()
  Real*8, Allocatable:: rDmp(:)
  Integer i
  Integer:: Len=5

  Call mma_allocate(rDmp,Len,Label='rDmp:RICD')

  rDmp(1) = DBLE(iRI_Type)
  i = 0
  If (LDF) i = 1
  rDmp(2) = i
  i = 0
  If (Do_RI) i = 1
  rDmp(3) = i
  i = 0
  If (Cholesky) i = 1
  rDmp(4) = i
  i = 0
  If (Do_acCD_Basis) i = 1
  rDmp(5) = i


  Call Put_dArray('RICD_Info',rDmp,Len)
  Call mma_deallocate(rDmp)
End Subroutine RICD_Info_Dmp


Subroutine RICD_Info_Get()
  Real*8, Allocatable:: rDmp(:)
  Integer:: Len=5

  Call mma_allocate(rDmp,Len,Label='rDmp:RICD')
  Call Get_dArray('RICD_Info',rDmp,Len)

  iRI_Type=NINT(rDmp(1))
  LDF = NINT(rDmp(2)).eq.1
  Do_RI = NINT(rDmp(3)).eq.1
  Cholesky = NINT(rDmp(4)).eq.1
  Do_acCD_Basis = NINT(rDmp(5)).eq.1

  Call mma_deallocate(rDmp)
End Subroutine RICD_Info_Get

End Module RICD_Info
