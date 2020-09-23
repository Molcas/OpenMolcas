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
Public :: iRI_Type, LDF, Do_RI, Cholesky, Do_acCD_Basis, Skip_High_AC, &
           Cho_OneCenter,  DiagCheck, LocalDF, Do_nacCD_Basis, &
          RICD_Info_Dmp, RICD_Info_Get
Integer :: iRI_Type=-1
Logical :: LDF=.False.
Logical :: Do_RI=.False.
Logical :: Cholesky=.False.
Logical :: Do_acCD_Basis=.True.
Logical :: Skip_High_AC=.False.
Logical :: Cho_OneCenter=.False.
Logical :: DiagCheck=.False.
Logical :: LocalDF=.False.
Logical :: Do_nacCD_Basis=.False.
#include "stdalloc.fh"

Contains

Subroutine RICD_Info_Dmp()
  Real*8, Allocatable:: rDmp(:)
  Integer i
  Integer:: Len=10

  Call mma_allocate(rDmp,Len,Label='rDmp:RICD')

  rDmp(1) = DBLE(iRI_Type)
  i = 0
  If (LDF) i = 1
  rDmp(2) = DBLE(i)
  i = 0
  If (Do_RI) i = 1
  rDmp(3) = DBLE(i)
  i = 0
  If (Cholesky) i = 1
  rDmp(4) = DBLE(i)
  i = 0
  If (Do_acCD_Basis) i = 1
  rDmp(5) = DBLE(i)
  i = 0
  If (Skip_High_AC) i = 1
  rDmp(6) = DBLE(i)
  i = 0
  If (Cho_OneCenter) i = 1
  rDmp(7) = DBLE(i)
  i = 0
  If (DiagCheck) i = 1
  rDmp(8) = DBLE(i)
  i = 0
  If (LocalDF) i = 1
  rDmp(9) = DBLE(i)
  i = 0
  If (Do_nacCD_Basis) i = 1
  rDmp(10) = DBLE(i)


  Call Put_dArray('RICD_Info',rDmp,Len)
  Call mma_deallocate(rDmp)
End Subroutine RICD_Info_Dmp


Subroutine RICD_Info_Get()
  Real*8, Allocatable:: rDmp(:)
  Integer:: Len=10

  Call mma_allocate(rDmp,Len,Label='rDmp:RICD')
  Call Get_dArray('RICD_Info',rDmp,Len)

  iRI_Type=NINT(rDmp(1))
  LDF = NINT(rDmp(2)).eq.1
  Do_RI = NINT(rDmp(3)).eq.1
  Cholesky = NINT(rDmp(4)).eq.1
  Do_acCD_Basis = NINT(rDmp(5)).eq.1
  Skip_High_AC = NINT(rDmp(6)).eq.1
  Cho_OneCenter = NINT(rDmp(7)).eq.1
  DiagCheck = NINT(rDmp(8)).eq.1
  LocalDF = NINT(rDmp(9)).eq.1
  Do_nacCD_Basis = NINT(rDmp(10)).eq.1

  Call mma_deallocate(rDmp)
End Subroutine RICD_Info_Get

End Module RICD_Info
