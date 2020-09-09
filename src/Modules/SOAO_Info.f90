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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************
#define _DEBUG_

Module SOAO_Info
Implicit None
Private
Public :: iSOInf, nSOInf, SOAO_Info_Init, SOAO_Info_Dmp, SOAO_Info_Get, SOAO_Info_Free
#include "stdalloc.fh"
Integer, Allocatable:: iSOInf(:,:)
Integer :: nSOInf=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
Interface
   Subroutine Abend()
   End Subroutine Abend
   Subroutine Put_iArray(Label,Data,nData)
   Character*(*) Label
   Integer       nData
   Integer       Data(nData)
   End Subroutine Put_iArray
   Subroutine Get_iArray(Label,Data,nData)
   Character*(*) Label
   Integer       nData
   Integer       Data(nData)
   End Subroutine Get_iArray
   Subroutine Qpg_iArray(Label,Found,nData)
   Character*(*) Label
   Logical       Found
   Integer       nData
   End Subroutine Qpg_iArray
End Interface
!
!***********************************************************************
!***********************************************************************
!
Contains
!
!***********************************************************************
!***********************************************************************
!
!     This to make either the initial allocation of dbsc and Shells according to the default sizes
!     as defined by the parameters in Molcas.fh or according to the actual sizes as recorded on the
!     run file.
!
Subroutine SOAO_Info_Init(nSize)
Implicit None
Integer nSize
If (Allocated(iSOInf)) Call SOAO_Info_Free()
nSOInf=nSize
Call mma_allocate(iSOInf,3,nSOInf,Label='iSOInf')
End Subroutine SOAO_Info_Init

Subroutine SOAO_Info_Dmp()
Implicit None
Call Put_iArray('iSOInf',iSOInf,3*nSOInf)
End Subroutine SOAO_Info_Dmp

Subroutine SOAO_Info_Get()
Implicit None
Logical Found
Call Qpg_iArray('iSOInf',Found,nSOInf)
If (.NOT.Found) Then
   Write (6,*) 'SOAO_Info_Get: iSOInf not found.'
   Call Abend()
End If
nSOInf=nSOInf/3
If (Allocated(iSOInf)) Call SOAO_Info_Free()
Call mma_allocate(iSOInf,3,nSOInf,Label='iSOInf')
Call Get_iArray('iSOInf',iSOInf,3*nSOInf)
End Subroutine SOAO_Info_Get

Subroutine SOAO_Info_Free()
Implicit None
If (.Not.Allocated(iSOInf)) Return
Call mma_deallocate(iSOInf)
nSOInf=0
End Subroutine SOAO_Info_Free
!
!***********************************************************************
!***********************************************************************
!
End Module SOAO_Info
