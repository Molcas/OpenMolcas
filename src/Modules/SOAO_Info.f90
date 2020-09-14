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
!#define _DEBUG_

Module SOAO_Info
Implicit None
Private
Public :: iSOInf, iAOtSO, nSOInf, SOAO_Info_Init, SOAO_Info_Dmp, SOAO_Info_Get, SOAO_Info_Free
#include "stdalloc.fh"
Integer, Allocatable:: iSOInf(:,:)
Integer, Allocatable:: iAOtSO(:,:)
Integer :: nSOInf=0
Integer :: nIrrep=0
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

!***********************************************************************

Subroutine SOAO_Info_Init(nSize,mIrrep)
Implicit None
#ifdef _DEBUG_
#endif
Integer nSize,mIrrep
If (Allocated(iSOInf).or.Allocated(iAOtSO)) Call SOAO_Info_Free()
nSOInf=nSize
nIrrep=mIrrep
Call mma_allocate(iSOInf,3,nSOInf,Label='iSOInf')
iSOInf(:,:)=-99999999
! Do not change this value, since it also explicitly signals symmetry information.
! This is explicitly used in some routines.
Call mma_allocate(iAOtSO,[1,nSOInf],[0,nIrrep-1],Label='iAOtSO')
iAOtSO(:,:)=-99999999        ! Dito
End Subroutine SOAO_Info_Init

!***********************************************************************

Subroutine SOAO_Info_Dmp()
Implicit None
#ifdef _DEBUG_
Write (6,*)
Write (6,*) 'Enter SOAO_Info_Dmp'
Write (6,*)
Block
Integer i
Write (6,*)
Write (6,*) 'iAOtSO:'
Do i = 0, nIrrep-1
   Write (6,'(8I9)') iAOtSO(1:nSOInf,i)
End Do
End Block
#endif
Call Put_iArray('iSOInf',iSOInf,3*nSOInf)
Call Put_iArray('iAOtSO',iAOtSO,nSOInf*nIrrep)
End Subroutine SOAO_Info_Dmp

!***********************************************************************

Subroutine SOAO_Info_Get()
Implicit None
Logical Found
#ifdef _DEBUG_
Write (6,*)
Write (6,*) 'Enter SOAO_Info_Get'
Write (6,*)
#endif
If (Allocated(iSOInf).or.Allocated(iAOtSO)) Call SOAO_Info_Free()
Call Qpg_iArray('iSOInf',Found,nSOInf)
If (.NOT.Found) Then
   Write (6,*) 'SOAO_Info_Get: iSOInf not found.'
   Call Abend()
End If
nSOInf=nSOInf/3
Call mma_allocate(iSOInf,3,nSOInf,Label='iSOInf')
Call Get_iArray('iSOInf',iSOInf,3*nSOInf)
!
Call Qpg_iArray('iAOtSO',Found,nIrrep)
If (.NOT.Found) Then
   Write (6,*) 'SOAO_Info_Get: iAOtSO not found.'
   Call Abend()
End If
nIrrep=nIrrep/nSOInf
Call mma_allocate(iAOtSO,[1,nSOInf],[0,nIrrep-1],Label='iAOtSO')
Call Get_iArray('iAOtSO',iAOtSO,nSOInf*nIrrep)
#ifdef _DEBUG_
Block
Integer i
Write (6,*)
Write (6,*) 'iAOtSO:'
Do i = 0, nIrrep-1
   Write (6,'(8I9)') iAOtSO(1:nSOInf,i)
End Do
Write (6,*)
Write (6,*) 'Exit SOAO_Info_Get'
Write (6,*)
End Block
#endif
End Subroutine SOAO_Info_Get

!***********************************************************************

Subroutine SOAO_Info_Free()
Implicit None
If (Allocated(iSOInf)) Call mma_deallocate(iSOInf)
If (Allocated(iAOtSO)) Call mma_deallocate(iAOtSO)
nSOInf=0
nIrrep=0
End Subroutine SOAO_Info_Free
!
!***********************************************************************
!***********************************************************************
!
End Module SOAO_Info
