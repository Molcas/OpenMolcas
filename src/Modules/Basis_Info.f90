!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************
!#define _DEBUG_
      Module Basis_Info
      Implicit None
      Private
      Public :: Basis_Info_Dmp, Basis_Info_Get, Basis_Info_Free, &
                Distinct_Basis_set_Centers, dbsc
#include "stdalloc.fh"
#include "Molcas.fh"
      Integer, Parameter :: Mxdbsc=MxAtom
!     Work in progress
!
      Type Distinct_Basis_set_centers
          Real*8, Allocatable:: Coor(:,:)
          Integer:: nCntr=0
      End Type Distinct_Basis_set_centers
!
      Type (Distinct_Basis_set_centers) :: dbsc(Mxdbsc)
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
         Subroutine Put_dArray(Label,Data,nData)
         Character*(*) Label
         Integer       nData
         Real*8        Data(nData)
         End Subroutine Put_dArray
         Subroutine Get_dArray(Label,Data,nData)
         Character*(*) Label
         Integer       nData
         Real*8        Data(nData)
         End Subroutine Get_dArray
         Subroutine Qpg_dArray(Label,Found,nData)
         Character*(*) Label
         Logical       Found
         Integer       nData
         End Subroutine Qpg_dArray
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
      Subroutine Basis_Info_Dmp()
!
      Integer i, j, nCnttp, nAtoms
      Integer, Allocatable:: iDmp(:)
      Real*8, Allocatable:: rDmp(:,:)
!     Write (*,*) 'Basis_Info_Dmp()'
!
!     Temporary code until nCnttp has been move over to the Module
!
      i = 0
      Do
          i=i+1
          If (i.gt.Mxdbsc .or. dbsc(i)%nCntr.eq.0) Exit
      End Do
      nCnttp=i-1
#ifdef _DEBUG_
      Write (6,*) 'Basis_Info_Dmp'
      Do i = 1, nCnttp
         Do j = 1, dbsc(i)%nCntr
            Write (6,*) dbsc(i)%Coor(:,j)
         End Do
      End Do
#endif
      Call mma_Allocate(iDmp,nCnttp,Label='iDmp')
      nAtoms=0
      Do i = 1, nCnttp
         iDmp(i) = dbsc(i)%nCntr
         nAtoms=nAtoms+dbsc(i)%nCntr
      End Do
      Call Put_iArray('iDmp',iDmp,nCnttp)
      Call mma_deallocate(iDmp)
!
      Call mma_allocate(rDmp,3,nAtoms,Label='rDmp')
      nAtoms = 0
      Do i = 1, nCnttp
!        Call RecPrt('dbsc(i)%Coor',' ',dbsc(i)%Coor(1,1),3,dbsc(i)%nCntr)
         Do j = 1, dbsc(i)%nCntr
            nAtoms=nAtoms+1
            rDmp(1:3,nAtoms)=dbsc(i)%Coor(1:3,j)
         End Do
      End Do
      Call Put_dArray('rDmp',rDmp,3*nAtoms)
      Call mma_deallocate(rDmp)
      Return
      End Subroutine Basis_Info_Dmp
!
!***********************************************************************
!
      Subroutine Basis_Info_Get()
!
      Integer, Allocatable:: iDmp(:)
      Real*8, Allocatable:: rDmp(:,:)
      Logical Found
      Integer Len, i, j, nCnttp, nAtoms
!     Write (*,*) 'Basis_Info_Get()'
!
      Call qpg_iArray('iDmp',Found,Len)
      nCnttp=Len
      Call mma_Allocate(iDmp,nCnttp,Label='iDmp')
      If (Found) Call Get_iArray('iDmp',iDmp,nCnttp)
      Do i = 1, nCnttp
         dbsc(i)%nCntr  = iDmp(i)
      End Do
      Call mma_deallocate(iDmp)
!
      Call qpg_dArray('rDmp',Found,Len)
      If (.Not.Found) Then
         Write (6,*) 'rDMP not found on the run file.'
         Call Abend()
      End If
      nAtoms=Len/3
      Call mma_allocate(rDmp,3,nAtoms,Label='rDmp')
      Call Get_dArray('rDmp',rDmp,3*nAtoms)
      nAtoms = 0
      Do i = 1, nCnttp
         If (.Not.Allocated(dbsc(i)%Coor)) Then
            Call mma_Allocate(dbsc(i)%Coor,3,dbsc(i)%nCntr,Label='dbsc:C')
         End If
         Do j = 1, dbsc(i)%nCntr
            nAtoms=nAtoms+1
            dbsc(i)%Coor(1:3,j)=rDmp(1:3,nAtoms)
         End Do
      End Do
      Call mma_deallocate(rDmp)
#ifdef _DEBUG_
      Write (6,*) 'Basis_Info_Get'
      Do i = 1, nCnttp
         Do j = 1, dbsc(i)%nCntr
            Write (6,*) dbsc(i)%Coor(:,j)
         End Do
      End Do
#endif
      Return
      End Subroutine Basis_Info_Get
!
!***********************************************************************
!
      Subroutine Basis_Info_Free()
      Integer i
!     Write (*,*) 'Basis_Info_Free()'
!
!     Deallocate all allocatable parts of dbsc.
!
      i = 0
      Do
         i=i+1
         If (i.gt.Mxdbsc .or. dbsc(i)%nCntr.eq.0) Exit
!
         If (allocated(dbsc(i)%Coor)) Call mma_deallocate(dbsc(i)%Coor)
         dbsc(i)%nCntr=-1
      End Do
!
      Return
      End Subroutine Basis_Info_Free
!
!***********************************************************************
!
      End Module Basis_Info
