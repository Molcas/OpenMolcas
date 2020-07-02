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
!
!     Work in progress
!
!     nCntr : number of centers associated with a dbsc
!     Coor  : the coordinates of a dbsc
!     nM1   : number of ECP M1 type terms on the i''th unique center
!     M1xp  : ECP M1-type exponents for i''th unq center
!     M1cf  : ECP M1 type coefficients for i''th unq cntr
!     nM2   : number of ECP M2 type terms on the i''th unique center
!     M2xp  : ECP M2-type exponents for i''th unq center
!     M2cf  : ECP M2 type coefficients for i''th unq cntr

      Type Distinct_Basis_set_centers
          Sequence
          Real*8, Allocatable:: Coor(:,:)
          Integer:: nCntr=0
          Integer:: nM1=0
          Real*8, Allocatable:: M1xp(:), M1cf(:)
          Integer:: nM2=0
          Real*8, Allocatable:: M2xp(:), M2cf(:)
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
      Integer i, j, nCnttp, nAtoms, nAux, nM1, nM2
      Integer, Allocatable:: iDmp(:,:)
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
      Call mma_Allocate(iDmp,3,nCnttp,Label='iDmp')
      nAtoms=0
      nAux   = 0
      Do i = 1, nCnttp
         iDmp(1,i) = dbsc(i)%nCntr
         iDmp(2,i) = dbsc(i)%nM1
         iDmp(3,i) = dbsc(i)%nM2
         nAtoms=nAtoms+dbsc(i)%nCntr
         nAux = nAux + dbsc(i)%nM1 + dbsc(i)%nM2
      End Do
      Call Put_iArray('iDmp',iDmp,3*nCnttp)
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
!
      If (nAux.ne.0) Then
      Call mma_allocate(rDmp,nAux,1,Label='rDmp')
      nAux=1
      Do i = 1, nCnttp
         nM1 = dbsc(i)%nM1
         If (nM1.gt.0) Then
            rDmp(nAux:nAux+nM1,1) = dbsc(i)%M1xp(:)
            nAux = nAux + nM1
            rDmp(nAux:nAux+nM1,1) = dbsc(i)%M1cf(:)
            nAux = nAux + nM1
         End If
         nM2 = dbsc(i)%nM2
         If (nM2.gt.0) Then
            rDmp(nAux:nAux+nM2,1) = dbsc(i)%M2xp(:)
            nAux = nAux + nM2
            rDmp(nAux:nAux+nM2,1) = dbsc(i)%M2cf(:)
            nAux = nAux + nM2
         End If
      End Do
      Call Put_dArray('rDmp:A',rDmp,nAux)
      Call mma_deallocate(rDmp)
      End If
!
      Return
      End Subroutine Basis_Info_Dmp
!
!***********************************************************************
!
      Subroutine Basis_Info_Get()
!
      Integer, Allocatable:: iDmp(:,:)
      Real*8, Allocatable:: rDmp(:,:)
      Logical Found
      Integer Len, i, j, nCnttp, nAtoms, nAux, nM1, nM2
!     Write (*,*) 'Basis_Info_Get()'
!
      Call qpg_iArray('iDmp',Found,Len)
      nCnttp=Len/3
      Call mma_Allocate(iDmp,3,nCnttp,Label='iDmp')
      If (Found) Call Get_iArray('iDmp',iDmp,3*nCnttp)
      nAux = 0
      Do i = 1, nCnttp
         dbsc(i)%nCntr  = iDmp(1,i)
         dbsc(i)%nM1    = iDmp(2,i)
         dbsc(i)%nM2    = iDmp(3,i)
         nAux = nAux + iDmp(2,i) + iDmp(3,i)
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
!
      If (nAux.gt.0) Then
      Call qpg_dArray('rDmp:A',Found,Len)
      Call mma_allocate(rDmp,nAux,1,Label='rDmp')
      nAux=1
      Do i = 1, nCnttp
         nM1 = dbsc(i)%nM1
         If (nM1.gt.0) Then
            If (.Not.Allocated(dbsc(i)%M1xp)) Call mma_allocate(dbsc(i)%M1xp,nM1,Label='dbsc:M1xp')
            dbsc(i)%M1xp(:)=rDmp(nAux:nAux+nM1,1)
            nAux=nAux+nM1
            If (.Not.Allocated(dbsc(i)%M1cf)) Call mma_allocate(dbsc(i)%M1cf,nM1,Label='dbsc:M1cf')
            dbsc(i)%M1cf(:)=rDmp(nAux:nAux+nM1,1)
            nAux=nAux+nM1
         End If
         nM2 = dbsc(i)%nM2
         If (nM2.gt.0) Then
            If (.Not.Allocated(dbsc(i)%M2xp)) Call mma_allocate(dbsc(i)%M2xp,nM2,Label='dbsc:M2xp')
            dbsc(i)%M2xp(:)=rDmp(nAux:nAux+nM2,1)
            nAux=nAux+nM2
            If (.Not.Allocated(dbsc(i)%M2cf)) Call mma_allocate(dbsc(i)%M2cf,nM2,Label='dbsc:M2cf')
            dbsc(i)%M2cf(:)=rDmp(nAux:nAux+nM2,1)
            nAux=nAux+nM2
         End If
      End Do
      Call mma_deallocate(rDmp)
      End If
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
         If (allocated(dbsc(i)%M1xp)) Call mma_deallocate(dbsc(i)%M1xp)
         If (allocated(dbsc(i)%M1cf)) Call mma_deallocate(dbsc(i)%M1cf)
         dbsc(i)%nM1=0
         If (allocated(dbsc(i)%M2xp)) Call mma_deallocate(dbsc(i)%M2xp)
         If (allocated(dbsc(i)%M2cf)) Call mma_deallocate(dbsc(i)%M2cf)
         dbsc(i)%nM2=0
      End Do
!
      Return
      End Subroutine Basis_Info_Free
!
!***********************************************************************
!
      End Module Basis_Info
