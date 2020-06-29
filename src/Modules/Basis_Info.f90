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
      Public :: Basis_Info_Dmp, Basis_Info_Get, &
                Distinct_Basis_set_Centers, dbsc
#include "stdalloc.fh"
#include "Molcas.fh"
      Integer, Parameter :: Mxdbsc=MxAtom
!     Work in progress
      Type Distinct_Basis_set_centers
          Integer:: ipCntr
          Integer:: nCntr=-1
      End Type Distinct_Basis_set_centers
!
      Type (Distinct_Basis_set_centers) :: dbsc(Mxdbsc)

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
      Integer nDmp, i, nCnttp
      Real*8, Allocatable:: iDmp(:,:)
!
!     Temporary code until nCnttp has been move over to the Module
!
      i = 0
      Do
          i=i+1
          If (i.gt.Mxdbsc .or. dbsc(i)%nCntr.eq.-1) Exit
      End Do
      nCnttp=i-1
#ifdef _DEBUG_
      Write (6,*) 'Basis_Info_Dmp'
      Do i = 1, nCnttp
         Write (6,*) dbsc(i)%ipCntr, dbsc(i)%nCntr
      End Do
#endif
      Call mma_Allocate(iDmp,2,nCnttp,Label='iDmp')
      Do i = 1, nCnttp
         iDmp(1,i) = dbsc(i)%ipCntr
         iDmp(2,i) = dbsc(i)%nCntr
      End Do
      Call Put_iArray('iDmp',iDmp,2*nCnttp)
      Call mma_deallocate(iDmp)
      End Subroutine Basis_Info_Dmp
!
!***********************************************************************
!
      Subroutine Basis_Info_Get()
      Real*8, Allocatable:: iDmp(:,:)
      Logical Found
      Integer Len, i, nCnttp
      Call qpg_iArray('iDmp',Found,Len)
      nCnttp=Len/2
      Call mma_Allocate(iDmp,2,nCnttp,Label='iDmp')
      If (Found) Call Get_iArray('iDmp',iDmp,2*nCnttp)
      Do i = 1, nCnttp
         dbsc(i)%ipCntr = iDmp(1,i)
         dbsc(i)%nCntr  = iDmp(2,i)
      End Do
      Call mma_deallocate(iDmp)
#ifdef _DEBUG_
      Write (6,*) 'Basis_Info_Get'
      Do i = 1, nCnttp
         Write (6,*) dbsc(i)%ipCntr, dbsc(i)%nCntr
      End Do
#endif
      End Subroutine Basis_Info_Get
!
!***********************************************************************
!
      End Module Basis_Info
