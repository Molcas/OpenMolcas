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
! Copyright (C) 2021, Roland Lindh                                     *
!***********************************************************************
!***********************************************************************
!         MODULE for TRANSFORMATION of CHOLESKY VECTORS                *
!              and GENERATION of TWO-ELECTRONS INTEGRALS               *
!***********************************************************************
Module Data_Structures
Private
Public:: CMO_Type, Deallocate_CMO, Allocate_CMO, Map_to_CMO
#include "stdalloc.fh"

Type V2
   Real*8, Pointer:: A(:,:)=>Null()
End Type V2

Type CMO_Type
 Integer:: nSym=0
 Real*8, Allocatable :: CMO_Full(:)
 Type (V2):: pA(8)
End Type CMO_Type

Contains
  Subroutine Deallocate_CMO(Adam)
  Implicit None
  Type (CMO_Type) Adam
  Integer iSym

  Do iSym = 1, Adam%nSym
     Adam%pA(iSym)%A => Null()
  End Do
  Call mma_deallocate(Adam%CMO_Full)
  Adam%nSym=0

  End Subroutine Deallocate_CMO

  Subroutine Allocate_CMO(Adam,n,m,nSym)
  Implicit None
  Type (CMO_Type),Target:: Adam
  Integer nSym
  Integer n(nSym), m(nSym)
  Integer iE, iS, iSym, MemTot

  Adam%nSym=nSym

  MemTot=0
  Do iSym = 1, nSym
     MemTot = MemTot + n(iSym)*m(iSym)
  End Do
  Call mma_allocate(Adam%CMO_Full,MemTot,Label='%CMO_Full')

  iE = 0
  Do iSym = 1, nSym
    iS = iE + 1
    iE = iE + n(iSym) * m(iSym)

    Adam%pA(iSym)%A(1:n(iSym),1:m(iSym)) => Adam%CMO_Full(iS:iE)
  End Do
  End Subroutine Allocate_CMO

  Subroutine Map_to_CMO(Adam,ipAdam)
  Implicit None
  Type (CMO_Type):: Adam
  Integer ipAdam(*)
  Integer, External:: ip_of_Work
  Integer iSym

  Do iSym=1, Adam%nSym
     ipAdam(iSym) = ip_of_Work(Adam%pA(iSym)%A(1,1))
  End Do

  End Subroutine Map_to_CMO

End Module Data_Structures
