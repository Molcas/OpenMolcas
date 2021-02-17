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
Public:: CMO_Type, Deallocate_CMO
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

End Module Data_Structures
