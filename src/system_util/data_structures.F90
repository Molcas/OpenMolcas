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
Public:: CMO_Type, Allocate_CMO, Deallocate_CMO, Map_to_CMO
Public:: Laq_Type, Allocate_Laq, Deallocate_Laq, Map_to_Laq
#include "stdalloc.fh"

Type V2
  Real*8, Pointer:: A(:,:)=>Null()
End Type V2

Type CMO_Type
  Integer:: nSym=0
  Real*8, Allocatable :: CMO_Full(:)
  Type (V2):: pA(8)
End Type CMO_Type

Type V3
  Real*8, Pointer:: A(:,:,:)=>Null()
End Type V3

Type Laq_type
  Integer:: iSwap=0
  Integer:: iSym=0
  Integer:: nSym=0
  Real*8, Allocatable :: Laq_Full(:)
  Type (V3):: pA(8)
End Type Laq_type


Contains


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


  Subroutine Allocate_Laq(Adam,n,m,NUMV,iSym,nSym,iSwap)
  Implicit None
  Type (Laq_Type),Target:: Adam
  Integer NUMV
  Integer iSym
  Integer nSym
  Integer n(nSym), m(nSym)
  Integer iSwap
  Integer iE, iS, iSyma, iSymb, MemTot

  Integer i, j, MulD2h
  MulD2h(i,j) = iEOR(i-1,j-1) + 1

  Adam%iSym=iSym
  Adam%nSym=nSym
  Adam%iSwap=iSwap

  MemTot=0

  Select Case (iSwap)
    Case(0)
       Do iSyma = 1, nSym
          iSymb = MulD2h(iSym,iSyma)
          MemTot = MemTot + n(iSyma)*m(iSymb)*NUMV
       End Do
    Case(1)
       Do iSyma = 1, nSym
          iSymb = MulD2h(iSym,iSyma)
          MemTot = MemTot + m(iSyma)*n(iSymb)*NUMV
       End Do
    Case(2)
       Do iSyma = 1, nSym
          iSymb = MulD2h(iSym,iSyma)
          MemTot = MemTot + n(iSyma)*NUMV*m(iSymb)
       End Do
    Case(3)
       Do iSyma = 1, nSym
          iSymb = MulD2h(iSym,iSyma)
          MemTot = MemTot + m(iSyma)*NUMV*n(iSymb)
       End Do
    Case Default
       Write (6,*) "Allocate_Laq: Illegal case."
       Call Abend()
  End Select

  Call mma_allocate(Adam%Laq_Full,MemTot,Label='%Laq_Full')


  iE = 0

  Select Case (iSwap)
    Case(0)
       Do iSyma = 1, nSym
          iSymb = MulD2h(iSym,iSyma)
          iS = iE + 1
          iE = iE + n(iSyma)*m(iSymb)*NUMV
          Adam%pA(iSyma)%A(1:n(iSyma),1:m(iSymb),1:NUMV) => Adam%Laq_Full(iS:iE)
       End Do
    Case(1)
       Do iSyma = 1, nSym
          iSymb = MulD2h(iSym,iSyma)
          iS = iE + 1
          iE = iE + m(iSyma)*n(iSymb)*NUMV
          Adam%pA(iSyma)%A(1:m(iSyma),1:n(iSymb),1:NUMV) => Adam%Laq_Full(iS:iE)
       End Do
    Case(2)
       Do iSyma = 1, nSym
          iSymb = MulD2h(iSym,iSyma)
          iS = iE + 1
          iE = iE + n(iSyma)*NUMV*m(iSymb)
          Adam%pA(iSyma)%A(1:n(iSyma),1:NUMV,1:m(iSymb)) => Adam%Laq_Full(iS:iE)
       End Do
    Case(3)
       Do iSyma = 1, nSym
          iSymb = MulD2h(iSym,iSyma)
          iS = iE + 1
          iE = iE + m(iSyma)*NUMV*n(iSymb)
          Adam%pA(iSyma)%A(1:m(iSyma),1:NUMV,1:n(iSymb)) => Adam%Laq_Full(iS:iE)
       End Do
    Case Default
       Write (6,*) "Allocate_Laq: Illegal case."
       Call Abend()
  End Select
  End Subroutine Allocate_Laq


  Subroutine Deallocate_Laq(Adam)
  Implicit None
  Type (Laq_Type) Adam
  Integer iSym

  Do iSym = 1, Adam%nSym
     Adam%pA(iSym)%A => Null()
  End Do
  Call mma_deallocate(Adam%Laq_Full)
  Adam%iSwap=0
  Adam%iSym=0
  Adam%nSym=0

  End Subroutine Deallocate_Laq

  Subroutine Map_to_Laq(Adam,ipAdam)
  Implicit None
  Type (Laq_Type):: Adam
  Integer ipAdam(*)
  Integer, External:: ip_of_Work
  Integer iSym

  Do iSym=1, Adam%nSym
     ipAdam(iSym) = ip_of_Work(Adam%pA(iSym)%A(1,1,1))
  End Do

  End Subroutine Map_to_Laq

End Module Data_Structures
