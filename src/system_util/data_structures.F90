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
Public:: twxy_Type, Allocate_twxy, Deallocate_twxy, Map_to_twxy
#include "stdalloc.fh"
#include "real.fh"

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
  Type (V2):: pA2(8)
End Type Laq_type

Type twxy_type
  Integer :: iCase=0
  Integer :: JSYM=0
  Integer :: nSym=0
  Real*8, Allocatable:: twxy_full(:)
  Type (V2):: pA(8,8)
End Type twxy_type


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
  Integer iE, iS, iSyma, iSymb, MemTot, n2Dim

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
    Case(4)
       Do iSyma = 1, nSym
          If (n(iSyma)/=m(iSyma)) Then
             Write (6,*) 'Allocate_Laq: iSwap=4 only valid if n(:)=m(:).'
             Call abend()
          End If
          iSymb = MulD2h(iSym,iSyma)
          If (iSyma==iSymb) Then
            n2Dim = n(iSyma)*(n(iSyma)+1)/2
          Else
            n2Dim = n(iSyma)*n(iSymb)
          End If
          MemTot = MemTot + n2dim*NUMV
       End Do
    Case(5)
       Do iSyma = 1, nSym
          If (n(iSyma)/=m(iSyma)) Then
             Write (6,*) 'Allocate_Laq: iSwap=5 only valid if n(:)=m(:).'
             Call abend()
          End If
          iSymb = MulD2h(iSym,iSyma)
          If (iSyma==iSymb) Then
            n2Dim = n(iSyma)*(n(iSyma)+1)/2
          Else If (iSymb>iSyma) Then
            n2Dim = n(iSyma)*n(iSymb)
          Else
            n2Dim = 0
          End If
          MemTot = MemTot + n2dim*NUMV
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
    Case(4)
       Do iSyma = 1, nSym
          iSymb = MulD2h(iSym,iSyma)
          iS = iE + 1
          If (iSyma==iSymb) Then
            n2Dim = n(iSyma)*(n(iSyma)+1)/2
          Else
            n2Dim = n(iSyma)*n(iSymb)
          End If
          iE = iE + n2Dim*NUMV
          Adam%pA2(iSymb)%A(1:n2Dim,1:NUMV) => Adam%Laq_Full(iS:iE)
       End Do
    Case(5)
       Do iSyma = 1, nSym
          iSymb = MulD2h(iSym,iSyma)
          If (iSymb>iSyma) Cycle
          iS = iE + 1
          If (iSyma==iSymb) Then
            n2Dim = n(iSyma)*(n(iSyma)+1)/2
          Else
            n2Dim = n(iSyma)*n(iSymb)
          End If
          iE = iE + n2Dim*NUMV
          Adam%pA2(iSymb)%A(1:n2Dim,1:NUMV) => Adam%Laq_Full(iS:iE)
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

  If (Adam%iSwap==4) Then
     Do iSym = 1, Adam%nSym
        Adam%pA2(iSym)%A => Null()
     End Do
  Else
     Do iSym = 1, Adam%nSym
        Adam%pA(iSym)%A => Null()
     End Do
  End If
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
  Integer iSym,jSym

  Integer i, j, MulD2h
  MulD2h(i,j) = iEOR(i-1,j-1) + 1

  If (Adam%iSwap<4) Then
     Do iSym=1, Adam%nSym
        ipAdam(iSym) = ip_of_Work(Adam%pA(iSym)%A(1,1,1))
     End Do
  Else
     Do iSym=1, Adam%nSym
        jsym=MulD2h(iSym,Adam%iSym)
        ipAdam(jSym) = ip_of_Work(Adam%pA2(jSym)%A(1,1))
     End Do
  End If

  End Subroutine Map_to_Laq



Subroutine Allocate_twxy(twxy,n,m,JSYM,nSym,iCase)
Implicit None
Type (twxy_type), Target:: twxy
Integer JSYM, nSym, iCase
Integer n(nSym), m(nSym)
Integer iSymx, iSymy, iSymt, iSymw
Integer iSyma
Integer mtwxy
Integer iS, iE, n1, n2

Integer i, j, MulD2h
MulD2h(i,j) = iEOR(i-1,j-1) + 1

twxy%iCase=iCase
twxy%JSYM=JSYM
twxy%nSym=nSym
! *** memory for the (tw|xy) integrals --- temporary array
mtwxy = 0
Select Case (iCase)
Case (0)  ! twxy
  Do iSymy=1,nSym
     If (n(iSymy)/=m(iSymy)) Then
        Write (6,*) 'Allocate_twxy: iCase=0 only valid if n(:)=m(:).'
        Call abend()
     End If
     iSymx=MulD2h(iSymy,JSYM)
     n2 = n(iSymx)*n(iSymy)
     Do iSymw=iSymy,nSym    ! iSymw.ge.iSymy (particle symmetry)
        iSymt=MulD2h(isymw,JSYM)
        n1 = n(iSymt)*n(iSymw)
        mtwxy = mtwxy + n1 * n2
     End Do
  End Do
Case (1) ! waxy
  Do iSymy=1,nSym
     iSymx=MulD2h(iSymy,JSYM)
     n2=n(iSymx)*n(iSymy)
     If (iSymx==iSymy) n2=n(iSymx)*(n(iSymx)+1)/2
     If (iSymx.le.iSymy) then
        Do iSyma=1,nSym
           iSymw=MulD2h(iSyma,JSYM)
           n1=n(iSymw)*m(iSyma)
           mtwxy = mtwxy + n1 * n2
        End Do
     End If
  End Do
Case (2)  ! twxy
  Do iSymy=1,nSym
     If (n(iSymy)/=m(iSymy)) Then
        Write (6,*) 'Allocate_twxy: iCase=2 only valid if n(:)=m(:).'
        Call abend()
     End If
     iSymx=MulD2h(iSymy,JSYM)
     n2=n(iSymx)*n(iSymy)
     If (iSymx==iSymy) n2=n(iSymx)*(n(iSymx)+1)/2
     If (iSymx.ge.iSymy) then
        Do iSymw=iSymy,nSym
           iSymt=MulD2h(iSymw,JSYM)
           If (iSymt.ge.iSymw) Then
              n1=n(iSymt)*n(iSymw)
              If (iSymt==iSymw) n1=n(iSymt)*(n(iSymt)+1)/2
              mtwxy = mtwxy + n1 * n2
           End If
        End Do
     End If
  End Do
Case Default
  Write (6,*) "Allocate_twxy: Illegal case."
  Call Abend()
End Select

Call mma_allocate(twxy%twxy_full,mtwxy,Label='twxy')
twxy%twxy_full(:)=Zero

! *** setup pointers to the symmetry blocks of (tw|xy)

iE = 0
Select Case (iCase)
Case (0)
  Do iSymy=1,nSym
     iSymx=MulD2h(iSymy,JSYM)
     n2 = n(iSymx)*n(iSymy)
     Do iSymw=iSymy,nSym   ! iSymw.ge.iSymy (particle symmetry)
        iSymt=MulD2h(isymw,JSYM)
        n1 = n(iSymt)*n(iSymw)
        iS = iE + 1
        iE = iE + n1*n2
        twxy%pA(iSymw,iSymy)%A(1:n1,1:n2) => twxy%twxy_full(iS:iE)
     End Do
  End Do
Case (1)
  Do iSymy=1,nSym
     iSymx=MulD2h(iSymy,JSYM)
     n2=n(iSymx)*n(iSymy)
     If (iSymx==iSymy) n2=n(iSymx)*(n(iSymx)+1)/2
     If (iSymx.le.iSymy) then
        Do iSyma=1,nSym
           iSymw=MulD2h(iSyma,JSYM)
           n1 =  n(iSymw)*m(iSyma)
           iS = iE + 1
           iE = iE + n1*n2
           twxy%pA(iSymw,iSymx)%A(1:n1,1:n2) => twxy%twxy_full(iS:iE)
        End Do
     End If
  End Do
Case (2) ! twxy
  Do iSymy=1,nSym
     iSymx=MulD2h(iSymy,JSYM)
     n2=n(iSymx)*n(iSymy)
     If (iSymx==iSymy) n2=n(iSymx)*(n(iSymx)+1)/2
     If (iSymx.ge.iSymy) then
        Do iSymw=iSymy,nSym ! iSymw.ge.iSymy
           iSymt=MulD2h(iSymw,JSYM)
           If (iSymt.ge.iSymw) Then
              n1=n(iSymt)*n(iSymw)
              If (iSymt==iSymw) n1=n(iSymt)*(n(iSymt)+1)/2
              iS = iE + 1
              iE = iE + n1 * n2
              twxy%pA(iSymw,iSymy)%A(1:n1,1:n2) => twxy%twxy_full(iS:iE)
              twxy%pA(iSymy,iSymw)%A(1:n1,1:n2) => twxy%twxy_full(iS:iE) ! symmetrization
           End If
        End Do
     End If
  End Do
End Select

End Subroutine Allocate_twxy

Subroutine Deallocate_twxy(twxy)
Implicit None
Type (twxy_type) twxy
Integer iSymy, iSymw

Integer i, j, MulD2h
MulD2h(i,j) = iEOR(i-1,j-1) + 1

Call mma_deallocate(twxy%twxy_full)

! *** setup pointers to the symmetry blocks of (tw|xy)

Do iSymy=1,8
   Do iSymw=1,8
      twxy%pA(iSymw,iSymy)%A => Null()
   End Do
End Do

End Subroutine Deallocate_twxy

Subroutine Map_to_twxy(Adam,ipAdam)
Implicit None
Type (twxy_type):: Adam
Integer ipAdam(8,8)
Integer, External:: ip_of_Work
Integer iSymx, iSymy, iSymt, iSymw
Integer iSyma

Integer i, j, MulD2h
MulD2h(i,j) = iEOR(i-1,j-1) + 1

ipAdam(:,:)=0
Select Case (Adam%iCase)
Case (0)
  Do iSymy=1,Adam%nSym
     iSymx=MulD2h(iSymy,Adam%JSYM)
     Do iSymw=iSymy,Adam%nSym   ! iSymw.ge.iSymy (particle symmetry)
        iSymt=MulD2h(isymw,Adam%JSYM)
        ipAdam(iSymw,iSymy) = ip_of_Work(Adam%pA(iSymw,iSymy)%A(1,1))
     End Do
  End Do
Case (1)
  Do iSymy=1,Adam%nSym
     iSymx=MulD2h(iSymy,Adam%JSYM)
     If (iSymx.le.iSymy) then
        Do iSyma=1,Adam%nSym
           iSymw=MulD2h(iSyma,Adam%JSYM)
           ipAdam(iSymw,iSymx) = ip_of_Work(Adam%pA(iSymw,iSymx)%A(1,1))
        End Do
     End If
  End Do
Case (2)
  Do iSymy=1,Adam%nSym
     iSymx=MulD2h(iSymy,Adam%JSYM)
     If (iSymx.ge.iSymy) then
        Do iSymw=iSymy,Adam%nSym ! iSymw.ge.iSymy
           iSymt=MulD2h(iSymw,Adam%JSYM)
           If (iSymt.ge.iSymw) Then
              ipAdam(iSymw,iSymy) = ip_of_Work(Adam%pA(iSymw,iSymy)%A(1,1))
              ipAdam(iSymy,iSymw) = ip_of_Work(Adam%pA(iSymw,iSymy)%A(1,1))
           End If
        End Do
     End If
  End Do
End Select

End Subroutine Map_to_twxy

End Module Data_Structures
