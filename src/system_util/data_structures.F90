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
Public:: DSBA_Type,   Allocate_DSBA,   Deallocate_DSBA, Map_to_DSBA
Public:: SBA_Type,    Allocate_SBA,    Deallocate_SBA,  Map_to_SBA
Public:: twxy_Type,   Allocate_twxy,   Deallocate_twxy, Map_to_twxy
Public:: NDSBA_Type,  Allocate_NDSBA,  Deallocate_NDSBA
Public:: G2_Type,     Allocate_G2,     Deallocate_G2
Public:: L_Full_Type, Allocate_L_Full, Deallocate_L_Full
Public:: Lab_Type, Allocate_Lab, Deallocate_Lab
#include "stdalloc.fh"
#include "real.fh"

Type SB_Type
  Real*8, Pointer:: A3(:,:,:)=>Null()
  Real*8, Pointer:: A2(:,:)=>Null()
  Real*8, Pointer:: A1(:)=>Null()
End Type  SB_Type

Type DSB_Type
  Real*8, Pointer:: A2(:,:)=>Null()
  Real*8, Pointer:: A1(:)=>Null()
End Type  DSB_Type

Type V1
  Real*8, Pointer:: A(:)=>Null()
End Type V1

Type V2
  Real*8, Pointer:: A(:,:)=>Null()
End Type V2

Type G2_pointers
  Real*8, Pointer:: A4(:,:,:,:)=>Null()
  Real*8, Pointer:: A2(:,:)=>Null()
End Type G2_pointers

Type L_Full_Pointers
  Real*8, Pointer :: A3(:,:,:)=>Null()
  Real*8, Pointer :: A21(:,:)=>Null()
  Real*8, Pointer :: A12(:,:)=>Null()
End Type L_Full_Pointers


Type NDSBA_Type
  Integer:: iCase=0
  Integer:: nSym=0
  Real*8, Allocatable :: A0(:)
  Type (DSB_Type):: SB(8,8)
End Type NDSBA_Type


Type DSBA_Type
  Integer:: iCase=0
  Integer:: nSym=0
  Logical:: Fake=.False.
  Logical:: Active=.False.
  Real*8, Pointer :: A0(:)
  Type (DSB_Type):: SB(8)
End Type DSBA_Type

Type SBA_Type
  Integer:: iCase=0
  Integer:: iSym=0
  Integer:: nSym=0
  Real*8, Allocatable :: A0(:)
  Type (SB_Type):: SB(8)
End Type SBA_Type

Type twxy_type
  Integer :: iCase=0
  Integer :: JSYM=0
  Integer :: nSym=0
  Real*8, Allocatable:: twxy_full(:)
  Type (V2):: SB(8,8)
End Type twxy_type

Type G2_type
  Integer :: iCase=0
  Integer :: nSym=0
  Real*8, Allocatable:: A0(:)
  Type (G2_Pointers):: SB(8,8,8)
End Type G2_type

Type L_Full_Type
  Integer :: iCase=0
  Integer :: iSym=0
  Integer :: nSym=0
  Integer :: nShell=0
  Real*8, Allocatable:: A0(:)
  Type (L_Full_Pointers), Allocatable :: SPB(:,:,:)
End Type L_Full_Type

Type Lab_Type
  Integer :: nSym=0
  Integer :: nDen=0
  Integer :: nShell=0
  Real*8, Allocatable:: A0(:)
  Logical, Allocatable:: Keep(:,:)
  Type (V1), Allocatable :: SB(:,:,:)
End Type Lab_Type

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     !
!                N D S B A - T Y P E   S E C T I O N                  !
!                                                                     !
!                Non-Diagonal Symmetry Blocked Arrays                 !
!                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine Allocate_NDSBA(Adam,n,m,nSym)
  Implicit None
  Type (NDSBA_Type),Target, Intent(Out) :: Adam
  Integer, Intent(In) :: nSym
  Integer, Intent(In) :: n(nSym), m(nSym)

  Integer iE, iS, iSym, jSym, MemTot

  Adam%iCase=1
  Adam%nSym=nSym

  MemTot=0
  Do jSym = 1, nSym
     Do iSym = jSym, nSym
        MemTot = MemTot + n(iSym)*m(jSym)
     End Do
  End Do
  Call mma_allocate(Adam%A0,MemTot,Label='%A0')

  iE = 0
  Do jSym = 1, nSym
     Do iSym = jSym, nSym
        iS = iE + 1
        iE = iE + n(iSym) * m(jSym)

        Adam%SB(jSym,iSym)%A2(1:n(iSym),1:m(jSym)) => Adam%A0(iS:iE)
        Adam%SB(jSym,iSym)%A1(1:n(iSym)*m(jSym))   => Adam%A0(iS:iE)
        Adam%SB(iSym,jSym)%A2(1:n(iSym),1:m(jSym)) => Adam%A0(iS:iE)
        Adam%SB(iSym,jSym)%A1(1:n(iSym)*m(jSym))   => Adam%A0(iS:iE)
     End Do
  End Do
  End Subroutine Allocate_NDSBA


  Subroutine Deallocate_NDSBA(Adam)
  Implicit None
  Type (NDSBA_Type) Adam
  Integer iSym, jSym

    Do iSym = 1, Adam%nSym
       Do jSym = 1, Adam%nSym
          Adam%SB(iSym,jSym)%A2=> Null()
          Adam%SB(iSym,jSym)%A1=> Null()
       End Do
    End Do
  Call mma_deallocate(Adam%A0)
  Adam%nSym=0
  Adam%iCase=0

  End Subroutine Deallocate_NDSBA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     !
!                D S B A - T Y P E   S E C T I O N                    !
!                                                                     !
!                Diagonal Symmetry Blocked Arrays                     !
!                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine Allocate_DSBA(Adam,n,m,nSym,Case,Ref)
  Implicit None
  Type (DSBA_Type),Target, Intent(Out) :: Adam
  Integer, Intent(In) :: nSym
  Integer, Intent(In) :: n(nSym), m(nSym)
  Character(LEN=3), Intent(In), Optional :: Case
  Real*8, Target, Optional :: Ref(*)

  Integer iE, iS, iSym, MemTot
  Integer :: iCase=0

  If (Present(Case)) Then
     Select Case (Case)
      Case ('TRI')
        iCase=2
        Do iSym = 1, nSym
           If (n(iSym)/=m(iSym)) Then
             Write (6,*) 'Allocate_DSBA: n(iSym)/=m(iSym), illegal if CASE="TRI".'
             Call Abend()
           End If
        End Do
      Case ('REC')
        iCase=1
      Case Default
        Write (6,*) 'Allocate_DSBA: Illegal Case parameter, Case=',Case
        Write (6,*) 'Allowed value are "TRI" and "REC".'
        Call Abend()
     End Select
  Else
    iCase=1
  End If
  Adam%iCase=iCase
  Adam%nSym=nSym

  MemTot=0
  If (iCase==1) Then
    Do iSym = 1, nSym
       MemTot = MemTot + n(iSym)*m(iSym)
    End Do
  Else
    Do iSym = 1, nSym
       MemTot = MemTot + n(iSym)*(n(iSym)+1)/2
    End Do
  End If
  If (Present(Ref)) Then
    Adam%Fake=.True.
    Adam%A0(1:MemTot) => Ref(1:MemTot)
  Else
    Adam%A0=>Null()
    Call mma_allocate(Adam%A0,MemTot,Label='%A0')
  End If

  Adam%Active=.True.
  iE = 0
  If (iCase==1) Then
    Do iSym = 1, nSym
      iS = iE + 1
      iE = iE + n(iSym) * m(iSym)

      Adam%SB(iSym)%A2(1:n(iSym),1:m(iSym)) => Adam%A0(iS:iE)
      Adam%SB(iSym)%A1(1:n(iSym)*m(iSym))   => Adam%A0(iS:iE)
    End Do
  Else
    Do iSym = 1, nSym
      iS = iE + 1
      iE = iE + n(iSym) * (n(iSym)+1)/2

      Adam%SB(iSym)%A1(1:n(iSym)*(n(iSym)+1)/2)   => Adam%A0(iS:iE)
    End Do
  End If
End Subroutine Allocate_DSBA



  Subroutine Deallocate_DSBA(Adam)
  Implicit None
  Type (DSBA_Type) Adam
  Integer iSym

  If (.NOT.Adam%Active) Return

  Adam%Active=.False.
  If (Adam%iCase==1) Then
    Do iSym = 1, Adam%nSym
       Adam%SB(iSym)%A2=> Null()
       Adam%SB(iSym)%A1=> Null()
    End Do
  Else
    Do iSym = 1, Adam%nSym
       Adam%SB(iSym)%A1=> Null()
    End Do
  End If
  If (Adam%Fake) Then
    Adam%A0=>Null()
    Adam%Fake=.False.
  Else
    Call mma_deallocate(Adam%A0)
  End If
  Adam%nSym=0
  Adam%iCase=0

  End Subroutine Deallocate_DSBA

  Subroutine Map_to_DSBA(Adam,ipAdam)
  Implicit None
  Type (DSBA_Type):: Adam
  Integer ipAdam(*)
  Integer, External:: ip_of_Work
  Integer iSym

  Do iSym=1, Adam%nSym
     ipAdam(iSym) = ip_of_Work(Adam%SB(iSym)%A1(1))
  End Do

  End Subroutine Map_to_DSBA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     !
!                  S B A - T Y P E   S E C T I O N                    !
!                                                                     !
!                    Symmetry Blocked Arrays                          !
!                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Subroutine Allocate_SBA(Adam,n,m,NUMV,iSym,nSym,iCase,Memory)
  Implicit None
  Type (SBA_Type),Target:: Adam
  Integer NUMV
  Integer iSym
  Integer nSym
  Integer n(nSym), m(nSym)
  Integer iCase
  Integer, Optional :: Memory


  Integer iE, iS, iSyma, iSymb, MemTot, n2Dim, n3Dim

  Integer i, j, MulD2h
  MulD2h(i,j) = iEOR(i-1,j-1) + 1

  MemTot=0

  Select Case (iCase)
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
             Write (6,*) 'Allocate_SBA: iCase=4 only valid if n(:)=m(:).'
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
             Write (6,*) 'Allocate_SBA: iCase=5 only valid if n(:)=m(:).'
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
    Case(6)
       Do iSyma = 1, nSym
          If (n(iSyma)/=m(iSyma)) Then
             Write (6,*) 'Allocate_SBA: iCase=5 only valid if n(:)=m(:).'
             Call abend()
          End If
          iSymb = MulD2h(iSym,iSyma)
          If (iSyma>=iSymb) Then
            n2Dim = n(iSyma)*n(iSymb)
          Else
            n2Dim = 0
          End If
          MemTot = MemTot + n2dim*NUMV
       End Do
    Case Default
       Write (6,*) "Allocate_SBA: Illegal case."
       Call Abend()
  End Select

  If (Present(Memory)) Then
     Memory=MemTot
     Return
  End If

  Adam%iSym=iSym
  Adam%nSym=nSym
  Adam%iCase=iCase


  Call mma_allocate(Adam%A0,MemTot,Label='%A0')


  iE = 0

  Select Case (iCase)
    Case(0)
       Do iSyma = 1, nSym
          iSymb = MulD2h(iSym,iSyma)
          iS = iE + 1
          iE = iE + n(iSyma)*m(iSymb)*NUMV
          n2Dim = n(iSyma)*m(iSymb)
          n3Dim = n2Dim*NUMV
          Adam%SB(iSyma)%A1(1:n3Dim) => Adam%A0(iS:iE)
          Adam%SB(iSyma)%A2(1:n2Dim,1:NUMV) => Adam%A0(iS:iE)
          Adam%SB(iSyma)%A3(1:n(iSyma),1:m(iSymb),1:NUMV) => Adam%A0(iS:iE)
       End Do
    Case(1)
       Do iSyma = 1, nSym
          iSymb = MulD2h(iSym,iSyma)
          iS = iE + 1
          iE = iE + m(iSyma)*n(iSymb)*NUMV
          n2Dim = m(iSyma)*n(iSymb)
          n3Dim = n2Dim*NUMV
          Adam%SB(iSyma)%A1(1:n3Dim) => Adam%A0(iS:iE)
          Adam%SB(iSyma)%A2(1:n2Dim,1:NUMV) => Adam%A0(iS:iE)
          Adam%SB(iSyma)%A3(1:m(iSyma),1:n(iSymb),1:NUMV) => Adam%A0(iS:iE)
       End Do
    Case(2)
       Do iSyma = 1, nSym
          iSymb = MulD2h(iSym,iSyma)
          iS = iE + 1
          iE = iE + n(iSyma)*NUMV*m(iSymb)
          Adam%SB(iSyma)%A3(1:n(iSyma),1:NUMV,1:m(iSymb)) => Adam%A0(iS:iE)
       End Do
    Case(3)
       Do iSyma = 1, nSym
          iSymb = MulD2h(iSym,iSyma)
          iS = iE + 1
          iE = iE + m(iSyma)*NUMV*n(iSymb)
          Adam%SB(iSyma)%A3(1:m(iSyma),1:NUMV,1:n(iSymb)) => Adam%A0(iS:iE)
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
          Adam%SB(iSymb)%A2(1:n2Dim,1:NUMV) => Adam%A0(iS:iE)
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
          Adam%SB(iSymb)%A2(1:n2Dim,1:NUMV) => Adam%A0(iS:iE)
       End Do
    Case(6)
       Do iSyma = 1, nSym
          iSymb = MulD2h(iSym,iSyma)
          If (iSymb>iSyma) Cycle
          iS = iE + 1
          n2Dim = n(iSyma)*n(iSymb)
          iE = iE + n2Dim*NUMV
          Adam%SB(iSymb)%A2(1:n2Dim,1:NUMV) => Adam%A0(iS:iE)
       End Do
    Case Default
       Write (6,*) "Allocate_SBA: Illegal case."
       Call Abend()
  End Select
  End Subroutine Allocate_SBA


  Subroutine Deallocate_SBA(Adam)
  Implicit None
  Type (SBA_Type) Adam
  Integer iSym

  Do iSym = 1, Adam%nSym
      Adam%SB(iSym)%A1 => Null()
      Adam%SB(iSym)%A2 => Null()
      Adam%SB(iSym)%A3 => Null()
  End Do
  Call mma_deallocate(Adam%A0)
  Adam%iCase=0
  Adam%iSym=0
  Adam%nSym=0

  End Subroutine Deallocate_SBA

  Subroutine Map_to_SBA(Adam,ipAdam,Tweak)
  Implicit None
  Type (SBA_Type):: Adam
  Integer ipAdam(*)
  Logical, Optional :: Tweak


  Integer, External:: ip_of_Work
  Integer iSym,jSym
  Logical :: Swap=.False.

  Integer i, j, MulD2h
  MulD2h(i,j) = iEOR(i-1,j-1) + 1

  If (Adam%iCase<4) Then
     Do iSym=1, Adam%nSym
        ipAdam(iSym) = ip_of_Work(Adam%SB(iSym)%A3(1,1,1))
     End Do
  Else
     If (Present(Tweak)) Swap=Tweak
     If (Swap) Then
        Do iSym=1, Adam%nSym
           jsym=MulD2h(iSym,Adam%iSym)
           If (.NOT.Associated(Adam%SB(jSym)%A2)) Cycle

           ipAdam(iSym) = ip_of_Work(Adam%SB(jSym)%A2(1,1))

        End Do
     Else
        Do iSym=1, Adam%nSym
           If (.NOT.Associated(Adam%SB(iSym)%A2)) Cycle

           ipAdam(iSym) = ip_of_Work(Adam%SB(iSym)%A2(1,1))
        End Do
     End If
  End If

  End Subroutine Map_to_SBA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     !
!                  T W X Y - T Y P E   S E C T I O N                  !
!                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
        twxy%SB(iSymw,iSymy)%A(1:n1,1:n2) => twxy%twxy_full(iS:iE)
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
           twxy%SB(iSymw,iSymx)%A(1:n1,1:n2) => twxy%twxy_full(iS:iE)
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
              twxy%SB(iSymw,iSymy)%A(1:n1,1:n2) => twxy%twxy_full(iS:iE)
              twxy%SB(iSymy,iSymw)%A(1:n1,1:n2) => twxy%twxy_full(iS:iE) ! symmetrization
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

Call mma_deallocate(twxy%twxy_full)

! *** setup pointers to the symmetry blocks of (tw|xy)

Do iSymy=1,8
   Do iSymw=1,8
      twxy%SB(iSymw,iSymy)%A => Null()
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
        ipAdam(iSymw,iSymy) = ip_of_Work(Adam%SB(iSymw,iSymy)%A(1,1))
     End Do
  End Do
Case (1)
  Do iSymy=1,Adam%nSym
     iSymx=MulD2h(iSymy,Adam%JSYM)
     If (iSymx.le.iSymy) then
        Do iSyma=1,Adam%nSym
           iSymw=MulD2h(iSyma,Adam%JSYM)
           ipAdam(iSymw,iSymx) = ip_of_Work(Adam%SB(iSymw,iSymx)%A(1,1))
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
              ipAdam(iSymw,iSymy) = ip_of_Work(Adam%SB(iSymw,iSymy)%A(1,1))
              ipAdam(iSymy,iSymw) = ip_of_Work(Adam%SB(iSymw,iSymy)%A(1,1))
           End If
        End Do
     End If
  End Do
End Select

End Subroutine Map_to_twxy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     !
!                      G 2 - T Y P E   S E C T I O N                  !
!                                                                     !
!              Symmetry block 2-particle-like arrays                  !
!                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine Allocate_G2(Adam,n,nSym,iCase)
Implicit None
Type (G2_Type), Target:: Adam
Integer nSym, iCase
Integer n(nSym)

Integer MemTot, ijSym, iSym, jSym, kSym, lSym, iE, iS, n1, n2, n3, n4, n12, n34

Adam%nSym=nSym
Adam%iCase=iCase

MemTot=0
Select Case (iCase)

Case(1)

Do ijsym=1,nsym
   Do isym=1,nsym
      jsym=iEOR(isym-1,ijsym-1)+1
      n12 = n(iSym)*n(jSym)
      Do kSym=1,nSym
         lSym=iEOR(kSym-1,ijSym-1)+1
         n34=n(kSym)*n(lSym)
         MemTot=MemTot+n12*n34
      End Do
   End Do
End Do

Case Default

Write (6,*) 'Allocate_G2: illegal case valeu=',iCase
Call Abend()

End Select

Call mma_allocate(Adam%A0,MemTot,Label='G2%A0')

iE = 0
Select Case (iCase)

Case(1)

Do ijsym=1,nsym
   Do isym=1,nsym
      jsym=iEOR(isym-1,ijsym-1)+1
      n1=n(iSym)
      n2=n(jSym)
      n12 = n1*n2
      Do kSym=1,nSym
         lSym=iEOR(kSym-1,ijSym-1)+1
         n3=n(kSym)
         n4=n(lSym)
         n34=n3*n4
         iS = iE + 1
         iE = iE + n12*n34
         Adam%SB(iSym,jSym,kSym)%A4(1:n1,1:n2,1:n3,1:n4) => Adam%A0(iS:iE)
         Adam%SB(iSym,jSym,kSym)%A2(1:n12,1:n34) => Adam%A0(iS:iE)
      End Do
   End Do
End Do

Case Default

Write (6,*) 'What?'
Call Abend()

End Select

End Subroutine Allocate_G2

Subroutine Deallocate_G2(Adam)
Implicit None
Type (G2_Type) Adam

Integer iSym, jSym, kSym

Adam%iCase=0

Call mma_deallocate(Adam%A0)

Do iSym=1,Adam%nSym
   Do jSym=1,Adam%nSym
      Do kSym=1,Adam%nSym
         Adam%SB(iSym,jSym,kSym)%A4=>Null()
         Adam%SB(iSym,jSym,kSym)%A2=>Null()
      End Do
   End Do
End Do
Adam%nSym=0
End Subroutine Deallocate_G2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     !
!                      L F u l l - T Y P E   S E C T I O N            !
!                                                                     !
!                      L  full storage shell-pair blocked             !
!                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine Allocate_L_Full(Adam,nShell,iShp_rs,JNUM,JSYM,nSym, Memory)
use ChoArr, only: nBasSh
use ChoSwp, only: nnBstRSh
Implicit None
Type (L_Full_Type), Target:: Adam
Integer nShell
Integer iShp_rs( nShell*(nShell+2)/2 )
Integer JNUM, JSYM, nSym
Integer, Optional, Intent(Out):: Memory

Integer iaSh, ibSh, iShp
Integer iSyma, iSymb
Integer LFULL
Integer iS, iE
Integer n1, n2

Integer i, j, MulD2h
MulD2h(i,j) = iEOR(i-1,j-1) + 1

LFULL=0
Do iaSh=1,nShell
   Do ibSh=1,iaSh
      iShp = iaSh*(iaSh-1)/2 + ibSh

      If (iShp_rs(iShp)<=0) Cycle

      If (nnBstRSh(Jsym,iShp_rs(iShp),1)<=0) Cycle

      Do iSymb=1,nSym
         iSyma=MulD2h(iSymb,Jsym)
         If (iSyma<iSymb) Cycle

         LFULL = LFULL + nBasSh(iSyma,iaSh)*nBasSh(iSymb,ibSh)
         If (iaSh==ibSh) Cycle

         LFULL = LFULL + nBasSh(iSyma,ibSh)*nBasSh(iSymb,iaSh)

      End Do

   End Do
End Do
LFULL=LFULL*JNUM

If (Present(Memory)) Then
   Memory=LFULL
   Return
End If

Adam%iCase=1
Adam%nSym=nSym
Adam%iSym=JSYM
Adam%nShell=nShell

Call mma_Allocate(Adam%A0,LFULL,Label='Adam%A0')

Allocate(Adam%SPB(nSym,nShell*(nShell+1)/2,2))

iE=0
Do iaSh=1,nShell
   Do ibSh=1,iaSh
      iShp = iaSh*(iaSh-1)/2 + ibSh

      If (iShp_rs(iShp)<=0) Cycle

      If (nnBstRSh(Jsym,iShp_rs(iShp),1)<=0) Cycle

      Do iSymb=1,nSym
         iSyma=MulD2h(iSymb,Jsym)
         If (iSyma<iSymb) Cycle

         iS = iE + 1

         n1 = nBasSh(iSyma,iaSh)
         n2 = nBasSh(iSymb,ibSh)

         iE = iE + n1*JNUM*n2

         Adam%SPB(iSyma,iShp_rs(iShp),1)%A3(1:n1,1:JNUM,1:n2) => Adam%A0(iS:iE)
         Adam%SPB(iSyma,iShp_rs(iShp),1)%A21(1:n1*JNUM,1:n2) => Adam%A0(iS:iE)
         Adam%SPB(iSyma,iShp_rs(iShp),1)%A12(1:n1,1:JNUM*n2) => Adam%A0(iS:iE)

         If (iaSh==ibSh) Cycle

         iS = iE + 1

         n1 = nBasSh(iSyma,ibSh)
         n2 = nBasSh(iSymb,iaSh)

         iE = iE + n1*JNUM*n2

         Adam%SPB(iSyma,iShp_rs(iShp),2)%A3(1:n1,1:JNUM,1:n2) => Adam%A0(iS:iE)
         Adam%SPB(iSyma,iShp_rs(iShp),2)%A21(1:n1*JNUM,1:n2) => Adam%A0(iS:iE)
         Adam%SPB(iSyma,iShp_rs(iShp),2)%A12(1:n1,1:JNUM*n2) => Adam%A0(iS:iE)

      End Do

   End Do
End Do

End Subroutine Allocate_L_Full


Subroutine deallocate_L_Full(Adam)
Implicit None
Type (L_Full_Type):: Adam

Integer iaSh, ibSh, iShp, iSyma

Do iaSh=1,Adam%nShell
   Do ibSh=1,iaSh
      iShp = iaSh*(iaSh-1)/2 + ibSh

      Do iSyma=1, Adam%nSym

         Adam%SPB(iSyma,iShp,1)%A3 => Null()
         Adam%SPB(iSyma,iShp,1)%A21=> Null()
         Adam%SPB(iSyma,iShp,1)%A12=> Null()
         Adam%SPB(iSyma,iShp,2)%A3 => Null()
         Adam%SPB(iSyma,iShp,2)%A21=> Null()
         Adam%SPB(iSyma,iShp,2)%A12=> Null()

      End Do

   End Do
End Do

deallocate(Adam%SPB)
call mma_deallocate(Adam%A0)
Adam%iCase=0
Adam%nSym=0
Adam%iSym=0
Adam%nShell=0

End Subroutine deallocate_L_Full


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     !
!                      L a b - T Y P E   S E C T I O N                !
!                                                                     !
!                      Lab storaged                                   !
!                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine Allocate_Lab(Lab,JNUM,nBasSh,nBas,nShell,nSym,nDen,Memory)
Implicit None

Type (Lab_Type), Target:: Lab
Integer JNUM, nShell, nSym, nDen
Integer nBasSh(nSym,nShell), nBas(nSym)
Integer, Optional :: Memory

Integer iSym, iDen, Lab_Memory
Integer iE, iS, iSh

Lab_Memory=0
Do iSym = 1, nSym
   Lab_Memory=Max(nBas(iSym),Lab_Memory)
End Do
Lab_Memory = Lab_Memory * JNUM * nDen

If (Present(Memory)) Then
   Memory=Lab_Memory
   Return
End If

Lab%nSym=nSym
Lab%nDen=nDen
Lab%nShell=nShell
Call mma_allocate(Lab%A0,Lab_Memory,Label='Lab%A0')
Call mma_allocate(Lab%Keep,nShell,nDen,Label='Lab%Keep')
Allocate(Lab%SB(nShell,nSym,nDen))

Do iSym = 1, nSym
   iE = 0
   Do iDen = 1, nDen
      Do iSh = 1, nShell

         iS = iE + 1
         iE = iE + nBasSh(iSym,iSh) * JNUM

         Lab%SB(iSh,iSym,iDen)%A(1:nBasSh(iSym,iSh)*JNUM) => Lab%A0(iS:iE)

      End Do
   End Do
End Do

End Subroutine Allocate_Lab

Subroutine Deallocate_Lab(Lab)
Implicit None
Type (Lab_Type) Lab

Integer iSym, iDen, iSh

Do iSym = 1, Lab%nSym
   Do iDen = 1, Lab%nDen
      Do iSh = 1, Lab%nShell

         Lab%SB(iSh,iSym,iDen)%A=>Null()

      End Do
   End Do
End Do

Lab%nSym=0
Lab%nDen=0
Lab%nShell=0
Call mma_deallocate(Lab%A0)
Call mma_deallocate(Lab%Keep)
Deallocate(Lab%SB)

End Subroutine Deallocate_Lab

End Module Data_Structures
