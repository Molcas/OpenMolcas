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
Public:: DSBA_Type, Allocate_DSBA, Deallocate_DSBA, Map_to_DSBA
Public:: SBA_Type, Allocate_SBA, Deallocate_SBA, Map_to_SBA
Public:: twxy_Type, Allocate_twxy, Deallocate_twxy, Map_to_twxy
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

Type V2
  Real*8, Pointer:: A(:,:)=>Null()
End Type V2

Type DSBA_Type
  Integer:: iCase=0
  Integer:: nSym=0
  Real*8, Allocatable :: A0(:)
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


Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     !
!                D S B A - T Y P E   S E C T I O N                    !
!                                                                     !
!                Diagonal Symmetry Blocked Arrays                     !
!                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Subroutine Allocate_DSBA(Adam,n,m,nSym,Case)
  Implicit None
  Type (DSBA_Type),Target, Intent(Out) :: Adam
  Integer, Intent(In) :: nSym
  Integer, Intent(In) :: n(nSym), m(nSym)
  Character(LEN=3), Intent(In), Optional :: Case

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
  Call mma_allocate(Adam%A0,MemTot,Label='%A0')

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
  Call mma_deallocate(Adam%A0)
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

  Subroutine Allocate_SBA(Adam,n,m,NUMV,iSym,nSym,iCase)
  Implicit None
  Type (SBA_Type),Target:: Adam
  Integer NUMV
  Integer iSym
  Integer nSym
  Integer n(nSym), m(nSym)
  Integer iCase
  Integer iE, iS, iSyma, iSymb, MemTot, n2Dim, n3Dim

  Integer i, j, MulD2h
  MulD2h(i,j) = iEOR(i-1,j-1) + 1

  Adam%iSym=iSym
  Adam%nSym=nSym
  Adam%iCase=iCase

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

End Module Data_Structures
