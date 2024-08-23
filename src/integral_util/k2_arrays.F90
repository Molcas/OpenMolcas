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

Module k2_arrays

!
!     DeDe is an array with desymmetrized 1-particle densities used
!     in the integral direct construction of the Fock matrix. In this
!     the array contain the subblocks with their associatated base
!     pointer.
!     ipDeDe: the set of densities to be used sorted according to the
!      integral shell pairs they contribute. Pointers to the individual
!      matrices are stored in ipOffD.
!     ipD00: some of the pointers in ipOffD point to an "empty" slot
!      and this pointer is to this part of DeDe.
!     ipDijS: points to an auxiliary pice of memory which is used in
!      case that used a subset of the elements of a matrix is used. In
!      this case picky_ will extract those elements and put them into
!      this part of DeDe on the fly.
!
Integer, Allocatable :: ipOffD(:,:)
Real*8, Allocatable:: FT(:), DeDe(:)
Integer  ipDeDe, ipD00, ipDijS
Integer  nDeDe, nDeDe_DFT, MaxDe, MxDij, MxFT, nFT
Logical  :: DoGrad_=.False., DoHess_=.False.
Real*8, Target, Allocatable:: Fq(:), Dq(:)
Real*8, Pointer:: pFq(:)=>Null(), pDq(:)=>Null()
Real*8, Allocatable:: Aux(:)
Integer, Allocatable:: iSOSym(:,:)
Logical :: XMem=.False.
Real*8, Allocatable, Target:: Sew_Scr(:)

Type BraKet_Type
  Real*8, Pointer:: Zeta(:)
  Real*8, Pointer:: ZInv(:)
  Real*8, Pointer:: KappaAB(:)
  Real*8, Pointer:: P(:,:)
  Real*8, Pointer:: xA(:)
  Real*8, Pointer:: xB(:)
  Real*8, Pointer:: Eta(:)
  Real*8, Pointer:: EInv(:)
  Real*8, Pointer:: KappaCD(:)
  Real*8, Pointer:: Q(:,:)
  Real*8, Pointer:: xG(:)
  Real*8, Pointer:: xD(:)
  Real*8, Pointer:: xPre(:)
  Integer, Pointer:: IndZet(:)
  Integer, Pointer:: IndEta(:)
End Type BraKet_Type

Type(BraKet_Type) BraKet
Real*8, Allocatable, Target:: BraKet_Base_R(:)
Integer, Allocatable, Target:: BraKet_Base_I(:)







Contains

Subroutine Create_BraKet_Base(nZeta)
use stdalloc, only: mma_allocate
Implicit None
Integer nZeta

Integer Mem


Mem =nZeta*16
If (DoHess_) Mem = Mem + nZeta**2
Call mma_allocate(BraKet_Base_R,Mem,Label='Base_R')
Mem =(nZeta+1)*2
Call mma_allocate(BraKet_Base_I,Mem ,Label='Base_I')

End Subroutine Create_BraKet_Base

Subroutine Destroy_BraKet_Base()
use stdalloc, only: mma_deallocate
Implicit None

If (Allocated(BraKet_Base_R)) Call mma_deallocate(BraKet_Base_R)
If (Allocated(BraKet_Base_I)) Call mma_deallocate(BraKet_Base_I)

End Subroutine Destroy_BraKet_Base

Subroutine Create_BraKet(nZeta,nEta)
Implicit None
Integer nZeta, nEta

Integer iS, iE

If (.Not.Allocated(BraKet_base_R) .or. .Not.Allocated(BraKet_base_I)) Then
   Write(6,*) 'Braket_Base not allocated!'
   Call Abend()
End If

If (nZeta*nEta==0) Return

iE=0

If (nZeta/=0) Then
iS=iE+1
iE=iE+nZeta
Braket%Zeta(1:nZeta) => BraKet_Base_R(iS:iE)
iS=iE+1
iE=iE+nZeta
Braket%ZInv(1:nZeta) => BraKet_Base_R(iS:iE)
iS=iE+1
iE=iE+nZeta
Braket%KappaAB(1:nZeta) => BraKet_Base_R(iS:iE)
iS=iE+1
iE=iE+3*nZeta
Braket%P(1:nZeta,1:3) => BraKet_Base_R(iS:iE)
iS=iE+1
iE=iE+nZeta
Braket%xA(1:nZeta) => BraKet_Base_R(iS:iE)
iS=iE+1
iE=iE+nZeta
Braket%xB(1:nZeta) => BraKet_Base_R(iS:iE)
End if

If (nEta/=0) Then
iS=iE+1
iE=iE+nEta
Braket%Eta(1:nEta) => BraKet_Base_R(iS:iE)
iS=iE+1
iE=iE+nEta
Braket%EInv(1:nEta) => BraKet_Base_R(iS:iE)
iS=iE+1
iE=iE+nEta
Braket%KappaCD(1:nEta) => BraKet_Base_R(iS:iE)
iS=iE+1
iE=iE+3*nEta
Braket%Q(1:nEta,1:3) => BraKet_Base_R(iS:iE)
iS=iE+1
iE=iE+nEta
Braket%xG(1:nEta) => BraKet_Base_R(iS:iE)
iS=iE+1
iE=iE+nEta
Braket%xD(1:nEta) => BraKet_Base_R(iS:iE)
End If

If (nZeta*nEta/=0 .and. DoHess_) Then
iS=iE+1
iE=iE+nZeta*nEta
Braket%xPre(1:nEta) => BraKet_Base_R(iS:iE)
End If

iE=0

If (nZeta/=0) Then
iS=iE+1
iE=iE+nZeta+1
Braket%IndZet(1:nZeta+1) => BraKet_Base_I(iS:iE)
End If

If (nEta/=0) Then
iS=iE+1
iE=iE+nEta+1
Braket%IndEta(1:nEta+1) => BraKet_Base_I(iS:iE)
End If

End Subroutine Create_BraKet

Subroutine Destroy_BraKet()
Implicit None

Nullify(Braket%Zeta,Braket%ZInv,Braket%KappaAB,Braket%P,Braket%xA,Braket%xB,  &
         Braket%Eta, Braket%EInv,Braket%KappaCD,Braket%Q,Braket%xG,Braket%xD, &
         Braket%IndZet,Braket%IndEta,BraKet%xPre)
End Subroutine Destroy_BraKet

End Module k2_arrays
