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
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************
!#define _DEBUGPRINT_

Subroutine sg_d2mat(SGS,CIS,EXS,Psi,nCSFs,PsiSym,D2MAT,nD2MAT)

use Index_functions, only: iTri
use sguga, only: CIStruct, EXStruct, SGStruct
use definitions, only: iwp, wp

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp), intent(in) :: PsiSym, nCSFs, nD2MAT
real(kind=wp), intent(in) :: Psi(nCSFs)
real(kind=wp), intent(out) :: D2MAT(nD2MAT)

!If (lRas==0) Then
   Call Mk_D2MAT_CASSCF(SGS,CIS,EXS,D2MAT,nD2MAT,Psi,nCSFs,PsiSym)
!Else
!   Call Mk_D2MAT_RASSCF(SGS,CIS,EXS,D2MAT,nD2MAT,Psi,nCSFs,PsiSym)
!End If

Contains

Subroutine Mk_D2MAT_CASSCF(SGS,CIS,EXS,D2MAT,nD2MAT,Psi,nCSFs,PsiSym)

use stdalloc, only: mma_allocate, mma_deallocate
use sguga, only: CIStruct, EXStruct, SGStruct
use Constants, only: Zero, One
use Definitions, only: iwp, wp

Implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp), intent(in) :: PsiSym, nCSFs, nD2MAT
real(kind=wp), intent(in) :: Psi(nCSFs)
real(kind=wp), intent(inout) :: D2MAT(nD2MAT)

real(kind=wp), Allocatable, Target :: Eij_Psi_X(:), Elk_Psi_X(:)
real(kind=wp), Pointer :: Eij_Psi(:)=>Null(), Elk_Psi(:)=>Null()
real(kind=wp), parameter :: CPQ=One
integer(kind=iwp) :: iOrb, jOrb, kOrb, lOrb

integer(kind=iwp) :: ijOrb, klOrb, klijOrb
integer(kind=iwp) :: ikOrb, ljOrb, ikljOrb, MaxDim, mCSFs
integer(kind=iwp) :: iSym, jSym, kSym, lSym, ijSym, klSym,lOrb_Min
real(kind=wp) :: D_ij=Zero, P_klij=Zero

MaxDim=MaxVal(CIS%nCSF(:))
Call mma_allocate(Eij_Psi_X,MaxDim,Label='Eij_Psi_X')
Call mma_allocate(Elk_Psi_X,MaxDim,Label='Elk_Psi_X')

D2MAT(:)=Zero

!
!   Compute P_lk,ij = <Psi|e_lkij|Psi> + <Psi|e_lkji|Psi>
!
!   Where e_lkij = E_lk E_ij - delta_ki E_lj
!
!   Hence,
!
!   P_lk,ij = sum_m <Psi|E_lk|m>(<m|Eij|Psi>+<m|E_ji|Psi>) - delta_ki D_lj -delta_kj D_li
!
!   In the symmetrization set we use the fact that for none zero elements we have that
!   <m|E_ij|m'>=<m'|E_ji|m>
!   This is trivial in the case of no symmetry. However, in the case of symmetry we have to be careful.
!
!   Let say that |m> belongs to the set of CSFs that are included in |Psi>, while |m'> doesn't. Then
!   there is no symmetrization to be achieved. We need to keep track on this in Mk_Eij_Psi.
!
!
! The size of the outer loop is nOrb*(nOrb+1)/2, i.e. nOrb=10 => 55 tasks to parallelize over.

Do iOrb =1, SGS%nLev
   iSym=SGS%ISM(iOrb)
Do jOrb =1, iOrb
   ijOrb=iTri(iOrb,jOrb)
   jSym=SGS%ISM(jOrb)
   ijSym=iEOR(iSym-1,jSym-1)+1
   mCSFs = CIS%nCSF(iEOR(PsiSym-1,ijSym-1)+1)

   Eij_Psi(1:mCSFs)=>Eij_Psi_X
   Elk_Psi(1:mCSFs)=>Elk_Psi_X

!  Operate with E_ij on |Psi> and produce E_ij|Psi>
!
!  Compute Dij to be distributed below
!
   Eij_Psi(:)=Zero
   Call SG_Epq_Psi(SGS,CIS,EXS,iOrb,jOrb,CPQ,PsiSym,Psi,Eij_Psi)
   If (iOrb==jOrb) Then
      D_ij=Dot_Product(Psi,Eij_Psi)
   Else
      If (iSym==jSym) D_ij=Dot_Product(Psi,Eij_Psi)
      Call SG_Epq_Psi(SGS,CIS,EXS,jOrb,iOrb,CPQ,PsiSym,Psi,Eij_Psi)
   End If

   Do kOrb=1,iOrb
      kSym=SGS%ISM(kOrb)
      lOrb=kOrb
      lSym=SGS%ISM(lOrb)
      klSym=iEOR(kSym-1,lSym-1)+1
      If (ijSym/=klSym) Cycle

!     Add the -d_il E_jk term, here in the form of a -d_kl E_ij term of Piklj
      ikOrb=iTri(iOrb,kOrb)
      ljOrb=iTri(lOrb,jOrb)
      ikljOrb=iTri(ikOrb,ljOrb)
      D2MAT(ikljOrb)=D2MAT(ikljOrb) - D_ij
   End Do

   Do kOrb=iOrb,SGS%nLev
      kSym=SGS%ISM(kOrb)
      lOrb_Min=1
      If (kOrb==iOrb) lOrb_Min=jOrb
   Do lOrb=lOrb_Min,kOrb
      klOrb=iTri(kOrb,lOrb)
      lSym=SGS%ISM(lOrb)
      klSym=iEOR(kSym-1,lSym-1)+1
      If (ijSym/=klSym) Cycle

!     Operate with E_lk on |Psi> and produce E_lk|Psi>

      Elk_Psi(1:mCSFs)=Zero
     ! <Psi|E_kl is computed as E_lk|Psi>
      Call SG_Epq_Psi(SGS,CIS,EXS,lOrb,kOrb,CPQ,PsiSym,Psi,Elk_Psi)
!     cap with <Psi|E_kl on E_ij|Psi>, contribution to Pklij
      P_klij = Dot_Product(Elk_Psi,Eij_Psi)

      klijOrb=iTri(klOrb,ijOrb)
      D2MAT(klijOrb)=D2MAT(klijOrb) + P_klij

   End Do
   End Do
   Eij_Psi=>Null()
   Elk_Psi=>Null()

End Do
End Do

Call mma_deallocate(Elk_Psi_X)
Call mma_deallocate(Eij_Psi_X)

End Subroutine Mk_D2MAT_CASSCF

End Subroutine sg_d2mat
