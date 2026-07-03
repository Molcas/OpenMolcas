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
! Copyright (C) 2025, Roland Lindh                                     *
!***********************************************************************

Subroutine sg_h_psi(SGS,CIS,EXS,Psi,nCSFs,PsiSym,Sigma,TUVX_Tri,nTUVX_Tri,TU_Tri,nTU_Tri)

use Index_functions, only: iTri
use symmetry_info, only: MUL
use sguga, only: SGStruct, CIStruct, EXStruct, sg_epq_psi
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: iwp, wp

Implicit none
type (SGStruct), intent(in)    :: SGS
type (CIStruct), intent(in)    :: CIS
type (EXStruct), intent(inout) :: EXS
integer(kind=iwp), intent(in) ::nCSFs, PsiSym, nTUVX_Tri, nTU_Tri
real(kind=wp), intent(in) :: Psi(nCSFs)
real(kind=wp), intent(out) :: Sigma(nCSFs)
real(kind=wp), intent(in) :: TUVX_Tri(nTUVX_Tri), TU_Tri(nTU_Tri)

real(kind=wp), Allocatable, Target :: Eij_Psi_X(:), Ekl_Eij_Psi(:,:)
real(kind=wp), Pointer :: Eij_Psi(:)=>Null()
integer(kind=iwp) :: iOrb, jOrb, kOrb, lOrb, nOrb
integer(kind=iwp) :: ijOrb, klOrb, klijOrb
integer(kind=iwp) :: ikOrb, ljOrb, ikljOrb, MaxDim, mCSFs
integer(kind=iwp) :: iSym, jSym, kSym, lSym, ijSym, klSym
real(kind=wp) :: OneInt, TwoInt, Fact
integer(kind=iwp) :: lOrb_Max
integer(kind=iwp) :: SigmaSym

integer(kind=iwp), Parameter:: nBuff=10
real(kind=wp) :: TUVX(nBUff)
integer(kind=iwp) :: iBuff=0, i
real(kind=wp), parameter :: Alpha=One, Beta=One
integer(kind=iwp), parameter :: incx=1, incy=1
real(kind=wp) :: CPQ

! SGUGA driven algorithm for H|Psi>, DGEMV version
! i>=j Symmetrize: E_ij + E_ji, k>=l Symmetrize E_kl + E_lj, ij>=kl"

Sigma(:)=Zero

MaxDim=MaxVal(CIS%nCSF(:))
Call mma_allocate(Eij_Psi_X,MaxDim,Label='Eij_Psi_X')
Call mma_allocate(Ekl_Eij_Psi,nCSFs,nBuff,Label='Ekl_Eij_Psi')

OneInt=Zero  ! Cardholder variable
TwoInt=Zero  ! Cardholder variable

!   The E_kl E_ij part

CPQ=One
nOrb=SGS%nLev

Do iOrb =1, nOrb
   iSym=SGS%ISM(iOrb)
Do jOrb =1, iOrb
   ijOrb=iTri(iOrb,jOrb)

   jSym=SGS%ISM(jOrb)
   ijSym=MUL(iSym,jSym)
   SigmaSym=MUL(PsiSym,ijSym)
   mCSFs = CIS%nCSF(SigmaSym)

   Eij_Psi(1:mCSFs)=>Eij_Psi_X

!  Operate with E_ij on |Psi> and produce E_ij|Psi>
   Eij_Psi(:)=Zero
   Call SG_Epq_Psi(SGS,CIS,EXS,iOrb,jOrb,CPQ,PsiSym,Psi,Eij_Psi)
   If (iOrb/=jOrb) Then
      Call SG_Epq_Psi(SGS,CIS,EXS,jOrb,iOrb,CPQ,PsiSym,Psi,Eij_Psi)
   End If

   If (ijSym==0) Then
!     ijSym==0, do the E_ij part of the sigma vector
      TwoInt=Zero
      Do kOrb=1, nOrb
         lOrb=kOrb
         ikOrb=iTri(iOrb,kOrb)
         ljOrb=iTri(lOrb,jOrb)
         ikljOrb=iTri(ikOrb,ljOrb)
         TwoInt = TwoInt +  TUVX_Tri(ikljOrb)
      End Do
      OneInt=TU_Tri(ijOrb) - Half * TwoInt
!     Call DaXpY_(nCSFs,OneInt,Eij_Psi(:),1,Sigma(:),1)
      Sigma(:)=Sigma(:)+OneInt*Eij_Psi(:)
   End If

   iBuff=0
   Do kOrb=1, iOrb
      kSym=SGS%ISM(kOrb)
      lOrb_Max=kOrb
      If (kOrb==iOrb) lOrb_Max=jOrb
   Do lOrb=1,lOrb_Max
      klOrb=iTri(kOrb,lOrb)

      lSym=SGS%ISM(lOrb)
      klSym=MUL(kSym,lSym)

      If (ijSym/=klSym) Cycle

      Fact=One
      If (klOrb/=ijOrb) Fact = Two

!     Operate with E_kl on E_ij_|Psi> and produce E_kl_E_ij_|Psi>
      iBuff=iBuff+1
      Ekl_Eij_Psi(1:nCSFs,iBuff)=Zero

      Call SG_Epq_Psi(SGS,CIS,EXS,kOrb,lOrb,CPQ,SigmaSym,Eij_Psi,Ekl_Eij_Psi)

      If (kOrb/=lOrb) Then
         Call SG_Epq_Psi(SGS,CIS,EXS,lOrb,kOrb,CPQ,SigmaSym,Eij_Psi,Ekl_Eij_Psi)
      End If

      klijOrb=iTri(klOrb,ijOrb)
      TUVX(iBuff)= Fact * Half*TUVX_Tri(klijOrb)

      If (iBuff==nBuff) Then
!         Call DGEMV_('N',nCSFs,nBuff,Alpha,Ekl_Eij_Psi(1:nCSFs,1:nBuff),nCSFs, &
!                     TUVX(1:nBuff),incx,Beta,Sigma,incy)
!         Sigma(:)=Sigma(:) + Mat_Mul(Ekl_Eij_Psi(1:nCSFs,1:nBuff),TUVX(1:nBuff))
          Do i = 1, nBuff
             Sigma(:)=Sigma(:) + Ekl_Eij_Psi(:,i)*TUVX(i)
          End Do
          iBuff=0
          Ekl_Eij_Psi(:,:)=Zero
      End If

   End Do
   End Do

   If (iBuff/=0) Then
!     Call DGEMV_('N',nCSFs,iBuff,Alpha,Ekl_Eij_Psi(1:nCSFs,1:iBuff),nCSFs, &
!                 TUVX(1:iBuff),incx,Beta,Sigma,incy)
!     Sigma(:)=Sigma(:) + Mat_Mul(Ekl_Eij_Psi(1:nCSFs,1:nBuff),TUVX(1:nBuff))
      Do i = 1, iBuff
         Sigma(:)=Sigma(:) + Ekl_Eij_Psi(:,i)*TUVX(i)
      End Do
      iBuff=0
   End If

   Eij_Psi=>Null()

End Do
End Do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Call mma_deallocate(Ekl_Eij_Psi)
Call mma_deallocate(Eij_Psi_X)

End Subroutine sg_h_psi
