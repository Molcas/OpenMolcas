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
use Definitions, only: iwp, wp, u6

Implicit none
type (SGStruct), intent(in)    :: SGS
type (CIStruct), intent(in)    :: CIS
type (EXStruct), intent(inout) :: EXS
integer(kind=iwp), intent(in) ::nCSFs, PsiSym, nTUVX_Tri, nTU_Tri
real(kind=wp), intent(in) :: Psi(nCSFs), TUVX_Tri(nTUVX_Tri), TU_Tri(nTU_Tri)
real(kind=wp), intent(out) :: Sigma(nCSFs)

real(kind=wp), Allocatable, Target :: Eij_Psi_X(:), Ekl_Eij_Psi(:,:)
real(kind=wp), Pointer :: Eij_Psi(:)=>Null()
integer(kind=iwp) :: iOrb, jOrb, kOrb, lOrb, nOrb
integer(kind=iwp) :: ijOrb, klOrb, klijOrb
integer(kind=iwp) :: ikOrb, ljOrb, ikljOrb, MaxDim, mCSFs
integer(kind=iwp) :: iSym, jSym, kSym, lSym, ijSym, klSym
real(kind=wp) :: OneInt, TwoInt
integer(kind=iwp) :: SigmaSym

integer(kind=iwp), Parameter:: nBuff=10
real(kind=wp) :: TUVX(nBUff)
integer(kind=iwp) :: iBuff=0
real(kind=wp), parameter :: Alpha=One, Beta=One
integer(kind=iwp), parameter :: incx=1, incy=1
real(kind=wp) :: CPQ

! SGUGA driven algorithm for H|Psi>, DGEMV version
! i>=j Symmetrize: E_ij + E_ji, k>=l Symmetrize E_kl + E_lj.
Write (u6,*) 'SGUGA driven algorithm for H|Psi>, DGEMV version'
Write (u6,*) 'i>=j Symmetrize: E_ij + E_ji, k>=l Symmetrize E_kl + E_lj.'

Sigma(:)=Zero

MaxDim = maxval(CIS%nCSF(:))
call mma_allocate(Eij_Psi_X,MaxDim,Label='Eij_Psi_X')
call mma_allocate(Ekl_Eij_Psi,nCSFs,nBuff,Label='Ekl_Eij_Psi')

OneInt=Zero  ! Cardholder variable
TwoInt=Zero  ! Cardholder variable

!   The E_kl E_ij part

CPQ=One
nOrb=SGS%nLev

do iOrb=1,nOrb
   iSym=SGS%ISM(iOrb)
  do jOrb=1,iOrb
   ijOrb=iTri(iOrb,jOrb)

   jSym=SGS%ISM(jOrb)
   ijSym=MUL(iSym,jSym)

   SigmaSym=MUL(PsiSym,ijSym)
   mCSFs = CIS%nCSF(SigmaSym)

   Eij_Psi(1:mCSFs)=>Eij_Psi_X

!  Operate with E_ij on |Psi> and produce E_ij|Psi>
   Eij_Psi(:)=Zero
    call SG_Epq_Psi(SGS,CIS,EXS,iOrb,jOrb,CPQ,PsiSym,Psi,Eij_Psi)
    if (iOrb /= jOrb) call SG_Epq_Psi(SGS,CIS,EXS,jOrb,iOrb,CPQ,PsiSym,Psi,Eij_Psi)

   if (ijSym==1) then

      TwoInt=Zero
      do kOrb=1,nOrb
         lOrb=kOrb
         ikOrb=iTri(iOrb,kOrb)
         ljOrb=iTri(lOrb,jOrb)
         ikljOrb=iTri(ikOrb,ljOrb)
         TwoInt = TwoInt +  TUVX_Tri(ikljOrb)
      end do
      OneInt=TU_Tri(ijOrb) - Half * TwoInt

!     Operate with E_ij on |Psi> and produce E_ij|Psi>
      !call DaXpY_(nCSFs,OneInt,Eij_Psi(:),1,Sigma(:),1)
      Sigma(:)=Sigma(:)+OneInt*Eij_Psi(:)
    end if

   iBuff=0
   do kOrb=1,nOrb
      kSym=SGS%ISM(kOrb)
   do lOrb=1,kOrb
      klOrb=iTri(kOrb,lOrb)

      lSym=SGS%ISM(lOrb)
      klSym=MUL(kSym,lSym)

      if (ijSym /= klSym) cycle
      iBuff=iBuff+1

!     Operate with E_kl on E_ij_|Psi> and produce E_kl_E_ij_|Psi>
      Ekl_Eij_Psi(1:nCSFs,iBuff)=Zero

      call SG_Epq_Psi(SGS,CIS,EXS,kOrb,lOrb,CPQ,SigmaSym,Eij_Psi,Ekl_Eij_Psi(1:nCSFs,iBuff))
      if (kOrb /= lOrb) Call SG_Epq_Psi(SGS,CIS,EXS,lOrb,kOrb,CPQ,SigmaSym,Eij_Psi,Ekl_Eij_Psi(1:nCSFs,iBuff))

      klijOrb=iTri(klOrb,ijOrb)
      TUVX(iBuff) = Half*TUVX_Tri(klijOrb)

      if (iBuff == nBuff) then
          Call DGEMV_('N',nCSFs,nBuff,Alpha,Ekl_Eij_Psi(1:nCSFs,1:nBuff),nCSFs, &
                      TUVX(1:nBuff),incx,Beta,Sigma,incy)
!         Sigma(:)=Sigma(:) + Mat_Mul(Ekl_Eij_Psi(1:nCSFs,1:nBuff),TUVX(1:nBuff))
!         Do i = 1, nBuff
!            Sigma(:)=Sigma(:) + Ekl_Eij_Psi(:,i)*TUVX(i)
!         End Do
          iBuff=0
!         Ekl_Eij_Psi(:,:)=Zero
      End If

      end do
    end do

   if (iBuff /= 0) Then
      Call DGEMV_('N',nCSFs,iBuff,Alpha,Ekl_Eij_Psi(1:nCSFs,1:iBuff),nCSFs,TUVX(1:iBuff),incx,Beta,Sigma,incy)
!     Sigma(:)=Sigma(:) + Mat_Mul(Ekl_Eij_Psi(1:nCSFs,1:nBuff),TUVX(1:nBuff))
!     do i = 1,iBuff
!        Sigma(:)=Sigma(:) + Ekl_Eij_Psi(:,i)*TUVX(i)
!     end do
      iBuff=0
   end if

    nullify(Eij_Psi)

  end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call mma_deallocate(Ekl_Eij_Psi)
call mma_deallocate(Eij_Psi_X)

end subroutine sg_h_psi
