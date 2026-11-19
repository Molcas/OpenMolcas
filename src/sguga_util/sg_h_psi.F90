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
use Constants, only: Zero, One, Half
use Definitions, only: iwp, wp

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
!Write (u6,*) 'SGUGA driven algorithm for H|Psi>, DGEMV version'
!Write (u6,*) 'i>=j Symmetrize: E_ij + E_ji, k>=l Symmetrize E_kl + E_lj.'

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

!> H|Psi> in a restricted (RAS) CSF space without intermediate states outside
!> that space (J. Olsen, B. O. Roos, P. Jorgensen, H. J. Aa. Jensen,
!> J. Chem. Phys. 89, 2185 (1988), doi:10.1063/1.455063, Eqs. 9a, 28, 29; as
!> Mk_Sigma_Vector_loop_driven5 of ASTA).
!>
!> With the ordered orbital pairs P = (i,j) indexed (i-1)*n+j,
!>
!>   H = sum_ab h'_ab E_ab + sum_{P<=Q} f_PQ (Q|P) E_Q E_P,  f_PP = 1/2, else 1,
!>   h'_ab = h_ab - sum_{j<a} (aj|jb) - [a>b] (aa|ab) - 1/2 delta_ab (aa|aa).
!>
!> If the orbitals are ordered RAS1 < RAS2 < RAS3, an intermediate E_P|Psi>
!> that leaves the RAS space cannot be brought back by E_Q with Q >= P, so
!> E_P|Psi> only needs its components in the RAS space, which is what
!> SG_Epq_Psi produces. The pair order therefore runs over the graph levels
!> (RAS1, RAS2, RAS3, each over all symmetries), while TU and TUVX are indexed
!> by the active orbitals in symmetry-blocked order (L2Orb). For CAS both
!> orders coincide and the result equals that of sg_h_psi.
subroutine sg_h_psi_ras(SGS,CIS,EXS,Psi,nCSFs,PsiSym,Sigma,TUVX_Tri,nTUVX_Tri,TU_Tri,nTU_Tri)

use Index_functions, only: iTri
use symmetry_info, only: MUL
use sguga, only: SGStruct, CIStruct, EXStruct, sg_epq_psi, SG_Level_Map
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: iwp, wp

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp), intent(in) :: nCSFs, PsiSym, nTUVX_Tri, nTU_Tri
real(kind=wp), intent(in) :: Psi(nCSFs), TUVX_Tri(nTUVX_Tri), TU_Tri(nTU_Tri)
real(kind=wp), intent(out) :: Sigma(nCSFs)
integer(kind=iwp), parameter :: nBuff = 10
integer(kind=iwp) :: iBuff, iA, iB, iLev, iOrb, ij, ijSym, jLev, jOrb, kl, klSym, kLev, kOrb, lLev, lOrb, mCSFs, mLev, nLev, &
                     SigmaSym
real(kind=wp) :: Fact, OneInt, TUVX(nBuff)
integer(kind=iwp), allocatable :: L2Orb(:)
real(kind=wp), allocatable :: Ekl_Eij_Psi(:,:), Eij_Psi(:)

nLev = SGS%nLev

! level -> active orbital (symmetry-blocked order of TU/TUVX)
call mma_allocate(L2Orb,nLev,Label='L2Orb')
call SG_Level_Map(SGS,L2Orb)

call mma_allocate(Eij_Psi,maxval(CIS%nCSF(:)),Label='Eij_Psi')
call mma_allocate(Ekl_Eij_Psi,nCSFs,nBuff,Label='Ekl_Eij_Psi')

Sigma(:) = Zero

do iLev=1,nLev
  iOrb = L2Orb(iLev)
  do jLev=1,nLev
    jOrb = L2Orb(jLev)
    ijSym = MUL(SGS%ISM(iLev),SGS%ISM(jLev))
    SigmaSym = MUL(PsiSym,ijSym)
    mCSFs = CIS%nCSF(SigmaSym)
    if (mCSFs == 0) cycle
    ij = (iLev-1)*nLev+jLev

    ! E_ij|Psi>, restricted to the CSF space
    Eij_Psi(1:mCSFs) = Zero
    call SG_Epq_Psi(SGS,CIS,EXS,iLev,jLev,One,PsiSym,Psi,Eij_Psi)

    ! one-electron part, h'_ij E_ij
    if (ijSym == 1) then
      OneInt = TU_Tri(iTri(iOrb,jOrb))
      do mLev=1,iLev-1
        iA = L2Orb(mLev)
        OneInt = OneInt-TUVX_Tri(iTri(iTri(iOrb,iA),iTri(iA,jOrb)))
      end do
      if (iLev > jLev) then
        OneInt = OneInt-TUVX_Tri(iTri(iTri(iOrb,iOrb),iTri(iOrb,jOrb)))
      else if (iLev == jLev) then
        OneInt = OneInt-Half*TUVX_Tri(iTri(iTri(iOrb,iOrb),iTri(iOrb,iOrb)))
      end if
      Sigma(:) = Sigma(:)+OneInt*Eij_Psi(1:nCSFs)
    end if

    ! two-electron part, f (kl|ij) E_kl E_ij for kl >= ij
    iBuff = 0
    do kLev=iLev,nLev
      kOrb = L2Orb(kLev)
      do lLev=1,nLev
        kl = (kLev-1)*nLev+lLev
        if (kl < ij) cycle
        klSym = MUL(SGS%ISM(kLev),SGS%ISM(lLev))
        if (klSym /= ijSym) cycle
        lOrb = L2Orb(lLev)
        iBuff = iBuff+1
        Ekl_Eij_Psi(:,iBuff) = Zero
        call SG_Epq_Psi(SGS,CIS,EXS,kLev,lLev,One,SigmaSym,Eij_Psi,Ekl_Eij_Psi(:,iBuff))
        Fact = One
        if (kl == ij) Fact = Half
        iA = iTri(kOrb,lOrb)
        iB = iTri(iOrb,jOrb)
        TUVX(iBuff) = Fact*TUVX_Tri(iTri(iA,iB))
        if (iBuff == nBuff) then
          call DGEMV_('N',nCSFs,nBuff,One,Ekl_Eij_Psi,nCSFs,TUVX,1,One,Sigma,1)
          iBuff = 0
        end if
      end do
    end do
    if (iBuff /= 0) call DGEMV_('N',nCSFs,iBuff,One,Ekl_Eij_Psi,nCSFs,TUVX,1,One,Sigma,1)

  end do
end do

call mma_deallocate(Ekl_Eij_Psi)
call mma_deallocate(Eij_Psi)
call mma_deallocate(L2Orb)

end subroutine sg_h_psi_ras

!> sg_h_psi for the module state iState (SGS/CIS/EXS of module SGUGA); lets
!> procedures of SGUGA call sg_h_psi with intrinsic arguments only.
subroutine sg_h_psi_state(iState,Psi,nCSFs,PsiSym,Sigma,TUVX_Tri,nTUVX_Tri,TU_Tri,nTU_Tri)

use sguga, only: SGS, CIS, EXS
use Definitions, only: iwp, wp

implicit none
integer(kind=iwp), intent(in) :: iState, nCSFs, PsiSym, nTUVX_Tri, nTU_Tri
real(kind=wp), intent(in) :: Psi(nCSFs), TUVX_Tri(nTUVX_Tri), TU_Tri(nTU_Tri)
real(kind=wp), intent(out) :: Sigma(nCSFs)

call sg_h_psi(SGS(iState),CIS(iState),EXS(iState),Psi,nCSFs,PsiSym,Sigma,TUVX_Tri,nTUVX_Tri,TU_Tri,nTU_Tri)

end subroutine sg_h_psi_state

!> sg_h_psi_ras for the module state iState (SGS/CIS/EXS of module SGUGA),
!> with intrinsic arguments only (see sg_h_psi_state).
subroutine sg_h_psi_ras_state(iState,Psi,nCSFs,PsiSym,Sigma,TUVX_Tri,nTUVX_Tri,TU_Tri,nTU_Tri)

use sguga, only: SGS, CIS, EXS
use Definitions, only: iwp, wp

implicit none
integer(kind=iwp), intent(in) :: iState, nCSFs, PsiSym, nTUVX_Tri, nTU_Tri
real(kind=wp), intent(in) :: Psi(nCSFs), TUVX_Tri(nTUVX_Tri), TU_Tri(nTU_Tri)
real(kind=wp), intent(out) :: Sigma(nCSFs)

call sg_h_psi_ras(SGS(iState),CIS(iState),EXS(iState),Psi,nCSFs,PsiSym,Sigma,TUVX_Tri,nTUVX_Tri,TU_Tri,nTU_Tri)

end subroutine sg_h_psi_ras_state
