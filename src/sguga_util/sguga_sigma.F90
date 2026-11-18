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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!               2025,2026, Roland Lindh                                *
!***********************************************************************

!> SGUGA: the Hamiltonian acting on CI vectors of a split graph (sigma
!> vectors H|Psi> of CAS and RAS graphs), its exact and fast diagonals, the
!> full CI matrix, and the diagonal of a one-electron occupation operator.
module SGUGA_Sigma

implicit none
private

public :: sg_ci_full, sg_h_psi, sg_h_psi_ras, sg_hdiag_fast, sg_hdiag_ref, SG_Occ_Diag

contains

subroutine sg_h_psi(SGS,CIS,EXS,Psi,nCSFs,PsiSym,Sigma,TUVX_Tri,nTUVX_Tri,TU_Tri,nTU_Tri)
  ! H|Psi> in the CSF space of a split graph (CAS; intermediate states in the
  ! graph), with S_ij = E_ij + E_ji (i > j), S_ii = E_ii:
  !   H|Psi> = sum_{i>=j} h'_ij S_ij|Psi> + 1/2 sum_{i>=j,k>=l} (kl|ij) S_kl S_ij|Psi>,
  !   h'_ij = h_ij - 1/2 sum_k (ik|kj).
  ! The two-electron part is resolved per symmetry block of the orbital
  ! pairs: D(:,p) = S_p|Psi>, W = D G with G(q,p) = 1/2 (q|p), and
  ! H|Psi> += sum_q S_q W(:,q): 2 n_pair applications of S per block (and one
  ! DGEMM) instead of n_pair**2. A block whose D does not fit in memory is
  ! done directly (S_kl S_ij|Psi> for every pair of pairs).

  use Index_functions, only: iTri
  use symmetry_info, only: MUL
  use sguga, only: SGStruct, CIStruct, EXStruct, sg_epq_psi
  use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
  use Constants, only: Zero, One, Half
  use Definitions, only: iwp, wp

  implicit none
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(EXStruct), intent(inout) :: EXS
  integer(kind=iwp), intent(in) :: nCSFs, PsiSym, nTUVX_Tri, nTU_Tri
  real(kind=wp), intent(in) :: Psi(nCSFs), TUVX_Tri(nTUVX_Tri), TU_Tri(nTU_Tri)
  real(kind=wp), intent(out) :: Sigma(nCSFs)
  integer(kind=iwp) :: nOrb, iOrb, jOrb, kOrb, ijSym, SigmaSym, mCSFs, nP, p, q, q0, nB, MaxMem
  integer(kind=iwp), allocatable :: PI(:), PJ(:)
  real(kind=wp) :: TwoInt
  real(kind=wp), allocatable :: D(:,:), W(:,:), G(:,:), X(:)

  Sigma(:) = Zero
  nOrb = SGS%nLev
  call mma_allocate(PI,max(1,nOrb*(nOrb+1)/2),Label='PI')
  call mma_allocate(PJ,max(1,nOrb*(nOrb+1)/2),Label='PJ')
  do ijSym=1,SGS%nSym
    ! the orbital pairs i >= j of symmetry ijSym
    nP = 0
    do iOrb=1,nOrb
      do jOrb=1,iOrb
        if (MUL(SGS%ISM(iOrb),SGS%ISM(jOrb)) /= ijSym) cycle
        nP = nP+1
        PI(nP) = iOrb
        PJ(nP) = jOrb
      end do
    end do
    SigmaSym = MUL(PsiSym,ijSym)
    mCSFs = CIS%nCSF(SigmaSym)
    if ((nP == 0) .or. (mCSFs == 0)) cycle
    call mma_maxDBLE(MaxMem)
    if (real(mCSFs,kind=wp)*real(nP+1,kind=wp)+real(nP,kind=wp)**2 > 0.8_wp*real(MaxMem,kind=wp)) then
      call Block_Direct()
      cycle
    end if
    ! D(:,p) = S_p|Psi>, and the one-electron part (totally symmetric pairs)
    call mma_allocate(D,mCSFs,nP,Label='D_HPsi')
    D(:,:) = Zero
    do p=1,nP
      call SG_Epq_Psi(SGS,CIS,EXS,PI(p),PJ(p),One,PsiSym,Psi,D(:,p))
      if (PI(p) /= PJ(p)) call SG_Epq_Psi(SGS,CIS,EXS,PJ(p),PI(p),One,PsiSym,Psi,D(:,p))
      if (ijSym == 1) then
        TwoInt = Zero
        do kOrb=1,nOrb
          TwoInt = TwoInt+TUVX_Tri(iTri(iTri(PI(p),kOrb),iTri(kOrb,PJ(p))))
        end do
        Sigma(:) = Sigma(:)+(TU_Tri(iTri(PI(p),PJ(p)))-Half*TwoInt)*D(:,p)
      end if
    end do
    ! G(q,p) = 1/2 (q|p)
    call mma_allocate(G,nP,nP,Label='G_HPsi')
    do p=1,nP
      do q=1,nP
        G(q,p) = Half*TUVX_Tri(iTri(iTri(PI(q),PJ(q)),iTri(PI(p),PJ(p))))
      end do
    end do
    ! W = D G in batches of columns, then H|Psi> += S_q W(:,q)
    call mma_maxDBLE(MaxMem)
    nB = int(min(real(nP,kind=wp),max(One,0.8_wp*real(MaxMem,kind=wp)/real(mCSFs,kind=wp))),kind=iwp)
    call mma_allocate(W,mCSFs,nB,Label='W_HPsi')
    do q0=1,nP,nB
      nB = min(nB,nP-q0+1)
      call DGEMM_('N','T',mCSFs,nB,nP,One,D,mCSFs,G(q0,1),nP,Zero,W,mCSFs)
      do q=q0,q0+nB-1
        call SG_Epq_Psi(SGS,CIS,EXS,PI(q),PJ(q),One,SigmaSym,W(:,q-q0+1),Sigma)
        if (PI(q) /= PJ(q)) call SG_Epq_Psi(SGS,CIS,EXS,PJ(q),PI(q),One,SigmaSym,W(:,q-q0+1),Sigma)
      end do
    end do
    call mma_deallocate(W)
    call mma_deallocate(G)
    call mma_deallocate(D)
  end do
  call mma_deallocate(PJ)
  call mma_deallocate(PI)

contains

  ! the block ijSym directly: S_kl S_ij|Psi> for every pair of pairs
  subroutine Block_Direct()
    integer(kind=iwp) :: r
    real(kind=wp), allocatable :: Y(:)
    call mma_allocate(X,mCSFs,Label='X_HPsi')
    call mma_allocate(Y,nCSFs,Label='Y_HPsi')
    do p=1,nP
      X(:) = Zero
      call SG_Epq_Psi(SGS,CIS,EXS,PI(p),PJ(p),One,PsiSym,Psi,X)
      if (PI(p) /= PJ(p)) call SG_Epq_Psi(SGS,CIS,EXS,PJ(p),PI(p),One,PsiSym,Psi,X)
      if (ijSym == 1) then
        TwoInt = Zero
        do kOrb=1,nOrb
          TwoInt = TwoInt+TUVX_Tri(iTri(iTri(PI(p),kOrb),iTri(kOrb,PJ(p))))
        end do
        Sigma(:) = Sigma(:)+(TU_Tri(iTri(PI(p),PJ(p)))-Half*TwoInt)*X(:)
      end if
      do r=1,nP
        Y(:) = Zero
        call SG_Epq_Psi(SGS,CIS,EXS,PI(r),PJ(r),One,SigmaSym,X,Y)
        if (PI(r) /= PJ(r)) call SG_Epq_Psi(SGS,CIS,EXS,PJ(r),PI(r),One,SigmaSym,X,Y)
        Sigma(:) = Sigma(:)+Half*TUVX_Tri(iTri(iTri(PI(r),PJ(r)),iTri(PI(p),PJ(p))))*Y(:)
      end do
    end do
    call mma_deallocate(Y)
    call mma_deallocate(X)
  end subroutine Block_Direct

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
  ! The same H|Psi> as sg_h_psi_ras_direct (operator order of Olsen et al.:
  ! f (kl|ij) E_kl E_ij|Psi> for kl >= ij, k >= i, so that the intermediate
  ! states E_ij|Psi> need no CSFs outside the space), resolved per symmetry
  ! block of the ordered level pairs: D(:,p) = E_ij|Psi> (restricted to the
  ! space, as there), W = D G^T with G(kl,ij) = f (kl|ij) for the allowed
  ! pairs and zero otherwise, H|Psi> += E_kl W(:,kl): 2 n_pair applications
  ! of E per block instead of n_pair**2/2. Falls back to
  ! sg_h_psi_ras_direct if the vectors of a block do not fit in memory.
  use Index_functions, only: iTri
  use symmetry_info, only: MUL
  use sguga, only: SGStruct, CIStruct, EXStruct, sg_epq_psi, SG_Level_Map
  use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
  use Constants, only: Zero, One, Half
  use Definitions, only: iwp, wp
  implicit none
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(EXStruct), intent(inout) :: EXS
  integer(kind=iwp), intent(in) :: nCSFs, PsiSym, nTUVX_Tri, nTU_Tri
  real(kind=wp), intent(in) :: Psi(nCSFs), TUVX_Tri(nTUVX_Tri), TU_Tri(nTU_Tri)
  real(kind=wp), intent(out) :: Sigma(nCSFs)
  integer(kind=iwp) :: iA, iLev, iOrb, ijSym, jLev, jOrb, kLev, lLev, mCSFs, mLev, nLev, SigmaSym, nP, ip, iq, q0, nB, &
                       MaxMem, MaxNeed
  real(kind=wp) :: Fact, OneInt
  integer(kind=iwp), allocatable :: L2Orb(:), PI(:), PJ(:)
  real(kind=wp), allocatable :: D(:,:), W(:,:), G(:,:)
  logical(kind=iwp), allocatable :: Used(:)

  nLev = SGS%nLev
  ! memory: the vectors of the largest block
  MaxNeed = 0
  do ijSym=1,SGS%nSym
    nP = 0
    do iLev=1,nLev
      do jLev=1,nLev
        if (MUL(SGS%ISM(iLev),SGS%ISM(jLev)) == ijSym) nP = nP+1
      end do
    end do
    MaxNeed = max(MaxNeed,(CIS%nCSF(MUL(PsiSym,ijSym))+1)*nP+nP**2)
  end do
  call mma_maxDBLE(MaxMem)
  if (real(MaxNeed,kind=wp) > 0.8_wp*real(MaxMem,kind=wp)) then
    call sg_h_psi_ras_direct(SGS,CIS,EXS,Psi,nCSFs,PsiSym,Sigma,TUVX_Tri,nTUVX_Tri,TU_Tri,nTU_Tri)
    return
  end if

  ! level -> active orbital (symmetry-blocked order of TU/TUVX)
  call mma_allocate(L2Orb,nLev,Label='L2Orb')
  call SG_Level_Map(SGS,L2Orb)
  call mma_allocate(PI,nLev**2,Label='PI')
  call mma_allocate(PJ,nLev**2,Label='PJ')
  Sigma(:) = Zero
  do ijSym=1,SGS%nSym
    ! the ordered level pairs (i,j) of symmetry ijSym, in the order of ij = (i-1)*nLev+j
    nP = 0
    do iLev=1,nLev
      do jLev=1,nLev
        if (MUL(SGS%ISM(iLev),SGS%ISM(jLev)) /= ijSym) cycle
        nP = nP+1
        PI(nP) = iLev
        PJ(nP) = jLev
      end do
    end do
    SigmaSym = MUL(PsiSym,ijSym)
    mCSFs = CIS%nCSF(SigmaSym)
    if ((nP == 0) .or. (mCSFs == 0)) cycle
    ! D(:,p) = E_ij|Psi> (restricted to the CSF space), and the one-electron part h'_ij E_ij
    call mma_allocate(D,mCSFs,nP,Label='D_HPsi')
    D(:,:) = Zero
    do ip=1,nP
      iLev = PI(ip)
      jLev = PJ(ip)
      call SG_Epq_Psi(SGS,CIS,EXS,iLev,jLev,One,PsiSym,Psi,D(:,ip))
      if (ijSym == 1) then
        iOrb = L2Orb(iLev)
        jOrb = L2Orb(jLev)
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
        Sigma(:) = Sigma(:)+OneInt*D(:,ip)
      end if
    end do
    ! G(kl,ij) = f (kl|ij) for kl >= ij and k >= i (f = 1/2 for kl = ij)
    call mma_allocate(G,nP,nP,Label='G_HPsi')
    call mma_allocate(Used,nP,Label='Used')
    G(:,:) = Zero
    Used(:) = .false.
    do ip=1,nP
      iLev = PI(ip)
      jLev = PJ(ip)
      do iq=ip,nP
        kLev = PI(iq)
        lLev = PJ(iq)
        if (kLev < iLev) cycle
        Fact = One
        if (iq == ip) Fact = Half
        G(iq,ip) = Fact*TUVX_Tri(iTri(iTri(L2Orb(kLev),L2Orb(lLev)),iTri(L2Orb(iLev),L2Orb(jLev))))
        Used(iq) = .true.
      end do
    end do
    ! W = D G^T in batches of columns, then H|Psi> += E_kl W(:,kl)
    call mma_maxDBLE(MaxMem)
    nB = int(min(real(nP,kind=wp),max(One,0.8_wp*real(MaxMem,kind=wp)/real(mCSFs,kind=wp))),kind=iwp)
    call mma_allocate(W,mCSFs,nB,Label='W_HPsi')
    do q0=1,nP,nB
      nB = min(nB,nP-q0+1)
      call DGEMM_('N','T',mCSFs,nB,nP,One,D,mCSFs,G(q0,1),nP,Zero,W,mCSFs)
      do iq=q0,q0+nB-1
        if (Used(iq)) call SG_Epq_Psi(SGS,CIS,EXS,PI(iq),PJ(iq),One,SigmaSym,W(:,iq-q0+1),Sigma)
      end do
    end do
    call mma_deallocate(W)
    call mma_deallocate(Used)
    call mma_deallocate(G)
    call mma_deallocate(D)
  end do
  call mma_deallocate(PJ)
  call mma_deallocate(PI)
  call mma_deallocate(L2Orb)
end subroutine sg_h_psi_ras

subroutine sg_h_psi_ras_direct(SGS,CIS,EXS,Psi,nCSFs,PsiSym,Sigma,TUVX_Tri,nTUVX_Tri,TU_Tri,nTU_Tri)

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

end subroutine sg_h_psi_ras_direct

subroutine sg_hdiag_ref(SGS,CIS,EXS,nCSFs,CSFSym, &
                        TUVX_Tri,nTUVX_Tri, &
                        TU_Tri,nTU_Tri, &
                        Diag)

  !#define _TEST_
#ifdef _TEST_
    use sguga, only: sg_Epq_Psi, sg_EpqEqp_diag
    use Definitions, only: u6
#endif
    use sguga, only: SGStruct, CIStruct, EXStruct
    use stdalloc, only: mma_allocate, mma_deallocate
    use Constants, only: Zero, One
    use Definitions, only: wp, iwp

    implicit none

    type(SGStruct), intent(in)    :: SGS
    type(CIStruct), intent(in)    :: CIS
    type(EXStruct), intent(inout) :: EXS

    integer(kind=iwp), intent(in) :: nCSFs
    integer(kind=iwp), intent(in) :: CSFSym

    integer(kind=iwp), intent(in) :: nTUVX_Tri
    integer(kind=iwp), intent(in) :: nTU_Tri

    real(kind=wp), intent(in) :: TUVX_Tri(nTUVX_Tri)
    real(kind=wp), intent(in) :: TU_Tri(nTU_Tri)

    real(kind=wp), intent(out) :: Diag(nCSFs)

    real(kind=wp), allocatable :: Psi(:)
    real(kind=wp), allocatable :: Sigma(:)

    integer(kind=iwp) :: iCSF, mCSFs

#ifdef _TEST_
    logical(kind=iwp) :: ReUse_Save
    integer(kind=iwp) :: p, q
    real(kind=wp), allocatable :: Epq(:), Eqp(:)
#endif

  ! Note, the sigma vector might be of a different length than Psi. For safty measurers allocate the largest possible size.
    mCSFs=MaxVal(CIS%nCSF)
    call mma_allocate(Psi,nCSFs,Label='Psi')
    call mma_allocate(Sigma,mCSFs,Label='Sigma')

#ifdef _TEST_
  ! turn off reuse of SGTMP
  ReUse_Save=EXS%Reuse_SGTMP
  EXS%Reuse_SGTMP=.False.
  call mma_allocate(Epq,nCSFs,Label='Epq')
  call mma_allocate(Eqp,nCSFs,Label='Eqp')

  do p=1,SGS%nLev
  do q=1,SGS%nLev
  !Write (u6,*) 'p,q=',p,q

  do iCSF=1,nCSFs
    Psi(:) = Zero
    Psi(iCSF) = One
    Sigma(:) = Zero
    call SG_Epq_Psi(SGS,CIS,EXS, &
                    p,q, &
                    One, &
                    CSFSym, &
                    Psi, &
                    Sigma)
  ! Call RecPrt('Sigma',' ',Sigma,1,nCSFs)
    Epq(iCSF)=Dot_product(Sigma,Sigma)
  end do
  !Call RecPrt('EqpEpq(exact)',' ',Epq,1,nCSFs)

  call SG_EpqEqp_Diag(SGS,CIS,EXS, &
                      p,q, &
                      One, &
                      CSFSym, &
                      Psi, &
                      Sigma)
  !Call RecPrt('EqpEpq(Direct)',' ',Sigma,1,nCSFs)

  Do iCSF=1,nCSFs
  If (Abs(Epq(iCSF)-Sigma(iCSF))>1.0e-10_wp) Then
     Write (u6,*) 'p, q, iCSF=',p, q, iCSF
     Write (u6,*) 'MidLev=',SGS%MidLev
     Write (u6,*) 'Level=',SGS%Level(1:SGS%nLev)
     Write (u6,*) 'L2Act=',SGS%L2Act(1:SGS%nLev)
     Call RecPrt('EqpEpq(exact)',' ',Epq,1,nCSFs)
     Call RecPrt('EqpEpq(Direct)',' ',Sigma,1,nCSFs)
     Call Abend()
  End If
  End Do

  end do
  end do

  call mma_deallocate(Eqp)
  call mma_deallocate(Epq)

  EXS%Reuse_SGTMP=ReUse_Save
  Write (u6,*)
  Write (u6,*) 'Test passed!'
  Write (u6,*)
  !Call Abend()
#endif

    do iCSF=1,nCSFs

       Psi(:)   = Zero
       Sigma(:) = Zero

       Psi(iCSF) = One

       ! sg_h_psi needs intermediate states outside a RAS space; there the
       ! operator-ordered sg_h_psi_ras is used (both read TU/TUVX in the
       ! symmetry-blocked orbital order)
       if (SGS%nRsPrt == 1) then
         call sg_h_psi(SGS,CIS,EXS, &
                       Psi,nCSFs,CSFSym, &
                       Sigma, &
                       TUVX_Tri,nTUVX_Tri, &
                       TU_Tri,nTU_Tri)
       else
         call sg_h_psi_ras(SGS,CIS,EXS,Psi,nCSFs,CSFSym,Sigma,TUVX_Tri,nTUVX_Tri,TU_Tri,nTU_Tri)
       end if

       Diag(iCSF) = Sigma(iCSF)

    end do
  !  call RecPrt('Diag',' ',Diag,1,nCSFs)

    call mma_deallocate(Sigma)
    call mma_deallocate(Psi)

end subroutine sg_hdiag_ref

subroutine sg_hdiag_fast(SGS,CIS,EXS,nCSFs,CSFSym, &
                         TUVX_Tri,nTUVX_Tri, &
                         TU_Tri,nTU_Tri, &
                         Diag)

  use sguga, only: SGStruct, CIStruct, EXStruct, SG_Epq_Psi, SG_EpqEqp_Diag, SG_Level_Ints
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One, Half, Two
  use Definitions, only: wp, iwp
  use Index_Functions, only: iTri

  implicit none

  type(SGStruct), intent(in)    :: SGS
  type(CIStruct), intent(in)    :: CIS
  type(EXStruct), intent(inout) :: EXS

  integer(kind=iwp), intent(in) :: nCSFs
  integer(kind=iwp), intent(in) :: CSFSym

  integer(kind=iwp), intent(in) :: nTUVX_Tri
  integer(kind=iwp), intent(in) :: nTU_Tri

  real(kind=wp), intent(in)  :: TUVX_Tri(nTUVX_Tri)
  real(kind=wp), intent(in)  :: TU_Tri(nTU_Tri)

  real(kind=wp), intent(out) :: Diag(nCSFs)

  real(kind=wp), allocatable :: Psi(:)
  real(kind=wp), allocatable :: Epp(:)
  real(kind=wp), allocatable :: Eqq(:)
  real(kind=wp), allocatable :: EqqEpp(:)
  real(kind=wp), allocatable :: EpqEqp(:)
  ! the integrals indexed by graph level (TU/TUVX: symmetry-blocked orbitals)
  real(kind=wp), allocatable :: TU_L(:), TUVX_L(:)

  integer(kind=iwp) :: p,q
  integer(kind=iwp) :: ppTri
  integer(kind=iwp) :: ppqqTri
  integer(kind=iwp) :: pqqpTri

  real(kind=wp) :: OneInt
  real(kind=wp) :: Coul
  real(kind=wp) :: Exch

  Diag(:) = Zero

  call mma_allocate(TU_L,nTU_Tri,Label='TU_L')
  call mma_allocate(TUVX_L,nTUVX_Tri,Label='TUVX_L')
  call SG_Level_Ints(SGS,TU_Tri,nTU_Tri,TUVX_Tri,nTUVX_Tri,TU_L,TUVX_L)

  call mma_allocate(Psi,nCSFs,Label='Psi')
  call mma_allocate(Epp,nCSFs,Label='Epp')
  call mma_allocate(Eqq,nCSFs,Label='Eqq')
  call mma_allocate(EqqEpp,nCSFs,Label='EqqEpp')
  call mma_allocate(EpqEqp,nCSFs,Label='EpqEqp')

  !
  ! Trick:
  !
  ! Psi(:)=1 allows SG_Epq_Psi to return the diagonal
  ! expectation values for ALL CSFs simultaneously.
  !
  Psi(:) = One

do p=1,SGS%nLev

   Epp(:) = Zero

   call SG_Epq_Psi(SGS,CIS,EXS, &
                   p,p, &
                   One, &
                   CSFSym, &
                   Psi, &
                   Epp)

!
! One-electron contribution
!
   ppTri = iTri(p,p)

   OneInt = TU_L(ppTri)

   Diag(:) = Diag(:) + OneInt*Epp(:)

!
!-------------------------------------------------------
! PPPP contribution
!
! 1/2 (pp|pp) [Epp*Epp - Epp]
!-------------------------------------------------------
!
   ppqqTri = iTri(iTri(p,p),iTri(p,p))
   Coul    = Half*TUVX_L(ppqqTri)

   Diag(:)=Diag(:) + Coul*(Epp(:)**2 - Epp(:))

! PPQQ and PQQP contributions

   do q=p+1,SGS%nLev

      EqqEpp(:)=Zero
      call SG_Epq_Psi(SGS,CIS,EXS, &
                      q,q, &
                      One, &
                      CSFSym, &
                      Epp, &
                      EqqEpp)
!
!-------------------------------------------------------
! PPQQ contribution
!
! (pp|qq) Epp Eqq
!-------------------------------------------------------
!
      ppqqTri = iTri(iTri(p,p),iTri(q,q))

      Diag(:)=Diag(:) + TUVX_L(ppqqTri)*EqqEpp(:)

!
!-------------------------------------------------------
! PQQP contribution
!
! 1/2 (pq|qp)
! [ <EpqEqp> + <EqpEpq> - Epp - Eqq ]
!-------------------------------------------------------
!
      Eqq(:)=Zero
      call SG_Epq_Psi(SGS,CIS,EXS, &
                      q,q, &
                      One, &
                      CSFSym, &
                      Psi, &
                      Eqq)

      EpqEqp(:)=Zero
      call SG_EpqEqp_Diag(SGS,CIS,EXS, &
                          p,q, &
                          One, &
                          CSFSym, &
                          Psi, &
                          EpqEqp)

      ! EpqEqp = <E_qp E_pq> (p < q), the walk |E_pq e_I|**2, is exact also
      ! in a RAS graph; <E_pq E_qp> = <E_qp E_pq> + n_p - n_q (its walk
      ! |E_qp e_I|**2 would need intermediate states outside a RAS space),
      ! so <E_pq E_qp> + <E_qp E_pq> - n_p - n_q = 2 (<E_qp E_pq> - n_q).
      pqqpTri = iTri(iTri(q,p),iTri(q,p))
      Exch    = Half*TUVX_L(pqqpTri)

      Diag(:)=Diag(:) + Exch * Two * ( EpqEqp(:) - Eqq(:) )
   end do
end do

  call mma_deallocate(EpqEqp)
  call mma_deallocate(TUVX_L)
  call mma_deallocate(TU_L)
  call mma_deallocate(EqqEpp)
  call mma_deallocate(Eqq)
  call mma_deallocate(Epp)
  call mma_deallocate(Psi)

end subroutine sg_hdiag_fast

subroutine sg_ci_full(SGS,CIS,EXS,nCSFs,CSFSym, &
                      TUVX_Tri,nTUVX_Tri, &
                      TU_Tri,nTU_Tri, &
                      Full)

  use sguga, only: SGStruct, CIStruct, EXStruct
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One
  use Definitions, only: wp, iwp

  implicit none

  type(SGStruct), intent(in)    :: SGS
  type(CIStruct), intent(in)    :: CIS
  type(EXStruct), intent(inout) :: EXS

  integer(kind=iwp), intent(in) :: nCSFs
  integer(kind=iwp), intent(in) :: CSFSym

  integer(kind=iwp), intent(in) :: nTUVX_Tri
  integer(kind=iwp), intent(in) :: nTU_Tri

  real(kind=wp), intent(in) :: TUVX_Tri(nTUVX_Tri)
  real(kind=wp), intent(in) :: TU_Tri(nTU_Tri)

  real(kind=wp), intent(out) :: Full(nCSFs,nCSFs)

  real(kind=wp), allocatable :: Psi(:)
  real(kind=wp), allocatable :: Sigma(:)

  integer(kind=iwp) :: iCSF, jCSF

  call mma_allocate(Sigma,nCSFs,Label='Sigma')
  call mma_allocate(Psi,nCSFs,Label='Psi')

  do iCSF=1,nCSFs

     Psi(:)   = Zero
     Sigma(:) = Zero

     Psi(iCSF) = One

     ! Construct H|I> = Sum_J <J|H|I> |J>

     call sg_h_psi(SGS,CIS,EXS, &
                   Psi,nCSFs,CSFSym, &
                   Sigma, &
                   TUVX_Tri,nTUVX_Tri, &
                   TU_Tri,nTU_Tri)

     Do jCSF=1,nCSFs
        Full(jCSF,iCSF) = Sigma(jCSF)
     End Do

  end do
  call RecPrt('CI matrix in CSF basis',' ',Full,nCSFs,nCSFs)

  call mma_deallocate(Sigma)
  call mma_deallocate(Psi)

end subroutine sg_ci_full

!> DIAG(I) = sum_lev n_lev(I) ETA(lev) for the CSFs I of symmetry ISYCI of
!> state iState (n the occupation numbers of the levels); from H0DIAG and
!> DIELMV of CASPT2: the upper and lower half walks of each (midvertex,
!> upper symmetry) block separately.
subroutine SG_Occ_Diag(iState,ISYCI,ETA,DIAG,nDiag)

  use Symmetry_Info, only: Mul
  use sguga, only: CIS, SGS
  use CSF_Graph, only: nPack
  use Constants, only: Zero
  use Definitions, only: wp, iwp

  implicit none
  integer(kind=iwp), intent(in) :: iState, ISYCI, nDiag
  real(kind=wp), intent(in) :: ETA(*)
  real(kind=wp), intent(out) :: DIAG(nDiag)
  integer(kind=iwp) :: ICS, IEMU, ISYDWN, ISYUP, JCS, MV, NDWN, NUP

  DIAG(:) = Zero
  IEMU = 1
  do MV=1,CIS(iState)%nMidV
    do ISYUP=1,SGS(iState)%nSym
      NUP = CIS(iState)%NOW(1,ISYUP,MV)
      if (NUP == 0) cycle
      ISYDWN = Mul(ISYUP,ISYCI)
      NDWN = CIS(iState)%NOW(2,ISYDWN,MV)
      if (NDWN == 0) cycle
      ICS = 1+CIS(iState)%IOW(1,ISYUP,MV)
      JCS = 1+CIS(iState)%IOW(2,ISYDWN,MV)
      call Block_Diag(ICS,JCS,NUP,NDWN,DIAG(IEMU))
      IEMU = IEMU+NUP*NDWN
    end do
  end do

  contains

    subroutine Block_Diag(ICS,JCS,NUP,NDWN,EMU)
      integer(kind=iwp), intent(in) :: ICS, JCS, NUP, NDWN
      real(kind=wp), intent(inout) :: EMU(NUP,NDWN)
      integer(kind=iwp) :: I, IC, IC1, II, IOC, ISTEP, LEV, LV1, nIpWlk, nLev, MidLev
      real(kind=wp) :: rSum

      nIpWlk = CIS(iState)%nIpWlk
      nLev = SGS(iState)%nLev
      MidLev = SGS(iState)%MidLev
      ! the upper half walks
      do I=1,NUP
        II = ICS-1+NIPWLK*(I-1)
        rSum = Zero
        do LV1=MidLev+1,nLev,nPack
          II = II+1
          IC = CIS(iState)%ICase(II)
          do LEV=LV1,min(LV1+(nPack-1),nLev)
            IC1 = IC/4
            ISTEP = IC-4*IC1
            IOC = (ISTEP+1)/2
            rSum = rSum+real(IOC,kind=wp)*ETA(LEV)
            IC = IC1
          end do
        end do
        EMU(I,:) = EMU(I,:)+rSum
      end do
      ! the lower half walks
      do I=1,NDWN
        II = JCS-1+NIPWLK*(I-1)
        rSum = Zero
        do LV1=1,MidLev,nPack
          II = II+1
          IC = CIS(iState)%ICase(II)
          do LEV=LV1,min(LV1+(nPack-1),MidLev)
            IC1 = IC/4
            ISTEP = IC-4*IC1
            IOC = (ISTEP+1)/2
            rSum = rSum+real(IOC,kind=wp)*ETA(LEV)
            IC = IC1
          end do
        end do
        EMU(:,I) = EMU(:,I)+rSum
      end do
    end subroutine Block_Diag

end subroutine SG_Occ_Diag

end module SGUGA_Sigma
