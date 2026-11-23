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

use Definitions, only: wp, iwp

implicit none
private

! work arrays of the routines below, kept between the calls (no new memory
! pages in every call for large CI spaces); freed by SG_Sigma_Free
real(kind=wp), allocatable, target :: Wrk_S1(:), Wrk_S2(:)

public :: SG_Check_HPsi_Op, SG_Sigma_Free, SG_Work_Need, sg_ci_full, sg_h_psi, sg_h_psi_op, sg_h_psi_ras, sg_h_psi_sq, &
          sg_hdiag_fast, sg_hdiag_ref, SG_Occ_Diag

contains

!> sg_h_op_core for a Hermitian operator of symmetry iOp given by
!> h_ij = TU_Tri(iTri(i,j)) and (ij|kl) = TUVX_Tri(iTri(iTri(i,j),iTri(k,l)))
!> over the active orbitals in symmetry-blocked order, as for sg_h_psi and
!> sg_h_psi_ras (mapped to the levels here; zero unless of symmetry iOp).
subroutine sg_h_psi_op(SGS,CIS,EXS,nPsi,Psi,PsiSym,nSig,Sigma,iOp,TUVX_Tri,nTUVX_Tri,TU_Tri,nTU_Tri)

  use Index_functions, only: iTri, nTri_Elem
  use sguga, only: SGStruct, CIStruct, EXStruct, SG_Level_Ints
  use stdalloc, only: mma_allocate, mma_deallocate
  use Definitions, only: iwp, wp

  implicit none
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(EXStruct), intent(inout) :: EXS
  integer(kind=iwp), intent(in) :: nPsi, PsiSym, nSig, iOp, nTUVX_Tri, nTU_Tri
  real(kind=wp), intent(in) :: Psi(nPsi), TUVX_Tri(nTUVX_Tri), TU_Tri(nTU_Tri)
  real(kind=wp), intent(out) :: Sigma(nSig)
  integer(kind=iwp) :: nL, n2, i, j, k, l
  real(kind=wp), allocatable :: TU_L(:), TUVX_L(:), H1(:,:), G2(:)

  nL = SGS%nLev
  n2 = nTri_Elem(nL*nL)
  call mma_allocate(TU_L,nTU_Tri,Label='HOp_TU_L')
  call mma_allocate(TUVX_L,nTUVX_Tri,Label='HOp_TUVX_L')
  call SG_Level_Ints(SGS,TU_Tri,nTU_Tri,TUVX_Tri,nTUVX_Tri,TU_L,TUVX_L)
  call mma_allocate(H1,nL,nL,Label='HOp_H1')
  call mma_allocate(G2,n2,Label='HOp_G2')
  do j=1,nL
    do i=1,nL
      H1(i,j) = TU_L(iTri(i,j))
      do l=1,nL
        do k=1,nL
          G2(iTri((j-1)*nL+i,(l-1)*nL+k)) = TUVX_L(iTri(iTri(i,j),iTri(k,l)))
        end do
      end do
    end do
  end do
  call mma_deallocate(TUVX_L)
  call mma_deallocate(TU_L)
  call sg_h_op_core(SGS,CIS,EXS,nPsi,Psi,PsiSym,nSig,Sigma,iOp,nL,H1,n2,G2)
  call mma_deallocate(G2)
  call mma_deallocate(H1)

end subroutine sg_h_psi_op

!> sg_h_op_core for a general (not necessarily Hermitian) operator of
!> symmetry iOp: H1(i,j) = h_ij (square) and G2(iTri(ij,kl)) = g(ij,kl)
!> with the square pair index ij = (j-1)*nA+i (g symmetric only under
!> (ij) <-> (kl)), over the nA active orbitals in symmetry-blocked order
!> (mapped to the levels here). For the time-dependent response of MCLR.
subroutine sg_h_psi_sq(SGS,CIS,EXS,nPsi,Psi,PsiSym,nSig,Sigma,iOp,nA,H1,n2,G2)

  use Index_functions, only: iTri
  use sguga, only: SGStruct, CIStruct, EXStruct, SG_Level_Map
  use stdalloc, only: mma_allocate, mma_deallocate
  use Definitions, only: iwp, wp, u6

  implicit none
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(EXStruct), intent(inout) :: EXS
  integer(kind=iwp), intent(in) :: nPsi, PsiSym, nSig, iOp, nA, n2
  real(kind=wp), intent(in) :: Psi(nPsi), H1(nA,nA), G2(n2)
  real(kind=wp), intent(out) :: Sigma(nSig)
  integer(kind=iwp) :: nL, i, j, k, l
  integer(kind=iwp), allocatable :: L2O(:)
  real(kind=wp), allocatable :: H1_L(:,:), G2_L(:)

  nL = SGS%nLev
  if ((nA /= nL) .or. (n2 /= nL*nL*(nL*nL+1)/2)) then
    write(u6,*) 'sg_h_psi_sq: wrong dimensions',nA,nL,n2
    call Abend()
  end if
  call mma_allocate(L2O,nL,Label='HSq_L2O')
  call SG_Level_Map(SGS,L2O)
  call mma_allocate(H1_L,nL,nL,Label='HSq_H1')
  call mma_allocate(G2_L,n2,Label='HSq_G2')
  do j=1,nL
    do i=1,nL
      H1_L(i,j) = H1(L2O(i),L2O(j))
      do l=1,nL
        do k=1,nL
          G2_L(iTri((j-1)*nL+i,(l-1)*nL+k)) = G2(iTri((L2O(j)-1)*nL+L2O(i),(L2O(l)-1)*nL+L2O(k)))
        end do
      end do
    end do
  end do
  call mma_deallocate(L2O)
  call sg_h_op_core(SGS,CIS,EXS,nPsi,Psi,PsiSym,nSig,Sigma,iOp,nL,H1_L,n2,G2_L)
  call mma_deallocate(G2_L)
  call mma_deallocate(H1_L)

end subroutine sg_h_psi_sq

!> Sigma = H Psi for a spin-free operator of symmetry iOp (Psi of symmetry
!> PsiSym, Sigma of symmetry Mul(PsiSym,iOp)), in a CAS or RAS split graph:
!>   H = sum_ij h_ij E_ij + 1/2 sum_ijkl g(ij,kl) (E_ij E_kl - delta_jk E_il),
!> with H1(i,j) = h_ij and G2(iTri(ij,kl)) = g(ij,kl) (square pair index
!> ij = (j-1)*nL+i), all in level order; h need not be symmetric, g only
!> under (ij) <-> (kl). The products E_q E_p (p applied first) are summed
!> over the CSFs of the graph; that is exact if p = (k,l) with k <= l or
!> q = (i,j) with i >= j (always in CAS). Otherwise E_q E_p = E_p E_q +
!> [E_q,E_p], which is exact then, and the commutator joins the one-body
!> part. Per symmetry block of the first pairs: D(:,p) = E_p Psi, W = D G,
!> Sigma += sum_q E_q W(:,q), in chunks of pairs if D and W do not fit.
subroutine sg_h_op_core(SGS,CIS,EXS,nPsi,Psi,PsiSym,nSig,Sigma,iOp,nL,H1,n2,G2)

  use Index_functions, only: iTri
  use symmetry_info, only: MUL
  use sguga, only: SGStruct, CIStruct, EXStruct, sg_epq_psi, SG_New_Tag, SG_No_Reuse
  use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
  use Constants, only: Zero, One, Half
  use Definitions, only: iwp, wp, u6

  implicit none
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(EXStruct), intent(inout) :: EXS
  integer(kind=iwp), intent(in) :: nPsi, PsiSym, nSig, iOp, nL, n2
  real(kind=wp), intent(in) :: Psi(nPsi), H1(nL,nL), G2(n2)
  real(kind=wp), intent(out) :: Sigma(nSig)
  integer(kind=iwp) :: SigSym, iT, iS, MidSym, nK, nP, nQ, nC, MaxMem, iP0, iQ0, cP, cQ, ip, iq, i, j, k, l, tPsi
  integer(kind=iwp), allocatable :: PI(:), PJ(:), QI(:), QJ(:)
  real(kind=wp) :: V
  real(kind=wp), allocatable :: F(:,:), D(:,:), W(:,:), G(:,:)

  SigSym = MUL(PsiSym,iOp)
  if ((nPsi < CIS%nCSF(PsiSym)) .or. (nSig < CIS%nCSF(SigSym))) then
    write(u6,*) 'sg_h_op_core: vectors shorter than the CSF spaces',nPsi,CIS%nCSF(PsiSym),nSig,CIS%nCSF(SigSym)
    call Abend()
  end if
  Sigma(:) = Zero
  if ((CIS%nCSF(PsiSym) == 0) .or. (CIS%nCSF(SigSym) == 0)) return
  tPsi = SG_New_Tag()

  ! the one-body part: h - 1/2 sum_j (ij|jl), plus the commutators of the
  ! products taken in the other order
  call mma_allocate(F,nL,nL,Label='HOp_F')
  F(:,:) = Zero
  do l=1,nL
    do i=1,nL
      if (MUL(SGS%ISm(i),SGS%ISm(l)) /= iOp) cycle
      F(i,l) = H1(i,l)
      do j=1,nL
        F(i,l) = F(i,l)-Half*Gv(i,j,j,l)
      end do
    end do
  end do
  do l=1,nL
    do k=1,nL
      do j=1,nL
        do i=1,nL
          if (Exact(i,j,k,l)) cycle
          ! (only the elements of the operator symmetry are read: the others
          ! may be undefined)
          if (MUL(MUL(SGS%ISm(i),SGS%ISm(j)),MUL(SGS%ISm(k),SGS%ISm(l))) /= iOp) cycle
          V = Half*Gv(i,j,k,l)
          if (V == Zero) cycle
          ! [E_ij,E_kl] = delta_jk E_il - delta_il E_kj
          if (j == k) F(i,l) = F(i,l)+V
          if (i == l) F(k,j) = F(k,j)-V
        end do
      end do
    end do
  end do
  do l=1,nL
    do i=1,nL
      if (MUL(SGS%ISm(i),SGS%ISm(l)) /= iOp) cycle
      if (F(i,l) == Zero) cycle
      call SG_Epq_Psi(SGS,CIS,EXS,i,l,F(i,l),PsiSym,Psi,Sigma,CI_Tag=tPsi)
    end do
  end do
  call mma_deallocate(F)

  ! the two-body part
  call mma_allocate(PI,nL*nL,Label='HOp_PI')
  call mma_allocate(PJ,nL*nL,Label='HOp_PJ')
  call mma_allocate(QI,nL*nL,Label='HOp_QI')
  call mma_allocate(QJ,nL*nL,Label='HOp_QJ')
  do iT=1,SGS%nSym
    ! first pairs p = (k,l) of symmetry iT, second pairs q = (i,j) of iS
    iS = MUL(iT,iOp)
    MidSym = MUL(PsiSym,iT)
    nK = CIS%nCSF(MidSym)
    if (nK == 0) cycle
    nP = 0
    nQ = 0
    do j=1,nL
      do i=1,nL
        if (MUL(SGS%ISm(i),SGS%ISm(j)) == iT) then
          nP = nP+1
          PI(nP) = i
          PJ(nP) = j
        end if
        if (MUL(SGS%ISm(i),SGS%ISm(j)) == iS) then
          nQ = nQ+1
          QI(nQ) = i
          QJ(nQ) = j
        end if
      end do
    end do
    if ((nP == 0) .or. (nQ == 0)) cycle
    call mma_maxDBLE(MaxMem)
    nC = int(0.4_wp*real(MaxMem,kind=wp)/real(2*nK+nL*nL,kind=wp),kind=iwp)
    nC = max(1,min(nC,max(nP,nQ)))
    call mma_allocate(D,nK,min(nC,nP),Label='HOp_D')
    call mma_allocate(W,nK,min(nC,nQ),Label='HOp_W')
    call mma_allocate(G,min(nC,nP),min(nC,nQ),Label='HOp_G')
    do iQ0=0,nQ-1,nC
      cQ = min(nC,nQ-iQ0)
      W(:,1:cQ) = Zero
      do iP0=0,nP-1,nC
        cP = min(nC,nP-iP0)
        do ip=1,cP
          D(:,ip) = Zero
          call SG_Epq_Psi(SGS,CIS,EXS,PI(iP0+ip),PJ(iP0+ip),One,PsiSym,Psi,D(:,ip),CI_Tag=tPsi)
        end do
        ! G(p,q): the coefficient of E_q E_p
        do iq=1,cQ
          i = QI(iQ0+iq)
          j = QJ(iQ0+iq)
          do ip=1,cP
            k = PI(iP0+ip)
            l = PJ(iP0+ip)
            G(ip,iq) = Zero
            if (Exact(i,j,k,l)) G(ip,iq) = G(ip,iq)+Half*Gv(i,j,k,l)
            if (.not. Exact(k,l,i,j)) G(ip,iq) = G(ip,iq)+Half*Gv(k,l,i,j)
          end do
        end do
        call DGEMM_('N','N',nK,cQ,cP,One,D,nK,G,size(G,1),One,W,nK)
      end do
      do iq=1,cQ
        call SG_Epq_Psi(SGS,CIS,EXS,QI(iQ0+iq),QJ(iQ0+iq),One,MidSym,W(:,iq),Sigma,CI_Tag=SG_No_Reuse)
      end do
    end do
    call mma_deallocate(G)
    call mma_deallocate(W)
    call mma_deallocate(D)
  end do
  call mma_deallocate(QJ)
  call mma_deallocate(QI)
  call mma_deallocate(PJ)
  call mma_deallocate(PI)

contains

  ! g(ij,kl), square pair indices
  pure function Gv(i,j,k,l)
    integer(kind=iwp), intent(in) :: i, j, k, l
    real(kind=wp) :: Gv
    Gv = G2(iTri((j-1)*nL+i,(l-1)*nL+k))
  end function Gv

  ! E_ij E_kl (E_kl first) summed over the graph is exact
  pure function Exact(i,j,k,l)
    integer(kind=iwp), intent(in) :: i, j, k, l
    logical(kind=iwp) :: Exact
    Exact = (k <= l) .or. (i >= j)
  end function Exact

end subroutine sg_h_op_core

!> sg_h_psi_op for arbitrary Hermitian integrals of every operator
!> symmetry iOp and every symmetry of Psi: Bra.(H Ket) against
!> sum_il h_il Rho1(i,l) + 1/2 sum (ij|kl) Rho2(ij,kl) with the transition
!> densities of sg_trans_two_pdm (checked against determinants), for two
!> arbitrary vectors; and for iOp = 1 element by element against sg_h_psi
!> (CAS) or sg_h_psi_ras (RAS).
subroutine SG_Check_HPsi_Op(SGS,CIS,EXS,nFail)

  use Index_functions, only: iTri, nTri_Elem
  use symmetry_info, only: MUL
  use sguga, only: SGStruct, CIStruct, EXStruct, SG_Level_Map
  use SGUGA_Densities, only: sg_trans_two_pdm
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One, Half
  use Definitions, only: iwp, wp, u6

  implicit none
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(EXStruct), intent(inout) :: EXS
  integer(kind=iwp), intent(out) :: nFail
  integer(kind=iwp), parameter :: MxCSF = 20000
  real(kind=wp), parameter :: Tol = 1.0e-10_wp
  integer(kind=iwp) :: nL, nL2, nTU, nTUVX, n2, iOp, PsiSym, SigSym, nPsi, nSig, i, j, k, l, ij, kl, nTest, ii, jj
  real(kind=wp) :: E1, E2, dMaxE, dMaxS, x
  integer(kind=iwp), allocatable :: L2O(:)
  real(kind=wp), allocatable :: TU(:), TUVX(:), TU_O(:), TUVX_O(:), Psi(:), Sig(:), Ref(:), Bra(:), Rho1(:,:), Rho2(:), &
                                H1L(:,:), H1O(:,:), G2L(:), G2O(:)
  real(kind=wp) :: dMaxQ

  nFail = 0
  nL = SGS%nLev
  if ((nL == 0) .or. (maxval(CIS%nCSF(1:SGS%nSym)) > MxCSF)) then
    write(u6,'(A)') ' sg_h_psi_op check skipped'
    return
  end if
  nL2 = nL*nL
  n2 = nTri_Elem(nL2)
  nTU = nTri_Elem(nL)
  nTUVX = nTri_Elem(nTU)
  call mma_allocate(TU,nTU,Label='CO_TU')
  call mma_allocate(TUVX,nTUVX,Label='CO_TUVX')
  call mma_allocate(TU_O,nTU,Label='CO_TU_O')
  call mma_allocate(TUVX_O,nTUVX,Label='CO_TUVX_O')
  call mma_allocate(L2O,nL,Label='CO_L2O')
  call SG_Level_Map(SGS,L2O)
  call mma_allocate(H1L,nL,nL,Label='CO_H1L')
  call mma_allocate(H1O,nL,nL,Label='CO_H1O')
  call mma_allocate(G2L,n2,Label='CO_G2L')
  call mma_allocate(G2O,n2,Label='CO_G2O')
  dMaxQ = Zero
  call mma_allocate(Rho1,nL,nL,Label='CO_Rho1')
  call mma_allocate(Rho2,n2,Label='CO_Rho2')
  nTest = 0
  dMaxE = Zero
  dMaxS = Zero
  do iOp=1,SGS%nSym
    ! arbitrary Hermitian integrals of symmetry iOp
    TU(:) = Zero
    TUVX(:) = Zero
    do i=1,nL
      do j=1,i
        if (MUL(SGS%ISm(i),SGS%ISm(j)) == iOp) TU(iTri(i,j)) = sin(0.37_wp*real(i+3*j,kind=wp)+real(iOp,kind=wp))
        do k=1,nL
          do l=1,k
            if (MUL(MUL(SGS%ISm(i),SGS%ISm(j)),MUL(SGS%ISm(k),SGS%ISm(l))) /= iOp) cycle
            x = real(iTri(iTri(i,j),iTri(k,l)),kind=wp)
            TUVX(iTri(iTri(i,j),iTri(k,l))) = Half*cos(0.113_wp*x+0.7_wp*real(iOp,kind=wp))
          end do
        end do
      end do
    end do
    ! (TU, TUVX in level order; TU_O, TUVX_O the same in orbital order)
    do i=1,nL
      do j=1,i
        TU_O(iTri(L2O(i),L2O(j))) = TU(iTri(i,j))
        do k=1,nL
          do l=1,k
            TUVX_O(iTri(iTri(L2O(i),L2O(j)),iTri(L2O(k),L2O(l)))) = TUVX(iTri(iTri(i,j),iTri(k,l)))
          end do
        end do
      end do
    end do
    ! arbitrary non-Hermitian integrals of symmetry iOp (sg_h_psi_sq): h
    ! general, g(ij,kl) symmetric only under (ij) <-> (kl); level order
    ! (H1L, G2L) and orbital order (H1O, G2O)
    H1L(:,:) = Zero
    G2L(:) = Zero
    do j=1,nL
      do i=1,nL
        ij = (j-1)*nL+i
        if (MUL(SGS%ISm(i),SGS%ISm(j)) == iOp) H1L(i,j) = sin(0.29_wp*real(i,kind=wp)+0.71_wp*real(j,kind=wp)+real(iOp,kind=wp))
        H1O(L2O(i),L2O(j)) = H1L(i,j)
        do l=1,nL
          do k=1,nL
            kl = (l-1)*nL+k
            if (kl > ij) cycle
            if (MUL(MUL(SGS%ISm(i),SGS%ISm(j)),MUL(SGS%ISm(k),SGS%ISm(l))) == iOp) &
              G2L(iTri(ij,kl)) = Half*cos(0.0731_wp*real(iTri(ij,kl),kind=wp)-0.3_wp*real(iOp,kind=wp))
            G2O(iTri((L2O(j)-1)*nL+L2O(i),(L2O(l)-1)*nL+L2O(k))) = G2L(iTri(ij,kl))
          end do
        end do
      end do
    end do
    do PsiSym=1,SGS%nSym
      SigSym = MUL(PsiSym,iOp)
      nPsi = CIS%nCSF(PsiSym)
      nSig = CIS%nCSF(SigSym)
      if ((nPsi == 0) .or. (nSig == 0)) cycle
      call mma_allocate(Psi,nPsi,Label='CO_Psi')
      call mma_allocate(Sig,nSig,Label='CO_Sig')
      call mma_allocate(Bra,nSig,Label='CO_Bra')
      do i=1,nPsi
        Psi(i) = cos(0.917_wp*real(i,kind=wp)+0.21_wp*real(PsiSym,kind=wp))
      end do
      do i=1,nSig
        Bra(i) = sin(0.613_wp*real(i,kind=wp)-0.33_wp*real(SigSym,kind=wp))
      end do
      call sg_h_psi_op(SGS,CIS,EXS,nPsi,Psi,PsiSym,nSig,Sig,iOp,TUVX_O,nTUVX,TU_O,nTU)
      E1 = dot_product(Bra,Sig)
      call sg_trans_two_pdm(SGS,CIS,EXS,nSig,Bra,SigSym,nPsi,Psi,PsiSym,nL,Rho1,n2,Rho2)
      E2 = Zero
      do l=1,nL
        do k=1,nL
          kl = (l-1)*nL+k
          E2 = E2+TU(iTri(k,l))*Rho1(k,l)
          do j=1,nL
            do i=1,nL
              ij = (j-1)*nL+i
              E2 = E2+Half*TUVX(iTri(iTri(i,j),iTri(k,l)))*Rho2(iTri(ij,kl))
            end do
          end do
        end do
      end do
      dMaxE = max(dMaxE,abs(E1-E2))
      if (abs(E1-E2) > Tol*max(One,abs(E2))) then
        nFail = nFail+1
        write(u6,'(A,2I3,2ES22.14)') ' sg_h_psi_op mismatch: iOp, PsiSym, Bra.H Ket, from densities =',iOp,PsiSym,E1,E2
      end if
      ! the general operator
      call sg_h_psi_sq(SGS,CIS,EXS,nPsi,Psi,PsiSym,nSig,Sig,iOp,nL,H1O,n2,G2O)
      E1 = dot_product(Bra,Sig)
      E2 = Zero
      do l=1,nL
        do k=1,nL
          kl = (l-1)*nL+k
          E2 = E2+H1L(k,l)*Rho1(k,l)
          do j=1,nL
            do i=1,nL
              ij = (j-1)*nL+i
              E2 = E2+Half*G2L(iTri(ij,kl))*Rho2(iTri(ij,kl))
            end do
          end do
        end do
      end do
      dMaxQ = max(dMaxQ,abs(E1-E2))
      if (abs(E1-E2) > Tol*max(One,abs(E2))) then
        nFail = nFail+1
        write(u6,'(A,2I3,2ES22.14)') ' sg_h_psi_sq mismatch: iOp, PsiSym, Bra.H Ket, from densities =',iOp,PsiSym,E1,E2
      end if
      ! the Hermitian operator through sg_h_psi_sq = sg_h_psi_op
      do jj=1,nL
        do ii=1,nL
          H1O(L2O(ii),L2O(jj)) = TU(iTri(ii,jj))
        end do
      end do
      call Sym_G2(G2O)
      call mma_allocate(Ref,nSig,Label='CO_Ref')
      call sg_h_psi_op(SGS,CIS,EXS,nPsi,Psi,PsiSym,nSig,Ref,iOp,TUVX_O,nTUVX,TU_O,nTU)
      call sg_h_psi_sq(SGS,CIS,EXS,nPsi,Psi,PsiSym,nSig,Sig,iOp,nL,H1O,n2,G2O)
      dMaxQ = max(dMaxQ,maxval(abs(Ref-Sig)))
      if (maxval(abs(Ref-Sig)) > Tol) then
        nFail = nFail+1
        write(u6,'(A,2I3,ES12.4)') ' sg_h_psi_sq vs sg_h_psi_op (Hermitian): iOp, PsiSym, max. difference',iOp,PsiSym, &
          maxval(abs(Ref-Sig))
      end if
      call mma_deallocate(Ref)
      ! restore the non-Hermitian operator in orbital order
      do jj=1,nL
        do ii=1,nL
          H1O(L2O(ii),L2O(jj)) = H1L(ii,jj)
        end do
      end do
      call Map_G2(G2L,G2O)
      if (iOp == 1) then
        call mma_allocate(Ref,nSig,Label='CO_Ref')
        if (SGS%nRsPrt > 1) then
          call sg_h_psi_ras(SGS,CIS,EXS,Psi,nPsi,PsiSym,Ref,TUVX_O,nTUVX,TU_O,nTU)
        else
          call sg_h_psi(SGS,CIS,EXS,Psi,nPsi,PsiSym,Ref,TUVX_O,nTUVX,TU_O,nTU)
        end if
        dMaxS = max(dMaxS,maxval(abs(Ref-Sig)))
        if (maxval(abs(Ref-Sig)) > Tol) then
          nFail = nFail+1
          write(u6,'(A,I3,ES12.4)') ' sg_h_psi_op vs sg_h_psi(_ras): PsiSym, max. difference',PsiSym,maxval(abs(Ref-Sig))
        end if
        call mma_deallocate(Ref)
      end if
      nTest = nTest+1
      call mma_deallocate(Bra)
      call mma_deallocate(Sig)
      call mma_deallocate(Psi)
    end do
  end do
  write(u6,'(/,A,I0,A)') ' sg_h_psi_op (sigma for operators of any symmetry): ',nTest,' (operator, CI) symmetry pairs'
  write(u6,'(A,ES10.3)') '   Bra.H Ket vs transition densities,  max. abs. difference ',dMaxE
  write(u6,'(A,ES10.3)') '   totally symmetric vs sg_h_psi(_ras), max. abs. difference ',dMaxS
  write(u6,'(A,ES10.3)') '   non-Hermitian (sg_h_psi_sq) vs densities, Hermitian vs sg_h_psi_op ',dMaxQ
  write(u6,'(A,I0)') '   failures: ',nFail
  call mma_deallocate(Rho2)
  call mma_deallocate(Rho1)
  call mma_deallocate(G2O)
  call mma_deallocate(G2L)
  call mma_deallocate(H1O)
  call mma_deallocate(H1L)
  call mma_deallocate(L2O)
  call mma_deallocate(TUVX_O)
  call mma_deallocate(TU_O)
  call mma_deallocate(TUVX)
  call mma_deallocate(TU)

contains

  ! G(iTri(ij,kl)) = (ij|kl) of TUVX_O (orbital order, square pair index)
  subroutine Sym_G2(G)
    real(kind=wp), intent(out) :: G(n2)
    integer(kind=iwp) :: a, bb, c, d
    do bb=1,nL
      do a=1,nL
        do d=1,nL
          do c=1,nL
            G(iTri((bb-1)*nL+a,(d-1)*nL+c)) = TUVX_O(iTri(iTri(a,bb),iTri(c,d)))
          end do
        end do
      end do
    end do
  end subroutine Sym_G2

  ! G (orbital order) from GL (level order)
  subroutine Map_G2(GL,G)
    real(kind=wp), intent(in) :: GL(n2)
    real(kind=wp), intent(out) :: G(n2)
    integer(kind=iwp) :: a, bb, c, d
    do bb=1,nL
      do a=1,nL
        do d=1,nL
          do c=1,nL
            G(iTri((L2O(bb)-1)*nL+L2O(a),(L2O(d)-1)*nL+L2O(c))) = GL(iTri((bb-1)*nL+a,(d-1)*nL+c))
          end do
        end do
      end do
    end do
  end subroutine Map_G2

end subroutine SG_Check_HPsi_Op

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
  use sguga, only: SGStruct, CIStruct, EXStruct, sg_epq_psi, SG_New_Tag, SG_No_Reuse
  use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
  use Constants, only: Zero, One, Half
  use Definitions, only: iwp, wp, u6

  implicit none
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(EXStruct), intent(inout) :: EXS
  integer(kind=iwp), intent(in) :: nCSFs, PsiSym, nTUVX_Tri, nTU_Tri
  real(kind=wp), intent(in) :: Psi(nCSFs), TUVX_Tri(nTUVX_Tri), TU_Tri(nTU_Tri)
  real(kind=wp), intent(out) :: Sigma(nCSFs)
  integer(kind=iwp) :: nOrb, iOrb, jOrb, kOrb, ijSym, SigmaSym, mCSFs, nP, p, q, q0, nB, nC, MaxMem, tPsi
  integer(kind=iwp), allocatable :: PI(:), PJ(:)
  real(kind=wp) :: TwoInt
  real(kind=wp), allocatable :: G(:,:)
  real(kind=wp), pointer :: D(:,:), W(:,:)
  logical(kind=iwp), save :: Warned = .false.

  Sigma(:) = Zero
  ! all E_pq|Psi> share the tag of Psi (intermediates reused by SG_Epq_Psi)
  tPsi = SG_New_Tag()
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
    MaxMem = MaxMem+Wrk_Held(1)
    if (real(mCSFs,kind=wp)*real(nP+1,kind=wp)+real(nP,kind=wp)**2 > 0.8_wp*real(MaxMem,kind=wp)) then
      ! D and W in chunks of nC pairs (both buffers and G in 80% of the memory)
      MaxMem = MaxMem+Wrk_Held(2)
      nC = int(0.8_wp*real(MaxMem,kind=wp)/real(2*mCSFs+nP,kind=wp),kind=iwp)
      if (.not. Warned) write(u6,'(1X,A,I0,A,I0,A,ES10.3,A,ES10.3,A,I0,A)') 'sg_h_psi: the block of ',nP,' pairs x ',mCSFs, &
        ' CSFs does not fit (',real(mCSFs,kind=wp)*real(nP+1,kind=wp),' words needed, ',real(MaxMem,kind=wp), &
        ' available); done in chunks of ',max(1,nC),' pairs'
      Warned = .true.
      if (nC >= 1) then
        call Block_Chunked()
      else
        call Block_Direct()
      end if
      cycle
    end if
    ! D(:,p) = S_p|Psi>, and the one-electron part (totally symmetric pairs)
    call Get_Wrk(Wrk_S1,mCSFs*nP,'SG_Wrk_S1')
    D(1:mCSFs,1:nP) => Wrk_S1(1:mCSFs*nP)
    D(:,:) = Zero
    do p=1,nP
      call SG_Epq_Psi(SGS,CIS,EXS,PI(p),PJ(p),One,PsiSym,Psi,D(:,p),CI_Tag=tPsi)
      if (PI(p) /= PJ(p)) call SG_Epq_Psi(SGS,CIS,EXS,PJ(p),PI(p),One,PsiSym,Psi,D(:,p),CI_Tag=tPsi)
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
    ! (D is in use: only the W buffer can be reused)
    MaxMem = MaxMem+Wrk_Held(2)
    nB = int(min(real(nP,kind=wp),max(One,0.8_wp*real(MaxMem,kind=wp)/real(mCSFs,kind=wp))),kind=iwp)
    call Get_Wrk(Wrk_S2,mCSFs*nB,'SG_Wrk_S2')
    W(1:mCSFs,1:nB) => Wrk_S2(1:mCSFs*nB)
    do q0=1,nP,nB
      nB = min(nB,nP-q0+1)
      call DGEMM_('N','T',mCSFs,nB,nP,One,D,mCSFs,G(q0,1),nP,Zero,W,mCSFs)
      do q=q0,q0+nB-1
        call SG_Epq_Psi(SGS,CIS,EXS,PI(q),PJ(q),One,SigmaSym,W(:,q-q0+1),Sigma,CI_Tag=SG_No_Reuse)
        if (PI(q) /= PJ(q)) call SG_Epq_Psi(SGS,CIS,EXS,PJ(q),PI(q),One,SigmaSym,W(:,q-q0+1),Sigma,CI_Tag=SG_No_Reuse)
      end do
    end do
    nullify(W)
    call mma_deallocate(G)
    nullify(D)
  end do
  call mma_deallocate(PJ)
  call mma_deallocate(PI)

contains

  ! the block ijSym in chunks of nC pairs: D(:,p) = S_p|Psi> for the pairs p
  ! of a chunk, then for every chunk of pairs q W = D G(q,p)^T and
  ! H|Psi> += S_q W(:,q): 2 n_pair (1 + n_chunk) applications of S instead of
  ! the n_pair**2 of Block_Direct, the same result.
  subroutine Block_Chunked()
    ! (local loop indices and work arrays; the host's are not touched)
    integer(kind=iwp) :: iP0, iQ0, cP, cQ, ip, iq, kk
    real(kind=wp) :: T2
    real(kind=wp), allocatable :: GC(:,:)
    real(kind=wp), pointer :: DC(:,:), WC(:,:)
    nC = min(nC,nP)
    call Get_Wrk(Wrk_S1,mCSFs*nC,'SG_Wrk_S1')
    call Get_Wrk(Wrk_S2,mCSFs*nC,'SG_Wrk_S2')
    DC(1:mCSFs,1:nC) => Wrk_S1(1:mCSFs*nC)
    WC(1:mCSFs,1:nC) => Wrk_S2(1:mCSFs*nC)
    call mma_allocate(GC,nC,nC,Label='GC_HPsi')
    do iP0=0,nP-1,nC
      cP = min(nC,nP-iP0)
      DC(:,1:cP) = Zero
      do ip=1,cP
        p = iP0+ip
        call SG_Epq_Psi(SGS,CIS,EXS,PI(p),PJ(p),One,PsiSym,Psi,DC(:,ip),CI_Tag=tPsi)
        if (PI(p) /= PJ(p)) call SG_Epq_Psi(SGS,CIS,EXS,PJ(p),PI(p),One,PsiSym,Psi,DC(:,ip),CI_Tag=tPsi)
        if (ijSym == 1) then
          T2 = Zero
          do kk=1,nOrb
            T2 = T2+TUVX_Tri(iTri(iTri(PI(p),kk),iTri(kk,PJ(p))))
          end do
          Sigma(:) = Sigma(:)+(TU_Tri(iTri(PI(p),PJ(p)))-Half*T2)*DC(:,ip)
        end if
      end do
      do iQ0=0,nP-1,nC
        cQ = min(nC,nP-iQ0)
        ! GC(ip,iq) = 1/2 (q|p)
        do iq=1,cQ
          q = iQ0+iq
          do ip=1,cP
            p = iP0+ip
            GC(ip,iq) = Half*TUVX_Tri(iTri(iTri(PI(q),PJ(q)),iTri(PI(p),PJ(p))))
          end do
        end do
        call DGEMM_('N','N',mCSFs,cQ,cP,One,DC,mCSFs,GC,nC,Zero,WC,mCSFs)
        do iq=1,cQ
          q = iQ0+iq
          call SG_Epq_Psi(SGS,CIS,EXS,PI(q),PJ(q),One,SigmaSym,WC(:,iq),Sigma,CI_Tag=SG_No_Reuse)
          if (PI(q) /= PJ(q)) call SG_Epq_Psi(SGS,CIS,EXS,PJ(q),PI(q),One,SigmaSym,WC(:,iq),Sigma,CI_Tag=SG_No_Reuse)
        end do
      end do
    end do
    call mma_deallocate(GC)
    nullify(WC)
    nullify(DC)
  end subroutine Block_Chunked

  ! the block ijSym directly: S_kl S_ij|Psi> for every pair of pairs
  subroutine Block_Direct()
    ! (local loop indices and work arrays; the host's are not touched)
    integer(kind=iwp) :: ip, ir, kk
    real(kind=wp) :: T2
    real(kind=wp), allocatable :: X(:), Y(:)
    call mma_allocate(X,mCSFs,Label='X_HPsi')
    call mma_allocate(Y,nCSFs,Label='Y_HPsi')
    do ip=1,nP
      X(:) = Zero
      call SG_Epq_Psi(SGS,CIS,EXS,PI(ip),PJ(ip),One,PsiSym,Psi,X,CI_Tag=tPsi)
      if (PI(ip) /= PJ(ip)) call SG_Epq_Psi(SGS,CIS,EXS,PJ(ip),PI(ip),One,PsiSym,Psi,X,CI_Tag=tPsi)
      if (ijSym == 1) then
        T2 = Zero
        do kk=1,nOrb
          T2 = T2+TUVX_Tri(iTri(iTri(PI(ip),kk),iTri(kk,PJ(ip))))
        end do
        Sigma(:) = Sigma(:)+(TU_Tri(iTri(PI(ip),PJ(ip)))-Half*T2)*X(:)
      end if
      do ir=1,nP
        Y(:) = Zero
        call SG_Epq_Psi(SGS,CIS,EXS,PI(ir),PJ(ir),One,SigmaSym,X,Y,CI_Tag=SG_No_Reuse)
        if (PI(ir) /= PJ(ir)) call SG_Epq_Psi(SGS,CIS,EXS,PJ(ir),PI(ir),One,SigmaSym,X,Y,CI_Tag=SG_No_Reuse)
        Sigma(:) = Sigma(:)+Half*TUVX_Tri(iTri(iTri(PI(ir),PJ(ir)),iTri(PI(ip),PJ(ip))))*Y(:)
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
  ! of E per block instead of n_pair**2/2. A block whose vectors do not fit
  ! is done in chunks of pairs (Block_Chunked); sg_h_psi_ras_direct only if
  ! not even two vectors fit.
  use Index_functions, only: iTri
  use symmetry_info, only: MUL
  use sguga, only: SGStruct, CIStruct, EXStruct, sg_epq_psi, SG_Level_Map, SG_New_Tag, SG_No_Reuse
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
  integer(kind=iwp) :: tPsi
  integer(kind=iwp) :: iA, iLev, iOrb, ijSym, jLev, jOrb, kLev, lLev, mCSFs, mLev, nLev, SigmaSym, nP, ip, iq, q0, nB, &
                       MaxMem, MaxNeed
  real(kind=wp) :: Fact, OneInt
  integer(kind=iwp), allocatable :: L2Orb(:), PI(:), PJ(:)
  real(kind=wp), allocatable :: G(:,:)
  real(kind=wp), pointer :: D(:,:), W(:,:)
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
    ! (the least: two chunks of one pair, D and W)
    MaxNeed = max(MaxNeed,2*CIS%nCSF(MUL(PsiSym,ijSym))+nP**2)
  end do
  call mma_maxDBLE(MaxMem)
  MaxMem = MaxMem+Wrk_Held(1)+Wrk_Held(2)
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
  ! all E_ij|Psi> share the tag of Psi
  tPsi = SG_New_Tag()
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
    call mma_maxDBLE(MaxMem)
    MaxMem = MaxMem+Wrk_Held(1)
    if (real(mCSFs,kind=wp)*real(nP+1,kind=wp)+real(nP,kind=wp)**2 > 0.8_wp*real(MaxMem,kind=wp)) then
      call Block_Chunked()
      cycle
    end if
    ! D(:,p) = E_ij|Psi> (restricted to the CSF space), and the one-electron part h'_ij E_ij
    call Get_Wrk(Wrk_S1,mCSFs*nP,'SG_Wrk_S1')
    D(1:mCSFs,1:nP) => Wrk_S1(1:mCSFs*nP)
    D(:,:) = Zero
    do ip=1,nP
      iLev = PI(ip)
      jLev = PJ(ip)
      call SG_Epq_Psi(SGS,CIS,EXS,iLev,jLev,One,PsiSym,Psi,D(:,ip),CI_Tag=tPsi)
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
    ! (D is in use: only the W buffer can be reused)
    MaxMem = MaxMem+Wrk_Held(2)
    nB = int(min(real(nP,kind=wp),max(One,0.8_wp*real(MaxMem,kind=wp)/real(mCSFs,kind=wp))),kind=iwp)
    call Get_Wrk(Wrk_S2,mCSFs*nB,'SG_Wrk_S2')
    W(1:mCSFs,1:nB) => Wrk_S2(1:mCSFs*nB)
    do q0=1,nP,nB
      nB = min(nB,nP-q0+1)
      call DGEMM_('N','T',mCSFs,nB,nP,One,D,mCSFs,G(q0,1),nP,Zero,W,mCSFs)
      do iq=q0,q0+nB-1
        if (Used(iq)) call SG_Epq_Psi(SGS,CIS,EXS,PI(iq),PJ(iq),One,SigmaSym,W(:,iq-q0+1),Sigma,CI_Tag=SG_No_Reuse)
      end do
    end do
    nullify(W)
    call mma_deallocate(Used)
    call mma_deallocate(G)
    nullify(D)
  end do
  call mma_deallocate(PJ)
  call mma_deallocate(PI)
  call mma_deallocate(L2Orb)

contains

  ! the block ijSym in chunks of nC pairs (D and W do not fit): D(:,p) for
  ! the pairs p of a chunk, then for every chunk of pairs q W = D G(q,p)^T
  ! and H|Psi> += E_q W(:,q); the same result with n_pair (1 + n_chunk)
  ! applications of E.
  subroutine Block_Chunked()
    ! (local loop indices and work arrays; the host's are not touched)
    integer(kind=iwp) :: nC, iP0, iQ0, cP, cQ, jp, jq, i_, j_, k_, l_, m_
    real(kind=wp) :: F_, H1
    real(kind=wp), allocatable :: GC(:,:)
    real(kind=wp), pointer :: DC(:,:), WC(:,:)
    logical(kind=iwp) :: Any_
    logical(kind=iwp), allocatable :: UsedC(:)
    MaxMem = MaxMem+Wrk_Held(2)
    nC = int(0.8_wp*real(MaxMem,kind=wp)/real(2*mCSFs+nP,kind=wp),kind=iwp)
    nC = max(1,min(nC,nP))
    call Get_Wrk(Wrk_S1,mCSFs*nC,'SG_Wrk_S1')
    call Get_Wrk(Wrk_S2,mCSFs*nC,'SG_Wrk_S2')
    DC(1:mCSFs,1:nC) => Wrk_S1(1:mCSFs*nC)
    WC(1:mCSFs,1:nC) => Wrk_S2(1:mCSFs*nC)
    call mma_allocate(GC,nC,nC,Label='GC_HPsi')
    call mma_allocate(UsedC,nP,Label='UsedC')
    ! E_kl is applied only if some (kl|ij) term has kl >= ij and k >= i
    UsedC(:) = .false.
    do jp=1,nP
      do jq=jp,nP
        if (PI(jq) >= PI(jp)) UsedC(jq) = .true.
      end do
    end do
    do iP0=0,nP-1,nC
      cP = min(nC,nP-iP0)
      DC(:,1:cP) = Zero
      do jp=1,cP
        i_ = PI(iP0+jp)
        j_ = PJ(iP0+jp)
        call SG_Epq_Psi(SGS,CIS,EXS,i_,j_,One,PsiSym,Psi,DC(:,jp),CI_Tag=tPsi)
        if (ijSym == 1) then
          H1 = TU_Tri(iTri(L2Orb(i_),L2Orb(j_)))
          do m_=1,i_-1
            iA = L2Orb(m_)
            H1 = H1-TUVX_Tri(iTri(iTri(L2Orb(i_),iA),iTri(iA,L2Orb(j_))))
          end do
          if (i_ > j_) then
            H1 = H1-TUVX_Tri(iTri(iTri(L2Orb(i_),L2Orb(i_)),iTri(L2Orb(i_),L2Orb(j_))))
          else if (i_ == j_) then
            H1 = H1-Half*TUVX_Tri(iTri(iTri(L2Orb(i_),L2Orb(i_)),iTri(L2Orb(i_),L2Orb(i_))))
          end if
          Sigma(:) = Sigma(:)+H1*DC(:,jp)
        end if
      end do
      do iQ0=0,nP-1,nC
        cQ = min(nC,nP-iQ0)
        ! GC(jp,jq) = f (kl|ij) for kl >= ij and k >= i, zero otherwise
        GC(:,:) = Zero
        Any_ = .false.
        do jq=1,cQ
          k_ = PI(iQ0+jq)
          l_ = PJ(iQ0+jq)
          do jp=1,cP
            if (iQ0+jq < iP0+jp) cycle
            i_ = PI(iP0+jp)
            j_ = PJ(iP0+jp)
            if (k_ < i_) cycle
            F_ = One
            if (iQ0+jq == iP0+jp) F_ = Half
            GC(jp,jq) = F_*TUVX_Tri(iTri(iTri(L2Orb(k_),L2Orb(l_)),iTri(L2Orb(i_),L2Orb(j_))))
            Any_ = .true.
          end do
        end do
        if (.not. Any_) cycle
        call DGEMM_('N','N',mCSFs,cQ,cP,One,DC,mCSFs,GC,nC,Zero,WC,mCSFs)
        do jq=1,cQ
          if (UsedC(iQ0+jq)) call SG_Epq_Psi(SGS,CIS,EXS,PI(iQ0+jq),PJ(iQ0+jq),One,SigmaSym,WC(:,jq),Sigma,CI_Tag=SG_No_Reuse)
        end do
      end do
    end do
    call mma_deallocate(UsedC)
    call mma_deallocate(GC)
    nullify(WC)
    nullify(DC)
  end subroutine Block_Chunked

end subroutine sg_h_psi_ras

subroutine sg_h_psi_ras_direct(SGS,CIS,EXS,Psi,nCSFs,PsiSym,Sigma,TUVX_Tri,nTUVX_Tri,TU_Tri,nTU_Tri)

  use Index_functions, only: iTri
  use symmetry_info, only: MUL
  use sguga, only: SGStruct, CIStruct, EXStruct, sg_epq_psi, SG_Level_Map, SG_New_Tag, SG_No_Reuse
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
  integer(kind=iwp) :: tPsi
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
  ! all E_ij|Psi> share the tag of Psi
  tPsi = SG_New_Tag()

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
      call SG_Epq_Psi(SGS,CIS,EXS,iLev,jLev,One,PsiSym,Psi,Eij_Psi,CI_Tag=tPsi)

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
          call SG_Epq_Psi(SGS,CIS,EXS,kLev,lLev,One,SigmaSym,Eij_Psi,Ekl_Eij_Psi(:,iBuff),CI_Tag=SG_No_Reuse)
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

  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
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
  mCSFs = maxval(CIS%nCSF)
  call mma_allocate(Psi,nCSFs,Label='Psi')
  call mma_allocate(Sigma,mCSFs,Label='Sigma')

#ifdef _TEST_
  ! turn off reuse of SGTMP
  ReUse_Save = EXS%Reuse_SGTMP
  EXS%Reuse_SGTMP = .false.
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
        Epq(iCSF) = dot_product(Sigma,Sigma)
      end do
      !Call RecPrt('EqpEpq(exact)',' ',Epq,1,nCSFs)

      call SG_EpqEqp_Diag(SGS,CIS,EXS, &
                          p,q, &
                          One, &
                          CSFSym, &
                          Psi, &
                          Sigma)
      !Call RecPrt('EqpEpq(Direct)',' ',Sigma,1,nCSFs)

      do iCSF=1,nCSFs
        if (abs(Epq(iCSF)-Sigma(iCSF)) > 1.0e-10_wp) then
          write(u6,*) 'p, q, iCSF=',p,q,iCSF
          write(u6,*) 'MidLev=',SGS%MidLev
          write(u6,*) 'Level=',SGS%Level(1:SGS%nLev)
          write(u6,*) 'L2Act=',SGS%L2Act(1:SGS%nLev)
          call RecPrt('EqpEpq(exact)',' ',Epq,1,nCSFs)
          call RecPrt('EqpEpq(Direct)',' ',Sigma,1,nCSFs)
          call Abend()
        end if
      end do

    end do
  end do

  call mma_deallocate(Eqp)
  call mma_deallocate(Epq)

  EXS%Reuse_SGTMP = ReUse_Save
  write(u6,*)
  write(u6,*) 'Test passed!'
  write(u6,*)
  !Call Abend()
#endif

  do iCSF=1,nCSFs

    Psi(:) = Zero
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

  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(EXStruct), intent(inout) :: EXS

  integer(kind=iwp), intent(in) :: nCSFs
  integer(kind=iwp), intent(in) :: CSFSym

  integer(kind=iwp), intent(in) :: nTUVX_Tri
  integer(kind=iwp), intent(in) :: nTU_Tri

  real(kind=wp), intent(in) :: TUVX_Tri(nTUVX_Tri)
  real(kind=wp), intent(in) :: TU_Tri(nTU_Tri)

  real(kind=wp), intent(out) :: Diag(nCSFs)

  real(kind=wp), allocatable :: Psi(:)
  real(kind=wp), allocatable :: Epp(:)
  real(kind=wp), allocatable :: Eqq(:)
  real(kind=wp), allocatable :: EqqEpp(:)
  real(kind=wp), allocatable :: EpqEqp(:)
  ! the integrals indexed by graph level (TU/TUVX: symmetry-blocked orbitals)
  real(kind=wp), allocatable :: TU_L(:), TUVX_L(:)

  integer(kind=iwp) :: p, q
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

    Diag(:) = Diag(:)+OneInt*Epp(:)

!
!-------------------------------------------------------
! PPPP contribution
!
! 1/2 (pp|pp) [Epp*Epp - Epp]
!-------------------------------------------------------
!
    ppqqTri = iTri(iTri(p,p),iTri(p,p))
    Coul = Half*TUVX_L(ppqqTri)

    Diag(:) = Diag(:)+Coul*(Epp(:)**2-Epp(:))

! PPQQ and PQQP contributions

    do q=p+1,SGS%nLev

      EqqEpp(:) = Zero
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

      Diag(:) = Diag(:)+TUVX_L(ppqqTri)*EqqEpp(:)

!
!-------------------------------------------------------
! PQQP contribution
!
! 1/2 (pq|qp)
! [ <EpqEqp> + <EqpEpq> - Epp - Eqq ]
!-------------------------------------------------------
!
      Eqq(:) = Zero
      call SG_Epq_Psi(SGS,CIS,EXS, &
                      q,q, &
                      One, &
                      CSFSym, &
                      Psi, &
                      Eqq)

      EpqEqp(:) = Zero
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
      Exch = Half*TUVX_L(pqqpTri)

      Diag(:) = Diag(:)+Exch*Two*(EpqEqp(:)-Eqq(:))
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

  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
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

    Psi(:) = Zero
    Sigma(:) = Zero

    Psi(iCSF) = One

    ! Construct H|I> = Sum_J <J|H|I> |J>

    call sg_h_psi(SGS,CIS,EXS, &
                  Psi,nCSFs,CSFSym, &
                  Sigma, &
                  TUVX_Tri,nTUVX_Tri, &
                  TU_Tri,nTU_Tri)

    do jCSF=1,nCSFs
      Full(jCSF,iCSF) = Sigma(jCSF)
    end do

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

!> The memory (words) the sigma (sg_h_psi, or sg_h_psi_ras with Ordered)
!> and the 2-RDM (sg_two_pdm) of CI vectors of symmetry PsiSym need for
!> their largest block of orbital pairs (D and W, A and B, G), so that the
!> callers (the Davidson vector store) leave it free.
function SG_Work_Need(SGS,CIS,PsiSym,Ordered) result(nNeed)
  use sguga, only: SGStruct, CIStruct
  use symmetry_info, only: MUL
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  integer(kind=iwp), intent(in) :: PsiSym
  logical(kind=iwp), intent(in) :: Ordered
  integer(kind=iwp) :: nNeed
  integer(kind=iwp) :: ijSym, i, j, nL, nT, nO, m
  nNeed = 0
  nL = SGS%nLev
  do ijSym=1,SGS%nSym
    nT = 0
    nO = 0
    do i=1,nL
      do j=1,nL
        if (MUL(SGS%ISM(i),SGS%ISM(j)) /= ijSym) cycle
        nO = nO+1
        if (j <= i) nT = nT+1
      end do
    end do
    m = CIS%nCSF(MUL(PsiSym,ijSym))
    if (Ordered) then
      nNeed = max(nNeed,2*m*nO+nO**2)
    else
      nNeed = max(nNeed,2*m*nT+nT**2)
    end if
    ! sg_two_pdm: A, B and their four overlap matrices
    nNeed = max(nNeed,2*m*nT+4*nT**2)
  end do
  ! the routines use a block only if it takes at most 80 % of the free memory
  nNeed = int(real(nNeed,kind=wp)/0.8_wp,kind=iwp)+1
end function SG_Work_Need

!> Free the work arrays kept between the calls.
subroutine SG_Sigma_Free()
  use stdalloc, only: mma_deallocate
  call mma_deallocate(Wrk_S1,safe='*')
  call mma_deallocate(Wrk_S2,safe='*')
end subroutine SG_Sigma_Free

!> Wrk has at least n elements (grown if needed, contents not kept).
subroutine Get_Wrk(Wrk,n,Label)
  use stdalloc, only: mma_allocate, mma_deallocate
  real(kind=wp), allocatable, intent(inout) :: Wrk(:)
  integer(kind=iwp), intent(in) :: n
  character(len=*), intent(in) :: Label
  if (allocated(Wrk)) then
    if (size(Wrk) >= n) return
    call mma_deallocate(Wrk)
  end if
  call mma_allocate(Wrk,max(1,n),Label=Label)
end subroutine Get_Wrk

!> The memory held by work array i (in words); the D blocks of the sigma
!> routines use Wrk_S1, W uses Wrk_S2 (W is sized while D is in use).
function Wrk_Held(i) result(n)
  integer(kind=iwp), intent(in) :: i
  integer(kind=iwp) :: n
  n = 0
  if (i == 1) then
    if (allocated(Wrk_S1)) n = size(Wrk_S1)
  else
    if (allocated(Wrk_S2)) n = size(Wrk_S2)
  end if
end function Wrk_Held

end module SGUGA_Sigma
