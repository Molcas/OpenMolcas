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

!> SGUGA: reduced density matrices of CI vectors of a split graph (one- and
!> two-particle, spin and transition densities).
module SGUGA_Densities

implicit none
private

public :: sg_one_pdm, sg_spin_pdm, sg_trans_pdm, sg_two_pdm

contains

subroutine sg_one_pdm(SGS,CIS,EXS,CI,lCI,ISYCI,D1MAT,lD1MAT)

  use Index_Functions, only: iTri
  use sguga, only: CIStruct, EXStruct, sg_epq_psi, SGStruct
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One
  use Definitions, only: wp, iwp

  implicit none
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(EXStruct), intent(inout) :: EXS
  integer(kind=iwp), intent(in) :: lCI, ISYCI, lD1MAT
  real(kind=wp), intent(in) :: CI(lCI)
  real(kind=wp), intent(out) :: D1MAT(lD1MAT)
  integer(kind=iwp) iOrb, jOrb, iSym, jSym
  real(kind=wp), allocatable :: SGM(:)
  real(kind=wp), parameter :: CPQ = One

  call mma_allocate(SGM,lCI,Label='SGM')

  D1MAT(:) = Zero
  do iOrb=1,SGS%nLev
    iSym = SGS%ISM(iOrb)
    do jOrb=1,iOrb
      jSym = SGS%ISM(jOrb)
      if (jSym /= iSym) cycle

      SGM(:) = Zero
      call SG_Epq_Psi(SGS,CIS,EXS,iOrb,jOrb,CPQ,ISYCI,CI,SGM)

      D1MAT(iTri(iOrb,jOrb)) = dot_product(CI,SGM)

    end do
  end do

  call mma_deallocate(SGM)

end subroutine sg_one_pdm

!> One-particle spin density D^S_pq = <S,S| E^a_pq - E^b_pq |S,S> (levels,
!> triangular packing as D1MAT of sg_one_pdm) from the spin-free densities,
!>
!>   D^S_pq = 1/(S+1) [ (2 - N/2) D_pq - Sum_r Gamma_{pr,rq} ],
!>   Sum_r Gamma_{pr,rq} = Sum_r <E_pr E_rq> - n D_pq,
!>
!> N active electrons, n levels, D1MAT the spin-free density (sg_one_pdm).
!> <E_pr E_rq> = (E_rp Psi).(E_rq Psi) is exact in a RAS graph (levels RAS1
!> < RAS2 < RAS3) for r <= max(p,q); otherwise <E_rq E_pr> + D_pq -
!> delta_pq D_rr = (E_qr Psi).(E_pr Psi) + ... is used, which is.
subroutine sg_spin_pdm(SGS,CIS,EXS,CI,lCI,ISYCI,D1MAT,SDMAT,lD1MAT)

  use Index_Functions, only: iTri
  use Symmetry_Info, only: Mul
  use sguga, only: CIStruct, EXStruct, sg_epq_psi, SGStruct
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One, Two, Half
  use Definitions, only: wp, iwp

  implicit none
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(EXStruct), intent(inout) :: EXS
  integer(kind=iwp), intent(in) :: lCI, ISYCI, lD1MAT
  real(kind=wp), intent(in) :: CI(lCI), D1MAT(lD1MAT)
  real(kind=wp), intent(out) :: SDMAT(lD1MAT)
  integer(kind=iwp) :: nL, p, q, r, MaxDim, iSymA
  real(kind=wp) :: S, X, Fact
  real(kind=wp), allocatable :: A(:,:), B(:,:), G(:)

  SDMAT(:) = Zero
  if (SGS%iSpin == 1) return
  nL = SGS%nLev
  S = Half*real(SGS%iSpin-1,kind=wp)
  MaxDim = maxval(CIS%nCSF(:))
  call mma_allocate(A,MaxDim,nL,Label='SD_A')
  call mma_allocate(B,MaxDim,nL,Label='SD_B')
  ! G = Sum_r <E_pr E_rq>, triangular
  call mma_allocate(G,lD1MAT,Label='SD_G')
  G(:) = Zero
  do r=1,nL
    ! A(:,q) = E_rq Psi, B(:,q) = E_qr Psi
    do q=1,nL
      A(:,q) = Zero
      B(:,q) = Zero
      call SG_Epq_Psi(SGS,CIS,EXS,r,q,One,ISYCI,CI,A(:,q))
      call SG_Epq_Psi(SGS,CIS,EXS,q,r,One,ISYCI,CI,B(:,q))
    end do
    do p=1,nL
      do q=1,p
        if (SGS%ISm(p) /= SGS%ISm(q)) cycle
        ! E_rq Psi has the symmetry of E_rp Psi
        iSymA = Mul(ISYCI,Mul(SGS%ISm(r),SGS%ISm(q)))
        if (r <= max(p,q)) then
          X = dot_product(A(1:CIS%nCSF(iSymA),p),A(1:CIS%nCSF(iSymA),q))
        else
          X = dot_product(B(1:CIS%nCSF(iSymA),q),B(1:CIS%nCSF(iSymA),p))+D1MAT(iTri(p,q))
          if (p == q) X = X-D1MAT(iTri(r,r))
        end if
        G(iTri(p,q)) = G(iTri(p,q))+X
      end do
    end do
  end do
  Fact = One/(S+One)
  do p=1,nL
    do q=1,p
      if (SGS%ISm(p) /= SGS%ISm(q)) cycle
      SDMAT(iTri(p,q)) = Fact*((Two-Half*real(SGS%nActEl,kind=wp))*D1MAT(iTri(p,q))- &
                               (G(iTri(p,q))-real(nL,kind=wp)*D1MAT(iTri(p,q))))
    end do
  end do
  call mma_deallocate(G)
  call mma_deallocate(B)
  call mma_deallocate(A)

end subroutine sg_spin_pdm

!> Transition densities between the CSF vectors Bra and Ket of symmetry
!> ISYCI (two states of the same spin S, M_S = S), square over the levels,
!> T1DM(p,q) = <Bra|E_pq|Ket>, and, if Spin, the transition spin density
!>
!>   T1SDM(p,q) = 1/(S+1) [ (2 - N/2) T_pq - Sum_r <Bra|E_pr E_rq|Ket> + n T_pq ]
!>
!> (as in sg_spin_pdm: S_+ annihilates both states, so <Bra|S.T|Ket> =
!> (S+1) <Bra|T_z|Ket>). The products are taken in the order that is exact
!> in a RAS graph: (E_rp Bra).(E_rq Ket) for r <= max(p,q), otherwise
!> (E_qr Bra).(E_pr Ket) + T_pq - delta_pq T_rr.
subroutine sg_trans_pdm(SGS,CIS,EXS,Bra,Ket,lCI,ISYCI,nL,T1DM,T1SDM,Spin)

  use Symmetry_Info, only: Mul
  use sguga, only: CIStruct, EXStruct, sg_epq_psi, SGStruct
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One, Two, Half
  use Definitions, only: wp, iwp

  implicit none
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(EXStruct), intent(inout) :: EXS
  integer(kind=iwp), intent(in) :: lCI, ISYCI, nL
  real(kind=wp), intent(in) :: Bra(lCI), Ket(lCI)
  real(kind=wp), intent(out) :: T1DM(nL,nL), T1SDM(nL,nL)
  logical(kind=iwp), intent(in) :: Spin
  integer(kind=iwp) :: p, q, r, MaxDim, nS
  real(kind=wp) :: S, X
  real(kind=wp), allocatable :: Sgm(:), AB(:,:), AK(:,:), BB(:,:), BK(:,:), G(:,:)

  T1DM(:,:) = Zero
  T1SDM(:,:) = Zero
  MaxDim = maxval(CIS%nCSF(:))
  call mma_allocate(Sgm,MaxDim,Label='TD_Sgm')
  ! T_pq = Bra.(E_pq Ket), totally symmetric part only
  do p=1,nL
    do q=1,nL
      if (SGS%ISm(p) /= SGS%ISm(q)) cycle
      Sgm(:) = Zero
      call SG_Epq_Psi(SGS,CIS,EXS,p,q,One,ISYCI,Ket,Sgm)
      T1DM(p,q) = dot_product(Bra(1:lCI),Sgm(1:lCI))
    end do
  end do
  call mma_deallocate(Sgm)
  if ((.not. Spin) .or. (SGS%iSpin == 1)) return

  S = Half*real(SGS%iSpin-1,kind=wp)
  call mma_allocate(AB,MaxDim,nL,Label='TD_AB')
  call mma_allocate(AK,MaxDim,nL,Label='TD_AK')
  call mma_allocate(BB,MaxDim,nL,Label='TD_BB')
  call mma_allocate(BK,MaxDim,nL,Label='TD_BK')
  call mma_allocate(G,nL,nL,Label='TD_G')
  G(:,:) = Zero
  do r=1,nL
    ! AB(:,p) = E_rp Bra, AK(:,q) = E_rq Ket, BB(:,q) = E_qr Bra, BK(:,p) = E_pr Ket
    do q=1,nL
      AB(:,q) = Zero
      AK(:,q) = Zero
      BB(:,q) = Zero
      BK(:,q) = Zero
      call SG_Epq_Psi(SGS,CIS,EXS,r,q,One,ISYCI,Bra,AB(:,q))
      call SG_Epq_Psi(SGS,CIS,EXS,r,q,One,ISYCI,Ket,AK(:,q))
      call SG_Epq_Psi(SGS,CIS,EXS,q,r,One,ISYCI,Bra,BB(:,q))
      call SG_Epq_Psi(SGS,CIS,EXS,q,r,One,ISYCI,Ket,BK(:,q))
    end do
    do p=1,nL
      do q=1,nL
        if (SGS%ISm(p) /= SGS%ISm(q)) cycle
        nS = CIS%nCSF(Mul(ISYCI,Mul(SGS%ISm(r),SGS%ISm(q))))
        if (r <= max(p,q)) then
          X = dot_product(AB(1:nS,p),AK(1:nS,q))
        else
          X = dot_product(BB(1:nS,q),BK(1:nS,p))+T1DM(p,q)
          if (p == q) X = X-T1DM(r,r)
        end if
        G(p,q) = G(p,q)+X
      end do
    end do
  end do
  do p=1,nL
    do q=1,nL
      if (SGS%ISm(p) /= SGS%ISm(q)) cycle
      T1SDM(p,q) = ((Two-Half*real(SGS%nActEl,kind=wp))*T1DM(p,q)-G(p,q)+real(nL,kind=wp)*T1DM(p,q))/(S+One)
    end do
  end do
  call mma_deallocate(G)
  call mma_deallocate(BK)
  call mma_deallocate(BB)
  call mma_deallocate(AK)
  call mma_deallocate(AB)

end subroutine sg_trans_pdm

!#define _DEBUGPRINT_
subroutine sg_two_pdm(SGS,CIS,EXS,Psi,nCSFs,PsiSym,P,PA,nP)
  ! The symmetric and antisymmetric 2-RDM (P, PA) of Psi, as in
  ! sg_two_pdm_direct (see there for the formulas), with the vectors
  ! E_ij|Psi> and E_ji|Psi> (i >= j) of a symmetry block of pairs computed
  ! once (A, B) and all their overlaps from four DGEMMs, instead of E_lk|Psi>
  ! again for every pair ij. Falls back to sg_two_pdm_direct if the vectors
  ! of a block do not fit in memory.
  use Index_functions, only: iTri
  use Symmetry_Info, only: Mul
  use sguga, only: CIStruct, EXStruct, sg_epq_psi, SGStruct
  use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
  use Constants, only: Zero, One, Half
  use Definitions, only: wp, iwp
  implicit none
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(EXStruct), intent(inout) :: EXS
  integer(kind=iwp), intent(in) :: nCSFs, PsiSym, nP
  real(kind=wp), intent(in) :: Psi(nCSFs)
  real(kind=wp), intent(inout) :: P(nP), PA(nP)
  integer(kind=iwp) :: nL, ijSym, iOrb, jOrb, kOrb, lOrb, lOrb_Min, mCSFs, nPr, ip, iq, ijOrb, klOrb, klijOrb, ikljOrb, &
                       MaxMem, MaxNeed
  integer(kind=iwp), allocatable :: PI(:), PJ(:), Loc(:,:)
  real(kind=wp) :: D_ij, P_klij
  real(kind=wp), allocatable :: A(:,:), B(:,:), AA(:,:), AB(:,:), BA(:,:), BB(:,:)

  nL = SGS%nLev
  ! memory: the vectors of the largest block
  MaxNeed = 0
  do ijSym=1,SGS%nSym
    nPr = 0
    do iOrb=1,nL
      do jOrb=1,iOrb
        if (Mul(SGS%ISM(iOrb),SGS%ISM(jOrb)) == ijSym) nPr = nPr+1
      end do
    end do
    MaxNeed = max(MaxNeed,2*CIS%nCSF(Mul(PsiSym,ijSym))*nPr+4*nPr**2)
  end do
  call mma_maxDBLE(MaxMem)
  if (real(MaxNeed,kind=wp) > 0.8_wp*real(MaxMem,kind=wp)) then
    call sg_two_pdm_direct(SGS,CIS,EXS,Psi,nCSFs,PsiSym,P,PA,nP)
    return
  end if

  P(:) = Zero
  PA(:) = Zero
  call mma_allocate(PI,max(1,nL*(nL+1)/2),Label='PI')
  call mma_allocate(PJ,max(1,nL*(nL+1)/2),Label='PJ')
  call mma_allocate(Loc,nL,nL,Label='Loc')
  do ijSym=1,SGS%nSym
    nPr = 0
    Loc(:,:) = 0
    do iOrb=1,nL
      do jOrb=1,iOrb
        if (Mul(SGS%ISM(iOrb),SGS%ISM(jOrb)) /= ijSym) cycle
        nPr = nPr+1
        PI(nPr) = iOrb
        PJ(nPr) = jOrb
        Loc(iOrb,jOrb) = nPr
      end do
    end do
    mCSFs = CIS%nCSF(Mul(PsiSym,ijSym))
    if ((nPr == 0) .or. (mCSFs == 0)) cycle
    ! A(:,p) = E_ij|Psi>, B(:,p) = E_ji|Psi> (zero for i = j)
    call mma_allocate(A,mCSFs,nPr,Label='A_2PDM')
    call mma_allocate(B,mCSFs,nPr,Label='B_2PDM')
    A(:,:) = Zero
    B(:,:) = Zero
    do ip=1,nPr
      call SG_Epq_Psi(SGS,CIS,EXS,PI(ip),PJ(ip),One,PsiSym,Psi,A(:,ip))
      if (PI(ip) /= PJ(ip)) call SG_Epq_Psi(SGS,CIS,EXS,PJ(ip),PI(ip),One,PsiSym,Psi,B(:,ip))
    end do
    call mma_allocate(AA,nPr,nPr,Label='AA')
    call mma_allocate(AB,nPr,nPr,Label='AB')
    call mma_allocate(BA,nPr,nPr,Label='BA')
    call mma_allocate(BB,nPr,nPr,Label='BB')
    call DGEMM_('T','N',nPr,nPr,mCSFs,One,A,mCSFs,A,mCSFs,Zero,AA,nPr)
    call DGEMM_('T','N',nPr,nPr,mCSFs,One,A,mCSFs,B,mCSFs,Zero,AB,nPr)
    call DGEMM_('T','N',nPr,nPr,mCSFs,One,B,mCSFs,A,mCSFs,Zero,BA,nPr)
    call DGEMM_('T','N',nPr,nPr,mCSFs,One,B,mCSFs,B,mCSFs,Zero,BB,nPr)
    do ip=1,nPr
      iOrb = PI(ip)
      jOrb = PJ(ip)
      ijOrb = iTri(iOrb,jOrb)
      D_ij = Zero
      if (ijSym == 1) D_ij = dot_product(Psi,A(:,ip))
      ! the -d_il E_jk term, as a -d_kl E_ij term of P_iklj
      if (ijSym == 1) then
        do kOrb=1,iOrb
          lOrb = kOrb
          ikljOrb = iTri(iTri(iOrb,kOrb),iTri(lOrb,jOrb))
          P(ikljOrb) = P(ikljOrb)-D_ij
          if ((iOrb /= kOrb) .and. (lOrb > jOrb)) PA(ikljOrb) = PA(ikljOrb)-D_ij
          if ((iOrb /= kOrb) .and. (lOrb < jOrb)) PA(ikljOrb) = PA(ikljOrb)+D_ij
        end do
      end if
      ! k >= l, (k,l) >= (i,j): <Psi|E_kl = (E_lk|Psi>)^T, E_lk|Psi> = B(:,q) (k /= l) or A(:,q) (k = l)
      do kOrb=iOrb,nL
        lOrb_Min = 1
        if (kOrb == iOrb) lOrb_Min = jOrb
        do lOrb=lOrb_Min,kOrb
          iq = Loc(kOrb,lOrb)
          if (iq == 0) cycle
          klOrb = iTri(kOrb,lOrb)
          klijOrb = iTri(klOrb,ijOrb)
          if (kOrb /= lOrb) then
            P_klij = BA(iq,ip)
          else
            P_klij = AA(iq,ip)
          end if
          P(klijOrb) = P(klijOrb)+P_klij
          if ((iOrb /= jOrb) .and. (kOrb /= lOrb)) PA(klijOrb) = PA(klijOrb)+P_klij
          if (iOrb /= jOrb) then
            if (kOrb /= lOrb) then
              P_klij = BB(iq,ip)
            else
              P_klij = AB(iq,ip)
            end if
            P(klijOrb) = P(klijOrb)+P_klij
            if (kOrb /= lOrb) PA(klijOrb) = PA(klijOrb)-P_klij
          end if
        end do
      end do
    end do
    call mma_deallocate(BB)
    call mma_deallocate(BA)
    call mma_deallocate(AB)
    call mma_deallocate(AA)
    call mma_deallocate(B)
    call mma_deallocate(A)
  end do
  call mma_deallocate(Loc)
  call mma_deallocate(PJ)
  call mma_deallocate(PI)
  P(:) = Half*P(:)
  PA(:) = Half*PA(:)
end subroutine sg_two_pdm

subroutine sg_two_pdm_direct(SGS,CIS,EXS,Psi,nCSFs,PsiSym,P,PA,nP)

  use Index_functions, only: iTri
  use Symmetry_Info, only: Mul
  use sguga, only: CIStruct, EXStruct, sg_epq_psi, SGStruct
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One, Half
  use Definitions, only: wp, iwp

  implicit none
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(EXStruct), intent(inout) :: EXS
  integer(kind=iwp), intent(in) :: nCSFs, PsiSym, nP
  real(kind=wp), intent(in) :: Psi(nCSFs)
  real(kind=wp), intent(inout) :: P(nP), PA(nP)
  integer(kind=iwp) :: ijOrb, ijSym, ikljOrb, ikOrb, iOrb, iSym, jOrb, jSym, klijOrb, klOrb, klSym, kOrb, kSym, ljOrb, lOrb, &
                       lOrb_Min, lSym, MaxDim, mCSFs
  real(kind=wp) :: D_ij, P_klij
  real(kind=wp), pointer :: Eij_Psi(:), Eji_Psi(:), Elk_Psi(:)
  real(kind=wp), allocatable, target :: Eij_Psi_X(:), Elk_Psi_X(:), Eji_Psi_X(:)
  real(kind=wp), parameter :: CPQ = One

  MaxDim = maxval(CIS%nCSF(:))
  call mma_allocate(Eij_Psi_X,MaxDim,Label='Eij_Psi_X')
  call mma_allocate(Eji_Psi_X,MaxDim,Label='Eji_Psi_X')
  call mma_allocate(Elk_Psi_X,MaxDim,Label='Elk_Psi_X')

  P(:) = Zero
  PA(:) = Zero
  D_ij = Zero

  ! Compute P_lk,ij = <Psi|e_lkij|Psi> + <Psi|e_lkji|Psi>
  ! Where e_lkij = E_lk E_ij - delta_ki E_lj
  ! Hence,
  !
  ! P_lk,ij = sum_m <Psi|E_lk|m>(<m|Eij|Psi>+<m|E_ji|Psi>) - delta_ki D_lj -delta_kj D_li
  !
  ! In the symmetrization set we use the fact that for none zero elements we have that
  ! <m|E_ij|m'>=<m'|E_ji|m>
  ! This is trivial in the case of no symmetry. However, in the case of symmetry we have to be careful.
  !
  ! Let say that |m> belongs to the set of CSFs that are included in |Psi>, while |m'> doesn't. Then
  ! there is no symmetrization to be achieved. We need to keep track on this in Mk_Eij_Psi.
  !
  ! The size of the outer loop is nOrb*(nOrb+1)/2, i.e. nOrb=10 => 55 tasks to parallelize over.

  do iOrb=1,SGS%nLev
    iSym = SGS%ISM(iOrb)
    do jOrb=1,iOrb
      ijOrb = iTri(iOrb,jOrb)
      jSym = SGS%ISM(jOrb)
      ijSym = Mul(iSym,jSym)
      mCSFs = CIS%nCSF(Mul(PsiSym,ijSym))

      Eij_Psi(1:mCSFs) => Eij_Psi_X
      Eji_Psi(1:mCSFs) => Eji_Psi_X
      Elk_Psi(1:mCSFs) => Elk_Psi_X

      Eij_Psi(:) = Zero
      Eji_Psi(:) = Zero
      ! Operate with E_ij on |Psi> and produce E_ij|Psi>
      call SG_Epq_Psi(SGS,CIS,EXS,iOrb,jOrb,CPQ,PsiSym,Psi,Eij_Psi)
      ! Compute Dij to be distributed below
      if (iSym == jSym) D_ij = dot_product(Psi,Eij_Psi)
      if (iOrb /= jOrb) call SG_Epq_Psi(SGS,CIS,EXS,jOrb,iOrb,CPQ,PsiSym,Psi,Eji_Psi)

      ! Note, in the case of a RASSCF the resulting sigma vector is incomplete. For the case of E_ij
      ! the complete expansion of E_ij|Psi> would require CSFs that are outside the set of CSFs defining the
      ! RAS expansion. However, if we cap this sigma vector with vectors that have zero coefficient for these
      ! external CSFs we experience no errors.

      do kOrb=1,iOrb
        kSym = SGS%ISM(kOrb)
        lOrb = kOrb
        lSym = SGS%ISM(lOrb)
        klSym = Mul(kSym,lSym)
        if (ijSym /= klSym) cycle

        ! Add the -d_il E_jk term, here in the form of a -d_kl E_ij term of Piklj
        ikOrb = iTri(iOrb,kOrb)
        ljOrb = iTri(lOrb,jOrb)
        ikljOrb = iTri(ikOrb,ljOrb)
        P(ikljOrb) = P(ikljOrb)-D_ij
        if ((iOrb /= kOrb) .and. (lOrb > jOrb)) PA(ikljOrb) = PA(ikljOrb)-D_ij
        if ((iOrb /= kOrb) .and. (lOrb < jOrb)) PA(ikljOrb) = PA(ikljOrb)+D_ij
      end do

      ! kOrb>=lOrb
      do kOrb=iOrb,SGS%nLev
        kSym = SGS%ISM(kOrb)
        lOrb_Min = 1
        if (kOrb == iOrb) lOrb_Min = jOrb
        do lOrb=lOrb_Min,kOrb
          klOrb = iTri(kOrb,lOrb)
          lSym = SGS%ISM(lOrb)
          klSym = Mul(kSym,lSym)
          if (ijSym /= klSym) cycle

          ! Operate with E_lk on |Psi> and produce E_lk|Psi>. Note, since k>=l E_lk|Psi>
          ! is completely definded by the CSFs of the RASSCF space.

          Elk_Psi(1:mCSFs) = Zero
          ! <Psi|E_kl is computed as E_lk|Psi>
          call SG_Epq_Psi(SGS,CIS,EXS,lOrb,kOrb,CPQ,PsiSym,Psi,Elk_Psi)
          ! cap with <Psi|E_kl on E_ij|Psi>, contribution to Pklij
          P_klij = dot_product(Elk_Psi,Eij_Psi)

          klijOrb = iTri(klOrb,ijOrb)
          P(klijOrb) = P(klijOrb)+P_klij
          if ((iOrb /= jOrb) .and. (kOrb /= lOrb)) PA(klijOrb) = PA(klijOrb)+P_klij

          if (iOrb /= jOrb) then
            P_klij = dot_product(Elk_Psi,Eji_Psi)
            P(klijOrb) = P(klijOrb)+P_klij
            if (kOrb /= lOrb) PA(klijOrb) = PA(klijOrb)-P_klij
          end if

        end do
      end do
      nullify(Eij_Psi,Eji_Psi,Elk_Psi)

    end do
  end do

  call mma_deallocate(Elk_Psi_X)
  call mma_deallocate(Eij_Psi_X)
  call mma_deallocate(Eji_Psi_X)

  P = Half*P
  PA = Half*PA

end subroutine sg_two_pdm_direct

end module SGUGA_Densities
