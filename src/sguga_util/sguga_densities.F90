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

use Definitions, only: wp, iwp

implicit none
private

! work arrays of the routines below, kept between the calls (no new memory
! pages in every call for large CI spaces); freed by SG_Densities_Free
real(kind=wp), allocatable, target :: Wrk_D1(:), Wrk_D2(:)

public :: SG_Check_Trans_Two_PDM, SG_Densities_Free, sg_one_pdm, sg_spin_pdm, sg_trans_pdm, sg_trans_two_pdm, sg_two_pdm

contains

subroutine sg_one_pdm(SGS,CIS,EXS,CI,lCI,ISYCI,D1MAT,lD1MAT)

  use Index_Functions, only: iTri
  use sguga, only: CIStruct, EXStruct, sg_epq_psi, SG_New_Tag, SGStruct
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

  integer(kind=iwp) :: tCI
  tCI = SG_New_Tag()
  call mma_allocate(SGM,lCI,Label='SGM')

  D1MAT(:) = Zero
  do iOrb=1,SGS%nLev
    iSym = SGS%ISM(iOrb)
    do jOrb=1,iOrb
      jSym = SGS%ISM(jOrb)
      if (jSym /= iSym) cycle

      SGM(:) = Zero
      call SG_Epq_Psi(SGS,CIS,EXS,iOrb,jOrb,CPQ,ISYCI,CI,SGM,CI_Tag=tCI)

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
  use sguga, only: CIStruct, EXStruct, SG_New_Tag, SGStruct
  use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
  use Constants, only: Zero, One, Two, Half
  use Definitions, only: wp, iwp

  implicit none
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(EXStruct), intent(inout) :: EXS
  integer(kind=iwp), intent(in) :: lCI, ISYCI, lD1MAT
  real(kind=wp), intent(in) :: CI(lCI), D1MAT(lD1MAT)
  real(kind=wp), intent(out) :: SDMAT(lD1MAT)
  integer(kind=iwp) :: nL, p, q, r, iSymO, iSymA, m, nA, nC, jq, ip, nq, np, i, j, MaxMem
  real(kind=wp) :: S, X, Fact
  integer(kind=iwp), allocatable :: Lv(:)
  real(kind=wp), allocatable :: G(:)
  real(kind=wp), pointer :: AQ(:,:), BQ(:,:), AP(:,:), BP(:,:)

  integer(kind=iwp) :: tCI
  tCI = SG_New_Tag()
  SDMAT(:) = Zero
  if (SGS%iSpin == 1) return
  nL = SGS%nLev
  S = Half*real(SGS%iSpin-1,kind=wp)
  call mma_maxDBLE(MaxMem)
  MaxMem = (4*(MaxMem+Wrk_Held()))/5
  call mma_allocate(Lv,nL,Label='SD_Lv')
  ! G = Sum_r <E_pr E_rq>, triangular
  call mma_allocate(G,lD1MAT,Label='SD_G')
  G(:) = Zero
  do r=1,nL
    ! the levels q of one symmetry, in chunks of nC: A(:,k) = E_rq Psi,
    ! B(:,k) = E_qr Psi (all of the symmetry at once if they fit, otherwise
    ! pairs of chunks, the second one recomputed)
    do iSymO=1,SGS%nSym
      m = 0
      do p=1,nL
        if (SGS%ISm(p) /= iSymO) cycle
        m = m+1
        Lv(m) = p
      end do
      if (m == 0) cycle
      ! E_rq Psi has the symmetry of E_rp Psi
      iSymA = Mul(ISYCI,Mul(SGS%ISm(r),iSymO))
      nA = CIS%nCSF(iSymA)
      if (nA == 0) cycle
      if (2*nA*m <= MaxMem) then
        nC = m
      else
        nC = max(1,min(m,MaxMem/(4*nA)))
      end if
      call Get_Wrk(Wrk_D1,2*nA*nC,'SG_Wrk_D1')
      AQ(1:nA,1:nC) => Wrk_D1(1:nA*nC)
      BQ(1:nA,1:nC) => Wrk_D1(nA*nC+1:2*nA*nC)
      if (nC < m) then
        call Get_Wrk(Wrk_D2,2*nA*nC,'SG_Wrk_D2')
        AP(1:nA,1:nC) => Wrk_D2(1:nA*nC)
        BP(1:nA,1:nC) => Wrk_D2(nA*nC+1:2*nA*nC)
      end if
      do jq=1,m,nC
        nq = min(nC,m-jq+1)
        call SD_Fill(SGS,CIS,EXS,r,Lv(jq),nq,ISYCI,CI,nA,AQ,BQ,tCI)
        do ip=jq,m,nC
          np = min(nC,m-ip+1)
          if (ip == jq) then
            AP => AQ
            BP => BQ
          else
            AP(1:nA,1:nC) => Wrk_D2(1:nA*nC)
            BP(1:nA,1:nC) => Wrk_D2(nA*nC+1:2*nA*nC)
            call SD_Fill(SGS,CIS,EXS,r,Lv(ip),np,ISYCI,CI,nA,AP,BP,tCI)
          end if
          do i=1,np
            p = Lv(ip+i-1)
            do j=1,nq
              q = Lv(jq+j-1)
              if (q > p) exit
              if (r <= p) then
                X = dot_product(AP(:,i),AQ(:,j))
              else
                X = dot_product(BQ(:,j),BP(:,i))+D1MAT(iTri(p,q))
                if (p == q) X = X-D1MAT(iTri(r,r))
              end if
              G(iTri(p,q)) = G(iTri(p,q))+X
            end do
          end do
        end do
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
  call mma_deallocate(Lv)
  nullify(AQ,BQ,AP,BP)
  call Wrk_Trim()

end subroutine sg_spin_pdm

!> A(:,k) = E_rq Psi, B(:,k) = E_qr Psi for the n levels q = Lv(1:n)
!> (for sg_spin_pdm).
subroutine SD_Fill(SGS,CIS,EXS,r,Lv,n,ISYCI,CI,nA,A,B,tCI)

  use sguga, only: CIStruct, EXStruct, sg_epq_psi, SGStruct
  use Constants, only: Zero, One
  use Definitions, only: wp, iwp

  implicit none
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(EXStruct), intent(inout) :: EXS
  integer(kind=iwp), intent(in) :: r, n, Lv(n), ISYCI, nA, tCI
  real(kind=wp), intent(in) :: CI(*)
  real(kind=wp), intent(out) :: A(nA,*), B(nA,*)
  integer(kind=iwp) :: k

  do k=1,n
    A(:,k) = Zero
    B(:,k) = Zero
    call SG_Epq_Psi(SGS,CIS,EXS,r,Lv(k),One,ISYCI,CI,A(:,k),CI_Tag=tCI)
    call SG_Epq_Psi(SGS,CIS,EXS,Lv(k),r,One,ISYCI,CI,B(:,k),CI_Tag=tCI)
  end do

end subroutine SD_Fill

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
  use sguga, only: CIStruct, EXStruct, sg_epq_psi, SG_New_Tag, SGStruct
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

  integer(kind=iwp) :: tBra, tKet
  tBra = SG_New_Tag()
  tKet = SG_New_Tag()
  T1DM(:,:) = Zero
  T1SDM(:,:) = Zero
  MaxDim = maxval(CIS%nCSF(:))
  call mma_allocate(Sgm,MaxDim,Label='TD_Sgm')
  ! T_pq = Bra.(E_pq Ket), totally symmetric part only
  do p=1,nL
    do q=1,nL
      if (SGS%ISm(p) /= SGS%ISm(q)) cycle
      Sgm(:) = Zero
      call SG_Epq_Psi(SGS,CIS,EXS,p,q,One,ISYCI,Ket,Sgm,CI_Tag=tKet)
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
      call SG_Epq_Psi(SGS,CIS,EXS,r,q,One,ISYCI,Bra,AB(:,q),CI_Tag=tBra)
      call SG_Epq_Psi(SGS,CIS,EXS,r,q,One,ISYCI,Ket,AK(:,q),CI_Tag=tKet)
      call SG_Epq_Psi(SGS,CIS,EXS,q,r,One,ISYCI,Bra,BB(:,q),CI_Tag=tBra)
      call SG_Epq_Psi(SGS,CIS,EXS,q,r,One,ISYCI,Ket,BK(:,q),CI_Tag=tKet)
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

!> The spin-free 1- and 2-body transition densities between Bra (symmetry
!> BraSym) and Ket (symmetry KetSym), in level order:
!>   Rho1(i,l) = <Bra|E_il|Ket>,
!>   Rho2(iTri(ij,kl)) = <Bra|E_ij E_kl - delta_jk E_il|Ket>,
!> with the square pair index ij = (j-1)*nL+i and the packing over
!> ij >= kl (e_ijkl = e_klij). Neither is symmetric for Bra /= Ket. The
!> operator symmetry is Mul(BraSym,KetSym).
!> The products are overlaps M(ij,kl) = (E_ji Bra).(E_kl Ket), summed over
!> the CSFs of the graph; for each symmetry block of pairs the vectors are
!> computed once (in chunks if they do not fit) and all overlaps come from
!> one DGEMM. The sum over the graph is exact for <E_ij E_kl> if k <= l or
!> i >= j (CAS: always; RAS: levels RAS1 < RAS2 < RAS3). Otherwise
!> e_ijkl = e_klij = <E_kl E_ij> - delta_il <E_kj>, which is exact then.
!> Not for GAS graphs (the products are not exact there).
subroutine sg_trans_two_pdm(SGS,CIS,EXS,nBra,Bra,BraSym,nKet,Ket,KetSym,nL,Rho1,n2,Rho2)

  use Index_functions, only: iTri
  use Symmetry_Info, only: Mul
  use sguga, only: CIStruct, EXStruct, sg_epq_psi, SG_New_Tag, SGStruct
  use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
  use Constants, only: Zero, One
  use Definitions, only: wp, iwp, u6

  implicit none
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(EXStruct), intent(inout) :: EXS
  integer(kind=iwp), intent(in) :: nBra, BraSym, nKet, KetSym, nL, n2
  real(kind=wp), intent(in) :: Bra(nBra), Ket(nKet)
  real(kind=wp), intent(out) :: Rho1(nL,nL), Rho2(n2)
  integer(kind=iwp) :: nL2, OpSym, iS, iT, KSym, nK, nX, nY, nC, MaxMem, iX0, iY0, cX, cY, ip, iq, i, j, k, l, ij, kl, tB, tK, &
                       MidSym
  integer(kind=iwp), allocatable :: PX(:), PY(:)
  real(kind=wp), allocatable :: M(:,:), X(:,:), Y(:,:), G(:,:), Sgm(:)

  nL2 = nL*nL
  if (n2 /= nL2*(nL2+1)/2) then
    write(u6,*) 'sg_trans_two_pdm: wrong size of Rho2',n2,nL2*(nL2+1)/2
    call Abend()
  end if
  if ((nBra < CIS%nCSF(BraSym)) .or. (nKet < CIS%nCSF(KetSym))) then
    write(u6,*) 'sg_trans_two_pdm: vectors shorter than the CSF spaces',nBra,CIS%nCSF(BraSym),nKet,CIS%nCSF(KetSym)
    call Abend()
  end if
  tB = SG_New_Tag()
  tK = SG_New_Tag()
  OpSym = Mul(BraSym,KetSym)
  Rho1(:,:) = Zero
  Rho2(:) = Zero

  ! Rho1(i,l) = Bra.(E_il Ket), E_il Ket in symmetry BraSym
  call mma_allocate(Sgm,max(1,CIS%nCSF(BraSym)),Label='TT_Sgm')
  do l=1,nL
    do i=1,nL
      if (Mul(SGS%ISm(i),SGS%ISm(l)) /= OpSym) cycle
      Sgm(:) = Zero
      call SG_Epq_Psi(SGS,CIS,EXS,i,l,One,KetSym,Ket,Sgm,CI_Tag=tK)
      Rho1(i,l) = dot_product(Bra(1:CIS%nCSF(BraSym)),Sgm(1:CIS%nCSF(BraSym)))
    end do
  end do
  call mma_deallocate(Sgm)

  call mma_allocate(M,nL2,nL2,Label='TT_M')
  M(:,:) = Zero
  call mma_allocate(PX,nL2,Label='TT_PX')
  call mma_allocate(PY,nL2,Label='TT_PY')
  do iS=1,SGS%nSym
    ! X(:,ij) = E_ji Bra (pairs ij of symmetry iS), Y(:,kl) = E_kl Ket
    ! (symmetry iT); both in the symmetry MidSym of the intermediate CSFs
    iT = Mul(iS,OpSym)
    MidSym = Mul(BraSym,iS)
    KSym = MidSym
    nK = CIS%nCSF(KSym)
    if (nK == 0) cycle
    nX = 0
    nY = 0
    do j=1,nL
      do i=1,nL
        if (Mul(SGS%ISm(i),SGS%ISm(j)) == iS) then
          nX = nX+1
          PX(nX) = (j-1)*nL+i
        end if
        if (Mul(SGS%ISm(i),SGS%ISm(j)) == iT) then
          nY = nY+1
          PY(nY) = (j-1)*nL+i
        end if
      end do
    end do
    if ((nX == 0) .or. (nY == 0)) cycle
    ! chunks: X and Y of nC columns each in at most half of the memory
    call mma_maxDBLE(MaxMem)
    nC = int(0.5_wp*real(MaxMem,kind=wp)/real(2*nK+nL2,kind=wp),kind=iwp)
    nC = max(1,min(nC,max(nX,nY)))
    call mma_allocate(X,nK,min(nC,nX),Label='TT_X')
    call mma_allocate(Y,nK,min(nC,nY),Label='TT_Y')
    call mma_allocate(G,min(nC,nX),min(nC,nY),Label='TT_G')
    do iX0=0,nX-1,nC
      cX = min(nC,nX-iX0)
      do ip=1,cX
        ij = PX(iX0+ip)
        i = mod(ij-1,nL)+1
        j = (ij-1)/nL+1
        X(:,ip) = Zero
        call SG_Epq_Psi(SGS,CIS,EXS,j,i,One,BraSym,Bra,X(:,ip),CI_Tag=tB)
      end do
      do iY0=0,nY-1,nC
        cY = min(nC,nY-iY0)
        do iq=1,cY
          kl = PY(iY0+iq)
          k = mod(kl-1,nL)+1
          l = (kl-1)/nL+1
          Y(:,iq) = Zero
          call SG_Epq_Psi(SGS,CIS,EXS,k,l,One,KetSym,Ket,Y(:,iq),CI_Tag=tK)
        end do
        call DGEMM_('T','N',cX,cY,nK,One,X,nK,Y,nK,Zero,G,size(G,1))
        do iq=1,cY
          do ip=1,cX
            M(PX(iX0+ip),PY(iY0+iq)) = G(ip,iq)
          end do
        end do
      end do
    end do
    call mma_deallocate(G)
    call mma_deallocate(Y)
    call mma_deallocate(X)
  end do
  call mma_deallocate(PY)
  call mma_deallocate(PX)

  ! e_ijkl from the product order that is exact in the graph
  do ij=1,nL2
    i = mod(ij-1,nL)+1
    j = (ij-1)/nL+1
    do kl=1,ij
      k = mod(kl-1,nL)+1
      l = (kl-1)/nL+1
      if ((k <= l) .or. (i >= j)) then
        Rho2(iTri(ij,kl)) = M(ij,kl)
        if (j == k) Rho2(iTri(ij,kl)) = Rho2(iTri(ij,kl))-Rho1(i,l)
      else
        Rho2(iTri(ij,kl)) = M(kl,ij)
        if (i == l) Rho2(iTri(ij,kl)) = Rho2(iTri(ij,kl))-Rho1(k,j)
      end if
    end do
  end do
  call mma_deallocate(M)

end subroutine sg_trans_two_pdm

!> sg_trans_two_pdm against the determinant expansion of the CSFs
!> (CSF_Det_Expand), for two arbitrary vectors of every pair of symmetries
!> (bra /= ket, also of different symmetry): Rho1 and Rho2 from
!> sum_{st} a+_is a+_kt a_lt a_js in the canonical operator order of
!> CSF_Det_Expand. Only for at most 10 levels and small CSF spaces
!> (dense determinant tables).
subroutine SG_Check_Trans_Two_PDM(SGS,CIS,EXS,nFail)

  use Index_functions, only: iTri
  use sguga, only: CIStruct, CSF_Det_Expand, EXStruct, SG_CSF_StepVec, SGStruct
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One
  use Definitions, only: wp, iwp, u6

  implicit none
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(EXStruct), intent(inout) :: EXS
  integer(kind=iwp), intent(out) :: nFail
  integer(kind=iwp), parameter :: MxLev = 10, MxCSF = 2000
  real(kind=wp), parameter :: Tol = 1.0e-11_wp
  integer(kind=iwp) :: nL, nL2, n2, nS, BraSym, KetSym, nB, nK, nTest, i, j, k, l, ij, kl, ia, ib, a, b, sj, sl, s1, &
                       nSgn, nDim, a1, b1, n1, a2, b2, n2s, a3, b3, n3
  integer(kind=iwp), allocatable :: Step(:), DA(:), DB(:)
  real(kind=wp), allocatable :: Bra(:), Ket(:), Rho1(:,:), Rho2(:), R1(:,:), R2(:,:), CB(:,:), CK(:,:), DC(:)
  real(kind=wp) :: d1, d2, dMax1, dMax2

  nFail = 0
  nL = SGS%nLev
  nS = SGS%nSym
  if ((nL == 0) .or. (nL > MxLev) .or. (maxval(CIS%nCSF(1:nS)) > MxCSF)) then
    write(u6,'(A,I0,A,I0,A)') ' sg_trans_two_pdm check skipped (',nL,' levels, ',maxval(CIS%nCSF(1:nS)),' CSFs)'
    return
  end if
  nL2 = nL*nL
  n2 = nL2*(nL2+1)/2
  nDim = 2**nL
  call mma_allocate(Step,nL,Label='TC_Step')
  call mma_allocate(Rho1,nL,nL,Label='TC_Rho1')
  call mma_allocate(Rho2,n2,Label='TC_Rho2')
  call mma_allocate(R1,nL,nL,Label='TC_R1')
  call mma_allocate(R2,nL2,nL2,Label='TC_R2')
  call mma_allocate(CB,[0,nDim-1],[0,nDim-1],Label='TC_CB')
  call mma_allocate(CK,[0,nDim-1],[0,nDim-1],Label='TC_CK')
  nTest = 0
  dMax1 = Zero
  dMax2 = Zero
  do BraSym=1,nS
    nB = CIS%nCSF(BraSym)
    if (nB == 0) cycle
    do KetSym=1,nS
      nK = CIS%nCSF(KetSym)
      if (nK == 0) cycle
      ! two arbitrary vectors
      call mma_allocate(Bra,nB,Label='TC_Bra')
      call mma_allocate(Ket,nK,Label='TC_Ket')
      do i=1,nB
        Bra(i) = sin(0.731_wp*real(i,kind=wp)+0.37_wp*real(BraSym,kind=wp))
      end do
      do i=1,nK
        Ket(i) = cos(1.217_wp*real(i,kind=wp)-0.53_wp*real(KetSym,kind=wp))
      end do
      call sg_trans_two_pdm(SGS,CIS,EXS,nB,Bra,BraSym,nK,Ket,KetSym,nL,Rho1,n2,Rho2)
      ! determinant coefficients of both vectors
      call Det_Coef(BraSym,nB,Bra,CB)
      call Det_Coef(KetSym,nK,Ket,CK)
      ! reference: R1(i,j) = sum_s <B|a+_is a_js|K>,
      ! R2(ij,kl) = sum_st <B|a+_is a+_kt a_lt a_js|K>
      R1(:,:) = Zero
      R2(:,:) = Zero
      do ib=0,nDim-1
        do ia=0,nDim-1
          if (CK(ia,ib) == Zero) cycle
          do j=1,nL
            do sj=1,2
              ! a_js
              a = ia
              b = ib
              nSgn = 0
              if (.not. Apply(j,sj,.false.)) cycle
              a1 = a
              b1 = b
              n1 = nSgn
              do i=1,nL
                a = a1
                b = b1
                nSgn = n1
                if (.not. Apply(i,sj,.true.)) cycle
                R1(i,j) = R1(i,j)+Sgn_Of(nSgn)*CK(ia,ib)*CB(a,b)
              end do
              do l=1,nL
                do sl=1,2
                  ! a_lt a_js
                  a = a1
                  b = b1
                  nSgn = n1
                  if (.not. Apply(l,sl,.false.)) cycle
                  a2 = a
                  b2 = b
                  n2s = nSgn
                  do k=1,nL
                    ! a+_kt a_lt a_js
                    a = a2
                    b = b2
                    nSgn = n2s
                    if (.not. Apply(k,sl,.true.)) cycle
                    a3 = a
                    b3 = b
                    n3 = nSgn
                    kl = (l-1)*nL+k
                    do i=1,nL
                      a = a3
                      b = b3
                      nSgn = n3
                      if (.not. Apply(i,sj,.true.)) cycle
                      ij = (j-1)*nL+i
                      R2(ij,kl) = R2(ij,kl)+Sgn_Of(nSgn)*CK(ia,ib)*CB(a,b)
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
      do ij=1,nL2
        do kl=1,ij
          d2 = abs(Rho2(iTri(ij,kl))-R2(ij,kl))
          dMax2 = max(dMax2,d2)
          if (d2 > Tol) then
            nFail = nFail+1
            if (nFail <= 20) write(u6,'(A,4I4,2ES20.12)') ' sg_trans_two_pdm mismatch: BraSym,KetSym,ij,kl, value, reference =', &
              BraSym,KetSym,ij,kl,Rho2(iTri(ij,kl)),R2(ij,kl)
          end if
        end do
      end do
      d1 = maxval(abs(Rho1(:,:)-R1(:,:)))
      dMax1 = max(dMax1,d1)
      if (d1 > Tol) nFail = nFail+1
      nTest = nTest+1
      call mma_deallocate(Bra)
      call mma_deallocate(Ket)
    end do
  end do
  write(u6,'(/,A,I0,A)') ' sg_trans_two_pdm (2-body transition densities) vs determinants: ',nTest,' pairs of symmetries'
  write(u6,'(A,2ES10.3)') '   max. abs. difference Rho1, Rho2 ',dMax1,dMax2
  write(u6,'(A,I0)') '   failures: ',nFail
  call mma_deallocate(Step)
  call mma_deallocate(Rho1)
  call mma_deallocate(Rho2)
  call mma_deallocate(R1)
  call mma_deallocate(R2)
  call mma_deallocate(CB)
  call mma_deallocate(CK)

contains

  ! the coefficients of the determinants (alpha string, beta string) of
  ! the vector V of symmetry iSym
  subroutine Det_Coef(iSym,n,V,C)
    integer(kind=iwp), intent(in) :: iSym, n
    real(kind=wp), intent(in) :: V(n)
    real(kind=wp), intent(out) :: C(0:nDim-1,0:nDim-1)
    integer(kind=iwp) :: I, iD, nD
    C(:,:) = Zero
    do I=1,n
      call SG_CSF_StepVec(SGS,CIS,iSym,I,Step)
      call CSF_Det_Expand(nL,Step,nD,DA,DB,DC)
      do iD=1,nD
        C(DA(iD),DB(iD)) = C(DA(iD),DB(iD))+V(I)*DC(iD)
      end do
      call mma_deallocate(DA)
      call mma_deallocate(DB)
      call mma_deallocate(DC)
    end do
  end subroutine Det_Coef

  ! a+ (Create) or a of level p, spin s (1 alpha, 2 beta) on (a,b); the
  ! operators left of it are those of the levels above p and, for alpha,
  ! the beta of level p (canonical order of CSF_Det_Expand)
  function Apply(p,s,Create) result(Ok)
    integer(kind=iwp), intent(in) :: p, s
    logical(kind=iwp), intent(in) :: Create
    logical(kind=iwp) :: Ok
    logical(kind=iwp) :: Occ
    if (s == 1) then
      Occ = btest(a,p-1)
    else
      Occ = btest(b,p-1)
    end if
    Ok = (Occ .neqv. Create)
    if (.not. Ok) return
    s1 = popcnt(ishft(a,-p))+popcnt(ishft(b,-p))
    if ((s == 1) .and. btest(b,p-1)) s1 = s1+1
    nSgn = nSgn+s1
    if (s == 1) then
      if (Create) then
        a = ibset(a,p-1)
      else
        a = ibclr(a,p-1)
      end if
    else
      if (Create) then
        b = ibset(b,p-1)
      else
        b = ibclr(b,p-1)
      end if
    end if
  end function Apply

  function Sgn_Of(n) result(f)
    integer(kind=iwp), intent(in) :: n
    real(kind=wp) :: f
    f = One
    if (mod(n,2) == 1) f = -One
  end function Sgn_Of

end subroutine SG_Check_Trans_Two_PDM

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
  use sguga, only: CIStruct, EXStruct, sg_epq_psi, SG_New_Tag, SGStruct
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
  real(kind=wp), allocatable :: AA(:,:), AB(:,:), BA(:,:), BB(:,:)
  real(kind=wp), pointer :: A(:,:), B(:,:)

  integer(kind=iwp) :: tCI
  tCI = SG_New_Tag()
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
  MaxMem = MaxMem+Wrk_Held()
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
    call Get_Wrk(Wrk_D1,mCSFs*nPr,'SG_Wrk_D1')
    call Get_Wrk(Wrk_D2,mCSFs*nPr,'SG_Wrk_D2')
    A(1:mCSFs,1:nPr) => Wrk_D1(1:mCSFs*nPr)
    B(1:mCSFs,1:nPr) => Wrk_D2(1:mCSFs*nPr)
    A(:,:) = Zero
    B(:,:) = Zero
    do ip=1,nPr
      call SG_Epq_Psi(SGS,CIS,EXS,PI(ip),PJ(ip),One,PsiSym,Psi,A(:,ip),CI_Tag=tCI)
      if (PI(ip) /= PJ(ip)) call SG_Epq_Psi(SGS,CIS,EXS,PJ(ip),PI(ip),One,PsiSym,Psi,B(:,ip),CI_Tag=tCI)
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
    nullify(B)
    nullify(A)
  end do
  call mma_deallocate(Loc)
  call mma_deallocate(PJ)
  call mma_deallocate(PI)
  P(:) = Half*P(:)
  PA(:) = Half*PA(:)
  call Wrk_Trim()
end subroutine sg_two_pdm

subroutine sg_two_pdm_direct(SGS,CIS,EXS,Psi,nCSFs,PsiSym,P,PA,nP)

  use Index_functions, only: iTri
  use Symmetry_Info, only: Mul
  use sguga, only: CIStruct, EXStruct, sg_epq_psi, SG_New_Tag, SGStruct
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

  integer(kind=iwp) :: tCI
  tCI = SG_New_Tag()
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
      call SG_Epq_Psi(SGS,CIS,EXS,iOrb,jOrb,CPQ,PsiSym,Psi,Eij_Psi,CI_Tag=tCI)
      ! Compute Dij to be distributed below
      if (iSym == jSym) D_ij = dot_product(Psi,Eij_Psi)
      if (iOrb /= jOrb) call SG_Epq_Psi(SGS,CIS,EXS,jOrb,iOrb,CPQ,PsiSym,Psi,Eji_Psi,CI_Tag=tCI)

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
          call SG_Epq_Psi(SGS,CIS,EXS,lOrb,kOrb,CPQ,PsiSym,Psi,Elk_Psi,CI_Tag=tCI)
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

!> Free the work arrays kept between the calls.
subroutine SG_Densities_Free()
  use stdalloc, only: mma_deallocate
  call mma_deallocate(Wrk_D1,safe='*')
  call mma_deallocate(Wrk_D2,safe='*')
end subroutine SG_Densities_Free

!> Keep the work arrays for the next call only if they hold at most a quarter
!> of the memory (held + free); otherwise free them.
subroutine Wrk_Trim()
  use stdalloc, only: mma_maxDBLE
  integer(kind=iwp) :: nFree, nHeld
  call mma_maxDBLE(nFree)
  nHeld = Wrk_Held()
  if (4*nHeld > nHeld+nFree) call SG_Densities_Free()
end subroutine Wrk_Trim

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

!> The memory held by the work arrays (in words).
function Wrk_Held() result(n)
  integer(kind=iwp) :: n
  n = 0
  if (allocated(Wrk_D1)) n = n+size(Wrk_D1)
  if (allocated(Wrk_D2)) n = n+size(Wrk_D2)
end function Wrk_Held

end module SGUGA_Densities
