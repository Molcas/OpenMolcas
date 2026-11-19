
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
    SDMAT(iTri(p,q)) = Fact*((Two-Half*real(SGS%nActEl,kind=wp))*D1MAT(iTri(p,q))-(G(iTri(p,q))-real(nL,kind=wp)*D1MAT(iTri(p,q))))
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
