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
!               2026, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
! The ASTA triplet operator T_pq = a+_pa a_qa - a+_pb a_qb: matrix
! elements, split-graph coupling tables, T_pq|Psi>, determinant references.
! (Default public, as SGUGA_Core; the interface is the module SGUGA.)
module SGUGA_Triplet

use SGUGA_Core

implicit none

contains

!> <I|E_pq|J> (Triplet = .false.) or <I|a+_pa a_qa - a+_pb a_qb|J>
!> (Triplet = .true.) of the CSFs with step vectors SI (bra) and SJ (ket)
!> of the same graph, as the product of ASTA segment values.
pure function ASTA_MatEl(nLev,SI,SJ,p,q,Triplet) result(V)
  integer(kind=iwp), intent(in) :: nLev, SI(nLev), SJ(nLev), p, q
  logical(kind=iwp), intent(in) :: Triplet
  real(kind=wp) :: V
  integer(kind=iwp) :: U(nLev), T(nLev), BU(0:nLev), BT(0:nLev), i, j, k
  if (p <= q) then
    U(:) = SI(:)
    T(:) = SJ(:)
    i = p
    j = q
  else
    U(:) = SJ(:)
    T(:) = SI(:)
    i = q
    j = p
  end if
  V = Zero
  do k=1,i-1
    if (U(k) /= T(k)) return
  end do
  call Walk_B(nLev,U,BU)
  call Walk_B(nLev,T,BT)
  if (i == j) then
    ! weight segment
    if (Triplet) then
      V = ASTA_Weight_T(U(i),T(i),BT(i)-BU(i),BT(i))
    else if (U(i) == T(i)) then
      V = real((1+U(i))/2,kind=wp)
    end if
  else
    V = ASTA_Tail(U(i),T(i),BT(i))
    do k=i+1,j-1
      if (V == Zero) return
      V = V*ASTA_Mid(U(k),T(k),BT(k)-BU(k),BT(k))
    end do
    if (V == Zero) return
    if (Triplet) then
      V = V*ASTA_Head_T(U(j),T(j),BT(j)-BU(j),BT(j))
    else if (BT(j) == BU(j)) then
      V = V*ASTA_Head(U(j),T(j),BT(j))
    else
      V = Zero
    end if
  end if
  ! above the loop head
  do k=j+1,nLev
    if (V == Zero) return
    if (Triplet) then
      V = V*ASTA_Top_T(U(k),T(k),BT(k)-BU(k),BT(k))
    else if (U(k) /= T(k)) then
      V = Zero
    end if
  end do
  if (BT(nLev) /= BU(nLev)) V = Zero
end function ASTA_MatEl

!> <I|E_pq|J> or <I|T_pq|J> (Triplet) for CSF I of symmetry ISYI and J of
!> symmetry ISYJ of state iState (ASTA walk product).
function SG_MatEl_Tpq(iState,ISYI,I,ISYJ,J,p,q,Triplet) result(V)
  integer(kind=iwp), intent(in) :: iState, ISYI, I, ISYJ, J, p, q
  logical(kind=iwp), intent(in) :: Triplet
  real(kind=wp) :: V
  integer(kind=iwp), allocatable :: SI(:), SJ(:)
  V = Zero
  if (Mul(Mul(SGS(iState)%ISm(p),SGS(iState)%ISm(q)),ISYJ) /= ISYI) return
  call mma_allocate(SI,SGS(iState)%nLev,Label='SI')
  call mma_allocate(SJ,SGS(iState)%nLev,Label='SJ')
  call SG_CSF_StepVec(SGS(iState),CIS(iState),ISYI,I,SI)
  call SG_CSF_StepVec(SGS(iState),CIS(iState),ISYJ,J,SJ)
  V = ASTA_MatEl(SGS(iState)%nLev,SI,SJ,p,q,Triplet)
  call mma_deallocate(SI)
  call mma_deallocate(SJ)
end function SG_MatEl_Tpq

!> Determinant expansion (M_S = S) of the genealogical CSF with step vector
!> Step, from the ASTA definition (paper I): level by level from the top,
!> u: sum_s C(+1/2,s) a+_{ks}, d: sum_s C(-1/2,s) a+_{ks}, 2: (-1)^b
!> a+_{kb} a+_{ka}. The creation operators stay in the canonical order
!> (levels descending from the left, beta before alpha in a level), so the
!> terms carry no reordering signs. DA/DB: alpha/beta occupations (bit k-1
!> = level k), DC: coefficients.
subroutine CSF_Det_Expand(nLev,Step,nDet,DA,DB,DC)
  integer(kind=iwp), intent(in) :: nLev, Step(nLev)
  integer(kind=iwp), intent(out) :: nDet
  integer(kind=iwp), allocatable, intent(out) :: DA(:), DB(:)
  real(kind=wp), allocatable, intent(out) :: DC(:)
  integer(kind=iwp) :: B(0:nLev), k, n, m, iT, nMax, s2, m2n
  integer(kind=iwp), allocatable :: M2(:), NA(:), NB(:), NM(:)
  real(kind=wp), allocatable :: NC(:)
  real(kind=wp) :: c
  call Walk_B(nLev,Step,B)
  nMax = 2**count((Step(:) == 1) .or. (Step(:) == 2))
  call mma_allocate(DA,nMax,Label='DA')
  call mma_allocate(DB,nMax,Label='DB')
  call mma_allocate(DC,nMax,Label='DC')
  call mma_allocate(M2,nMax,Label='M2')
  call mma_allocate(NA,nMax,Label='NA')
  call mma_allocate(NB,nMax,Label='NB')
  call mma_allocate(NC,nMax,Label='NC')
  call mma_allocate(NM,nMax,Label='NM')
  n = 1
  DA(1) = 0
  DB(1) = 0
  DC(1) = One
  M2(1) = B(nLev)
  do k=nLev,1,-1
    m = 0
    do iT=1,n
      select case (Step(k))
        case (0)
          m = m+1
          NA(m) = DA(iT)
          NB(m) = DB(iT)
          NC(m) = DC(iT)
          NM(m) = M2(iT)
        case (3)
          m = m+1
          NA(m) = ibset(DA(iT),k-1)
          NB(m) = ibset(DB(iT),k-1)
          NC(m) = Sgn(B(k))*DC(iT)
          NM(m) = M2(iT)
        case (1,2)
          do s2=1,-1,-2
            ! 2*sigma = s2, new projection M - sigma
            m2n = M2(iT)-s2
            if (Step(k) == 1) then
              if (abs(m2n) > B(k)-1) cycle
              c = sqrt(real(B(k)+s2*M2(iT),kind=wp)/real(2*B(k),kind=wp))
            else
              if (abs(m2n) > B(k)+1) cycle
              c = Sgn(B(k)+1)*real(s2,kind=wp)*sqrt(real(B(k)+2-s2*M2(iT),kind=wp)/real(2*(B(k)+2),kind=wp))
            end if
            if (c == Zero) cycle
            m = m+1
            if (s2 == 1) then
              NA(m) = ibset(DA(iT),k-1)
              NB(m) = DB(iT)
            else
              NA(m) = DA(iT)
              NB(m) = ibset(DB(iT),k-1)
            end if
            NC(m) = c*DC(iT)
            NM(m) = m2n
          end do
      end select
    end do
    n = m
    DA(1:n) = NA(1:n)
    DB(1:n) = NB(1:n)
    DC(1:n) = NC(1:n)
    M2(1:n) = NM(1:n)
  end do
  nDet = n
  call mma_deallocate(M2)
  call mma_deallocate(NA)
  call mma_deallocate(NB)
  call mma_deallocate(NC)
  call mma_deallocate(NM)
end subroutine CSF_Det_Expand

!> <I|a+_pa a_qa + sb a+_pb a_qb|J> from the determinant expansions of I
!> and J (CSF_Det_Expand, canonical operator order); sb = 1 for E_pq, -1
!> for T_pq.
pure function Det_MatEl(nI,AI,BI,CI,nJ,AJ,BJ,CJ,p,q,sb) result(V)
  integer(kind=iwp), intent(in) :: nI, AI(nI), BI(nI), nJ, AJ(nJ), BJ(nJ), p, q
  real(kind=wp), intent(in) :: CI(nI), CJ(nJ), sb
  real(kind=wp) :: V
  integer(kind=iwp) :: iJ, iI, iSpin, a, bb, nLeft
  real(kind=wp) :: w
  V = Zero
  do iJ=1,nJ
    do iSpin=1,2
      a = AJ(iJ)
      bb = BJ(iJ)
      if (iSpin == 1) then
        ! alpha: annihilate a_qa, the operators left of it are all of the
        ! levels above q and the beta of level q
        if (.not. btest(a,q-1)) cycle
        nLeft = popcnt(ishft(a,-q))+popcnt(ishft(bb,-q))
        if (btest(bb,q-1)) nLeft = nLeft+1
        a = ibclr(a,q-1)
        if (btest(a,p-1)) cycle
        nLeft = nLeft+popcnt(ishft(a,-p))+popcnt(ishft(bb,-p))
        if (btest(bb,p-1)) nLeft = nLeft+1
        a = ibset(a,p-1)
        w = One
      else
        ! beta: the operators left of it are those of the levels above
        if (.not. btest(bb,q-1)) cycle
        nLeft = popcnt(ishft(a,-q))+popcnt(ishft(bb,-q))
        bb = ibclr(bb,q-1)
        if (btest(bb,p-1)) cycle
        nLeft = nLeft+popcnt(ishft(a,-p))+popcnt(ishft(bb,-p))
        bb = ibset(bb,p-1)
        w = sb
      end if
      w = w*CJ(iJ)
      if (mod(nLeft,2) == 1) w = -w
      do iI=1,nI
        if ((AI(iI) == a) .and. (BI(iI) == bb)) V = V+CI(iI)*w
      end do
    end do
  end do
end function Det_MatEl

!> <I|T_pq|J> from the tabulated segments (catalog ITVPT/IC1/IC2, ISGM,
!> VSGM, IVR of MKSEG), walking both CSFs down from the graph head: the
!> left (bra of the raising operator) and right (ket) vertices at the top
!> of each level give the class of the segment, which with the two steps
!> identifies its type; the segment must connect the actual vertices. For
!> the verification of the triplet catalog against ASTA_MatEl.
function Table_MatEl_T(SGS,CIS,SI,SJ,p,q) result(V)
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  integer(kind=iwp), intent(in) :: SI(SGS%nLev), SJ(SGS%nLev), p, q
  real(kind=wp) :: V
  integer(kind=iwp) :: U(SGS%nLev), T(SGS%nLev), i, j, k, vL, vR, db, iCls, iBeg, iEnd, ISGT, kS, vRB
  if (p <= q) then
    U(:) = SI(:)
    T(:) = SJ(:)
    i = p
    j = q
  else
    U(:) = SJ(:)
    T(:) = SI(:)
    i = q
    j = p
  end if
  V = One
  vL = 1
  vR = 1
  do k=SGS%nLev,1,-1
    db = SGS%DRT(vR,CSFG_BTab)-SGS%DRT(vL,CSFG_BTab)
    if (k > j) then
      ! top segments (triplet classes 4-6)
      iBeg = iSegTTopBeg
      iEnd = iSegTTopEnd
      iCls = 4+merge(1,0,db == -2)+merge(2,0,db == 2)
      if (abs(db) /= 2 .and. (db /= 0)) iCls = -1
    else if (k == j) then
      if (i == j) then
        iBeg = iSegTWBeg
        iEnd = iSegTWEnd
      else
        iBeg = iSegTHeadBeg
        iEnd = iSegTHeadEnd
      end if
      iCls = 4+merge(1,0,db == -2)+merge(2,0,db == 2)
      if (abs(db) /= 2 .and. (db /= 0)) iCls = -1
    else if (k > i) then
      ! intermediate segments (singlet)
      iBeg = 9
      iEnd = 18
      iCls = merge(1,2,db == -1)
      if (abs(db) /= 1) iCls = -1
    else if (k == i) then
      ! tail segments (singlet)
      iBeg = 19
      iEnd = 22
      iCls = merge(1,2,db == -1)
      if (abs(db) /= 1) iCls = -1
    else
      ! tail walk
      iBeg = 23
      iEnd = 26
      iCls = 3
      if (db /= 0) iCls = -1
    end if
    if (iCls < 0) then
      V = Zero
      return
    end if
    ISGT = 0
    do kS=iBeg,iEnd
      if ((ITVPT(kS) == iCls) .and. (IC1(kS) == U(k)) .and. (IC2(kS) == T(k))) ISGT = kS
    end do
    if (ISGT == 0) then
      V = Zero
      return
    end if
    if (CIS%ISGM(vL,ISGT) == 0) then
      V = Zero
      return
    end if
    ! the segment must start at the right vertex vR and end at the walks
    select case (iCls)
      case (1,2)
        if (CIS%IVR(vL,iCls) /= vR) call Abend()
      case (5)
        if (CIS%IVR(vL,3) /= vR) call Abend()
      case (6)
        if (CIS%IVR(vL,4) /= vR) call Abend()
      case default
        if (vL /= vR) call Abend()
    end select
    vRB = SGS%Down(vR,T(k))
    if (CIS%ISGM(vL,ISGT) /= SGS%Down(vL,U(k))) call Abend()
    V = V*CIS%VSGM(vL,ISGT)
    vL = SGS%Down(vL,U(k))
    vR = vRB
    if ((vL == 0) .or. (vR == 0)) call Abend()
  end do
end function Table_MatEl_T

!> The step vectors of all CSFs of state iState (all symmetries, CSF k of
!> symmetry ISY at column ASTA_Off(ISY)+k) and, per symmetry, the packed
!> step vectors sorted (ASTA_Key, with the CSF index ASTA_Idx) for the
!> lookup of a walk; built on first use, freed by ASTA_Free (SG_Free).
subroutine ASTA_Setup(iState)
  integer(kind=iwp), intent(in) :: iState
  integer(kind=iwp) :: ISY, I, k, n, nTot
  integer(kind=iwp), allocatable :: Ord(:)
  if (allocated(ASTA_Steps) .and. (ASTA_iState == iState)) return
  call ASTA_Free()
  ASTA_iState = iState
  ASTA_Off(:) = 0
  do ISY=2,SGS(iState)%nSym
    ASTA_Off(ISY) = ASTA_Off(ISY-1)+CIS(iState)%NCSF(ISY-1)
  end do
  nTot = ASTA_Off(SGS(iState)%nSym)+CIS(iState)%NCSF(SGS(iState)%nSym)
  call mma_allocate(ASTA_Steps,SGS(iState)%nLev,max(1,nTot),Label='ASTA_Steps')
  call mma_allocate(ASTA_Key,max(1,nTot),Label='ASTA_Key')
  call mma_allocate(ASTA_Idx,max(1,nTot),Label='ASTA_Idx')
  do ISY=1,SGS(iState)%nSym
    n = CIS(iState)%NCSF(ISY)
    if (n == 0) cycle
    do I=1,n
      call SG_CSF_StepVec(SGS(iState),CIS(iState),ISY,I,ASTA_Steps(:,ASTA_Off(ISY)+I))
      ASTA_Key(ASTA_Off(ISY)+I) = Walk_Key(SGS(iState)%nLev,ASTA_Steps(:,ASTA_Off(ISY)+I))
    end do
    call mma_allocate(Ord,n,Label='Ord')
    call Sort_Keys(ASTA_Key(ASTA_Off(ISY)+1:ASTA_Off(ISY)+n),1,n,Ord)
    ASTA_Idx(ASTA_Off(ISY)+1:ASTA_Off(ISY)+n) = Ord(:)
    ASTA_Key(ASTA_Off(ISY)+1:ASTA_Off(ISY)+n) = ASTA_Key(ASTA_Off(ISY)+Ord(:))
    do k=2,n
      if (ASTA_Key(ASTA_Off(ISY)+k) == ASTA_Key(ASTA_Off(ISY)+k-1)) then
        write(u6,*) 'ASTA_Setup: two CSFs with the same step vector',ISY,k
        call Abend()
      end if
    end do
    call mma_deallocate(Ord)
  end do
end subroutine ASTA_Setup

!> A step vector packed into one integer, 2 bits per level (level 1 lowest).
pure function Walk_Key(nLev,Step) result(Key)
  integer(kind=iwp), intent(in) :: nLev, Step(nLev)
  integer(kind=iwp) :: Key
  integer(kind=iwp) :: k
  Key = 0
  do k=nLev,1,-1
    Key = 4*Key+Step(k)
  end do
end function Walk_Key

!> The CSF of symmetry ISY with the packed step vector Key (0 if none).
function ASTA_Find(ISY,n,Key) result(I)
  integer(kind=iwp), intent(in) :: ISY, n, Key
  integer(kind=iwp) :: I
  integer(kind=iwp) :: iLo, iHi, iMid
  I = 0
  iLo = 1
  iHi = n
  do while (iLo <= iHi)
    iMid = (iLo+iHi)/2
    if (ASTA_Key(ASTA_Off(ISY)+iMid) == Key) then
      I = ASTA_Idx(ASTA_Off(ISY)+iMid)
      return
    else if (ASTA_Key(ASTA_Off(ISY)+iMid) < Key) then
      iLo = iMid+1
    else
      iHi = iMid-1
    end if
  end do
end function ASTA_Find

!> Sigma = Sigma + X T_pq CI (Triplet: T_pq = a+_pa a_qa - a+_pb a_qb), or
!> + X E_pq CI (Triplet = .false.; the ASTA reference of SG_Epq_Psi), for a
!> CI vector of symmetry ISYCI of state iState; Sigma has the symmetry
!> Gamma(p) x Gamma(q) x ISYCI. For each ket CSF J the partner walks with a
!> nonzero product of ASTA segment values are generated level by level from
!> the graph tail (ASTA_Gen) and located by their packed step vectors.
subroutine SG_Xpq_Psi(iState,p,q,X,ISYCI,CI,Sigma,Triplet)
  integer(kind=iwp), intent(in) :: iState, p, q, ISYCI
  real(kind=wp), intent(in) :: X, CI(*)
  real(kind=wp), intent(inout) :: Sigma(*)
  logical(kind=iwp), intent(in) :: Triplet
  integer(kind=iwp) :: ISYS, J, nLev, i_, j_
  integer(kind=iwp), allocatable :: P_(:), BF(:), BP(:)
  logical(kind=iwp) :: FixedIsKet
  call ASTA_Setup(iState)
  nLev = SGS(iState)%nLev
  ISYS = Mul(Mul(SGS(iState)%ISm(p),SGS(iState)%ISm(q)),ISYCI)
  if ((CIS(iState)%NCSF(ISYCI) == 0) .or. (CIS(iState)%NCSF(ISYS) == 0)) return
  ! raising (p < q): the fixed walk J is the ket of the raising operator;
  ! lowering: <I|X_pq|J> = <J|X_qp|I>, J is its bra
  FixedIsKet = (p <= q)
  i_ = min(p,q)
  j_ = max(p,q)
  call mma_allocate(P_,nLev,Label='ASTA_P')
  call mma_allocate(BF,[0,nLev],Label='ASTA_BF')
  call mma_allocate(BP,[0,nLev],Label='ASTA_BP')
  do J=1,CIS(iState)%NCSF(ISYCI)
    if (CI(J) == Zero) cycle
    call Walk_B(nLev,ASTA_Steps(:,ASTA_Off(ISYCI)+J),BF)
    BP(0) = 0
    call ASTA_Gen(SGS(iState),1,SGS(iState)%nVert,One,nLev,ASTA_Steps(:,ASTA_Off(ISYCI)+J),BF,P_,BP,i_,j_,FixedIsKet, &
                  Triplet,ISYS,CIS(iState)%NCSF(ISYS),X*CI(J),Sigma)
  end do
  call mma_deallocate(P_)
  call mma_deallocate(BF)
  call mma_deallocate(BP)
end subroutine SG_Xpq_Psi

!> Sigma = Sigma + X T_pq CI (T_pq = a+_pa a_qa - a+_pb a_qb), see SG_Xpq_Psi:
!> the reference of SG_Tpq_Psi (the split-graph tables).
subroutine SG_Tpq_Psi_Ref(iState,p,q,X,ISYCI,CI,Sigma)
  integer(kind=iwp), intent(in) :: iState, p, q, ISYCI
  real(kind=wp), intent(in) :: X, CI(*)
  real(kind=wp), intent(inout) :: Sigma(*)
  call SG_Xpq_Psi(iState,p,q,X,ISYCI,CI,Sigma,.true.)
end subroutine SG_Tpq_Psi_Ref

!***********************************************************************
! Coupling tables of the triplet operator T_pq in the split graph (TRP_*).
! A matrix element of a raising operator (tail level i <= head level j) is
! the product of an upper and a lower factor, each a pair of half walks
! (left = bra, right = ket) with the product of the catalog segment values
! over the levels of the half:
!   j <= MidLev:          upper: top segments only (category 1, class c =
!                         4, 5, 6 at the mid level, independent of i, j);
!                         lower: top segments, head (or weight, i = j) at j,
!                         intermediate and tail segments (category 4: c,i,j)
!   i <= MidLev < j:      upper: top segments, head at j, intermediate
!                         segments (category 2: class c = 1, 2, j); lower:
!                         intermediate and tail segments (category 5: c, i)
!   MidLev < i:           upper: the whole loop (category 3: i, j); the
!                         lower walks are identical
! The tuples (left walk, right walk, value) are generated by a depth-first
! walk over the segment catalog (TRP_Gen) and sorted by (category, left
! midvertex); half walks are numbered as in CIS (lookup TRP_HW_*). Built on
! first use (TRP_Setup), freed by TRP_Free (SG_Free).
!***********************************************************************

!> The tables of state iState.
subroutine TRP_Setup(iState)
  integer(kind=iwp), intent(in) :: iState
  integer(kind=iwp) :: IV, ICls, IVR_, n
  integer(kind=iwp), allocatable :: Ord(:), K2(:,:)
  if (allocated(TRP_Key) .and. (TRP_iState == iState)) return
  call TRP_Free()
  TRP_iState = iState
  call TRP_HW_Setup(iState)
  TRP_n = 0
  TRP_Cap = 1024
  call TRP_Alloc(TRP_Cap)
  call mma_allocate(TRP_SL,SGS(iState)%nLev,Label='TRP_SL')
  call mma_allocate(TRP_SR,SGS(iState)%nLev,Label='TRP_SR')
  ! upper half: from the graph head (db = 0, class 4)
  TRP_Half = 1
  call TRP_Gen(iState,SGS(iState)%nLev,1,1,4,1,0,0,One)
  ! lower half: from every midvertex pair of the classes 1, 2, 4, 5, 6
  TRP_Half = 2
  do IV=SGS(iState)%MVSta,SGS(iState)%MVEnd
    do ICls=1,6
      if (ICls == 3) cycle
      select case (ICls)
        case (1,2)
          IVR_ = CIS(iState)%IVR(IV,ICls)
        case (5)
          IVR_ = CIS(iState)%IVR(IV,3)
        case (6)
          IVR_ = CIS(iState)%IVR(IV,4)
        case default
          IVR_ = IV
      end select
      if (IVR_ == 0) cycle
      TRP_MVL0 = IV-SGS(iState)%MVSta+1
      TRP_Start = ICls
      call TRP_Gen(iState,SGS(iState)%MidLev,IV,IVR_,ICls,merge(2,1,ICls <= 2),0,0,One)
    end do
  end do
  call mma_deallocate(TRP_SL)
  call mma_deallocate(TRP_SR)
  ! sort by (category, left midvertex)
  n = TRP_n
  call mma_allocate(Ord,max(1,n),Label='TRP_Ord')
  call mma_allocate(K2,1,max(1,n),Label='TRP_K2')
  K2(1,1:n) = TRP_Key(1:n)
  call Sort_Keys(K2,1,n,Ord)
  call TRP_Permute(Ord,n)
  call mma_deallocate(K2)
  call mma_deallocate(Ord)
end subroutine TRP_Setup

subroutine TRP_Alloc(nCap)
  integer(kind=iwp), intent(in) :: nCap
  call mma_allocate(TRP_Key,nCap,Label='TRP_Key')
  call mma_allocate(TRP_IL,nCap,Label='TRP_IL')
  call mma_allocate(TRP_IR,nCap,Label='TRP_IR')
  call mma_allocate(TRP_SYL,nCap,Label='TRP_SYL')
  call mma_allocate(TRP_SYR,nCap,Label='TRP_SYR')
  call mma_allocate(TRP_MVR,nCap,Label='TRP_MVR')
  call mma_allocate(TRP_V,nCap,Label='TRP_V')
end subroutine TRP_Alloc

!> Reorder the tuples by Ord (also used to grow the arrays, Ord = identity).
subroutine TRP_Permute(Ord,n)
  integer(kind=iwp), intent(in) :: n, Ord(n)
  integer(kind=iwp), allocatable :: I1(:)
  real(kind=wp), allocatable :: R1(:)
  call mma_allocate(I1,max(1,n),Label='TRP_I1')
  call mma_allocate(R1,max(1,n),Label='TRP_R1')
  I1(1:n) = TRP_Key(Ord(:))
  TRP_Key(1:n) = I1(1:n)
  I1(1:n) = TRP_IL(Ord(:))
  TRP_IL(1:n) = I1(1:n)
  I1(1:n) = TRP_IR(Ord(:))
  TRP_IR(1:n) = I1(1:n)
  I1(1:n) = TRP_SYL(Ord(:))
  TRP_SYL(1:n) = I1(1:n)
  I1(1:n) = TRP_SYR(Ord(:))
  TRP_SYR(1:n) = I1(1:n)
  I1(1:n) = TRP_MVR(Ord(:))
  TRP_MVR(1:n) = I1(1:n)
  R1(1:n) = TRP_V(Ord(:))
  TRP_V(1:n) = R1(1:n)
  call mma_deallocate(I1)
  call mma_deallocate(R1)
end subroutine TRP_Permute

subroutine TRP_Push(Key,IL,IR,SYL,SYR,MVR,V)
  integer(kind=iwp), intent(in) :: Key, IL, IR, SYL, SYR, MVR
  real(kind=wp), intent(in) :: V
  integer(kind=iwp), allocatable :: K(:), A(:), B(:), C(:), D(:), E(:)
  real(kind=wp), allocatable :: R(:)
  integer(kind=iwp) :: n
  if (TRP_n == TRP_Cap) then
    n = TRP_n
    call mma_allocate(K,n,Label='TRP_t1')
    call mma_allocate(A,n,Label='TRP_t2')
    call mma_allocate(B,n,Label='TRP_t3')
    call mma_allocate(C,n,Label='TRP_t4')
    call mma_allocate(D,n,Label='TRP_t5')
    call mma_allocate(E,n,Label='TRP_t6')
    call mma_allocate(R,n,Label='TRP_t7')
    K(:) = TRP_Key(1:n)
    A(:) = TRP_IL(1:n)
    B(:) = TRP_IR(1:n)
    C(:) = TRP_SYL(1:n)
    D(:) = TRP_SYR(1:n)
    E(:) = TRP_MVR(1:n)
    R(:) = TRP_V(1:n)
    call mma_deallocate(TRP_Key)
    call mma_deallocate(TRP_IL)
    call mma_deallocate(TRP_IR)
    call mma_deallocate(TRP_SYL)
    call mma_deallocate(TRP_SYR)
    call mma_deallocate(TRP_MVR)
    call mma_deallocate(TRP_V)
    TRP_Cap = 2*TRP_Cap
    call TRP_Alloc(TRP_Cap)
    TRP_Key(1:n) = K(:)
    TRP_IL(1:n) = A(:)
    TRP_IR(1:n) = B(:)
    TRP_SYL(1:n) = C(:)
    TRP_SYR(1:n) = D(:)
    TRP_MVR(1:n) = E(:)
    TRP_V(1:n) = R(:)
    call mma_deallocate(K)
    call mma_deallocate(A)
    call mma_deallocate(B)
    call mma_deallocate(C)
    call mma_deallocate(D)
    call mma_deallocate(E)
    call mma_deallocate(R)
  end if
  TRP_n = TRP_n+1
  TRP_Key(TRP_n) = Key
  TRP_IL(TRP_n) = IL
  TRP_IR(TRP_n) = IR
  TRP_SYL(TRP_n) = SYL
  TRP_SYR(TRP_n) = SYR
  TRP_MVR(TRP_n) = MVR
  TRP_V(TRP_n) = V
end subroutine TRP_Push

!> Sort key of a tuple: category (type 1-5, class, i, j) and left midvertex.
pure function TRP_Cat(nLev,nMidV,iType,iCls,i,j,MV) result(Key)
  integer(kind=iwp), intent(in) :: nLev, nMidV, iType, iCls, i, j, MV
  integer(kind=iwp) :: Key
  Key = ((((iType*8+iCls)*(nLev+1)+i)*(nLev+1)+j)*(nMidV+1))+MV
end function TRP_Cat

!> The half walks of state iState: per half, the packed step vectors with
!> their midvertex, sorted (TRP_HW_K(1:2,:) = step key, MV), with symmetry
!> and index of the walk in its (symmetry, MV) block.
subroutine TRP_HW_Setup(iState)
  integer(kind=iwp), intent(in) :: iState
  integer(kind=iwp) :: iHalf, ISY, MV, w, n, k, iLev, lev0, nL, IOW0, IPSHFT, IC, Key
  integer(kind=iwp), allocatable :: Ord(:), KT(:,:), ST(:), IT(:)
  n = 0
  do iHalf=1,2
    n = n+sum(CIS(iState)%NOW(iHalf,1:SGS(iState)%nSym,1:CIS(iState)%nMidV))
  end do
  call mma_allocate(TRP_HW_K,3,max(1,n),Label='TRP_HW_K')
  call mma_allocate(TRP_HW_SY,max(1,n),Label='TRP_HW_SY')
  call mma_allocate(TRP_HW_ID,max(1,n),Label='TRP_HW_ID')
  k = 0
  do iHalf=1,2
    if (iHalf == 1) then
      lev0 = SGS(iState)%MidLev
      nL = SGS(iState)%nLev-SGS(iState)%MidLev
    else
      lev0 = 0
      nL = SGS(iState)%MidLev
    end if
    do MV=1,CIS(iState)%nMidV
      do ISY=1,SGS(iState)%nSym
        do w=1,CIS(iState)%NOW(iHalf,ISY,MV)
          Key = 0
          do iLev=lev0+nL,lev0+1,-1
            IOW0 = CIS(iState)%IOW(iHalf,ISY,MV)
            IPSHFT = 2*(iLev-1-lev0)
            IC = CIS(iState)%ICase(1+IOW0-CIS(iState)%nIpWlk+IPSHFT/(2*CSFG_nPack)+w*CIS(iState)%nIpWlk)
            IPSHFT = mod(IPSHFT,2*CSFG_nPack)
            Key = 4*Key+mod(IC/2**IPSHFT,4)
          end do
          k = k+1
          TRP_HW_K(1,k) = Key
          TRP_HW_K(2,k) = MV
          TRP_HW_K(3,k) = iHalf
          TRP_HW_SY(k) = ISY
          TRP_HW_ID(k) = w
        end do
      end do
    end do
  end do
  call mma_allocate(Ord,max(1,n),Label='TRP_HWOrd')
  call Sort_Keys(TRP_HW_K,3,n,Ord)
  call mma_allocate(KT,3,max(1,n),Label='TRP_KT')
  call mma_allocate(ST,max(1,n),Label='TRP_ST')
  call mma_allocate(IT,max(1,n),Label='TRP_IT')
  KT(:,1:n) = TRP_HW_K(:,Ord(1:n))
  ST(1:n) = TRP_HW_SY(Ord(1:n))
  IT(1:n) = TRP_HW_ID(Ord(1:n))
  TRP_HW_K(:,1:n) = KT(:,1:n)
  TRP_HW_SY(1:n) = ST(1:n)
  TRP_HW_ID(1:n) = IT(1:n)
  TRP_HW_n = n
  call mma_deallocate(KT)
  call mma_deallocate(ST)
  call mma_deallocate(IT)
  call mma_deallocate(Ord)
end subroutine TRP_HW_Setup

!> Symmetry and index of the half walk (iHalf) with packed steps Key ending
!> (upper) or starting (lower) at midvertex MV; ISY = 0 if none.
subroutine TRP_HW_Find(iHalf,MV,Key,ISY,ID)
  integer(kind=iwp), intent(in) :: iHalf, MV, Key
  integer(kind=iwp), intent(out) :: ISY, ID
  integer(kind=iwp) :: iLo, iHi, iMid, c
  ISY = 0
  ID = 0
  iLo = 1
  iHi = TRP_HW_n
  do while (iLo <= iHi)
    iMid = (iLo+iHi)/2
    ! compare (iHalf, MV, Key), the most significant first
    c = TRP_HW_K(3,iMid)-iHalf
    if (c == 0) c = TRP_HW_K(2,iMid)-MV
    if (c == 0) c = TRP_HW_K(1,iMid)-Key
    if (c == 0) then
      ISY = TRP_HW_SY(iMid)
      ID = TRP_HW_ID(iMid)
      return
    else if (c < 0) then
      iLo = iMid+1
    else
      iHi = iMid-1
    end if
  end do
end subroutine TRP_HW_Find

!> One level k of the generation of the tuples: vertices vL, vR at the top
!> of level k, class iCls, phase (1: above the loop head, top segments; 2:
!> inside the loop; 3: below the loop tail), tail level i and head level j
!> (0 if not yet passed), product V. TRP_SL/TRP_SR hold the steps.
recursive subroutine TRP_Gen(iState,k,vL,vR,iCls,iPhase,i,j,V)
  integer(kind=iwp), intent(in) :: iState, k, vL, vR, iCls, iPhase, i, j
  real(kind=wp), intent(in) :: V
  integer(kind=iwp) :: kBot, ISGT, iSet, iBeg, iEnd, vLB, vRB, nPh, ni, nj
  kBot = merge(SGS(iState)%MidLev+1,1,TRP_Half == 1)
  if (k < kBot) then
    call TRP_Leaf(iState,vL,vR,iCls,iPhase,i,j,V)
    return
  end if
  do iSet=1,3
    ! the segment sets allowed in this phase
    select case (iPhase)
      case (1)
        select case (iSet)
          case (1)
            iBeg = iSegTTopBeg
            iEnd = iSegTTopEnd
          case (2)
            iBeg = iSegTHeadBeg
            iEnd = iSegTHeadEnd
          case default
            iBeg = iSegTWBeg
            iEnd = iSegTWEnd
        end select
      case (2)
        if (iSet == 3) cycle
        iBeg = merge(9,19,iSet == 1)
        iEnd = merge(18,22,iSet == 1)
      case default
        if (iSet > 1) cycle
        iBeg = 23
        iEnd = 26
    end select
    do ISGT=iBeg,iEnd
      if (ITVPT(ISGT) /= iCls) cycle
      vLB = CIS(iState)%ISGM(vL,ISGT)
      if (vLB == 0) cycle
      vRB = SGS(iState)%Down(vR,IC2(ISGT))
      if (vRB == 0) cycle
      TRP_SL(k) = IC1(ISGT)
      TRP_SR(k) = IC2(ISGT)
      nPh = iPhase
      ni = i
      nj = j
      if ((iPhase == 1) .and. (iSet == 2)) then
        nPh = 2
        nj = k
      else if ((iPhase == 1) .and. (iSet == 3)) then
        nPh = 3
        ni = k
        nj = k
      else if ((iPhase == 2) .and. (iSet == 2)) then
        nPh = 3
        ni = k
      end if
      call TRP_Gen(iState,k-1,vLB,vRB,IBVPT(ISGT),nPh,ni,nj,V*CIS(iState)%VSGM(vL,ISGT))
    end do
  end do
end subroutine TRP_Gen

!> A complete pair of half walks: store the tuple.
subroutine TRP_Leaf(iState,vL,vR,iCls,iPhase,i,j,V)
  integer(kind=iwp), intent(in) :: iState, vL, vR, iCls, iPhase, i, j
  real(kind=wp), intent(in) :: V
  integer(kind=iwp) :: KeyL, KeyR, k, k1, k2, MVL, MVR, SYL, SYR, IL, IR, iType, nLev, nMidV
  nLev = SGS(iState)%nLev
  nMidV = CIS(iState)%nMidV
  if (TRP_Half == 1) then
    k1 = SGS(iState)%MidLev+1
    k2 = nLev
    MVL = vL-SGS(iState)%MVSta+1
    MVR = vR-SGS(iState)%MVSta+1
    select case (iPhase)
      case (1)
        iType = 1
      case (2)
        iType = 2
      case default
        iType = 3
    end select
  else
    ! a lower factor must end below the loop tail (at the graph tail)
    if (iPhase /= 3) return
    k1 = 1
    k2 = SGS(iState)%MidLev
    MVL = TRP_MVL0
    MVR = 0
    iType = merge(5,4,TRP_Start <= 2)
  end if
  KeyL = 0
  KeyR = 0
  do k=k2,k1,-1
    KeyL = 4*KeyL+TRP_SL(k)
    KeyR = 4*KeyR+TRP_SR(k)
  end do
  if (TRP_Half == 2) then
    ! the right midvertex of the lower factor: partner of the start class
    select case (TRP_Start)
      case (1,2)
        MVR = CIS(iState)%IVR(MVL+SGS(iState)%MVSta-1,TRP_Start)
      case (5)
        MVR = CIS(iState)%IVR(MVL+SGS(iState)%MVSta-1,3)
      case (6)
        MVR = CIS(iState)%IVR(MVL+SGS(iState)%MVSta-1,4)
      case default
        MVR = MVL+SGS(iState)%MVSta-1
    end select
    MVR = MVR-SGS(iState)%MVSta+1
  end if
  call TRP_HW_Find(TRP_Half,MVL,KeyL,SYL,IL)
  call TRP_HW_Find(TRP_Half,MVR,KeyR,SYR,IR)
  if ((SYL == 0) .or. (SYR == 0)) then
    write(u6,*) 'TRP_Leaf: half walk not found',TRP_Half,MVL,MVR,KeyL,KeyR
    call Abend()
  end if
  select case (iType)
    case (1)
      k = TRP_Cat(nLev,nMidV,1,iCls,0,0,MVL)
    case (2)
      k = TRP_Cat(nLev,nMidV,2,iCls,0,j,MVL)
    case (3)
      k = TRP_Cat(nLev,nMidV,3,0,i,j,MVL)
    case (4)
      k = TRP_Cat(nLev,nMidV,4,TRP_Start,i,j,MVL)
    case default
      k = TRP_Cat(nLev,nMidV,5,TRP_Start,i,0,MVL)
  end select
  if (iType == 1 .and. iCls < 4) return
  call TRP_Push(k,IL,IR,SYL,SYR,MVR,V)
end subroutine TRP_Leaf

!> First tuple with key Key and the number of them (binary search).
subroutine TRP_Range(Key,iFirst,nT)
  integer(kind=iwp), intent(in) :: Key
  integer(kind=iwp), intent(out) :: iFirst, nT
  integer(kind=iwp) :: iLo, iHi, iMid
  iLo = 1
  iHi = TRP_n+1
  do while (iLo < iHi)
    iMid = (iLo+iHi)/2
    if (TRP_Key(iMid) < Key) then
      iLo = iMid+1
    else
      iHi = iMid
    end if
  end do
  iFirst = iLo
  nT = 0
  do while (iFirst+nT <= TRP_n)
    if (TRP_Key(iFirst+nT) /= Key) exit
    nT = nT+1
  end do
end subroutine TRP_Range

!> Sigma = Sigma + X T_pq CI (T_pq = a+_pa a_qa - a+_pb a_qb) for a CI
!> vector of symmetry ISYCI of state iState, from the split-graph coupling
!> tables of the triplet operator (TRP_Setup). Raising operators (p <= q)
!> map the right (ket) walks of the tuples onto the left (bra) ones,
!> lowering operators the reverse (transpose).
subroutine SG_Tpq_Psi(iState,p,q,X,ISYCI,CI,Sigma)
  integer(kind=iwp), intent(in) :: iState, p, q, ISYCI
  real(kind=wp), intent(in) :: X, CI(*)
  real(kind=wp), intent(inout) :: Sigma(*)
  integer(kind=iwp) :: i, j, nLev, nMidV, MidLev, MV, iCls, iU, nU, iL, nL, kU, kL, ISYD, w, nW
  logical(kind=iwp) :: Raise
  call TRP_Setup(iState)
  nLev = SGS(iState)%nLev
  nMidV = CIS(iState)%nMidV
  MidLev = SGS(iState)%MidLev
  i = min(p,q)
  j = max(p,q)
  Raise = (p <= q)
  do MV=1,nMidV
    if (j <= MidLev) then
      do iCls=4,6
        call TRP_Range(TRP_Cat(nLev,nMidV,1,iCls,0,0,MV),iU,nU)
        if (nU == 0) cycle
        call TRP_Range(TRP_Cat(nLev,nMidV,4,iCls,i,j,MV),iL,nL)
        do kU=iU,iU+nU-1
          do kL=iL,iL+nL-1
            call TRP_Add(iState,Raise,ISYCI,X,CI,Sigma,MV,TRP_MVR(kU),TRP_SYL(kU),TRP_IL(kU),TRP_SYL(kL),TRP_IL(kL), &
                         TRP_SYR(kU),TRP_IR(kU),TRP_SYR(kL),TRP_IR(kL),TRP_V(kU)*TRP_V(kL))
          end do
        end do
      end do
    else if (i <= MidLev) then
      do iCls=1,2
        call TRP_Range(TRP_Cat(nLev,nMidV,2,iCls,0,j,MV),iU,nU)
        if (nU == 0) cycle
        call TRP_Range(TRP_Cat(nLev,nMidV,5,iCls,i,0,MV),iL,nL)
        do kU=iU,iU+nU-1
          do kL=iL,iL+nL-1
            call TRP_Add(iState,Raise,ISYCI,X,CI,Sigma,MV,TRP_MVR(kU),TRP_SYL(kU),TRP_IL(kU),TRP_SYL(kL),TRP_IL(kL), &
                         TRP_SYR(kU),TRP_IR(kU),TRP_SYR(kL),TRP_IR(kL),TRP_V(kU)*TRP_V(kL))
          end do
        end do
      end do
    else
      call TRP_Range(TRP_Cat(nLev,nMidV,3,0,i,j,MV),iU,nU)
      do kU=iU,iU+nU-1
        do ISYD=1,SGS(iState)%nSym
          nW = CIS(iState)%NOW(2,ISYD,MV)
          do w=1,nW
            call TRP_Add(iState,Raise,ISYCI,X,CI,Sigma,MV,TRP_MVR(kU),TRP_SYL(kU),TRP_IL(kU),ISYD,w,TRP_SYR(kU),TRP_IR(kU), &
                         ISYD,w,TRP_V(kU))
          end do
        end do
      end do
    end if
  end do

end subroutine SG_Tpq_Psi

!> SG_Tpq_Psi: the contribution of one pair of CSFs, left (bra of the
!> raising operator) at midvertex MV with upper walk (SUL,IUL) and lower
!> walk (SDL,IDL), right at midvertex MVR_ with (SUR,IUR), (SDR,IDR).
subroutine TRP_Add(iState,Raise,ISYCI,X,CI,Sigma,MV,MVR_,SUL,IUL,SDL,IDL,SUR,IUR,SDR,IDR,V)
  integer(kind=iwp), intent(in) :: iState, ISYCI, MV, MVR_, SUL, IUL, SDL, IDL, SUR, IUR, SDR, IDR
  logical(kind=iwp), intent(in) :: Raise
  real(kind=wp), intent(in) :: X, CI(*), V
  real(kind=wp), intent(inout) :: Sigma(*)
  integer(kind=iwp) :: ISYL, ISYR, IL_, IR_
  ISYL = Mul(SUL,SDL)
  ISYR = Mul(SUR,SDR)
  if (Raise) then
    if (ISYR /= ISYCI) return
  else
    if (ISYL /= ISYCI) return
  end if
  IL_ = CIS(iState)%IOCSF(SUL,MV,ISYL)+IUL+(IDL-1)*CIS(iState)%NOW(1,SUL,MV)
  IR_ = CIS(iState)%IOCSF(SUR,MVR_,ISYR)+IUR+(IDR-1)*CIS(iState)%NOW(1,SUR,MVR_)
  if (Raise) then
    Sigma(IL_) = Sigma(IL_)+X*V*CI(IR_)
  else
    Sigma(IR_) = Sigma(IR_)+X*V*CI(IL_)
  end if
end subroutine TRP_Add

!> One level k of the generation of the partner walks of SG_Xpq_Psi: the
!> partner steps P(1:k-1) and their b values BP(0:k-1) are set, the walk is
!> at vertex iVert (level k-1), with the product V of the segment values so
!> far. F/BF: the fixed walk. At the graph head Sigma(I) += Y*V.
recursive subroutine ASTA_Gen(SGS,k,iVert,V,nLev,F,BF,P,BP,i,j,FixedIsKet,Triplet,ISYS,nS,Y,Sigma)
  type(SGStruct), intent(in) :: SGS
  integer(kind=iwp), intent(in) :: k, iVert, nLev, F(nLev), BF(0:nLev), i, j, ISYS, nS
  real(kind=wp), intent(in) :: V, Y
  integer(kind=iwp), intent(inout) :: P(nLev), BP(0:nLev)
  logical(kind=iwp), intent(in) :: FixedIsKet, Triplet
  real(kind=wp), intent(inout) :: Sigma(*)
  integer(kind=iwp) :: s, iV, u, t, bU, bT, I_
  real(kind=wp) :: W
  if (k > nLev) then
    if (BP(nLev) /= BF(nLev)) return
    I_ = ASTA_Find(ISYS,nS,Walk_Key(nLev,P))
    if (I_ == 0) then
      write(u6,*) 'ASTA_Gen: generated walk not found in symmetry',ISYS
      call Abend()
    end if
    Sigma(I_) = Sigma(I_)+Y*V
    return
  end if
  do s=0,3
    if ((k < i) .and. (s /= F(k))) cycle
    if ((.not. Triplet) .and. (k > j) .and. (s /= F(k))) cycle
    iV = SGS%Up(iVert,s)
    if (iV == 0) cycle
    P(k) = s
    select case (s)
      case (1)
        BP(k) = BP(k-1)+1
      case (2)
        BP(k) = BP(k-1)-1
      case default
        BP(k) = BP(k-1)
    end select
    if (FixedIsKet) then
      u = s
      t = F(k)
      bU = BP(k)
      bT = BF(k)
    else
      u = F(k)
      t = s
      bU = BF(k)
      bT = BP(k)
    end if
    if (k < i) then
      W = One
    else if (k == i) then
      if (i == j) then
        if (Triplet) then
          W = ASTA_Weight_T(u,t,bT-bU,bT)
        else if (u == t) then
          W = real((1+u)/2,kind=wp)
        else
          W = Zero
        end if
      else
        W = ASTA_Tail(u,t,bT)
      end if
    else if (k < j) then
      W = ASTA_Mid(u,t,bT-bU,bT)
    else if (k == j) then
      if (Triplet) then
        W = ASTA_Head_T(u,t,bT-bU,bT)
      else if (bT == bU) then
        W = ASTA_Head(u,t,bT)
      else
        W = Zero
      end if
    else if (Triplet) then
      W = ASTA_Top_T(u,t,bT-bU,bT)
    else
      W = One
    end if
    if (W == Zero) cycle
    call ASTA_Gen(SGS,k+1,iV,V*W,nLev,F,BF,P,BP,i,j,FixedIsKet,Triplet,ISYS,nS,Y,Sigma)
  end do
end subroutine ASTA_Gen

end module SGUGA_Triplet
