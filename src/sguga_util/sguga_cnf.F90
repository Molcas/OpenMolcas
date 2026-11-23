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
! Configuration lists of a split graph and the explicit CI block
! (PHPCSF selection, block of H over the selected configurations).
! (Default public, as SGUGA_Core; the interface is the module SGUGA.)
module SGUGA_Cnf

use SGUGA_Core
use SGUGA_MatEl
use SGUGA_Tables

implicit none

contains

!> The configurations of the CSFs of symmetry LSM of state iState, from the
!> occupations of their walks (see SGCnfList); independent of any other CI
!> library. Also sets up the single-CSF matrix elements (SG_MatEl_Init).
subroutine SG_Mk_Cnf_List(iState,LSM,Cnf)
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, LSM
  type(SGCnfList), intent(inout) :: Cnf
  integer(kind=iwp), parameter :: nPerWord = (bit_size(0_iwp)-2)/2
  integer(kind=iwp) :: nC, nL, nW, I, k, iLev, iW
  integer(kind=iwp), allocatable :: Occ(:), KeyC(:,:), Ord(:)

  call SG_Free_Cnf_List(Cnf)
  call SG_MatEl_Init(iState,LSM)
  nC = CIS(iState)%NCSF(LSM)
  nL = SGS(iState)%nLev
  nW = max(1,(nL+nPerWord-1)/nPerWord)
  Cnf%LSM = LSM
  Cnf%nLev = nL
  Cnf%nWord = nW

  ! packed occupation of each CSF
  call mma_allocate(Occ,nL,Label='Cnf_Occ')
  call mma_allocate(KeyC,nW,max(1,nC),Label='Cnf_KeyC')
  KeyC(:,:) = 0
  do I=1,nC
    call CSF_Occupation(LSM,I,Occ)
    do iLev=1,nL
      iW = (iLev-1)/nPerWord+1
      KeyC(iW,I) = KeyC(iW,I)+Occ(iLev)*4_iwp**(iLev-1-(iW-1)*nPerWord)
    end do
  end do
  call mma_deallocate(Occ)

  ! CSFs sorted by key (ties by CSF index), configurations = runs of equal keys
  call mma_allocate(Ord,max(1,nC),Label='Cnf_Ord')
  call Sort_Keys(KeyC,nW,nC,Ord)
  Cnf%nCnf = 0
  do I=1,nC
    if (I == 1) then
      Cnf%nCnf = 1
    else if (any(KeyC(:,Ord(I)) /= KeyC(:,Ord(I-1)))) then
      Cnf%nCnf = Cnf%nCnf+1
    end if
  end do
  call mma_allocate(Cnf%iOff,Cnf%nCnf+1,Label='Cnf%iOff')
  call mma_allocate(Cnf%iCSF,max(1,nC),Label='Cnf%iCSF')
  call mma_allocate(Cnf%Key,nW,max(1,Cnf%nCnf),Label='Cnf%Key')
  k = 0
  do I=1,nC
    if (I == 1) then
      k = 1
      Cnf%iOff(1) = 0
      Cnf%Key(:,1) = KeyC(:,Ord(1))
    else if (any(KeyC(:,Ord(I)) /= KeyC(:,Ord(I-1)))) then
      k = k+1
      Cnf%iOff(k) = I-1
      Cnf%Key(:,k) = KeyC(:,Ord(I))
    end if
    Cnf%iCSF(I) = Ord(I)
  end do
  Cnf%iOff(Cnf%nCnf+1) = nC
  call mma_deallocate(Ord)
  call mma_deallocate(KeyC)

end subroutine SG_Mk_Cnf_List

subroutine SG_Free_Cnf_List(Cnf)
  use stdalloc, only: mma_deallocate
  type(SGCnfList), intent(inout) :: Cnf
  call mma_deallocate(Cnf%iOff,safe='*')
  call mma_deallocate(Cnf%iCSF,safe='*')
  call mma_deallocate(Cnf%Key,safe='*')
  Cnf%nCnf = 0
end subroutine SG_Free_Cnf_List

!> Occupation numbers (0, 1, 2 per level) of configuration k.
subroutine SG_Cnf_Occ(Cnf,k,Occ)
  type(SGCnfList), intent(in) :: Cnf
  integer(kind=iwp), intent(in) :: k
  integer(kind=iwp), intent(out) :: Occ(:)
  integer(kind=iwp), parameter :: nPerWord = (bit_size(0_iwp)-2)/2
  integer(kind=iwp) :: iLev, iW
  do iLev=1,Cnf%nLev
    iW = (iLev-1)/nPerWord+1
    Occ(iLev) = int(mod(Cnf%Key(iW,k)/4_iwp**(iLev-1-(iW-1)*nPerWord),4_iwp),kind=iwp)
  end do
end subroutine SG_Cnf_Occ

!> Explicit subspace for the Davidson start, with the rule of PHPCSF: whole
!> configurations in the order of increasing key (ties within Acc in the
!> order of the configurations) until the next one does not fit into MxDim
!> CSFs. The key of a configuration is the mean of the exact diagonal
!> elements Diag (split-graph order) of its CSFs, the trace of H over the
!> configuration divided by its number of CSFs, independent of the spin
!> coupling scheme. On output IPCNF(1:NPCNF) are the configurations and
!> IPCSF(1:NPCSF) their CSFs (split-graph indices, configuration by
!> configuration).
subroutine SG_Select_Explicit(Cnf,Diag,MxDim,NPCNF,IPCNF,NPCSF,IPCSF)
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: One
  type(SGCnfList), intent(in) :: Cnf
  real(kind=wp), intent(in) :: Diag(*)
  integer(kind=iwp), intent(in) :: MxDim
  integer(kind=iwp), intent(out) :: NPCNF, IPCNF(*), NPCSF, IPCSF(*)
  real(kind=wp), parameter :: Acc = 1.0e-13_wp
  integer(kind=iwp) :: k, kMin, n
  real(kind=wp) :: XMax, XMin
  real(kind=wp), allocatable :: Key(:)

  call mma_allocate(Key,max(1,Cnf%nCnf),Label='Cnf_Key')
  do k=1,Cnf%nCnf
    n = Cnf%iOff(k+1)-Cnf%iOff(k)
    Key(k) = sum(Diag(Cnf%iCSF(Cnf%iOff(k)+1:Cnf%iOff(k+1))))/real(n,kind=wp)
  end do
  XMax = maxval(Key(1:Cnf%nCnf))
  NPCNF = 0
  NPCSF = 0
  do while (NPCNF < Cnf%nCnf)
    XMin = XMax+One
    kMin = 0
    do k=1,Cnf%nCnf
      if (Key(k)+Acc < XMin) then
        XMin = Key(k)
        kMin = k
      end if
    end do
    if (kMin == 0) exit
    n = Cnf%iOff(kMin+1)-Cnf%iOff(kMin)
    if (NPCSF+n > MxDim) exit
    NPCNF = NPCNF+1
    IPCNF(NPCNF) = kMin
    IPCSF(NPCSF+1:NPCSF+n) = Cnf%iCSF(Cnf%iOff(kMin)+1:Cnf%iOff(kMin+1))
    NPCSF = NPCSF+n
    Key(kMin) = XMax+One
  end do
  call mma_deallocate(Key)

end subroutine SG_Select_Explicit

!> <I|H|J> over the split-graph CSFs iG(1:nG) of symmetry LSM, given
!> configuration by configuration (nGCnf(1:nGC) CSFs each, all CSFs of a
!> configuration sharing their occupations), packed as a lower triangle
!> H(nTri(nG)), without core energy. TU/TUVX in the symmetry-blocked order
!> of RASSCF (reordered to the levels here). Pairs of configurations more
!> than a double excitation apart are skipped. Requires SG_MatEl_Init.
subroutine SG_Block_H(LSM,nGC,nGCnf,nG,iG,TU,nTU,TUVX,nTUVX,H)
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: LSM, nGC, nGCnf(nGC), nG, iG(nG), nTU, nTUVX
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX)
  real(kind=wp), intent(out) :: H(nG*(nG+1)/2)
  integer(kind=iwp) :: nL, kL, kR, IILB, IIRB, IIL, IIR, IIRMAX, I, J, nX, t
  integer(kind=iwp), allocatable :: Occ(:,:)
  real(kind=wp), allocatable :: TU_L(:), TUVX_L(:)

  nL = SGS(ME_iState)%nLev
  call mma_allocate(TU_L,nTU,Label='TU_L')
  call mma_allocate(TUVX_L,nTUVX,Label='TUVX_L')
  call SG_Level_Ints(SGS(ME_iState),TU,nTU,TUVX,nTUVX,TU_L,TUVX_L)
  call mma_allocate(Occ,nL,max(1,nGC),Label='Blk_Occ')
  IILB = 1
  do kL=1,nGC
    call CSF_Occupation(LSM,iG(IILB),Occ(:,kL))
    IILB = IILB+nGCnf(kL)
  end do
  IILB = 1
  do kL=1,nGC
    IIRB = 1
    do kR=1,kL
      ! excitation level between the two configurations
      nX = 0
      do t=1,nL
        nX = nX+max(0_iwp,Occ(t,kL)-Occ(t,kR))
      end do
      do IIL=1,nGCnf(kL)
        I = IILB-1+IIL
        IIRMAX = merge(IIL,nGCnf(kR),IILB == IIRB)
        do IIR=1,IIRMAX
          J = IIRB-1+IIR
          if (nX > 2) then
            H(I*(I-1)/2+J) = Zero
          else
            H(I*(I-1)/2+J) = MatEl_H_Direct(iG(I),iG(J),Occ(:,kL),Occ(:,kR),TU_L,nTU,TUVX_L,nTUVX)
          end if
        end do
      end do
      IIRB = IIRB+nGCnf(kR)
    end do
    IILB = IILB+nGCnf(kL)
  end do
  call mma_deallocate(Occ)
  call mma_deallocate(TUVX_L)
  call mma_deallocate(TU_L)

end subroutine SG_Block_H

!> Explicit Hamiltonian over the CSFs of the configurations IPCNF(1:NPCNF)
!> of Cnf (SG_Select_Explicit), in the order of IPCSF there, packed as a
!> lower triangle, with ECORE on the diagonal; SGUGA only.
subroutine SG_Explicit_Block(iState,Cnf,NPCNF,IPCNF,ECORE,TU,nTU,TUVX,nTUVX,H,iRel,RelMap)
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, NPCNF, IPCNF(NPCNF), nTU, nTUVX
  type(SGCnfList), intent(in) :: Cnf
  real(kind=wp), intent(in) :: ECORE, TU(nTU), TUVX(nTUVX)
  real(kind=wp), intent(out) :: H(*)
  ! GAS: the matrix elements are evaluated in the relaxed space iRel (it
  ! holds every intermediate state of the products), RelMap(I) the index
  ! there of CSF I of iState
  integer(kind=iwp), intent(in), optional :: iRel, RelMap(:)
  integer(kind=iwp) :: k, n, nG, I
  integer(kind=iwp), allocatable :: nGCnf(:), iG(:)

  if (present(iRel)) then
    call SG_MatEl_Init(iRel,Cnf%LSM)
  else
    call SG_MatEl_Init(iState,Cnf%LSM)
  end if
  call mma_allocate(nGCnf,max(1,NPCNF),Label='nGCnf')
  nG = 0
  do k=1,NPCNF
    nGCnf(k) = Cnf%iOff(IPCNF(k)+1)-Cnf%iOff(IPCNF(k))
    nG = nG+nGCnf(k)
  end do
  call mma_allocate(iG,max(1,nG),Label='iG')
  nG = 0
  do k=1,NPCNF
    n = nGCnf(k)
    iG(nG+1:nG+n) = Cnf%iCSF(Cnf%iOff(IPCNF(k))+1:Cnf%iOff(IPCNF(k)+1))
    nG = nG+n
  end do
  if (present(RelMap)) iG(1:nG) = RelMap(iG(1:nG))
  call SG_Block_H(Cnf%LSM,NPCNF,nGCnf,nG,iG,TU,nTU,TUVX,nTUVX,H)
  do I=1,nG
    H(I*(I+1)/2) = H(I*(I+1)/2)+ECORE
  end do
  call mma_deallocate(iG)
  call mma_deallocate(nGCnf)

end subroutine SG_Explicit_Block

!> Configuration-average diagonal (the CI-solver diagonal of Lucia): for
!> each configuration of Cnf the average of the diagonal elements of its
!> determinants with M_S = S, assigned to all its CSFs (Diag, split-graph
!> order):
!>
!>   E_avg = Sum_i n_i h_ii + Sum_{i<j} n_i n_j (ii|jj) + Sum_{n_i=2} (ii|ii)
!>           - Sum_{i<j} c_ij (ij|ji),
!>
!> c_ij = 2 (both doubly occupied), 1 (one doubly, one singly), and for two
!> open shells the probability that they have the same spin,
!> [na(na-1) + nb(nb-1)]/[No(No-1)], with No open shells, na = (No+2S)/2,
!> nb = (No-2S)/2. TU/TUVX in the symmetry-blocked order of RASSCF; no core
!> energy.
subroutine SG_Cnf_Avg_Diag(iState,Cnf,TU,nTU,TUVX,nTUVX,Diag)
  use stdalloc, only: mma_allocate, mma_deallocate
  integer(kind=iwp), intent(in) :: iState, nTU, nTUVX
  type(SGCnfList), intent(in) :: Cnf
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX)
  real(kind=wp), intent(inout) :: Diag(*)
  integer(kind=iwp) :: k, i, j, nL, nO, na, nb
  real(kind=wp) :: E, PSame, cij
  integer(kind=iwp), allocatable :: Occ(:)
  real(kind=wp), allocatable :: TU_L(:), TUVX_L(:)

  nL = SGS(iState)%nLev
  call mma_allocate(TU_L,nTU,Label='TU_L')
  call mma_allocate(TUVX_L,nTUVX,Label='TUVX_L')
  call SG_Level_Ints(SGS(iState),TU,nTU,TUVX,nTUVX,TU_L,TUVX_L)
  call mma_allocate(Occ,nL,Label='Avg_Occ')
  do k=1,Cnf%nCnf
    call SG_Cnf_Occ(Cnf,k,Occ)
    nO = count(Occ(1:nL) == 1)
    na = (nO+SGS(iState)%iSpin-1)/2
    nb = (nO-SGS(iState)%iSpin+1)/2
    PSame = Zero
    if (nO > 1) PSame = real(na*(na-1)+nb*(nb-1),kind=wp)/real(nO*(nO-1),kind=wp)
    E = Zero
    do i=1,nL
      if (Occ(i) == 0) cycle
      E = E+Occ(i)*TU_L(Tri_SG(i,i))
      if (Occ(i) == 2) E = E+TUVX_L(Tri_SG(Tri_SG(i,i),Tri_SG(i,i)))
      do j=i+1,nL
        if (Occ(j) == 0) cycle
        E = E+Occ(i)*Occ(j)*TUVX_L(Tri_SG(Tri_SG(i,i),Tri_SG(j,j)))
        if ((Occ(i) == 2) .and. (Occ(j) == 2)) then
          cij = Two
        else if ((Occ(i) == 2) .or. (Occ(j) == 2)) then
          cij = One
        else
          cij = PSame
        end if
        E = E-cij*TUVX_L(Tri_SG(Tri_SG(i,j),Tri_SG(j,i)))
      end do
    end do
    Diag(Cnf%iCSF(Cnf%iOff(k)+1:Cnf%iOff(k+1))) = E
  end do
  call mma_deallocate(Occ)
  call mma_deallocate(TUVX_L)
  call mma_deallocate(TU_L)
end subroutine SG_Cnf_Avg_Diag

end module SGUGA_Cnf
