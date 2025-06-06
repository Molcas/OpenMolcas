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
! Copyright (C) 1990,1991,1992,2000,2007, Roland Lindh                 *
!               1990, IBM                                              *
!***********************************************************************

!#define _DEBUGPRINT_
!#define _CD_TIMING_
subroutine Drvg1_3Center_RI(Temp,nGrad,ij2,nij_Eff)
!***********************************************************************
!                                                                      *
!  Object: driver for two-electron integrals. The four outermost loops *
!          will control the type of the two-electron integral, eg.     *
!          (ss|ss), (sd|pp), etc. The next four loops will generate    *
!          list of symmetry distinct centers that do have basis func-  *
!          tions of the requested type.                                *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
! For RI-HF gradients read:                                            *
! "Analytical Gradients of Hartee-Fock Exchange with Density Fitting   *
! Approximations", J. Bostrom, F. Aquilante, T. B. Pedersen and R.     *
! Lindh, JCTC,  9:204-212 (2013).                                      *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified for k2 loop. August '91                         *
!             Modified for gradient calculation. January '92           *
!             Modified for SetUp_Ints. January '00                     *
!             Modified for 3-center RI gradients, March '07            *
!***********************************************************************

use setup, only: mSkal, MxPrm
use Index_Functions, only: iTri, nTri_Elem
use iSD_data, only: iSD
use pso_stuff, only: DMdiag, lPSO, lSA, n_Txy, nCalAO, nG1, nnP, nZ_p_k, Thpkl, Txy, Z_p_k, ReadBPT2
use k2_arrays, only: Sew_Scr
use Basis_Info, only: nBas, nBas_Aux
use Sizes_of_Seward, only: S
use Gateway_Info, only: CutInt
use RICD_Info, only: Do_RI, RI_3C
use Symmetry_Info, only: nIrrep
use RI_glob, only: BklK, BMP2, CijK, CilK, CMOi, DMLT, DoCholExch, iAdrCVec, iBDsh, iMP2prpt, iOff_Ymnij, LuCVector, MxChVInShl, &
                   nAdens, nAvec, nChOrb, nIJ1, nIJR, nJdens, nKdens, nKvec, NumAuxVec, nYmnij, Yij, Ymnij
use Cholesky, only: nSym, ThrCom
use Data_Structures, only: Deallocate_DT
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nGrad, nij_Eff, ij2(2,nij_Eff)
real(kind=wp), intent(out) :: Temp(nGrad)
integer(kind=iwp) :: i, iAdrC, iAng, ib, id, iiQ, ij, ijQ, ijS, iMOleft, iMOright, iOpt, iS, iS_, ish, iSO, iSym, itmp, j, jb, &
                     jjQ, jS, jS_, jsh, jSym, jSym_s, KAux, klS, klS_, kS, kSym, lB_mp2, lBklK, lCijK, lCilK, lMaxDens, lS, &
                     maxnAct, maxnnP, mij, mj, MumOrb, MxBasSh, MxInShl, nAct(0:7), ni, nij, nIJ1Max, nIJRMax, nIMax, nj, nK, &
                     nnSkal, nPrev, nSkal, nSkal2, nSkal2_, nSkal_Auxiliary, nSkal_Valence, nThpkl, nTMax, NumOrb, NumOrb_i, nXki
real(kind=wp) :: A_int, A_int_ij, A_int_kl, Dm_ij, ExFac, PZmnij, SDGmn, ThrAO, TMax_all, XDm_ii, XDm_ij, XDm_jj, XDm_max, xfk, &
                 Xik, Xil, Xjk, Xjl
logical(kind=iwp) :: DoFock, DoGrad, FlipFlop, Found, Indexation
character(len=80) :: KSDFT
character(len=8) :: Method
integer(kind=iwp), allocatable :: LBList(:), Shij(:,:), Shij2(:,:)
real(kind=wp), allocatable :: CVec(:,:), CVec2(:,:,:), MaxDens(:), SDG(:), Thhalf(:), TMax_Auxiliary(:), TMax_Valence(:,:), &
                              Tmp(:,:), Xmi(:,:,:,:)
integer(kind=iwp), external :: Cho_irange
real(kind=wp), external :: Get_ExFac
logical(kind=iwp), external :: Rsv_Tsk2

!                                                                      *
!***********************************************************************
!                                                                      *
Temp(:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
xfk = 1.0e-3_wp  ! changing this parameter tunes LK-type screening thr
!xfk = 1.0e-12_wp ! Debugging
!xfk = Zero       ! Debugging
!                                                                      *
!***********************************************************************
!                                                                      *
call Get_dScalar('Cholesky Threshold',ThrCom)
ThrCom = max(ThrCom,1.0e-6_wp) ! not to sacrify efficiency too much
!                                                                      *
!***********************************************************************
!                                                                      *
! Handle mixed basis set

if (Do_RI) then
  call Set_Basis_Mode('Auxiliary')
  call Nr_Shells(nSkal_Auxiliary)
  call Set_Basis_Mode('WithAuxiliary')
  RI_3C = .true.
else
  call Set_Basis_Mode('Valence')
  nSkal_Auxiliary = 0
end if
call SetUp_iSD()
!                                                                      *
!***********************************************************************
!                                                                      *
! Precompute k2 entities.

Indexation = .false.
DoFock = .false.
DoGrad = .true.
ThrAO = Zero
call SetUp_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
nSkal_Valence = nSkal-nSkal_Auxiliary
mSkal = nSkal
!                                                                      *
!***********************************************************************
!                                                                      *
MxPrm = 0
do iAng=0,S%iAngMx
  MxPrm = max(MxPrm,S%MaxPrm(iAng))
end do
!                                                                      *
!***********************************************************************
!                                                                      *
maxnAct = 0
if (lPSO) then
  call Get_iArray('nAsh',nAct,nIrrep)
  maxnnP = nnP(0)
  maxnAct = nAct(0)
  do i=1,nIrrep-1
    maxnnP = max(maxnnP,nnP(i))
    maxnAct = max(maxnAct,nAct(i))
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute entities for prescreening at shell level

call mma_allocate(TMax_Valence,nSkal_Valence,nSkal_Valence,Label='TMax_Valence')
nTMax = nSkal_Auxiliary
if (Do_RI) nTMax = max(1,nTMax-1)
call mma_allocate(TMax_Auxiliary,nTMax,Label='TMax_Auxiliary')

call mma_allocate(Tmp,nSkal,nSkal,Label='Tmp')
call Shell_MxSchwz(nSkal,Tmp)

TMax_all = Zero
do iS=1,nSkal_Valence
  do jS=1,iS
    TMax_Valence(iS,jS) = Tmp(iS,jS)
    TMax_Valence(jS,iS) = Tmp(iS,jS)
    TMax_all = max(TMax_all,Tmp(iS,jS))
  end do
end do
if (Do_RI) then
  do iS=1,nSkal_Auxiliary-1
    iS_ = iS+nSkal_Valence
    jS_ = nSkal_Valence+nSkal_Auxiliary
    TMax_Auxiliary(iS) = Tmp(iS_,jS_)
  end do
end if

call mma_deallocate(Tmp)

!                                                                      *
!***********************************************************************
!                                                                      *
MxInShl = 1
do i=1,nSkal_Valence
  MxInShl = max(MxInShl,iSD(3,i)*iSD(2,i))
end do

! Calculate maximum density value for each shellpair

lMaxDens = nTri_Elem(nSkal_Valence)
call mma_allocate(MaxDens,lMaxDens,Label='MaxDens')
MaxDens(:) = Zero

do iSym=0,nSym-1
  kS = 1+nSkal_Valence*iSym ! note diff wrt declaration of iBDsh
  do j=1,nBas(iSym)
    jsh = Cho_Irange(j,iBDsh(kS),nSkal_Valence,.true.)
    do i=1,j
      ish = Cho_Irange(i,iBDsh(kS),nSkal_Valence,.true.)
      ijS = nTri_Elem(jsh-1)+ish
      do iSO=1,nJDens
        if (.not. DMLT(iSO)%Active) cycle
        ij = iTri(j,i)
        Dm_ij = abs(DMLT(iSO)%SB(iSym+1)%A1(ij))
        MaxDens(ijS) = max(MaxDens(ijS),Dm_ij)
      end do
    end do
  end do
end do
call mma_deallocate(iBDsh)

do i=1,5
  if (DMLT(i)%Active) call Deallocate_DT(DMLT(i))
end do

! Create list of non-vanishing pairs

! 1) For the valence basis set

call mma_allocate(Shij,2,nTri_Elem(nSkal_Valence),Label='Shij')
nSkal2 = 0
do iS=1,nSkal_Valence
  iiQ = nTri_Elem(iS)
  XDm_ii = MaxDens(iiQ)
  do jS=1,iS
    jjQ = nTri_Elem(jS)
    XDm_jj = MaxDens(jjQ)
    ijQ = iTri(iS,jS)
    XDm_ij = MaxDens(ijQ)
    XDm_max = max(XDm_ij,XDm_ii,XDm_jj)
    A_int_ij = TMax_Valence(iS,jS)
    if (TMax_All*A_int_ij >= CutInt) then

      ! FAQ: more aggressive screening to discard shprs formed
      !      by AOs contributing mainly to the virtual MO space.

      if (A_int_ij*XDm_max >= CutInt) then
        nSkal2 = nSkal2+1
        Shij(1,nSkal2) = iS
        Shij(2,nSkal2) = jS
      end if

    end if
  end do
end do

! 2) For the auxiliary basis set

if (Do_RI) then
  mij = nSkal_Auxiliary
  call mma_allocate(Shij2,2,mij,Label='Shij2')
  nij = 0
  do jS=nSkal_Valence+1,nSkal-1
    if (TMax_All*TMax_Auxiliary(jS-nSkal_Valence) >= CutInt) then
      nij = nij+1
      Shij2(1,nij) = nSkal
      Shij2(2,nij) = jS
    end if
  end do
else
  mij = nij_Eff
  call mma_allocate(Shij2,2,mij,Label='Shij2')
  nij = 0
  do ij=1,nij_Eff
    iS = ij2(1,ij)
    jS = ij2(2,ij)
    if (TMax_All*TMax_Valence(iS,jS) >= CutInt) then
      nij = nij+1
      Shij2(1,nij) = iS
      Shij2(2,nij) = jS
    end if
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (DoCholExch) then

  ! Find the largest number of contractions in any given shell of
  ! auxiliary functions.

  MxChVInShl = 1
  if (Do_RI) then
    do i=nSkal_Valence+1,nSkal_Valence+nSkal_Auxiliary
      MxChVInShl = max(MxChVInShl,iSD(3,i))
    end do
  else
    write(u6,*) 'Not implemented for Cholesky yet!'
    call Abend()
  end if

  ! Find the largest set of ij. The i and j basis are due to the CD
  ! of the one-particle density matrix.

  ! nIJ1: diagonal blocks are triangularized.
  ! nIJR: diagonal blocks are square.
  ! nIMax: largest number of i basis in any irrep.

  nIJ1Max = 0
  nIJRMax = 0
  nIMax = 0
  do iSym=1,nIrrep
    do iSO=1,nKDens
      nIMax = max(nIMax,nChOrb(iSym-1,iSO))
    end do
    do jSym=1,nIrrep
      do iSO=1,nKVec
        nIJ1Max = max(nIJ1Max,nIJ1(iSym,jSym,iSO))
        nIJRMax = max(nIJRMax,nIJR(iSym,jSym,iSO))
      end do
    end do
  end do

  ! Allocate scratch memory for step 4 (Eq. 16). This is done
  ! in two steps.

  ! First step: Sum(j) X_lj C_ij^K = C_il^K; (l=valence basis)
  ! Second step: Sum(i) C_il^K X_ki = B_kl^K

  lCijK = nIJRMax*MxChVInShl
  lCilK = MxInShl*nIMax*MxChVInShl
  lCilK = max(lCilK,lCijK) ! it is used as scratch in pget
  call mma_allocate(CijK,lCilK,Label='CijK')
  if (lPSO) lCilK = max(lCilK,maxnAct) ! used as scratch
  call mma_allocate(CilK,lCilK,Label='CilK')
  lBklK = MxInShl**2*MxChVInShl
  call mma_allocate(BklK,lBklK,Label='BklK')

  if (iMp2prpt == 2) then
    lB_mp2 = mxChVInShl*nBas(0)*nBas(0)
    call mma_allocate(BMP2,lB_mp2,2,Label='BMP2')
  end if
  !                                                                    *
  !--------------------------------------------------------------------*
  !                                                                    *
  ! The C_ij^K vectors are stored in triangular form. We now
  ! change this to stricked rectangular/square form. Diagonal
  ! elements are rescaled. In case of symmetry this is only
  ! done for the blocks with isym=jSym, i.e. kSym=1

  kSym = 1
  nK = NumAuxVec(kSym)

  do iSO=1,nKVec
    if (lSA) cycle
    do iSym=1,nSym
      jSym = iSym
      !jSym = Mul(kSym,jSym)

      ! Read a whole block of C_ij^K

      iAdrC = iAdrCVec(kSym,iSym,iSO)
      call mma_allocate(CVec,nIJ1(iSym,jSym,iSO),nK,Label='CVec')
      call dDaFile(LuCVector(kSym,iSO),2,CVec,nIJ1(iSym,jSym,iSO)*nK,iAdrC)

      ni = nChOrb(iSym-1,iSO)
      call mma_allocate(CVec2,ni,ni,nK,Label='CVec2')
      !nj = ni

      do KAux=1,nK

        ij = 0
        do i=1,ni
          do j=1,i-1
            ij = ij+1
            CVec2(i,j,KAux) = CVec(ij,KAux)
            CVec2(j,i,KAux) = CVec(ij,KAux)
          end do
          ij = ij+1
          CVec2(i,i,KAux) = CVec(ij,KAux)*sqrt(Two)
        end do

      end do

      ! Write back to disk. Note that the file is prepared for
      ! rectangular/square storage so that one can safely write
      ! back the expanded set to disk without any overwrite problems.

      iAdrC = iAdrCVec(kSym,iSym,iSO)
      call dDaFile(LuCVector(kSym,iSO),1,CVec2,nIJR(iSym,jSym,iSO)*nK,iAdrC)

      call mma_deallocate(CVec2)
      call mma_deallocate(CVec)
    end do
  end do
  !                                                                    *
  !--------------------------------------------------------------------*
  !                                                                    *
  ! Stuff used in the prescreening!

  MumOrb = 0
  NumOrb = 0
  nPrev = 0
  do iSO=1,nKDens
    do jSym=1,nSym
      NumOrb = NumOrb+nChOrb(jSym-1,iSO)
      MumOrb = max(MumOrb,nChOrb(jSym-1,iSO))
    end do
    call mma_allocate(Ymnij(iSO)%A,NumOrb-nPrev,label='Ymnij%A')
    nPrev = NumOrb
  end do

  ! Scratch store the index of the MOs which finds the estimate
  ! according to Eq. 18 to be larger than the threshold.

  do i=1,nKDens
    Ymnij(i)%A(:) = 0
  end do

  ! Make a list for each shell-pair over the largest element
  ! SQRT(ABS( (mu,nu|mu,nu) ))

  nnSkal = nTri_Elem(nSkal_valence)
  call mma_allocate(SDG,nnSkal,Label='SDG')
  call get_maxDG(SDG,nnSkal,MxBasSh)

  ! Scratch for reduced lists of X_mi. Used in pget.

  nXki = MumOrb*MxBasSh*nSym
  call mma_allocate(Yij,nXki,2,nKDens,Label='Yij')

  ! Make a list the largest element X_mu,i for each valence shell
  ! and a fixed i. X_mu,i defined in Eq. 13.

  call mma_allocate(Xmi,MumOrb,nSkal_Valence,nIrrep,nKDens,Label='Xmi')

  do iSO=1,nKDens
    call get_mXOs(iSO,Xmi(:,:,:,iSO),MumOrb,nSkal_valence,nIrrep,nChOrb(0,iSO))
  end do

else

  nXki = 0
  NumOrb = 0
  MumOrb = 0

end if
!                                                                      *
!***********************************************************************
!                                                                      *
! For CASSCF process the active space contribution.

if (lPSO) then
  nBas_Aux(0) = nBas_Aux(0)-1
  call mma_allocate(Thhalf,maxnnP,Label='Thhalf')
  nThpkl = MxChVInShl*MxInShl**2
  call mma_allocate(Thpkl,nThpkl,Label='Thpkl')

  call contract_Zpk_Tpxy(Z_p_k,nZ_p_k,Txy,n_Txy,Thhalf,maxnnP,DMdiag,nG1,nnP,nBas_Aux,nADens,nAvec,nAct,nIrrep)

  call mma_deallocate(Thhalf)
  nBas_Aux(0) = nBas_Aux(0)+1
else
  nThpkl = 1
  call mma_allocate(Thpkl,nThpkl,Label='Thpkl')
end if
nCalAO = 0 !! for CASPT2
!                                                                      *
!***********************************************************************
!                                                                      *
! For a parallel implementation the iterations over shell-pairs
! are parallelized.

! If only Coulombic terms are to be processed use dynamic setup.
! Otherwise do batches exactly in the same order as Seward did the
! 2-center terms.

call Get_cArray('Relax Method',Method,8)
if (Method /= 'KS-DFT  ') then
  iOpt = 1
else
  call Get_cArray('DFT functional',KSDFT,80)
  ExFac = Get_ExFac(KSDFT)
  iOpt = 0
  if (ExFac /= Zero) iOpt = 1
end if
if (.not. Do_RI) iOpt = 0

if (iOpt == 1) then
  call qpg_iArray('LBList',Found,nSkal2_)
  if (Found) then
    call mma_allocate(LBList,nSkal2_,Label='LBList')
    call Get_iArray('LBList',LBList,nSkal2_)
  else
    call WarningMessage(2,'LBList not found!')
    call Abend()
  end if
else
  call mma_allocate(LBList,1,Label='LBList')
end if
call Init_Tsk2(id,nSkal2,iOpt,LBList)
call mma_deallocate(LBList)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
!do klS=1,nSkal2
do while (Rsv_Tsk2(id,klS))

  kS = Shij(1,klS)
  lS = Shij(2,klS)

  A_Int_kl = TMax_Valence(kS,lS)

  klS_ = iTri(kS,lS)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Prescreening stuff for exchange

  if (DoCholExch) then

    ! For the shell-pair, (kS,lS), pick up the largest element
    ! Sqrt(Abs(  (kappa,lambda|kappa,lambda) ))

    SDGmn = SDG(klS_)

    ! Loop over the MO basis, jb and ib and approximate Y_ij (Eq. 18)

    do iSO=1,nKVec
      FlipFlop = .true.
      iMOleft = iSO
      iMOright = iSO
      if (lSA) iMOright = iSO+2
      do
        nj = 0
        do jSym=1,nSym

          mj = 0
          do jb=1,nChOrb(jSym-1,iMOleft)
            Xjk = Xmi(jb,kS,jSym,iMOleft)
            Xjl = Xmi(jb,lS,jSym,iMOleft)

            jSym_s = jSym
            if ((ks /= ls) .or. (iMOright /= iMOleft)) jSym_s = 1
            loop1: do iSym=jsym_s,nSym
              NumOrb_i = nChOrb(iSym-1,iMOright)
              if ((iSym == jSym) .and. (ks == ls) .and. (iMOright == iMOleft)) NumOrb_i = jb

              do ib=NumOrb_i,1,-1
                Xik = Xmi(ib,kS,iSym,iMOright)
                Xil = Xmi(ib,lS,iSym,iMOright)

                ! Yij[mn] = (1+Pij) Xim * (mn|mn)^1/2 * Xjn

                PZmnij = (Xik*Xjl+Xil*Xjk)*SDGmn

                ! If larger than the threshold put j in the list and exit the loop.

                if (PZmnij >= xfk*ThrCom) then
                  ! orbital in the list
                  mj = mj+1
                  Ymnij(iMOleft)%A(mj+nj) = jb
                  exit loop1
                end if
              end do     ! ib
            end do loop1 ! iSym
          end do

          ! The first element is to keep track on how many elements that were saved.

          ! nOrbs in the list ==> dim(ij)=nOrbs**2
          nYmnij(jSym,iMOleft) = mj
          iOff_Ymnij(jSym,iMOleft) = nj
          nj = nj+mj

        end do ! jSym
        if (.not. (lSA .and. FlipFlop)) exit
        FlipFlop = .false.
        itmp = iMOleft
        iMOleft = iMOright
        iMOright = itmp
      end do
    end do

  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ReadBPT2 = .true.
  do ijS=1,nij
    iS = Shij2(1,ijS)
    jS = Shij2(2,ijS)

    if (Do_RI) then
      A_int = A_Int_kl*TMax_Auxiliary(jS-nSkal_Valence)
    else
      A_int = A_Int_kl*TMax_Valence(iS,jS)
    end if
    if (A_Int < CutInt) cycle

    call Eval_g1_ijkl(iS,jS,kS,lS,Temp,nGrad,A_Int)

  end do

end do
! End of big task loop
!                                                                      *
!***********************************************************************
!                                                                      *
RI_3C = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
!                         E P I L O G U E                              *
!                                                                      *
!***********************************************************************
!                                                                      *
! Deallocate scratch for exchange term

if (DoCholExch) then
  call mma_deallocate(Xmi)
  call mma_deallocate(Yij)
  call mma_deallocate(SDG)
  do i=1,nKDens
    call mma_deallocate(Ymnij(i)%A)
  end do
end if
call mma_deallocate(CijK,safe='*')
call mma_deallocate(CilK,safe='*')
call mma_deallocate(BklK,safe='*')
do i=1,nKDens
  call Deallocate_DT(CMOi(i))
end do
call mma_deallocate(MaxDens)

call mma_deallocate(BMP2,safe='*')
call mma_deallocate(Thpkl,safe='*')

call mma_deallocate(Sew_Scr,safe='*')
call Free_Tsk2(id)
call mma_deallocate(Shij2)
call mma_deallocate(Shij)
call mma_deallocate(TMax_Auxiliary)
call mma_deallocate(TMax_Valence)
!                                                                      *
!***********************************************************************
!                                                                      *
call Term_Ints()
!                                                                      *
!***********************************************************************
!                                                                      *
call Free_iSD()
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Drvg1_3Center_RI
