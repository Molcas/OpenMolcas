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
! Copyright (C) 2007, Francesco Aquilante                              *
!               2011, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine CHO_GET_GRAD(irc,nDen,DLT,DLT2,MSQ,Txy,nTxy,ipTxy,DoExchange,lSA,nChOrb_,AOrb,nAorb,DoCAS,Estimate,Update,V_k,nV_k,U_k, &
                        Z_p_k,nZ_p_k,nnP,npos)
!***********************************************************************
!  Author : F. Aquilante (visiting F. Illas group in Barcelona, Spain, *
!                                                    March-April 2007) *
!                                                                      *
!  Purpose:                                                            *
!         Computation of the relevant quantities for RI                *
!         (and Cholesky) gradient code                                 *
!                                                                      *
!         Coulomb term :  V_k = sum_gd  D_gd L_gd_k                    *
!                                                                      *
!         MP2 Coulomb term : U_k = sum_gd D(MP2)_gd L_gd_k             *
!                                                                      *
!         Active term :   Z_p_k = sum_xy T(xy,p) L_xy_k                *
!                                                                      *
!         Inact. Exchange term: the quantity returned on disk is       *
!                                                                      *
!                     L_ij_k = sum_gd L_gd_k C_gi C_dj                 *
!                                                                      *
!                                                                      *
!  Input:                                                              *
!                                                                      *
!         nDen : is equal to 2 iff Spin Unrestricted                   *
!                4 for SA-CASSCF, otherwise nDen=1                     *
!                                                                      *
!         DLT : the LT-packed and symm. blocked one-body Dmat.         *
!                 For spin unrestricted, Dmat = Dalpha + Dbeta         *
!                                                                      *
!         DLT2: pointer to the LT-packed and symm. blocked             *
!                 one body MP2 Dmat.                                   *
!                                                                      *
!         MSQ : Cholesky MOs stored as C(a,i) symm. blocked            *
!                 with blocks of dim. (nBas,nBas). These are           *
!                 obtained from CD of the 1-particle DMAT.             *
!                 (Two pointers iff alpha and beta spinorbitals)       *
!                                                                      *
!         ipTxy : array (8,8) of pointers to the symm. blocks          *
!                 of the Cholesky decomposed MO-basis (symmetrized)    *
!                 2-body density matrix                                *
!                 T(xy,p) : is stored by compound symmetry JSYM        *
!                            the indices {xy} are stored as PACKED     *
!                            (sym x <= sym y)                          *
!                                                                      *
!         DoExchange : logical to activate exchange grad. components   *
!                                                                      *
!         nChOrb_ : array of nr. of Cholesky orbitals in each irrep    *
!                                                                      *
!         nAorb : array with # of Active orbitals in each irrep        *
!                 (The same orbital basis                              *
!                 in which the 2-body Dmat is expressed)               *
!                                                                      *
!         DoCAS : logical to activate CASSCF grad. components          *
!                                                                      *
!         nScreen : See e.g. LK-screening docum. in SCF                *
!                   or CASSCF read-input routines. Default = 10        *
!                                                                      *
!         dmpK : damping for the LK-screening threshold. Def: 1.0      *
!                                                                      *
!         Estimate : logical for LK-screening. Default: .false.        *
!                                                                      *
!         Update : logical for LK-screening. Default: .true.           *
!                                                                      *
!         nnP : array of # of Cholesky vectors for the dec 2-body      *
!               density matrix in each compound symmetry               *
!                                                                      *
!                                                                      *
!  Output:                                                             *
!         irc : return code                                            *
!                                                                      *
!         V_k : array Real*8 for the Coulomb interm. Size=NumCho(1)    *
!                                                                      *
!         U_k : array Real*8 for the mp2 Coulomb interm. Size=NumCho(1)*
!                                                                      *
!         Z_p_k : array Real*8 for the active grad. components.        *
!                  Must be zeroed by the calling routine. Stored       *
!                  according to jSym and blocked after symm. blocks    *
!                  of the active orbitals (square storage).            *
!                                                                      *
!  Modifications:                                                      *
!    August 24, 2011, Thomas Bondo Pedersen:                           *
!       Allow zero vectors on a node.                                  *
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use Symmetry_Info, only: Mul
use Cholesky, only: iiBstR, IndRed, InfVec, MaxRed, nBas, nBasSh, nDimRS, nnBstR, nnBstRSh, nnBstRT, nnShl, nnShl_tot, nShell, &
                    nSym, NumCho, NumChT, timings, ThrCom
use Data_Structures, only: DSBA_Type, NDSBA_Type, SBA_Type, V2
use Cholesky_Structures, only: Allocate_DT, Deallocate_DT, L_Full_Type, Lab_Type
use RI_glob, only: CMOi, dmpK, iBDsh, iMP2prpt, nAdens, nIJ1, nIJR, nJdens, nKdens, nKvec, nScreen
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: nDen, nTxy, ipTxy(8,8,2), nChOrb_(8,5), nAorb(8), nV_k, nZ_p_k, nnP(8), npos(8,3)
type(DSBA_Type), intent(in) :: DLT(5), DLT2, MSQ(nDen), AOrb(*)
real(kind=wp), intent(in) :: Txy(nTxy)
logical(kind=iwp), intent(in) :: DoExchange, lSA, DoCAS, Estimate, Update
real(kind=wp), intent(_OUT_) :: V_k(nV_k,*), U_k(*)
real(kind=wp), intent(inout) :: Z_p_k(nZ_p_k,*)
#include "Molcas.fh"
#include "print.fh"
!#define _CD_TIMING_
#ifdef _CD_TIMING_
#include "temptime.fh"
#endif
integer(kind=iwp) :: i, iAdr, iaSh, iAvec, iBatch, ibcount, ibs, ibs_a, ibSh, iE, ij, ik, iLoc, iml, iMO1, iMO2, iMOleft, &
                     iMOright, ioff, iOffShb, iOffZp, iPrint, ipZp, ir, ired1, IREDC, iRout, iS, iSeed, ish, iShp, iSSa, iStart, &
                     iSwap, iSwap_lxy, ISYM, iSym1, iSym2, iSyma, iSymb, iSymv, iSymx, iSymy, it, itk, iTmp, iTxy, IVEC2, iVrs, j, &
                     jDen, jGam, jK, jK_a, jml, jmlmax, JNUM, JRED, JRED1, JRED2, jrs, jSym, jvc, JVEC, k, kMOs, kOff(8,5), krs, &
                     kscreen, kSym, l, l1, LFMAX, LKsh, LREAD, lSym, LuRVec(8,3), LWORK, MaxB, MaxRedT, mDen, mrs, MUSED, n1, n2, &
                     nAv, NAw, nBatch, nI2t, nIJMax, Nik, nInd, nIt(5), nItmx, nL_Full(2), nLab(2), nLaq, nLik, nLxy, nLxy0, nMat, &
                     nMOs, nnA(8,8), npos2, nQo, nQoT, nRik, nRS, NumCV, numSh, NUMV, NumVT, nVec, nVrs
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: myJRED1, NNBSTMX, ntv0
#endif
real(kind=wp) :: SkSh, tau, tcasg(2), TCC1, TCC2, tcoul(2), TCR1, TCR2, TCS1, TCS2, TCT1, TCT2, TCX1, TCX2, temp, thrv, tmotr(2), &
                 tmotr2(2), TOTCPU, TOTCPU1, TOTCPU2, TOTWALL, TOTWALL1, TOTWALL2, tread(2), tscrn(2), TWC1, TWC2, TWR1, TWR2, &
                 TWS1, TWS2, TWT1, TWT2, TWX1, TWX2, xtau, xTmp, YMax, YshMax
logical(kind=iwp) :: add, BatchWarn, DoScreen
character(len=50) :: CFmt
character(len=6) :: Fname
character :: mode
type(L_Full_Type) :: L_Full
type(Lab_Type) :: Lab
type(NDSBA_Type) :: DiaH
type(SBA_Type) :: Laq(1), Lxy
type(V2) :: SumClk(5)
integer(kind=iwp), allocatable :: Indik(:,:), Indx(:,:), iShp_rs(:), kOffSh(:,:)
real(kind=wp), allocatable :: AbsC(:), Diag(:), Drs(:,:), Drs2(:,:), Lrs(:,:), MLk(:), SvShp(:,:), Ylk(:,:)
#ifdef _MOLCAS_MPP_
real(kind=wp), allocatable :: DiagJ(:)
#endif
real(kind=wp), allocatable, target :: Aux(:), Aux0(:), Yik(:)
real(kind=wp), pointer :: Lik(:,:), pYik(:,:), Rik(:)
logical(kind=iwp), parameter :: DoRead = .false.
character(len=*), parameter :: SECNAM = 'CHO_GET_GRAD'
integer(kind=iwp), external :: IsFreeUnit
real(kind=wp), external :: ddot_

!                                                                      *
!***********************************************************************
!                                                                      *
!     General Initialization                                           *
!                                                                      *
!***********************************************************************
!                                                                      *

iRout = 9
iPrint = nPrint(iRout)

call CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

! 1 --> CPU   2 --> Wall
tread(:) = zero  !time read vectors
tcoul(:) = zero  !time for computing V_k
tcasg(:) = zero  !time for computing Z_p_k
tmotr(:) = zero  !time for the MO transf of vectors
tmotr2(:) = zero  !time for the 2nd MO transf of vectors
tscrn(:) = zero  !time for screening overhead

IREDC = -1  ! unknown reduced set in core

BatchWarn = .true.
nInd = 0

call set_nnA(nSym,nAorb,nnA)

! Various offsets

MaxB = nBas(1)
do ISYM=2,NSYM
  MaxB = max(MaxB,nBas(iSym))
end do

nI2t = 0
nItmx = 0
nIt(:) = 0
do jDen=nDen,1,-1
  kOff(1,jDen) = 0
  nIt(jDen) = nChOrb_(1,jDen)
  do i=2,nSym
    kOff(i,jDen) = nIt(jDen)
    nIt(jDen) = nIt(jDen)+nChOrb_(i,jDen)
  end do
  nI2t = nI2t+nIt(jDen)
  nItmx = max(nItmx,nIt(jDen))
end do

! Initialize pointers to avoid compiler warnings

thrv = Zero
xtau = Zero

! Construct iBDsh for later use

call mma_allocate(iBDsh,nShell*nSym,label='iBDsh')
do iSyma=1,nSym
  LKsh = 0
  do iaSh=1,nShell
    iSSa = nShell*(iSyma-1)+iaSh
    iBDsh(iSSa) = LKsh
    LKsh = LKsh+nBasSh(iSyma,iaSh)
  end do
end do

! iShp_rs
call mma_allocate(iShp_rs,nnShl_tot,Label='iShp_rs')

!***********************************************************************
!                                                                      *
!     Initialize a few things for ij-screening //Jonas B               *
!                                                                      *
!***********************************************************************
if (DoExchange) then

  ! Define the screening thresholds

  call Get_dScalar('Cholesky Threshold',ThrCom)

  tau = (ThrCom/real(max(1,nItmx),kind=wp))*dmpK

  MaxRedT = MaxRed
  call GAIGOP_SCAL(MaxRedT,'+')

  if (Estimate) tau = tau/real(MaxRedT,kind=wp)
  xtau = sqrt(tau)

  NumVT = NumChT
  call GAIGOP_SCAL(NumVT,'+')
  ! Vector MO transformation screening thresholds
  thrv = (sqrt(ThrCom/real(max(1,nItmx)*NumVT,kind=wp)))*dmpK

# ifdef _MOLCAS_MPP_
  if (Is_Real_Par() .and. Update) then
    NNBSTMX = 0
    do i=1,nSym
      NNBSTMX = max(NNBSTMX,NNBSTR(i,1))
    end do
    call mma_allocate(DiagJ,NNBSTMX,Label='DiagJ')
    DiagJ(:) = Zero
  end if
# endif

  ! Read the diagonal integrals (stored as 1st red set)

  call mma_allocate(DIAG,NNBSTRT(1),Label='Diag')
  if (Update) call CHO_IODIAG(DIAG,2) ! 2 means "read"


  ! Allocate memory

  ! sqrt(D(a,b)) stored in full (squared) dim
  call Allocate_DT(DiaH,nBas,nBas,nSym)
  DiaH%A0(:) = Zero

  call mma_allocate(AbsC,MaxB,Label='AbsC')

  call mma_allocate(Ylk,MaxB,nItmx,Label='Ylk')

  call mma_allocate(Yik,nItmx**2,Label='Yik') ! Yi[k] vectors

  ! used to be nShell*something
  ! ML[k] lists of largest elements in significant shells
  call mma_allocate(MLk,nShell,Label='MLk')

  ! list of S:= sum_l abs(C(l)[k])
  call mma_allocate(Aux0,nShell*nI2t,Label='Aux0')
  iE = 0
  do i=1,nDen
    iS = iE+1
    iE = iE+nShell*nIt(i)
    SumClk(i)%A(1:nShell,1:nIt(i)) => Aux0(iS:iE)
  end do

  ! Indx and Indik must be stored for each density, symmetry, etc.
  ! in case of a batched procedure
  do jDen=1,nKvec
    do kSym=1,nSym
      nInd = nInd+nChOrb_(kSym,jDen)
    end do
  end do

  ! Index array
  call mma_allocate(Indx,[0,nShell],[1,nInd],Label='Indx')

  ! Yi[k] Index array
  call mma_allocate(Indik,(nItmx+1)*nItmx+1,nInd,Label='Indik')

  ! kOffSh
  call mma_allocate(kOffSh,nShell,nSym,Label='kOffSh')

  ! shell-pair Frobenius norm of the vectors
  call mma_allocate(SvShp,nnShl,2,Label='SvShp')

  ! Jonas - June 2010:
  ! allocate memory for rearranged CMO-matrix

  do i=1,nDen
    call Allocate_DT(CMOi(i),nChOrb_(:,i),nBas,nSym)
  end do

  nQoT = 0

  ! Compute Shell Offsets ( MOs and transformed vectors)

  do iSyma=1,nSym
    LKsh = 0
    do iaSh=1,nShell    ! kOffSh(iSh,iSym)

      kOffSh(iaSh,iSyma) = LKsh

      LKsh = LKsh+nBasSh(iSyma,iaSh)
    end do
  end do

  ! Determine S:= sum_l C(l)[k]^2  in each shell of C(a,k)

  do jDen=1,nDen
    do kSym=1,nSym

      do jK=1,nChOrb_(kSym,jDen)
        jK_a = jK+kOff(kSym,jDen)

        do iaSh=1,nShell

          iS = kOffSh(iaSh,kSym)+1
          iE = kOffSh(iaSh,kSym)+nBasSh(kSym,iaSh)

          SKSh = Zero
          do ik=iS,iE
            SKsh = SKsh+MSQ(jDen)%SB(kSym)%A2(ik,jK)**2
          end do

          SumClk(jDen)%A(iaSh,jK_a) = SkSh

        end do
      end do
    end do
  end do

  ! Reorder CMO-matrix, Needed to construct B-matrix for exchange
  ! Jonas - June 2010

  do jDen=1,nKdens
    do kSym=1,nSym

      ! If the orbitals come from eigenvalue decomposition, change sign

      if (lSA .and. (jDen >= 3)) then
        npos2 = npos(ksym,jDen-2)
        do jK=1,nPos2
          do jGam=1,nBas(kSym)
            CMOi(jDen)%SB(kSym)%A2(jK,jGam) = MSQ(jDen)%SB(kSym)%A2(jGam,jK)
          end do
        end do
        do jK=npos2+1,nChOrb_(kSym,jDen)
          do jGam=1,nBas(kSym)
            CMOi(jDen)%SB(kSym)%A2(jK,jGam) = -MSQ(jDen)%SB(kSym)%A2(jGam,jK)
          end do
        end do
      else

        do jK=1,nChOrb_(kSym,jDen)
          do jGam=1,nBas(kSym)
            CMOi(jDen)%SB(kSym)%A2(jK,jGam) = MSQ(jDen)%SB(kSym)%A2(jGam,jK)
          end do
        end do
      end if
    end do
  end do
end if

! Mapping shell pairs from the full to the reduced set

call Mk_iShp_rs(iShp_rs,nShell)

!                                                                      *
!***********************************************************************
!                                                                      *
!     BIG LOOP OVER VECTORS SYMMETRY                                   *
!                                                                      *
!***********************************************************************
!                                                                      *
do jSym=1,nSym
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  NumCV = NumCho(jSym)
  call GAIGOP_SCAL(NumCV,'max')
  if (NumCV < 1) cycle

  ! offsets for active term

  iOffZp = 0
  do j=1,jSym-1
    iOffZp = iOffZp+nnP(j)*NumCho(j)
  end do

  ! Open some files to store exchange auxiliary vectors

  if (DoExchange) then
    iSeed = 7
    do i=1,nSym
      k = Mul(jSym,i)
      LuRVec(i,1) = IsFreeUnit(iSeed)
      write(Fname,'(A4,I1,I1)') 'CHTA',i,k
      call DANAME_MF_WA(LuRVec(i,1),Fname)
      iSeed = iSeed+1
      if (nKvec >= 2) then
        LuRVec(i,2) = IsFreeUnit(iSeed)
        write(Fname,'(A4,I1,I1)') 'CHTB',i,k
        call DANAME_MF_WA(LuRVec(i,2),Fname)
        iSeed = iSeed+1
      end if
    end do
  end if
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  !        M E M O R Y   M A N A G E M E N T   S E C T I O N           *
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !
  ! For one Cholesky vector, JNUM=1, compute the amount of memory
  ! needed for the various vectors.

  JNUM = 1

  ! L_Full
  call Allocate_DT(L_Full,nShell,iShp_rs,JNUM,JSYM,nSym,Memory=nL_Full)
  ! Lab
  mDen = 1
  call Allocate_DT(Lab,JNUM,nBasSh,nBas,nShell,nSym,mDen,Memory=nLab)
  if (DoCas) then
    iSwap = 0  ! Lvb,J are returned
    call Allocate_DT(Laq(1),nAorb,nBas,JNUM,JSYM,nSym,iSwap,Memory=nLaq)
    nLxy = 0
    do iMO1=1,nAdens
      iSwap_lxy = 5
      if (iMO1 == 2) iSwap_lxy = 6
      call Allocate_DT(Lxy,nAorb,nAorb,JNUM,JSYM,nSym,iSwap_lxy,Memory=nLxy0)
      nLxy = max(nLxy,nLxy0)
    end do
  else
    nLaq = 0
    nLxy = 0
  end if

  ! compute memory needed to store at least 1 vector of JSYM
  ! and do all the subsequent calculations

  nLik = 0
  nRik = 0
  do l=1,nSym
    k = Mul(l,JSYM)
    do jDen=1,nDen
      nRik = max(nRik,nChOrb_(l,jDen)*nChOrb_(k,jDen))
      if (nChOrb_(k,jDen) > 0) then
        nLik = max(nLik,nChOrb_(l,jDen))
      end if
    end do
  end do

  ! re-use memory for the active vec
  LFMAX = max(nLaq+nLxy,nL_Full(1)+nRik+nLik+nLab(1))
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
!                                                                      *

  iLoc = 3 ! use scratch location in reduced index arrays

  if (NumCho(jSym) < 1) then
    JRED1 = 1
    JRED2 = 1
  else
    JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
    JRED2 = InfVec(NumCho(jSym),2,jSym) ! red set of the last vec
  end if
# ifdef _MOLCAS_MPP_
  myJRED1 = JRED1 ! first red set present on this node
  ntv0 = 0
# endif

  ! entire red sets range for parallel run
  call GAIGOP_SCAL(JRED1,'min')
  call GAIGOP_SCAL(JRED2,'max')

  ! MGD does it need to be so?

  DoScreen = .true.
  kscreen = 1
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  do JRED=JRED1,JRED2
    !                                                                  *
    !*******************************************************************
    !                                                                  *

    if (NumCho(jSym) < 1) then
      iVrs = 0
      nVrs = 0
    else
      call Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)
    end if

    if (nVrs /= 0) then  ! some vector in that (jred,jsym)

      if (nVrs < 0) then
        write(u6,*) SECNAM//': Cho_X_nVecRS returned nVrs<0. STOP!'
        call Abend()
      end if

      ! set index arrays at iLoc
      call Cho_X_SetRed(irc,iLoc,JRED)
      if (irc /= 0) then
        write(u6,*) SECNAM,': cho_X_setred non-zero return code. rc= ',irc
        call Abend()
      end if

      IREDC = JRED

      nRS = nDimRS(JSYM,JRED)

      if (JSYM == 1) then
        call mma_allocate(Drs,nRS,nJdens,Label='Drs')
        Drs(:,:) = Zero
        if (iMp2prpt == 2) then
          call mma_allocate(Drs2,nRS,1,Label='Drs2')
        end if
      end if

      call mma_maxDBLE(LWORK)

      nVec = min((LWORK-nL_Full(2)-nLab(2))/(nRS+LFMAX),nVrs)

      if (nVec < 1) then
        write(u6,*) SECNAM//': Insufficient memory for batch'
        write(u6,*) ' LWORK= ',LWORK
        write(u6,*) ' min. mem. need= ',nRS+LFMAX
        write(u6,*) ' jsym= ',jsym
        write(u6,*) ' nRS = ',nRS
        write(u6,*) ' LFMAX = ',LFMAX
        write(u6,*)
        write(u6,*) ' nL_Full = ',nL_Full
        write(u6,*) ' nRik = ',nRik
        write(u6,*) ' nLik = ',nLik
        write(u6,*) ' nLab = ',nLab
        write(u6,*)
        write(u6,*) ' nLaq = ',nLaq
        write(u6,*) ' nLxy = ',nLxy
        irc = 33
        call Abend()
        nBatch = -9999  ! dummy assignment
      end if

      !                                                                *
      !*****************************************************************
      !                                                                *
      LREAD = nRS*nVec

      call mma_allocate(Lrs,nRS,nVec,Label='Lrs')

      if (JSYM == 1) then
        ! Transform the densities to reduced set storage
        add = .false.
        nMat = 1
        do jDen=1,nJdens
          call swap_full2rs(irc,iLoc,nRS,nMat,JSYM,DLT(jDen),Drs(:,jDen),add)
        end do
        if (iMp2prpt == 2) then
          call swap_full2rs(irc,iLoc,nRS,nMat,JSYM,[DLT2],Drs2(:,1),add)
        end if
      end if

      ! BATCH over the vectors

      nBatch = (nVrs-1)/nVec+1

      if (BatchWarn .and. (nBatch > 1)) then
        if (iPrint >= 6) then
          write(u6,'(A)') repeat('-',60)
          write(u6,*) ' Batch procedure used. Increase memory if possible!'
          write(u6,'(A)') repeat('-',60)
          write(u6,*)
          call XFlush(u6)
        end if
        BatchWarn = .false.
      end if

      !                                                                *
      !*****************************************************************
      !                                                                *
      do iBatch=1,nBatch
        !                                                              *
        !***************************************************************
        !                                                              *
        if (iBatch == nBatch) then
          JNUM = nVrs-nVec*(nBatch-1)
        else
          JNUM = nVec
        end if

        JVEC = nVec*(iBatch-1)+iVrs
        IVEC2 = JVEC-1+JNUM

        call CWTIME(TCR1,TWR1)

        call CHO_VECRD(Lrs,LREAD,JVEC,IVEC2,JSYM,NUMV,IREDC,MUSED)

        if ((NUMV <= 0) .or. (NUMV /= JNUM)) then
          irc = 77
          return
        end if

        call CWTIME(TCR2,TWR2)
        tread(1) = tread(1)+(TCR2-TCR1)
        tread(2) = tread(2)+(TWR2-TWR1)

        !***************************************************************
        !***************************************************************
        !**                                                           **
        !**      Coulomb term                                         **
        !**      V{#J} = sum_ab  L(ab,{#J}) * D(ab)                   **
        !**                                                           **
        !***************************************************************
        !***************************************************************
        if (JSYM == 1) then

          call CWTIME(TCC1,TWC1)

          ! Inactive Coulomb term

          do jden=1,nJdens
            call DGEMV_('T',nRS,JNUM,One,Lrs,nRS,Drs(1,jden),1,zero,V_k(jVec,jDen),1)
          end do

          ! MP2 Coulomb term

          if (iMp2prpt == 2) then
            call DGEMV_('T',nRS,JNUM,One,Lrs,nRS,Drs2(:,1),1,zero,U_k(jVec),1)
          end if

          call CWTIME(TCC2,TWC2)
          tcoul(1) = tcoul(1)+(TCC2-TCC1)
          tcoul(2) = tcoul(2)+(TWC2-TWC1)
        end if
        !***************************************************************
        !***************************************************************
        !**                                                           **
        !**      E X C H A N G E    T E R M                           **
        !**                                                           **
        !***************************************************************
        !***************************************************************

        if (DoExchange) then

          call CWTIME(TCS1,TWS1)
          !*************************************************************
          !                                                            *
          ! 1) Screening                                               *
          !                                                            *
          !    Select only important  ij pairs                         *
          !    For this, one computes the quantity                     *
          !       Yik = sum_mu_nu (mu nu | mu nu)^1/2 X_mu_i X_nu_k    *
          !    with (mu nu | mu nu) = sum_J (L_mu_nu,J)^2              *
          !                                                            *
          !                                                            *
          !    a) Estimate the diagonals:                              *
          !         D(mu,nu) = sum_J (L_mu_nu,J)^2                     *
          !                                                            *
          !*************************************************************
          if (Estimate) then

            Diag(iiBstR(jSym,1)+1:iiBstR(jSym,1)+NNBSTR(jSym,1)) = Zero

            do krs=1,nRS

              mrs = iiBstR(JSYM,iLoc)+krs
              jrs = IndRed(mrs,iLoc) ! address in 1st red set

              do jvc=1,JNUM

                Diag(jrs) = Diag(jrs)+Lrs(krs,jvc)**2

              end do

            end do

          end if

          call CWTIME(TCS2,TWS2)
          tscrn(1) = tscrn(1)+(TCS2-TCS1)
          tscrn(2) = tscrn(2)+(TWS2-TWS1)
          !                                                            *
          !*************************************************************
          !*************************************************************
          !*************************************************************
          !                                                            *
          call Allocate_DT(L_Full,nShell,iShp_rs,JNUM,JSYM,nSym)
          call mma_allocate(Aux,(nRik+nLik)*nVec,Label='Aux')
          call Allocate_DT(Lab,JNUM,nBasSh,nBas,nShell,nSym,mDen)

          call CWTIME(TCX1,TWX1)

          ! Reorder vectors to Full-dimensions
          !
          ! Vectors are returned in the storage LaJ,b with the restriction:
          !    Sym(a) >= Sym(b)
          ! and blocked in shell pairs

          call CHO_getShFull(Lrs,lread,JNUM,JSYM,IREDC,L_Full,SvShp,nnShl,iShp_rs,nnShl_tot)

          call CWTIME(TCX2,TWX2)
          tmotr(1) = tmotr(1)+(TCX2-TCX1)
          tmotr(2) = tmotr(2)+(TWX2-TWX1)

          !*************************************************************
          !                                                            *
          ! 1) Screening                                               *
          !                                                            *
          !    b) DH(mu,nu)=sqrt(D(mu,nu))                             *
          !       Only the symmetry blocks with compound symmetry JSYM *
          !       are computed                                         *
          !                                                            *
          !*************************************************************
          if (DoScreen) then

            call CWTIME(TCS1,TWS1)

            ired1 = 1 ! location of the 1st red set
            call swap_tosqrt(irc,ired1,NNBSTRT(1),JSYM,DIAH,DIAG)

            call CWTIME(TCS2,TWS2)
            tscrn(1) = tscrn(1)+(TCS2-TCS1)
            tscrn(2) = tscrn(2)+(TWS2-TWS1)

          end if

          !*************************************************************
          !                                                            *
          ! 1) Screening                                               *
          !                                                            *
          !    c) 1st MO transformation of DH(mu,nu)                   *
          !          Y(mu)[k] = sum_nu  DH(mu,nu) * |C(nu)[k]|         *
          !                                                            *
          !*************************************************************

          nInd = 1
          do jDen=1,nKvec

            ! Choose which MO sets on each side

            iMOleft = jDen
            iMOright = jDen

            n1 = nIt(iMOright)
            n2 = nItMx

            pYik(1:n1,1:n2) => Yik(1:n1*n2)

            if (DoCAS .and. lSA) iMOright = jDen+2

            do kSym=1,nSym

              lSym = Mul(JSYM,kSym)
              Nik = nChOrb_(kSym,iMOleft)*nChOrb_(lSym,iMOright)
              nIJR(kSym,lSym,jDen) = Nik
              nIJ1(kSym,lSym,jDen) = Nik
              if ((JSYM == 1) .and. (iMOleft == iMOright)) Nik = nTri_Elem(nChOrb_(kSym,iMOleft))
              nIJ1(kSym,lSym,jDen) = Nik

              if (Nik == 0) cycle

              iS = 1
              iE = nChOrb_(lSym,iMOright)*JNUM

              Lik(1:JNUM,1:nChOrb_(lSym,iMOright)) => Aux(iS:iE)

              iS = iE+1
              iE = iE+Nik*JNUM

              Rik(1:Nik*JNUM) => Aux(iS:iE)

              Rik(:) = Zero

              do jK=1,nChOrb_(kSym,iMOleft)
                jK_a = jK+kOff(kSym,iMOleft)

                Lik(:,:) = Zero
                Lab%A0(1:nBas(lSym)*JNUM) = Zero

                if (DoScreen .and. (iBatch == 1)) then
                  call CWTIME(TCS1,TWS1)
                  !-----------------------------------------------------
                  ! Setup the screening
                  !-----------------------------------------------------

                  AbsC(1:nBas(kSym)) = abs(MSQ(iMOleft)%SB(kSym)%A2(:,jK))

                  if (lSym >= kSym) then

                    mode = 'N'
                    n1 = nBas(lSym)
                    n2 = nBas(kSym)

                  else ! lSym<kSym

                    mode = 'T'
                    n1 = nBas(kSym)
                    n2 = nBas(lSym)

                  end if

                  if (n1 > 0) call DGEMV_(Mode,n1,n2,ONE,DiaH%SB(lSym,kSym)%A2,n1,AbsC,1,ZERO,Ylk(1,jK_a),1)

                  !*****************************************************
                  !                                                    *
                  ! 1) Screening                                       *
                  !                                                    *
                  !    d) 2nd MO transformation of DH(mu,nu)           *
                  !          Y(i)[k] = sum_mu  |C(mu)[i]| * Y(mu)[k]   *
                  !                                                    *
                  !*****************************************************

                  if ((kSym /= lSym) .or. (iMOleft /= iMOright)) then
                    iStart = 1
                  else
                    iStart = jK
                  end if

                  nQo = 0
                  do i=iStart,nChOrb_(lSym,iMOright)

                    AbsC(1:nBas(lSym)) = abs(MSQ(iMOright)%SB(lSym)%A2(:,i))

                    pYik(i,jK_a) = ddot_(nBas(lSym),AbsC,1,Ylk,1)

                    if (pYik(i,jK_a) >= xtau) then
                      nQo = nQo+1
                      if ((iBatch == 1) .and. (JRED == 1)) then
                        nQoT = nQoT+1
                        if ((lSym == kSym) .and. (i /= jK) .and. (iMOright == iMOleft)) nQoT = nQoT+1
                      end if
                      Indik(1+nQo,nInd) = i
                    end if

                  end do
                  Indik(1,nInd) = nQo
                  !*****************************************************
                  !                                                    *
                  ! 1) Screening                                       *
                  !                                                    *
                  !    e) List the shells present in Y(l)[k] by the    *
                  !       largest element and sort the list            *
                  !                                                    *
                  !*****************************************************

                  do ish=1,nShell
                    YshMax = zero
                    do ibs=1,nBasSh(lSym,ish)
                      ibs_a = koffSh(ish,lSym)+ibs
                      YshMax = max(YshMax,Ylk(ibs_a,1))
                    end do
                    MLk(ish) = YshMax
                  end do

                  do ish=1,nShell
                    Indx(ish,nInd) = ish
                  end do

                  !*****************************************************
                  !                                                    *
                  ! 1) Screening                                       *
                  !                                                    *
                  !    f) Screening                                    *
                  !                                                    *
                  ! Here we use a non-exact bound for the exchange     *
                  ! matrix to achieve linear scaling. The positive     *
                  ! definiteness of the exchange matrix combined with  *
                  ! the structure of the density matrix makes this     *
                  ! bound acceptable and likely to be almost exact for *
                  ! what concerns the exchange energy                  *
                  !                                                    *
                  ! The exact bounds (quadratic scaling of the MO      *
                  ! transformation) would be                           *
                  !    If (MLk(jml)*MLk(1) >= tau) then                *
                  !                                                    *
                  !*****************************************************

                  numSh = 0  ! # of significant shells
                  jml = 1
                  do while (jml <= nShell)

                    YMax = MLk(jml)
                    jmlmax = jml

                    do iml=jml+1,nShell  ! get the max
                      if (MLk(iml) > YMax) then
                        YMax = MLk(iml)
                        jmlmax = iml
                      end if
                    end do

                    if (jmlmax /= jml) then  ! swap positions
                      xTmp = MLk(jml)
                      iTmp = Indx(jml,nInd)
                      MLk(jml) = YMax
                      Indx(jml,nInd) = Indx(jmlmax,nInd)
                      MLk(jmlmax) = xTmp
                      Indx(jmlmax,nInd) = iTmp
                    end if

                    if (MLk(jml) >= xtau) then
                      numSh = numSh+1
                    else
                      jml = nShell  ! exit the loop
                    end if

                    jml = jml+1

                  end do

                  Indx(0,nInd) = numSh

                  call CWTIME(TCS2,TWS2)
                  tscrn(1) = tscrn(1)+(TCS2-TCS1)
                  tscrn(2) = tscrn(2)+(TWS2-TWS1)
                  !-----------------------------------------------------
                end if  ! Screening setup

                call CWTIME(TCT1,TWT1)

                !*******************************************************
                !                                                      *
                !             E X C H A N G E    T E R M               *
                !                                                      *
                ! 2) MO transformation                                 *
                !    a) 1st half transformation                        *
                !                                                      *
                !    Transform vectors for shells in the list ML[k]    *
                !                                                      *
                !    Screening based on the Frobenius norm:            *
                !       sqrt(sum_ij  A(i,j)^2)                         *
                !       || La,J[k] ||  <=  || Lab,J || * || Cb[k] ||   *
                !                                                      *
                !*******************************************************

                do iSh=1,Indx(0,nInd)

                  iaSh = Indx(iSh,nInd)

                  Lab%Keep(iaSh,1) = .true.

                  ibcount = 0

                  do ibSh=1,nShell

                    iOffShb = kOffSh(ibSh,kSym)

                    iShp = iTri(iaSh,ibSh)

                    if (iShp_rs(iShp) <= 0) cycle

                    if ((nnBstRSh(JSym,iShp_rs(iShp),iLoc)*nBasSh(lSym,iaSh)*nBasSh(kSym,ibSh) > 0) .and. &
                        (sqrt(abs(SumClk(iMOleft)%A(ibSh,jK_a)*SvShp(iShp_rs(iShp),1))) >= thrv)) then

                      ibcount = ibcount+1

                      if (lSym >= kSym) then

                        l1 = 1
                        if (iaSh < ibSh) l1 = 2

                        ! LaJ,[k] = sum_b  L(aJ,b) * C(b)[k]
                        !------------------------------------

                        Mode = 'N'
                        n1 = nBasSh(lSym,iaSh)*JNUM
                        n2 = nBasSh(kSym,ibSh)

                        call DGEMV_(Mode,n1,n2,One,L_Full%SPB(lSym,iShp_rs(iShp),l1)%A21,n1, &
                                    MSQ(iMOleft)%SB(kSym)%A2(iOffShb+1:,jK),1,ONE,Lab%SB(iaSh,lSym,1)%A,1)

                      else   ! lSym < kSym

                        l1 = 1
                        if (ibSh < iaSh) l1 = 2

                        ! LJa,[k] = sum_b  L(b,Ja) * C(b)[k]
                        !------------------------------------

                        Mode = 'T'
                        n1 = nBasSh(kSym,ibSh)
                        n2 = JNUM*nBasSh(lSym,iaSh)

                        call DGEMV_(Mode,n1,n2,One,L_Full%SPB(kSym,iShp_rs(iShp),l1)%A12,n1, &
                                    MSQ(iMOleft)%SB(kSym)%A2(iOffShb+1:,jK),1,ONE,Lab%SB(iaSh,lSym,1)%A,1)

                      end if

                    end if

                  end do

                  ! The following re-assignement is used later on to check if the
                  ! iaSh vector LaJ[k] can be neglected because identically zero

                  if (ibcount == 0) Lab%Keep(iaSh,1) = .false.

                end do

                !*******************************************************
                !                                                      *
                ! 2) MO transformation                                 *
                !    b) 2nd half transformation                        *
                !                                                      *
                !*******************************************************

                nQo = Indik(1,nInd)

                do ir=1,nQo

                  it = Indik(1+ir,nInd)

                  do iSh=1,Indx(0,nInd)

                    iaSh = Indx(iSh,nInd)

                    if (.not. Lab%Keep(iaSh,1)) cycle

                    iS = kOffsh(iaSh,lSym)+1

                    if (lSym >= kSym) then

                      ! LJi[k] = sum_a  LaJ[k] * Cai
                      !------------------------------

                      Mode = 'T'
                      n1 = nBasSh(lSym,iaSh)
                      n2 = JNUM

                    else   ! lSym < kSym

                      ! LJi[k] = sum_a  LJa[k] * Cai
                      !------------------------------

                      Mode = 'N'
                      n1 = JNUM
                      n2 = nBasSh(lSym,iaSh)

                    end if

                    call DGEMV_(Mode,n1,n2,One,Lab%SB(iaSh,lSym,1)%A,n1,MSQ(iMOright)%SB(lSym)%A2(iS:,it),1,one,Lik(:,it),1)

                  end do

                  ! Copy LJi[k] in the standard ordered matrix Lik,J

                  if ((jSym == 1) .and. (iMOright == iMOleft)) then
                    itk = nTri_Elem(it-1)+jK
                  else
                    itk = nChOrb_(lSym,iMOright)*(jK-1)+it
                  end if
                  call dcopy_(JNUM,Lik(:,it),1,Rik(itk:),Nik)

                end do

                nInd = nInd+1

                call CWTIME(TCT2,TWT2)
                tmotr2(1) = tmotr(1)+(TCT2-TCT1)
                tmotr2(2) = tmotr(2)+(TWT2-TWT1)

              end do  ! loop over k MOs

              call CWTIME(TCT1,TWT1)

              !*********************************************************
              !                                                        *
              ! 3) Put to disk                                         *
              !                                                        *
              !*********************************************************
              iAdr = Nik*(JVEC-1)
              call DDAFILE(LuRVec(lSym,jDen),1,Rik,Nik*JNUM,iAdr)

              call CWTIME(TCT2,TWT2)
              tmotr(1) = tmotr(1)+(TCT2-TCT1)
              tmotr(2) = tmotr(2)+(TWT2-TWT1)

            end do   ! loop over MOs symmetry

            nullify(pYik)

          end do   ! loop over densities

          call Deallocate_DT(Lab)
          call mma_deallocate(Aux)
          call Deallocate_DT(L_Full)
          !                                                            *
          !*************************************************************
          !*************************************************************
          !*************************************************************
          !                                                            *

          ! ************  END EXCHANGE CONTRIBUTION  ****************

          ! Diagonals updating. It only makes sense if nScreen > 0

          if (Update .and. (nScreen > 0)) then

            call CWTIME(TCS1,TWS1)
            !-------------------------------------------------------------
            ! update the diagonals :   D(a,b) = D(a,b) - sum_J (Lab,J)^2

            ! subtraction is done in the 1st reduced set
#           ifdef _MOLCAS_MPP_
            if (Is_Real_Par()) then

              do krs=1,nRS
                mrs = iiBstR(JSYM,iLoc)+krs
                jrs = IndRed(mrs,iLoc)-iiBstR(JSYM,1)
                do jvc=1,JNUM
                  DiagJ(jrs) = DiagJ(jrs)+Lrs(krs,jvc)**2
                end do
              end do

            else

              do krs=1,nRS
                mrs = iiBstR(JSYM,iLoc)+krs
                jrs = IndRed(mrs,iLoc) ! address in 1st red set
                do jvc=1,JNUM
                  Diag(jrs) = Diag(jrs)-Lrs(krs,jvc)**2
                end do
              end do

            end if

#           else
            do krs=1,nRS
              mrs = iiBstR(JSYM,iLoc)+krs
              jrs = IndRed(mrs,iLoc) ! address in 1st red set
              do jvc=1,JNUM
                Diag(jrs) = Diag(jrs)-Lrs(krs,jvc)**2
              end do
            end do
#           endif

            call CWTIME(TCS2,TWS2)
            tscrn(1) = tscrn(1)+(TCS2-TCS1)
            tscrn(2) = tscrn(2)+(TWS2-TWS1)

          end if

        end if ! DoExchange

        !***************************************************************
        !***************************************************************
        !**                                                           **
        !**    Active term                                            **
        !**                                                           **
        !***************************************************************
        !***************************************************************
        if (DoCAS) then

          call CWTIME(TCC1,TWC1)

          ! Set up the skipping flags
          ! The memory used before for the full-dimension AO-vectors
          !     is now re-used to store half and full transformed
          !     vectors in the active space

          iSwap = 0  ! Lvb,J are returned
          call Allocate_DT(Laq(1),nAorb,nBas,JNUM,JSYM,nSym,iSwap)

          iMO2 = 1
          do iMO1=1,nAdens

            ! iSwap_lxy=5 diagonal blocks are triangular
            ! iSwap_lxy=6 diagonal blocks are square
            iSwap_lxy = 5
            if (iMO1 == 2) iSwap_lxy = 6
            call Allocate_DT(Lxy,nAorb,nAorb,JNUM,JSYM,nSym,iSwap_lxy)

            !***********************************************************
            !                                                          *
            !  MO transformation of Cholesky vectors                   *
            !                                                          *
            !      1) Lvb,J = sum_a  C(v,a) * Lab,J                    *
            !                                                          *
            !***********************************************************

            kMOs = 1
            nMOs = 1  ! Active MOs (1st set)

            call CHO_X_getVtra(irc,Lrs,LREAD,jVEC,JNUM,JSYM,iSwap,IREDC,nMOs,kMOs,Aorb(iMO1),Laq(1),DoRead)

            if (irc /= 0) then
              return
            end if

            !***********************************************************
            !                                                          *
            !  MO transformation of Cholesky vectors                   *
            !                                                          *
            !      2) Lvw,J = sum_b  Lvb,J * C(w,b)                    *
            !                                                          *
            !***********************************************************
            if ((JSYM == 1) .and. (iMO1 == iMO2)) then

              do iSymb=1,nSym

                NAv = nAorb(iSymb)

                if (NAv < 1) cycle

                do JVC=1,JNUM
                  ! triangular blocks
                  call DGEMM_Tri('N','T',NAv,NAv,NBAS(iSymb),One,Laq(1)%SB(iSymb)%A3(:,:,JVC),NAv,Aorb(iMO2)%SB(iSymb)%A2,NAv, &
                                 Zero,Lxy%SB(iSymb)%A2(:,JVC),NAv)

                end do

              end do

            else

              do iSymb=1,nSym

                iSymv = Mul(JSYM,iSymb)
                NAv = nAorb(iSymv)
                NAw = nAorb(iSymb) ! iSymb=iSymw

                if ((NAv*NAw /= 0) .and. (iSymv <= iSymb)) then

                  do JVC=1,JNUM

                    ! square or rectangular blocks
                    call DGEMM_('N','T',NAv,NAw,NBAS(iSymb),One,Laq(1)%SB(iSymv)%A3(:,:,JVC),NAv,Aorb(iMO2)%SB(iSymb)%A2,NAw,Zero, &
                                Lxy%SB(iSymv)%A2(:,JVC),NAv)

                  end do

                end if

              end do

            end if
            !***********************************************************
            !                                                          *
            !  Evaluation of the Z_p_k                                 *
            !                                                          *
            !      Z(p){#J} = sum_xy  T(xy,p) * L(xy,{#J})             *
            !                                                          *
            !  T(xy,p) : is stored by compound symmetry JSYM           *
            !            the indices {xy} are stored as PACKED         *
            !            (sym x < = sym y)                             *
            !                                                          *
            !***********************************************************
            do iTxy=iMO1,nAdens
              iAvec = iMO1+iTxy-1
              do iSymy=1,nSym

                iSymx = Mul(iSymy,JSYM)

                if ((iSymx <= iSymy) .and. (nnA(iSymx,iSymy) /= 0)) then

                  ipZp = iOffZp+nnP(JSYM)*(JVEC-1)+1

                  if (iMO1 == iMO2) then

                    ! diagonal symmetry blocks are triangular
                    call DGEMM_('T','N',nnP(JSYM),JNUM,nnA(iSymx,iSymy),ONE,Txy(ipTxy(iSymx,iSymy,iTxy)),nnP(JSYM), &
                                Lxy%SB(iSymx)%A2,nnA(iSymx,iSymy),ONE,Z_p_k(ipZp,iAvec),nnP(JSYM))

                  else
                    !MGD may rearrange the loops

                    do i=1,nnP(JSYM)
                      ioff = ipTxy(iSymx,iSymy,iTxy)+nnA(iSymx,iSymy)*(i-1)

                      do j=1,JNUM

                        !MGD don't work with symmetry
                        temp = Zero
                        do k=0,nAOrb(iSymx)-1
                          do l=0,k
                            temp = temp+Half*Txy(ioff+iTri(k+1,l))*(Lxy%SB(iSymx)%A2(l+1+nAOrb(iSymx)*k,j)+ &
                                   Lxy%SB(iSymx)%A2(k+1+nAOrb(iSymx)*l,j))
                          end do
                        end do

                        ij = ipZp-1+i+nnP(JSYM)*(j-1)

                        Z_p_k(ij,iAvec) = Z_p_k(ij,iAvec)+temp

                      end do ! j
                    end do   ! i

                  end if

                end if

              end do
            end do

            call Deallocate_DT(Lxy)
          end do

          call CWTIME(TCC2,TWC2)
          tcasg(1) = tcasg(1)+(TCC2-TCC1)
          tcasg(2) = tcasg(2)+(TWC2-TWC1)

          call Deallocate_DT(Laq(1))

        end if  ! DoCAS

      end do  ! end batch loop
      !                                                                *
      !*****************************************************************
      !*****************************************************************
      !**                                                             **
      !**    Epilogue                                                 **
      !**                                                             **
      !*****************************************************************
      !*****************************************************************
      !                                                                *

      ! free memory
      call mma_deallocate(Lrs)

      if (JSYM == 1) then
        call mma_deallocate(Drs)
        if (iMp2prpt == 2) call mma_deallocate(Drs2)
      end if

    end if

    ! Screening control section
    DoScreen = kscreen == nScreen

    if (.not. DoScreen) then
      kscreen = kscreen+1
    else
      kscreen = 1
    end if

    if (DoExchange) then
#     ifdef _MOLCAS_MPP_
      if (Is_Real_Par() .and. Update .and. DoScreen) then
        call GaDsum(DiagJ,nnBSTR(JSYM,1))
        Diag(iiBstR(JSYM,1)+1:iiBstR(JSYM,1)+nnBstR(JSYM,1)) = Diag(iiBstR(JSYM,1)+1:iiBstR(JSYM,1)+nnBstR(JSYM,1))- &
                                                               DiagJ(1:nnBSTR(JSYM,1))
        DiagJ(1:nnBSTR(JSYM,1)) = Zero
      end if
      ! Need to activate the screening to setup the contributing shell
      ! indices the first time the loop is entered .OR. whenever other nodes
      ! have performed screening in the meanwhile
      if (Is_Real_Par() .and. (.not. DoScreen) .and. (nVrs == 0)) then
        ntv0 = ntv0+1
        DoScreen = ((JRED < myJRED1) .or. (ntv0 >= nScreen))
        if (DoScreen) ntv0 = 0
      end if
#     endif
    end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  end do   ! loop over red sets
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (DoExchange) then
    do jDen=1,nKvec
      do i=1,nSym
        call DACLOS(LuRVec(i,jDen))
      end do
    end do
  end if

  !                                                                    *
  !*********************************************************************
  !                                                                    *
end do  ! loop over JSYM
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate a field to be used by Compute_A_jk later
! since allocations cannot be made at that stage
!                                                                      *
!***********************************************************************
!                                                                      *
if (DoExchange) then
  nIJMax = 0
  do jDen=1,nKvec
    do iSym1=1,nSym
      do iSym2=1,nSym
        nIJMax = max(nIJMax,nIJR(iSym1,iSym2,jDen))
      end do
    end do
  end do
end if

call mma_deallocate(iShp_rs)
if (DoExchange) then
  call mma_deallocate(SvShp)
  call mma_deallocate(kOffSh)
  call mma_deallocate(Indik)
  call mma_deallocate(Indx)
  do i=1,nDen
    nullify(SumClk(i)%A)
  end do
  call mma_deallocate(Aux0)
  call mma_deallocate(MLk)
  call mma_deallocate(Yik)
  call mma_deallocate(Ylk)
  call mma_deallocate(AbsC)
  call Deallocate_DT(DiaH)
# ifdef _MOLCAS_MPP_
  if (Is_Real_Par() .and. Update) call mma_deallocate(DiagJ)
# endif
  call mma_deallocate(Diag)
end if

call CWTIME(TOTCPU2,TOTWALL2)
TOTCPU = TOTCPU2-TOTCPU1
TOTWALL = TOTWALL2-TOTWALL1
#ifdef _CD_TIMING_
ChoGet_CPU = TOTCPU
ChoGet_Wall = TOTWALL
#endif

! Write out timing information
if (timings) then

  CFmt = '(2x,A)'
  write(u6,*)
  write(u6,CFmt) 'Cholesky Gradients timing from '//SECNAM
  write(u6,CFmt) '----------------------------------------'
  write(u6,*)
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,CFmt) '                                CPU       WALL   '
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'

  write(u6,'(2x,A26,2f10.2)') 'READ VECTORS                              ',tread(1),tread(2)
  write(u6,'(2x,A26,2f10.2)') 'COULOMB CONTRIB.                          ',tcoul(1),tcoul(2)
  write(u6,'(2x,A26,2f10.2)') 'SCREENING OVERHEAD                        ',tscrn(1),tscrn(2)
  write(u6,'(2x,A26,2f10.2)') 'INACT MO-TRANSFORM VECTORS                ',tmotr(1),tmotr(2)
  write(u6,'(2x,A26,2f10.2)') 'INACT MO-TRANSFORM VECTORS 2              ',tmotr2(1),tmotr2(2)
  write(u6,'(2x,A26,2f10.2)') 'ACTIVE CONTRIB.                           ',tcasg(1),tcasg(2)
  write(u6,*)
  write(u6,'(2x,A26,2f10.2)') 'TOTAL                                     ',TOTCPU,TOTWALL
  write(u6,CFmt) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,*)

end if

irc = 0

return

end subroutine CHO_GET_GRAD
