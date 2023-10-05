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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_Drv_ParTwoStep(irc)
!
! Thomas Bondo Pedersen, April 2010.
!
! Purpose: Parallel two-step decomposition of two-electron
!          integrals.

use Index_Functions, only: iTri, nTri_Elem
use Cholesky, only: BlockSize, ChkOnly, Cho_1Center, Cho_DecAlg, Cho_DecAlg_Def, Cho_DiaChk, Cho_Fake_Par, Cho_IntChk, Cho_IOVec, &
                    Cho_MinChk, Cho_No2Center, Cho_PreScreen, Cho_ReOrd, Cho_SimP, Cho_SimRI, Cho_SScreen, Cho_TrcNeg, &
                    Cho_TstScreen, CHo_UseAbs, Diag, Diag_G, Diag_G_Hidden, Diag_Hidden, Did_DecDrv, Frac_ChvBuf, HaltIt, iAlQua, &
                    iAtomShl, Idle, IFCSew, INF_PASS, INF_TIMING, IPRINT, iShP2Q, iShP2RS, LuPri, MaxQual, MinQual, Mode_Screen, &
                    ModRst, Mx2Sh, MxShPr, N1_Qual, N1_VecRd, N2_Qual, N2_VecRd, N_Subtr, nDecom, nDGM_call, nInteg, nMisc, &
                    nSection, nShell, nSym, nSys_call, RstCho, RstDia, ScDiag, SSTau, tDecDrv, tDecom, Thr_PreScreen, Thr_SimRI, &
                    ThrCom, ThrDiag, TimSec, tInteg, tMisc, Tol_DiaChk, Trace_Idle
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp) :: BlockSize_Bak, Cho_DecAlg_Bak, Cho_DecAlg_Def_Bak, Cho_IOVec_Bak, i1, iAlQua_Bak, iBlock, ijBlock, &
                     iPrint_Bak, iSec, iSym, jBlock, l_Z, MaxQual_Bak, MinQual_Bak, Mode_Screen_Bak, ModRst_Bak, MxShPr_Bak, n, &
                     N1_Qual_Bak, N1_VecRD_Bak, N2_Qual_Bak, N2_VecRd_Bak, N_Subtr_Bak, nB, nB_Max, nDGM_Call_Bak, ni, nj, &
                     nnBlock, nSys_Call_Bak
real(kind=wp) :: Byte, C0, C1, Frac_ChVBuf_Bak, SSTau_Bak, tC0, tC1, tCPU0, tCPU1, tDecom_Bak(2,nDecom), Thr_PreScreen_Bak, &
                 Thr_SimRI_Bak, ThrDiag_Bak, TimSec_Bak(4,nSection), tInteg_Bak(2,nInteg), tMisc_Bak(2,nMisc), Tol_DiaChk_Bak, &
                 tW0, tW1, tWall0, tWall1, W0, W1
logical(kind=iwp) :: ChkOnly_Bak, Cho_1Center_Bak, Cho_DiaChk_Bak, Cho_Fake_Par_Bak, Cho_IntChk_Bak, Cho_MinChk_Bak, &
                     Cho_No2Center_Bak, Cho_PreScreen_Bak, Cho_ReOrd_Bak, Cho_SimP_Bak, Cho_SimRI_Bak, Cho_SScreen_Bak, &
                     Cho_TrcNeg_Bak, Cho_TstScreen_Bak, Cho_UseAbs_Bak, Did_DecDrv_Bak, Free_Z, HaltIt_Bak, lConv, RstCho_Bak, &
                     RstDia_Bak, ScDiag_Bak, Trace_Idle_Bak
character(len=2) :: Unt
integer(kind=iwp), allocatable :: iV1Block(:,:), nBlock(:), nVBlock(:,:), NVT(:), ZBlock(:,:)
real(kind=wp), allocatable :: Check(:), Err(:), Z(:)
real(kind=wp), parameter :: DumTst = 0.123456789_wp, DumTol = 1.0e-15_wp
character(len=*), parameter :: myName = 'DPTS', SecNam = 'Cho_Drv_ParTwoStep'

! Preliminaries.
! ==============

#ifdef _DEBUGPRINT_
call Cho_PrtMaxMem('Start of '//SecNam)
#endif

! Start overall timing
if (iPrint >= Inf_Timing) call CWTime(tCPU0,tWall0)

! Init return code
irc = 0

! make a dummy allocation
call mma_allocate(Check,1,Label='Check')
Check(1) = DumTst

lConv = .false.

! Initialization.
! ===============

iSec = 1
if (iPrint >= Inf_Timing) call CWTime(TimSec(1,iSec),TimSec(3,iSec))
call Cho_Init(.false.,.true.)
call GASync()
if (iPrint >= Inf_Timing) then
  call CWTime(TimSec(2,iSec),TimSec(4,iSec))
  call Cho_PrtTim('Cholesky initialization',TimSec(2,iSec),TimSec(1,iSec),TimSec(4,iSec),TimSec(3,iSec),1)
end if
#ifdef _DEBUGPRINT_
call Cho_PrtMaxMem(SecNam//': After init')
#endif

! Get diagonal.
! =============

iSec = 2
if (iPrint >= Inf_Timing) then
  call CWTime(TimSec(1,iSec),TimSec(3,iSec))
  write(LuPri,'(/,A)') '***** Starting Cholesky diagonal setup *****'
  call XFlush(LuPri)
end if
call Cho_GetDiag(lConv)
call GASync()
if (lConv) then
  ! restart is not possible, so it cannot be converged!!
  write(LuPri,'(A,A)') SecNam,': logical error: converged but not restart?!?!'
  call Cho_Quit('Error in '//SecNam,103)
end if
if (iPrint >= Inf_Timing) then
  call CWTime(TimSec(2,iSec),TimSec(4,iSec))
  call Cho_PrtTim('Cholesky diagonal setup',TimSec(2,iSec),TimSec(1,iSec),TimSec(4,iSec),TimSec(3,iSec),1)
end if
#ifdef _DEBUGPRINT_
call Cho_PrtMaxMem(SecNam//': After diagonal')
#endif

! Cholesky decomposition in two steps.
! ====================================

! Time the rest as "decomposition driver"
call CWTime(C0,W0)

iSec = 3
if (iPrint >= Inf_Timing) then
  call CWTime(TimSec(1,iSec),TimSec(3,iSec))
  write(LuPri,'(/,A)') '***** Starting Cholesky decomposition *****'
  call XFlush(LuPri)
end if

! Perform first CD step: get parent diagonals and Z vectors.
! Z vectors contain more elements than needed at this stage!
! (more elements than vectors). The final Z vectors are obtained
! below (Cho_GetZ).
call Cho_P_SetAddr()
call Cho_DecDrv(Diag)
call GASync()
if (iPrint >= Inf_Timing) then
  call CWTime(tC1,tW1)
  call Cho_PrtTim('Cholesky map generation',tC1,TimSec(1,iSec),tW1,TimSec(3,iSec),2)
end if
#ifdef _DEBUGPRINT_
call Cho_PrtMaxMem(SecNam//': After 1st step')
#endif

! Shut down and re-initialize.
! This is done to get rid of the vast amount of different data
! generated by the parallel decomposition in particular.
! In addition, it will make it easier to use the vector calculator
! externally (i.e. after seward execution)
if (iPrint >= Inf_Timing) call CWTime(tC0,tW0)
TimSec_Bak(:,:) = TimSec(:,:)
tInteg_Bak(:,:) = tInteg(:,:)
tDecom_Bak(:,:) = tDecom(:,:)
tMisc_Bak(:,:) = tMisc(:,:)
Trace_Idle_Bak = Trace_Idle
Cho_UseAbs_Bak = Cho_UseAbs
Cho_DiaChk_Bak = Cho_DiaChk
Cho_Fake_Par_Bak = Cho_Fake_Par
Cho_SimP_Bak = Cho_SimP
Cho_ReOrd_Bak = Cho_ReOrd
ChkOnly_Bak = ChkOnly
Cho_IntChk_Bak = Cho_IntChk
Cho_MinChk_Bak = Cho_MinChk
Cho_TrcNeg_Bak = Cho_TrcNeg
Cho_TstScreen_Bak = Cho_TstScreen
RstDia_Bak = RstDia
RstCho_Bak = RstCho
Mode_Screen_Bak = Mode_Screen
Cho_DecAlg_Def_Bak = Cho_DecAlg_Def
ModRst_Bak = ModRst
N_Subtr_Bak = N_Subtr
Did_DecDrv_Bak = Did_DecDrv
HaltIt_Bak = HaltIt
Tol_DiaChk_Bak = Tol_DiaChk
Cho_1Center_Bak = Cho_1Center
Cho_No2Center_Bak = Cho_No2Center
Cho_PreScreen_Bak = Cho_PreScreen
Thr_PreScreen_Bak = Thr_PreScreen
ThrDiag_Bak = ThrDiag
ScDiag_Bak = ScDiag
MinQual_Bak = MinQual
MaxQual_Bak = MaxQual
N1_Qual_Bak = N1_Qual
N2_Qual_Bak = N2_Qual
MxShPr_Bak = MxShPr
iAlQua_Bak = iAlQua
Cho_IOVec_Bak = Cho_IOVec
Frac_ChVBuf_Bak = Frac_ChVBuf
Cho_SScreen_Bak = Cho_SScreen
SStau_Bak = SStau
Cho_SimRI_Bak = Cho_SimRI
Thr_SimRI_Bak = Thr_SimRI
Cho_DecAlg_Bak = Cho_DecAlg
BlockSize_Bak = BlockSize
iPrint_Bak = iPrint
Cho_IOVec_Bak = Cho_IOVec
N1_VecRd_Bak = N1_VecRd
N2_VecRd_Bak = N2_VecRd
nSys_Call_Bak = nSys_Call
nDGM_Call_Bak = nDGM_Call
call Cho_P_WrDiag()
call Cho_Final(.true.)
call Cho_P_OpenVR(2)
call Cho_X_Init(irc,Zero)
if (irc /= 0) then
  write(LuPri,*) SecNam,': Cho_X_Init returned code ',irc
  irc = 1
  ! clear memory and return
  call Finish_this()
  return
end if
if (Cho_1Center) then
  if (.not. allocated(iAtomShl)) then
    call mma_allocate(iAtomShl,nShell,Label='iAtomShl')
    call Cho_SetAtomShl(irc,iAtomShl,size(iAtomShl))
    if (irc /= 0) then
      write(LuPri,'(A,A,I8)') SecNam,': Cho_SetAtomShl returned code',irc
      irc = 1
      ! clear memory and return
      call Finish_this()
      return
    end if
  end if
end if
TimSec(:,:) = TimSec_Bak(:,:)
tInteg(:,:) = tInteg_Bak(:,:)
tDecom(:,:) = tDecom_Bak(:,:)
tMisc(:,:) = tMisc_Bak(:,:)
Trace_Idle = Trace_Idle_Bak
Cho_UseAbs = Cho_UseAbs_Bak
Cho_DiaChk = Cho_DiaChk_Bak
Cho_Fake_Par = Cho_Fake_Par_Bak
Cho_SimP = Cho_SimP_Bak
Cho_ReOrd = Cho_ReOrd_Bak
ChkOnly = ChkOnly_Bak
Cho_IntChk = Cho_IntChk_Bak
Cho_MinChk = Cho_MinChk_Bak
Cho_TrcNeg = Cho_TrcNeg_Bak
Cho_TstScreen = Cho_TstScreen_Bak
RstDia = RstDia_Bak
RstCho = RstCho_Bak
Mode_Screen = Mode_Screen_Bak
Cho_DecAlg_Def = Cho_DecAlg_Def_Bak
ModRst = ModRst_Bak
N_Subtr = N_Subtr_Bak
Did_DecDrv = Did_DecDrv_Bak
HaltIt = HaltIt_Bak
Tol_DiaChk = Tol_DiaChk_Bak
Cho_1Center = Cho_1Center_Bak
Cho_No2Center = Cho_No2Center_Bak
Cho_PreScreen = Cho_PreScreen_Bak
Thr_PreScreen = Thr_PreScreen_Bak
ThrDiag = ThrDiag_Bak
ScDiag = ScDiag_Bak
MinQual = MinQual_Bak
MaxQual = MaxQual_Bak
N1_Qual = N1_Qual_Bak
N2_Qual = N2_Qual_Bak
N1_VecRd = N1_VecRd_Bak
N2_VecRd = N2_VecRd_Bak
MxShPr = MxShPr_Bak
iAlQua = iAlQua_Bak
Cho_IOVec = Cho_IOVec_Bak
Frac_ChVBuf = Frac_ChVBuf_Bak
Cho_SScreen = Cho_SScreen_Bak
SStau = SStau_Bak
Cho_SimRI = Cho_SimRI_Bak
Thr_SimRI = Thr_SimRI_Bak
Cho_DecAlg = Cho_DecAlg_Bak
BlockSize = BlockSize_Bak
iPrint = iPrint_Bak
Cho_IOVec = Cho_IOVec_Bak
N1_VecRd = N1_VecRd_Bak
N2_VecRd = N2_VecRd_Bak
nSys_Call = nSys_Call_Bak
nDGM_Call = nDGM_Call_Bak
! Allocate memory for extracting integral columns directly in
! reduced set from Seward. Set Seward interface to 3 (to treat
! columns shell pair-wise).
IFCSEW = 3
call mma_allocate(iShP2RS,2,Mx2Sh,Label='iShP2RS')
call mma_allocate(iShP2Q,2,Mx2Sh,Label='iShP2Q ')
if (iPrint >= Inf_Timing) then
  call CWTime(tC1,tW1)
  call Cho_PrtTim('Cholesky reinitialization',tC1,tC0,tW1,tW0,2)
end if
#ifdef _DEBUGPRINT_
call Cho_PrtMaxMem(SecNam//': After re-init')
#endif

! Get Z vectors in core
if (iPrint >= Inf_Timing) call CWTime(tC0,tW0)
call mma_allocate(NVT,nSym,Label='NVT')
call Cho_X_GetTotV(NVT,size(NVT))
call Cho_ZMem(irc,l_Z,NVT,size(NVT),iPrint >= Inf_Timing,.true.)
if (irc /= 0) then
  write(LuPri,'(A,A,I6)') SecNam,': Cho_ZMem returned code',irc
  if (irc == 999) then
    if (iPrint < Inf_Timing) call Cho_ZMem(irc,l_Z,NVT,size(NVT),.true.,.false.)
    call mma_maxDBLE(l_Z)
    call Cho_Word2Byte(l_Z,8,Byte,Unt)
    write(LuPri,'(A,I12,A,F7.3,1X,A,A)') 'Largest available memory block:',l_Z,' words (',Byte,Unt,')'
    write(LuPri,'(A)') '=> INSUFFICIENT MEMORY FOR STORING Z VECTORS!'
    write(LuPri,'(/,A,/,A)') 'You have the following options:','(a) Increase available memory (MOLCAS_MEM),'
    write(LuPri,'(A)') 'and/or'
    write(LuPri,'(A,/,A,/,A,/,A,/,A,/,A)') '(b) -SERIAL EXECUTION:', &
                                           '       Use the serial two-step algorithm by specifying the keywords', &
                                           '          ChoInput','          TwoStep','          EndChoInput', &
                                           '       in Seward input.'
    write(LuPri,'(A,/,A,/,A,/,A,/,A,/,A,/,A)') '    -PARALLEL EXECUTION:', &
                                               '       Use the parallel one-step algorithm by specifying the keywords', &
                                               '          ChoInput','          OneStep','          Parallel', &
                                               '          EndChoInput','       in Seward input.'
    call Cho_Quit(SecNam//': Insufficient memory for Z vectors',101)
  end if
  irc = 1
  ! clear memory and return
  call Finish_this()
  return
end if
call mma_allocate(Z,l_Z,Label='Z')
call mma_allocate(nBlock,nSym,Label='nBlock')
nB_Max = 0
do iSym=1,nSym
  nB = (NVT(iSym)-1)/BlockSize+1
  nB_Max = max(nB_Max,nB)
  nBlock(iSym) = nB
end do
call mma_allocate(nVBlock,nB_Max,nSym,Label='nVBlock')
call mma_allocate(iV1Block,nB_Max,nSym,Label='iV1Block')
nVBlock(:,:) = 0
iV1Block(:,:) = 0

do iSym=1,nSym
  i1 = 1
  do iBlock=1,nBlock(iSym)-1
    nVBlock(iBlock,iSym) = BlockSize
    iV1Block(iBlock,iSym) = i1
    i1 = i1+BlockSize
  end do
  nVBlock(nBlock(iSym),iSym) = NVT(iSym)-BlockSize*(nBlock(iSym)-1)
  iV1Block(nBlock(iSym),iSym) = i1
end do
nnBlock = nTri_Elem(nB_Max)
call mma_allocate(ZBlock,nnBlock,nSym,Label='ZBlock')
ZBlock(:,:) = 0

n = 1
do iSym=1,nSym
  do jBlock=1,nBlock(iSym)
    nj = nVBlock(jBlock,iSym)
    ijBlock = iTri(jBlock,jBlock)
    ZBlock(ijBlock,iSym) = n
    n = n+nTri_ELem(nj)
    do iBlock=jBlock+1,nBlock(iSym)
      ni = nVBlock(iBlock,iSym)
      ijBlock = iTri(iBlock,jBlock)
      ZBlock(ijBlock,iSym) = n
      n = n+ni*nj
    end do
  end do
end do
call Cho_GetZ(irc,NVT,size(NVT),nBlock,size(nBlock),nVBlock,nB_Max,nSym,iV1Block,nB_Max,nSym,ZBlock,nnBlock,nSym,Z,l_Z)
if (irc /= 0) then
  write(LuPri,*) SecNam,': Cho_GetZ returned code ',irc
  irc = 1
  ! clear memory and return
  call Finish_this()
  return
end if
if (iPrint >= Inf_Timing) then
  call CWTime(tC1,tW1)
  call Cho_PrtTim('Cholesky Z vector fetching',tC1,tC0,tW1,tW0,2)
end if
#ifdef _DEBUGPRINT_
call Cho_PrtMaxMem(SecNam//': After GetZ')
#endif

! Perform second step: calculation of Cholesky vectors from Z
! vectors and integrals.
if (iPrint >= Inf_Timing) then
  call CWTime(tC0,tW0)
end if
Free_Z = .true.
call Cho_X_CompVec(irc,NVT,size(NVT),nBlock,size(nBlock),nVBlock,nB_Max,nSym,iV1Block,nB_Max,nSym,ZBlock,nnBlock,nSym,Z,l_Z,Free_Z)
if (Free_Z) call mma_deallocate(Z)
if (irc /= 0) then
  write(LuPri,*) SecNam,': Cho_X_CompVec returned code ',irc
  irc = 1
  ! clear memory and return
  call Finish_this()
  return
end if
! Write restart files
call Cho_PTS_WrRst(irc,NVT,size(NVT))
if (irc /= 0) then
  write(LuPri,*) SecNam,': Cho_PTS_WrRst returned code ',irc
  irc = 1
  ! clear memory and return
  call Finish_this()
  return
end if
if (iPrint >= Inf_Timing) then
  call CWTime(tC1,tW1)
  call Cho_PrtTim('Cholesky vector generation',tC1,tC0,tW1,tW0,2)
end if

! Final timing of decomposition section
if (iPrint >= Inf_Timing) then
  call CWTime(TimSec(2,iSec),TimSec(4,iSec))
  call Cho_PrtTim('Cholesky decomposition',TimSec(2,iSec),TimSec(1,iSec),TimSec(4,iSec),TimSec(3,iSec),1)
end if

! time as "decomposition driver"
call CWTime(C1,W1)
tDecDrv(1) = C1-C0
tDecDrv(2) = W1-W0
#ifdef _DEBUGPRINT_
call Cho_PrtMaxMem(SecNam//': After 2nd step.')
#endif

! Check diagonal.
! ===============

iSec = 4
if (iPrint >= Inf_Timing) then
  call CWTime(TimSec(1,iSec),TimSec(3,iSec))
  write(LuPri,'(/,A)') '***** Starting Cholesky diagonal check *****'
  call XFlush(LuPri)
end if
call mma_allocate(Err,4,Label='Err')
iPrint_Bak = iPrint
if (iPrint < Inf_Pass) iPrint = -99999999 ! suppress printing in Cho_X_CheckDiag
call Cho_X_CheckDiag(irc,Err)
iPrint = iPrint_Bak
if (irc /= 0) then
  write(LuPri,*) SecNam,': Cho_X_CheckDiag returned code ',irc
  irc = 1
  ! release memory and return
  call Finish_this()
  return
end if
if (Err(2) > ThrCom) then
  write(LuPri,'(/,A)') 'Cholesky decomposition failed!'
  write(LuPri,'(3X,A,1P,D15.6)') 'Largest integral diagonal..',Err(2)
  write(LuPri,'(3X,A,1P,D15.6)') 'Decomposition threshold....',ThrCom
  irc = 1
  ! release memory and return
  call Finish_this
  return
end if
call mma_deallocate(Err)
if (iPrint >= Inf_Timing) then
  call CWTime(TimSec(2,iSec),TimSec(4,iSec))
  call Cho_PrtTim('Cholesky diagonal check',TimSec(2,iSec),TimSec(1,iSec),TimSec(4,iSec),TimSec(3,iSec),1)
end if
#ifdef _DEBUGPRINT_
call Cho_PrtMaxMem(SecNam//': After diagonal check.')
#endif

! Finalization.
! =============

iSec = 8
if (iPrint >= Inf_Timing) then
  call CWTime(TimSec(1,iSec),TimSec(3,iSec))
  write(LuPri,'(/,A)') '***** Starting Cholesky finalization *****'
  call XFlush(LuPri)
end if
if (allocated(Idle)) call mma_deallocate(Idle)
call Cho_PTS_Final(NVT,size(NVT))
if (iPrint >= Inf_Timing) then
  call CWTime(TimSec(2,iSec),TimSec(4,iSec))
  call Cho_PrtTim('Cholesky finalization',TimSec(2,iSec),TimSec(1,iSec),TimSec(4,iSec),TimSec(3,iSec),1)
end if
#ifdef _DEBUGPRINT_
call Cho_PrtMaxMem(SecNam//': After finalization')
#endif

! Statistics.
! ===========

if (iPrint >= 1) then
  iSec = 9
  if (iPrint >= Inf_Timing) then
    call CWTime(TimSec(1,iSec),TimSec(3,iSec))
    write(LuPri,'(/,A)') '***** Starting Cholesky statistics *****'
    call XFlush(LUPRI)
  end if
  call Cho_PTS_Stat()
  call GASync()
  if (iPrint >= Inf_Timing) then
    call CWTime(TimSec(2,iSec),TimSec(4,iSec))
    call Cho_PrtTim('Cholesky statistics',TimSec(2,iSec),TimSec(1,iSec),TimSec(4,iSec),TimSec(3,iSec),1)
  end if
end if
#ifdef _DEBUGPRINT_
call Cho_PrtMaxMem(SecNam//': After statistics')
#endif

! Wrap it up and return.
! ======================

call mma_deallocate(iV1Block)
call mma_deallocate(nVBlock)
call mma_deallocate(ZBlock)
call mma_deallocate(nBlock)
call mma_deallocate(NVT)

! Close vector and restart files
call Cho_OpenVR(2,2)

call Finish_this()

contains

! error termination point
subroutine Finish_this()

  ! check memory
  if (abs(DumTst-Check(1)) > DumTol) then
    write(LuPri,*) SecNam,': memory has been out of bounds [2]'
    irc = 2
  end if

  if (allocated(Diag_Hidden)) call mma_deallocate(Diag_Hidden)
  if (allocated(Diag_G_Hidden)) call mma_deallocate(Diag_G_Hidden)
  nullify(Diag)
  nullify(Diag_G)
  call mma_deallocate(Check)

  ! Print total timing
  if ((iPrint >= Inf_Timing) .and. (irc == 0)) then
    call CWTime(tCPU1,tWall1)
    call Cho_PrtTim('Cholesky Procedure',tCPU1,tCPU0,tWall1,tWall0,1)
  end if

  call XFlush(LuPri)
# ifdef _DEBUGPRINT_
  call Cho_PrtMaxMem('End of '//SecNam)
# endif

end subroutine Finish_this

end subroutine Cho_Drv_ParTwoStep
