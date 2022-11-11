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
! Copyright (C) 1990,1991,1993,1998,2006,2007, Roland Lindh            *
!               1990, IBM                                              *
!***********************************************************************

subroutine Drv2El_3Center_RI(Integral_WrOut,ThrAO)
!***********************************************************************
!                                                                      *
!  Object: driver for the 3 center integrals in the RI approach.       *
!                                                                      *
!          This code has three sections                                *
!          1) a 2-center section to generate the Q-vectors             *
!          2) a 3-center section to generate the R-vectors             *
!          3) a partial transpose section to generate the RI vectors   *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified for k2 loop. August '91                         *
!             Modified to minimize overhead for calculations with      *
!             small basis sets and large molecules. Sept. '93          *
!             Modified driver. Jan. '98                                *
!             Modified to 3-center ERIs for RI Jan '06                 *
!             Modified to out-of-core version Feb '07                  *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use RI_procedures, only: Drv2El_2Center_RI
use iSD_data, only: iSD
use Basis_Info, only: dbsc, nBas, nBas_Aux
use Gateway_global, only: force_out_of_core
use Gateway_Info, only: CutInt
use RICD_Info, only: LDF
use Symmetry_Info, only: nIrrep
use RI_glob, only: iShij, iSSOff, klS, Lu_Q, nBasSh, nChV, nSkal_Valence, nSO, ShlSO, SOShl
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
external :: Integral_WrOut
real(kind=wp), intent(in) :: ThrAO
#include "Molcas.fh"
#include "print.fh"
#include "lRI.fh"
#include "iTOffs.fh"
integer(kind=iwp) :: i, iAddr, iAddr_R(0:7), iAdr_AB, iCase, iCenter, iChoVec, id, iIrrep, iLB, iLO, iMax_R(2,0:7), IncVec, &
                     iOff_3C(3,0:7), iOff_Rv(0:7), ip_R, iPass, iPL, iPrint, irc, iRed, iRout, iS, iS_, iSeed, iSO_Aux, iSym, &
                     iTask, iTtmp(0:7), iVec, j_e, j_s, jCenter, jS, jS_, kCenter, kCnttp, klCenter, klS_, kQv, kS, lCenter, &
                     lCnttp, LenVec, LenVec_Red, lJ, lS, Lu_AB, Lu_R(0:7), m3C, MaxCntr, MaxMem, MemLow, MemSew, mMuNu, mQv, &
                     MuNu_e, MuNu_s, n3C, n3CMax, n_Rv, nB_Aux, nCase, nDiag, nMuNu, NoChoVec(0:7), nQv, nRv, nRvMax, nSkal, &
                     nSkal2, nSkal_Auxiliary, nTask, NumVec, NumVec_
real(kind=wp) :: A_int, A_int_kl, TC0, TC1, TCpu1, TCpu2, TMax_all, TW0, TW1, TWall1, Twall2
character(len=6) :: Name_R
logical(kind=iwp) :: DoFock, DoGrad, FreeK2, Indexation, Out_of_Core, Skip, Verbose
integer(kind=iwp), allocatable :: AB(:,:), Addr(:), iRv(:), LBList(:), NuMu(:,:), SO2C(:), TmpList(:)
real(kind=wp), allocatable :: A_Diag(:), Arr_3C(:), Diag(:), Local_A(:,:), Qv(:), Rv(:), TMax_Auxiliary(:), TMax_Valence(:,:), &
                              Tmp(:,:)
integer(kind=iwp), external :: iPrintLevel, IsFreeUnit, nSize_3C, nSize_Rv
logical(kind=iwp), external :: Reduce_Prt, Rsv_Tsk
interface
  subroutine Post_2Center_LDF(A_Diag,AB,MaxCntr,Lu_AB,Local_A,SO2C,nSO_Aux)
    import :: wp, iwp
    real(kind=wp) :: A_Diag(*)
    real(kind=wp), allocatable :: Local_A(:,:)
    integer(kind=iwp), allocatable :: SO2C(:), AB(:,:)
    integer(kind=iwp) :: MaxCntr, Lu_AB, nSO_Aux
  end subroutine Post_2Center_LDF
end interface

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 9

! Get global print level

iPL = iPrintLevel(-1)
if (iPL == 2) then
  iPL = 5
else if (iPL == 3) then
  iPL = 6
else if (iPL == 4) then
  iPL = 99
else if (iPL == 5) then
  iPL = 99
end if
nPrint(iRout) = iPL

! Reduce print level if iterating

if (Reduce_Prt() .and. (iPL <= 5)) then
  nPrint(iRout) = 4
end if
iPrint = nPrint(iRout)

if (iPrint >= 6) call CWTime(TC0,TW0)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
!     2 - C E N T E R   S E C T I O N                                  *
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Compute the two-center integrals over the auxiliary basis

call Drv2El_2Center_RI(ThrAO,A_Diag,nSO_Aux,MaxCntr,SO2C)

! Post processing to generate the Q-vectors.

if (LDF) then

  ! Local RI

  call Post_2Center_LDF(A_Diag,AB,MaxCntr,Lu_AB,Local_A,SO2C,nSO_Aux)

else

  ! Standard RI

  call Post_2Center_RI(A_Diag)

end if

call mma_deallocate(A_Diag)

call Set_Basis_Mode('Auxiliary')
call Nr_Shells(nSkal_Auxiliary)

if (iPrint >= 6) then
  write(u6,'(A)') ' 2-center integrals:'
  call CWTime(TC1,TW1)
  write(u6,'(A,F8.2,A,/,A,F8.2,A)') '      CPU time :',TC1-TC0,' sec.','      Wall time:',TW1-TW0,' sec.'
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
!     3 - C E N T E R   S E C T I O N                                  *
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *

call StatusLine(' Seward:',' Computing 3-center RI integrals')

! Handle both the valence and the auxiliary basis set

call Set_Basis_Mode('WithAuxiliary')
call SetUp_iSD()
!                                                                      *
!***********************************************************************
!                                                                      *
! Initialize for 2-electron integral evaluation. Do not generate
! tables for indexation.

Indexation = .false.
DoGrad = .false.
DoFock = .false.
call Setup_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
nSkal_Valence = nSkal-nSkal_Auxiliary
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute entities for prescreening at shell level

call mma_allocate(TMax_Valence,nSkal_Valence,nSkal_Valence,Label='TMax_Valence')
call mma_allocate(TMax_Auxiliary,nSkal_Auxiliary,Label='TMax_Auxiliary')

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
do iS=1,nSkal_Auxiliary-1
  iS_ = iS+nSkal_Valence
  jS_ = nSkal_Valence+nSkal_Auxiliary
  TMax_Auxiliary(iS) = Tmp(jS_,iS_)
  TMax_all = max(TMax_all,Tmp(jS_,iS_))
end do

call mma_deallocate(Tmp)
!                                                                      *
!***********************************************************************
!                                                                      *
! Set up indexation for Gaussian pairs.

! Generate some offsets and dimensions for the J12 matrix and
! the RI vectors.

call Setup_Aux(nIrrep,nBas,nSkal_Valence,nSkal_Auxiliary,nSO,TMax_Valence,CutInt,nSkal2,nBas_Aux,nChV,iTOffs)

call mma_Allocate(iRv,nSkal2,Label='iRv')
iRv(:) = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! Let us now decide on the memory partitioning

!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Preallocate some core for Seward!

call mma_maxDBLE(MemSew)
MemLow = min(MemSew/2,1024*128)
MemSew = max(MemSew/10,MemLow)
call xSetMem_Ints(MemSew)

!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! During this phase we will have three memory sections
!
! 1) the three center integrals for a fixed {kl}
! 2) a similar block for the R-vectors
! 3) a buffer to contain subsets of the Q-vectors

! Compute the max size of 1 and 2

n3CMax = 0
nRvMax = 0
iMax_R(:,0:nIrrep-1) = 0
do klS_=1,nSkal2
  kS = iShij(1,klS_)
  lS = iShij(2,klS_)
  nRv = nSize_Rv(kS,lS,nBasSh,nSkal-1,nIrrep,iOff_Rv,nChV)
  nRvMax = max(nRvMax,nRv)
  n3C = nSize_3C(kS,lS,nBasSh,nSkal-1,nIrrep,iOff_3C,nBas_Aux)
  n3CMax = max(n3CMax,n3C)
  iMax_R(1,0:nIrrep-1) = max(iMax_R(1,0:nIrrep-1),iOff_3C(1,0:nIrrep-1))
  iMax_R(2,0:nIrrep-1) = iMax_R(2,0:nIrrep-1)+iOff_3C(1,0:nIrrep-1)
end do

call mma_allocate(Arr_3C,n3CMax,Label='Arr_3C')
call mma_allocate(Rv,nRvMax,Label='Rv')

call mma_maxDBLE(MaxMem)
nQv = 0
do iIrrep=0,nIrrep-1
  lJ = nBas_Aux(iIrrep)
  if (iIrrep == 0) lJ = lJ-1 ! remove dummy basis function
  nQv = nQv+lJ*nChV(iIrrep)
end do

! The Q-vectors can be read in a single whole block or in chunks.

if (Force_Out_of_Core) MaxMem = (8*nQv)/10
Out_of_Core = nQv > MaxMem
nQv = min(nQv,MaxMem)  ! note that nQv is effectively reset here
call mma_allocate(Qv,nQv,Label='Qv')
!                                                                      *
!***********************************************************************
!                                                                      *
! In case of in-core mode read Q-vectors only once!

if (.not. Out_of_Core) then
  mQv = 1
  do iIrrep=0,nIrrep-1
    lJ = nBas_Aux(iIrrep)
    if (iIrrep == 0) lJ = lJ-1 ! remove dummy basis function

    if (lJ > 0) then
      iAddr = 0
      kQv = lJ*nChV(iIrrep)
      call dDaFile(Lu_Q(iIrrep),2,Qv(mQv),kQv,iAddr)
      mQv = mQv+kQv
    end if
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Open files for the R-vectors.

do iIrrep=0,nIrrep-1
  nB_Aux = nBas_Aux(iIrrep)
  if (iIrrep == 0) nB_Aux = nB_Aux-1
  if (nB_Aux /= 0) then
    iSeed = 55+iIrrep
    Lu_R(iIrrep) = IsFreeUnit(iSeed)
    write(Name_R,'(A4,I2.2)') 'RVec',iIrrep
    call DaName_MF_WA(Lu_R(iIrrep),Name_R)
  end if
  iAddr_R(iIrrep) = 0
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TCpu1,TWall1)

kCenter = 0  ! dummy initialize
lCenter = 0  ! dummy initialize
iS = nSkal ! point to dummy shell
! Save this field for the time being!
iTtmp(0:nIrrep-1) = iTOffs(3:3*nIrrep:3)

call Init_Tsk(id,nSkal2)

klS = 0
iTask = 0
!do klS=1,nSkal2
do while (Rsv_Tsk(id,klS))
  !write(u6,*) 'Processing shell-pair:',klS
  iTask = iTask+1

  iRv(iTask) = klS
  kS = iShij(1,klS)
  lS = iShij(2,klS)

  ! Logic to avoid integrals with mixed muonic and electronic basis.

  kCnttp = iSD(13,kS)
  lCnttp = iSD(13,lS)

  if (LDF) then

    ! Pick up the corresponding (K|L)^{-1} block

    kCenter = iSD(10,kS)
    lCenter = iSD(10,lS)
    !write(u6,*) 'kCenter, lCenter=',kCenter, lCenter
    klCenter = nTri_Elem(kCenter-1)+lCenter
    iAdr_AB = AB(1,klCenter)
    nAB = AB(2,klCenter)
    call dDaFile(Lu_AB,2,Local_A(:,2),nAB**2,iAdr_AB)
    !call RecPrt('A^-1',' ',Local_A,nAB,nAB)

    ! Now I need some lookup tables to be used below. I need to
    ! go from SO index to lO index and from a given lO index
    ! back to the SO index.

    ISO2LO(:,:) = 0
    iLO = 0
    nCase = 1
    if (kCenter /= lCenter) nCase = 2
    do iCase=1,nCase
      if (iCase == 1) then
        jCenter = kCenter
      else
        jCenter = lCenter
      end if
      do iSO_Aux=1,nSO_Aux
        iCenter = SO2C(iSO_Aux)
        !write(u6,*) 'iCenter=',iCenter
        if (iCenter == jCenter) then
          iLO = iLO+1
          !write(u6,*) 'iLO,iSO_Aux=',iLO,iSO_Aux
          iSO2LO(1,iSO_Aux) = iLO
          iSO2LO(2,iLO) = iSO_Aux
        end if
      end do
    end do
  end if

  A_int_kl = TMax_Valence(kS,lS)
  if (dbsc(kCnttp)%fMass /= dbsc(lCnttp)%fMass) A_int_kl = Zero

  nRv = nSize_Rv(kS,lS,nBasSh,nSkal-1,nIrrep,iOff_Rv,nChV)
  n3C = nSize_3C(kS,lS,nBasSh,nSkal-1,nIrrep,iOff_3C,nBas_Aux)
  Arr_3C(1:n3C) = Zero
  Rv(1:nRv) = Zero

  iTOffs(3:3*nIrrep:3) = iOff_3C(1,0:nIrrep-1)

  ! Loop over the auxiliary basis set

  do jS=nSkal_Valence+1,nSkal-1
    !write(u6,*) 'jS,kS,lS=',jS,kS,lS
    Skip = .false.
    if (LDF) then
      jCenter = iSD(10,jS)
      if ((jCenter /= kCenter) .and. (jCenter /= lCenter)) Skip = .true.
      !write(u6,*) 'jCenter=',jCenter
    end if

    if (.not. Skip) then
      A_int = A_int_kl*TMax_Auxiliary(jS-nSkal_Valence)

#     ifdef _DEBUGPRINT_
      write(u6,*)
      write(u6,*) 'iS,jS,kS,lS=',iS,jS,kS,lS
      write(u6,*) 'A_Int,CutInt=',A_Int,CutInt
      write(u6,*)
#     endif
      if (A_Int >= CutInt) call Eval_IJKL(iS,jS,kS,lS,Arr_3C,n3C,Integral_WrOut)
    end if

  end do    ! jS
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Multiply the 3-center integrals with the Q-vectors

  ! Compute HQ

  call Mult_3C_Qv_S(Arr_3C,n3C,Qv,nQv,Rv,nRv,nChV,iOff_3C,nIrrep,Out_of_Core,Lu_Q,'N')
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Write the R-vectors to disk. These will be retrieved and sort
  ! afterwards in step 3.

  do iIrrep=0,nIrrep-1
    ip_R = 1+iOff_Rv(iIrrep)
    nRv = iOff_3C(1,iIrrep)*nChV(iIrrep)
    !write(u6,*) 'iAddr_R(iIrrep)=',iAddr_R(iIrrep)
    if (nRv > 0) call dDaFile(Lu_R(iIrrep),1,Rv(ip_R),nRv,iAddr_R(iIrrep))
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !end do    ! klS

end do
nTask = iTask

! Restore iTOffs(3,*)
iTOffs(3:3*nIrrep:3) = iTtmp(0:nIrrep-1)

call CWTime(TCpu2,TWall2)
!                                                                      *
!***********************************************************************
!                                                                      *
!                         E P I L O G U E                              *
!                                                                      *
!***********************************************************************
!                                                                      *
call Free_Tsk(id)

! Set up array to store the load balance if this might be needed in
! a gradient calculation.

call mma_allocate(TmpList,nSkal2,Label='TmpList')
TmpList(:) = 0
call mma_allocate(LBList,nSkal2,Label='LBList')
LBList(:) = -1
do iTask=1,nTask
  klS_ = iRv(iTask)
  TmpList(klS_) = 1
end do

iLB = 1
do klS_=1,nSkal2
  if (TmpList(klS_) == 1) then
    LBList(iLB) = klS_
    iLB = iLB+1
  end if
end do

call Put_iArray('LBList',LBList,nSkal2)

call mma_deallocate(LBList)
call mma_deallocate(TmpList)

call mma_deallocate(Rv)
call mma_deallocate(Arr_3C)
call mma_deallocate(Qv)
call xRlsMem_Ints()
call mma_deallocate(TMax_Auxiliary)
call mma_deallocate(TMax_Valence)
if (LDF) then
  call mma_deallocate(SO2C)
  call mma_deallocate(AB)
  call mma_deallocate(Local_A)
  call DaClos(Lu_AB)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Each node does now have an incomplete set of R-vectors!
!                                                                      *
!***********************************************************************
!                                                                      *
! Terminate integral environment.

Verbose = .false.
FreeK2 = .true.
call Term_Ints(Verbose,FreeK2)

call mma_deallocate(iSSOff)
call mma_deallocate(ShlSO)
call mma_deallocate(SOShl)
call Free_iSD()

! Let go off the Q-vectors for now!

do iIrrep=0,nIrrep-1
  nB_Aux = nBas_Aux(iIrrep)
  if (iIrrep == 0) nB_Aux = nB_Aux-1
  if (nB_Aux /= 0) call DaClos(Lu_Q(iIrrep))
end do

if (iPrint >= 6) then
  write(u6,'(A)') ' 3-center integrals:'
  call CWTime(TC0,TW0)
  write(u6,'(A,F8.2,A,/,A,F8.2,A)') '      CPU time :',TC0-TC1,' sec.','      Wall time:',TW0-TW1,' sec.'
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
!     P A R T I A L   T R A N S P O S E   S E C T I O N                *
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! For the interface to work fix the tables of Seward

call Set_Basis_Mode('Valence')
call SetUp_iSD()
!                                                                      *
!***********************************************************************
!                                                                      *
! Initialize for 2-electron integral evaluation. Do generate
! tables for indexation.

Indexation = .true.
call Setup_Ints(nSkal_Valence,Indexation,ThrAO,DoFock,DoGrad)

! Initiate stuff for Cholesky style storage.

call IniCho_RI(nSkal_Valence,nChV,nIrrep,iTOffs,iShij,nSkal2)

call mma_allocate(Addr,nSkal2,Label='Addr') ! addr for read
call mma_allocate(NuMu,2,nSkal2,Label='NuMu')
!                                                                      *
!***********************************************************************
!                                                                      *
! Write out the RI vectors in Cholesky format

! Here we will read one chuck from the R-vector file, while we will
! store an as large part of the RI vectors in Cholesky format.

LenVec = 0
do iIrrep=0,nIrrep-1
  iChoVec = 0

  nB_Aux = nBas_Aux(iIrrep)
  if (iIrrep == 0) nB_Aux = nB_Aux-1
  if (nB_Aux /= 0) then

    iSym = iIrrep+1

    ! NumVec: is no longer equal to the # of auxiliary functions

    NumVec = iTOffs(3*iIrrep+1)
    if (NumVec /= 0) then

      Addr(1) = 0
      do i=2,nTask  ! init the addr for reading vectors
        klS_ = iRv(i-1)
        kS = iShij(1,klS_)
        lS = iShij(2,klS_)
        n3C = nSize_3C(kS,lS,nBasSh,nSkal-1,nIrrep,iOff_3C,nBas_Aux)
        nMuNu = iOff_3C(1,iIrrep)
        Addr(i) = Addr(i-1)+nMuNu*NumVec
      end do

      LenVec_Red = iMax_R(1,iIrrep)
      n_Rv = NumVec*LenVec_Red
      call mma_allocate(Rv,n_Rv,Label='Rv')

      ! LenVec: # of valence Gaussian products in this irrep

      LenVec = iMax_R(2,iIrrep)
      call Create_Chunk(LenVec,NumVec,IncVec)

      do iVec=1,NumVec,IncVec
        NumVec_ = min(NumVec-iVec+1,IncVec)
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Read now the R-vectors for a fixed shell-pair and
        ! irrep, but for all auxiliary functions.

        mMuNu = 0
        do klS_=1,nSkal2
          kS = iShij(1,klS_)
          lS = iShij(2,klS_)
          n3C = nSize_3C(kS,lS,nBasSh,nSkal-1,nIrrep,iOff_3C,nBas_Aux)
          nMuNu = iOff_3C(1,iIrrep)
          m3C = nMuNu*NumVec_

          if (m3C > 0) then
            MuNu_s = mMuNu+1
            MuNu_e = mMuNu+nMuNu

            NuMu(1,klS_) = MuNu_s
            NuMu(2,klS_) = MuNu_e
          end if

          mMuNu = mMuNu+nMuNu
        end do

        do i=1,nTask
          klS_ = iRv(i)
          kS = iShij(1,klS_)
          lS = iShij(2,klS_)

          n3C = nSize_3C(kS,lS,nBasSh,nSkal-1,nIrrep,iOff_3C,nBas_Aux)
          nMuNu = iOff_3C(1,iIrrep)
          m3C = nMuNu*NumVec_

          if (m3C <= 0) cycle

          call dDaFile(Lu_R(iIrrep),2,Rv,m3C,Addr(i))

          ! Copy the appropriate section into the RI vectors in Cholesky format.

          MuNu_s = NuMu(1,klS_)
          MuNu_e = NuMu(2,klS_)
          j_s = 1
          j_e = NumVec_
          call Put_Chunk(MuNu_s,MuNu_e,j_s,j_e,Rv,nMuNu,LenVec)

        end do
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Now transfer the RI vectors to disk

        call Get_Chunk(LenVec,NumVec_,iChoVec,iSym,iVec)

      end do   ! iVec = 1, NumVec, IncVec

      call Destroy_Chunk()
      call mma_deallocate(Rv)

    end if

    ! Let go of the R-vectors for good!

    call DaClos(Lu_R(iIrrep))
  end if
  NoChoVec(iIrrep) = iChoVec

end do    ! iIrrep
call mma_deallocate(NuMu)
call mma_deallocate(Addr)
call mma_deallocate(iRv)
call mma_deallocate(nBasSh)
call mma_deallocate(iShij)
!                                                                      *
!***********************************************************************
!                                                                      *
iPass = 1
iRed = 1
call Cho_RI_PutInfo(iPass,iRed)
!                                                                      *
!***********************************************************************
!                                                                      *
! Terminate integral environment.

Verbose = .false.
FreeK2 = .true.
call Term_Ints(Verbose,FreeK2)
!
if (iPrint >= 6) then
  write(u6,'(A)') ' Block-transpose:'
  call CWTime(TC1,TW1)
  write(u6,'(A,F8.2,A,/,A,F8.2,A)') '      CPU time :',TC1-TC0,' sec.','      Wall time:',TW1-TW0,' sec.'
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
!     D I A G O N A L   S E C T I O N                                  *
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
nDiag = 0
do iIrrep=0,nIrrep-1
  nDiag = nDiag+nBas(iIrrep)
end do
nDiag = nTri_Elem(nDiag)
call mma_allocate(Diag,nDiag,Label='Diag')
Diag(:) = Zero

call Drv2El_RI_Diag(ThrAO,Diag,nDiag)

! Write the diagonal to disk

call Cho_IODiag(Diag,1)

call mma_deallocate(Diag)

if (iPrint >= 6) then
  write(u6,*) 'Diagonal vector:'
  call CWTime(TC0,TW0)
  write(u6,'(A,F8.2,A,/,A,F8.2,A)') '      CPU time :',TC0-TC1,' sec.','      Wall time:',TW0-TW1,' sec.'
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Terminate Cholesky stuff here.

irc = 0
call TermCho_RI(irc,NoChoVec,8)
if (irc /= 0) then
  write(u6,*) 'TermCho_RI returned ',irc
  call SysAbendMsg('Drv2El_3Center_RI','Cholesky termination failed!',' ')
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Drv2El_3Center_RI
