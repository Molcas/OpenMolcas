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
! Copyright (C) Jonas Bostrom                                          *
!***********************************************************************

subroutine Mult_with_Q_CASPT2(nBas_aux,nBas,nV_t,nIrrep,SubAux)
!***********************************************************************
! Author: Jonas Bostrom
!
! Purpose: Multiply CASPT2 A~_sep and B~_sep with inverse cholesky factors.
!
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use Cholesky, only: nSym, NumCho
use pso_stuff, only: A_PT2, LuGamma2
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp
#ifdef _MOLCAS_MPP_
use para_info, only: is_real_par
#endif
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nIrrep, nBas_Aux(nIrrep), nBas(nIrrep), nV_t(nIrrep)
logical(kind=iwp), intent(in) :: SubAux
integer(kind=iwp) :: i, iAdrQ, iAO, iAOlast, iAOstart, id, iOffQ1, iOpt, iost, ip_B, ip_B2, iSym, j, jAO, jSym, jVec, kSym, kVec, &
                     l_A_ht, l_A_t, l_B_t, l_Q, lRealName, Lu_Q, LUAPT2, LUGAMMA, lVec, MaxMem, myRank, nBas2, nBasTri, nCalAO, &
                     nCalAO_tot, nLR, nLRb(8), NPROCS, nSkal2_, NumAux, NumCV, NumCVt, NumVecJ, NumVecK, nVec
real(kind=wp) :: aaa, Fac, TotCPU0, TotCPU1, TotWall0, TotWall1
logical(kind=iwp) :: is_error, Found
character(len=4096) :: RealName
character(len=6) :: Name_Q
integer(kind=iwp), allocatable :: AOList(:,:), IWRK(:,:), LBList(:), nList_AO(:), nList_Shell(:)
real(kind=wp), allocatable :: A_ht(:), A_t(:), B_t(:), QVec(:)
character(len=*), parameter :: SECNAM = 'Mult_with_Q_CASPT2'
integer(kind=iwp), external :: IsFreeUnit
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
logical(kind=iwp) :: bStat
integer(kind=iwp) :: iHi1, iHi2, iLo1, iLo2, jHi1, jHi2, jLo1, jLo2, lg_V1, lg_V2, mV1, mV2, nDim
integer(kind=iwp), allocatable :: nList(:)
#endif

#ifdef _MOLCAS_MPP_
if (is_real_par()) then
  NPROCS = GA_nNodes()
  myRank = GA_NodeID()
else
#endif
  NPROCS = 1
  myRank = 0
#ifdef _MOLCAS_MPP_
end if
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TotCPU0,TotWall0)
!                                                                      *
!***********************************************************************
!                                                                      *
do iSym=1,nSym
  NumCV = NumCho(iSym)
  NumAux = nBas_Aux(iSym)-1
  if (SubAux) NumAux = nBas_Aux(iSym)-1
  nLR = 0
  do jSym=1,nSym
    kSym = Mul(iSym,jSym)
    nLR = nLR+nBas(jSym)*nBas(kSym)
  end do
  nLRb(iSym) = nLR
end do
!                                                                      *
!***********************************************************************
!                                                                      *
do iSym=1,nSym

  nBas2 = nLRb(iSym)
  !write(u6,*) 'nBas2 = ',nBas2
  NumCV = NumCho(iSym) ! local # of vecs in parallel run
  NumCVt = nV_t(iSym)  ! total # of vecs
  NumAux = nBas_Aux(iSym)-1
  if (SubAux) NumAux = nBas_Aux(iSym)-1

  ! Get Q-vectors from disk
  ! -----------------------

  l_Q = NumCVt*NumAux
  call mma_allocate(QVec,l_Q,Label='Q_Vector')

  Lu_Q = IsFreeUnit(7)
  write(Name_Q,'(A4,I2.2)') 'QVEC',iSym-1
  call DaName_MF_WA(Lu_Q,Name_Q)

  iOpt = 2
  iAdrQ = 0
  call dDaFile(Lu_Q,iOpt,QVec,l_Q,iAdrQ)

  ! Get MP2 A-tilde vectors from disk
  ! ---------------------------------

  l_A_t = NumCVt*NumCVt
  l_A_ht = NumAux*NumCVt
  call mma_allocate(A_t,l_A_t,Label='A_t')
  call mma_allocate(A_ht,l_A_ht,Label='A_ht')

  ! Read A_PT2 from LUAPT2
  LuAPT2 = isFreeUnit(68)
  call daname_mf_wa(LUAPT2,'A_PT2')
  id = 0
  call ddafile(LUAPT2,2,A_t,l_A_t,id)

  ! Symmetrized A_PT2
  do i=1,NumCVt
    do j=1,i
      aaa = Half*(A_t((j-1)*NumCVt+i)+A_t((i-1)*NumCVt+j))
      A_t((j-1)*NumCVt+i) = aaa
      A_t((i-1)*NumCVt+j) = aaa
    end do
  end do

# ifdef _DEBUGPRINT_
  write(u6,*) 'Q-vectors'
  do i=1,l_Q
    write(u6,*) QVec(i)
  end do

  write(u6,*) 'A-vectors'
  do i=1,l_A_t
    write(u6,*) A_t(i)
  end do
# endif

  ! Make first halftransformation to cholesky-base
  ! ----------------------------------------------

  call dGemm_('N','N',NumAux,NumCVt,NumCVt,One,QVec,NumAux,A_t,NumCVt,Zero,A_ht,NumAux)

  call dGemm_('N','T',NumAux,NumAux,NumCVt,One,A_ht,NumAux,QVec,NumAux,Zero,A_PT2,NumAux)

  ! Leave the vectors in A_PT2

  call mma_deallocate(A_t)
  call mma_deallocate(A_ht)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Construct AO list for all nodes (should be ouside the symmetry loop)

  call mma_allocate(nList_Shell,NPROCS,Label='nList_Shell')
  call mma_allocate(nList_AO,NPROCS,Label='nList_AO')
  call mma_allocate(AOList,2,nBas(1)**2,Label='AOList')
  nList_Shell = 0
  nList_AO = 0
  AOList = 0

  ! On exit, we will have the offset array of Shell and AO.
  ! nCalAO is the number of AOs of the present rank,
  ! and AOList is the AO index to be processed in the present rank
  nCalAO = 0
  call construct_AOList()

  call mma_allocate(IWRK,2,nCalAO,Label='IWRK')
  IWRK(:,:) = AOList(:,1:nCalAO)
  AOList(:,:) = 0
  AOList(:,iAOstart:iAOlast) = IWRK(:,:)
  call mma_deallocate(IWRK)
# ifdef _MOLCAS_MPP_
  if (is_real_par()) call GAIGOP(AOList,2*nCalAO_tot,'+')
# endif
  ! nCalAO will be the number of AOs to be computed in all nodes
  nCalAO = nCalAO_tot
  call mma_maxDBLE(MaxMem)
  MaxMem = 9*(MaxMem/10)
  call mma_allocate(B_t,MaxMem,Label='B_t')

  ! Applicable to C1 only
  nBasTri = max(nTri_Elem(nBas(1)),nCalAO)
  nVec = MaxMem/(2*nBasTri+1) ! MaxMem/(2*nLRb(iSym)+1)
  nVec = min(max(NumCVt,NumAux),nVec)
  if (nVec < 1) call SysAbendMsg(SecNam,'nVec is non-positive','[1]')

  l_B_t = nBasTri*nVec ! nLRb(iSym)*nVec
  ip_B = 1+l_B_t
  ip_B2 = ip_B+l_B_t

# ifdef _MOLCAS_MPP_
  if (is_real_par()) then
    ! If GA is used, we do not use much replicated memory
    call mma_deallocate(B_t)
    call mma_allocate(B_t,max(nBas2,NumAux),Label='B_t')

    ! How to inquire the available GA memory?
    ! original vector is vertial stripe
    call mma_allocate(nList,NPROCS,Label='nList')
    nList = 0
    nList(myRank+1) = NumCho(iSym)
    call GAIGOP(nList,NPROCS,'+')
    kVec = 1
    do i=1,NPROCS
      lVec = nList(i)
      nList(i) = kVec
      kVec = kVec+lVec
    end do
    bStat = GA_CREATE_IRREG(MT_DBL,nCalAO,NumCVt,'WRK1',[1],1,NLIST,NPROCS,LG_V1)
    ! resultant vector is horizontal stripe
    ! This global array is constructed so that each local memory corresponds to the AOs processed on each node
    bStat = bStat .and. GA_CREATE_IRREG(MT_DBL,nCalAO,NumCVt,'WRK2',nList_AO,NPROCS,[1],1,LG_V2)

    if (.not. bStat) then
      bStat = GA_destroy(lg_V1)
      bStat = GA_destroy(lg_V2)
      call mma_deallocate(B_t)
      call mma_allocate(B_t,MaxMem,Label='B_t')
    end if
    call mma_deallocate(nList)
  end if
# endif

  LuGAMMA = IsFreeUnit(65)
  call PrgmTranslate('GAMMA',RealName,lRealName)
  call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),'DIRECT','UNFORMATTED',iost,.true.,nBas2*8,'OLD',is_error)

  LuGamma2 = IsFreeUnit(67)
  call PrgmTranslate('GAMMA2',RealName,lRealName)
  call MOLCAS_Open_Ext2(LuGamma2,RealName(1:lRealName),'DIRECT','UNFORMATTED',iost,.true.,NumAux*8,'REPLACE',is_error)

  ! The B-vectors should be read one batch at the time
  ! --------------------------------------------------

# ifdef _MOLCAS_MPP_
  ! Cholesky -> Auxiliary transformation using distributed arrays
  if (is_real_par()) then
    ! with sufficient (how to know?) distributed memory, use GA
    ! lg_V1 is used as a working space, but probably inappropriate
    call GA_Distribution(lg_V1,myRank,iLo1,iHi1,jLo1,jHi1)
    call GA_Access(lg_V1,iLo1,iHi1,jLo1,jHi1,mV1,nDim)
    do lVec=1,NumCV
      read(LuGAMMA,rec=lVec) B_t(1:nBas2)
      ! symmetrize (mu nu | P)
      do i=1,nCalAO
        iAO = AOList(1,i)
        jAO = AOList(2,i)
        aaa = B_t(iAO+nBas(1)*(jAO-1))+B_t(jAO+nBas(1)*(iAO-1))
        DBL_MB(mV1+i-1+nCalAO*(lVec-1)) = aaa
      end do
    end do

    ! Put in lg_V2 (horizontal)
    call GA_Put(lg_V2,iLo1,iHi1,jLo1,jHi1,DBL_MB(mV1),nDim)
    ! Destroy the vertical GA
    call GA_Release(lg_V1,iLo1,iHi1,jLo1,jHi1)
    bStat = GA_destroy(lg_V1)
    ! Remake a horizontal GA, it is used as a resultant array
    bStat = GA_CREATE_IRREG(MT_DBL,nCalAO,NumAux,'WRK1',nList_AO,NPROCS,[1],1,lg_V1)

    ! DGEMM can be done locally
    call GA_Distribution(lg_V1,myRank,iLo1,iHi1,jLo1,jHi1)
    call GA_Access(lg_V1,iLo1,iHi1,jLo1,jHi1,mV1,nDim)
    call GA_Distribution(lg_V2,myRank,iLo2,iHi2,jLo2,jHi2)
    call GA_Access(lg_V2,iLo2,iHi2,jLo2,jHi2,mV2,nDim)
    ndim = iAOlast-iAOstart+1 ! Number of AOs processed in this node
    call DGEMM_('N','T',ndim,NumAux,NumCVt,Half,DBL_MB(mV2),ndim,QVec,NumAux,Zero,DBL_MB(mV1),ndim)
    call GA_Release(lg_V2,iLo2,iHi2,jLo2,jHi2)
    bStat = GA_destroy(lg_V2)

    ! Get the local chunk and write to disk
    ! It is of course possible to avoid using disk
    do jVec=1,ndim
      call DCopy_(NumAux,DBL_MB(mV1+jVec-1),nDim,B_t,1)
      write(LuGAMMA2,rec=jVec) B_t(1:NumAux)
    end do
    call GA_Release(lg_V1,iLo1,iHi1,jLo1,jHi1)
    bStat = GA_destroy(lg_V1)
  else
# endif
    do kVec=1,NumAux,nVec
      NumVecK = min(nVec,NumAux-(kVec-1))

      do jVec=1,NumCV,nVec
        NumVecJ = min(nVec,NumCV-(jVec-1))

        do lVec=1,NumVecJ
          read(LuGAMMA,rec=jVec+lVec-1) B_t(ip_B2:ip_B2+nBas2-1)
          ! symmetrize (mu nu | P)
          ! only the lower triangle part is used
          do i=1,nCalAO
            iAO = AOList(1,i)
            jAO = AOList(2,i)
            aaa = B_t(ip_B2+iAO-1+nBas(1)*(jAO-1))+B_t(ip_B2+jAO-1+nBas(1)*(iAO-1))
            B_t(i+nCalAO*(lVec-1)) = aaa
          end do
        end do

        Fac = Zero
        if (jVec /= 1) Fac = One
        iOffQ1 = kVec+NumAux*(jVec-1)
        call dGemm_('N','T',nCalAO,NumVecK,NumVecJ,Half,B_t,nCalAO,QVec(iOffQ1),NumAux,Fac,B_t(ip_B),nCalAO)
      end do

      ! (mu nu | P) --> (P | mu nu)
      ! For parallel, write only the density used in their nodes
      if (max(NumCV,NumAux) == nVec) then
        do jVec=iAOstart,iAOlast
          call DCopy_(NumAux,B_t(ip_B+jVec-1),nCalAO,B_t,1)
          write(LuGAMMA2,rec=jVec-iAOstart+1) B_t(1:NumAux)
        end do
      else
        do jVec=iAOstart,iAOlast
          if (kVec /= 1) read(LuGAMMA2,rec=jVec-iAOstart+1) B_t(1:kVec-1)
          call DCopy_(NumVecK,B_t(ip_B+jVec-1),nCalAO,B_t(kVec),1)
          write(LuGAMMA2,rec=jVec-iAOstart+1) B_t(1:kVec+NumVecK-1)
        end do
      end if
    end do
# ifdef _MOLCAS_MPP_
  end if
# endif

  call mma_deallocate(nList_Shell)
  call mma_deallocate(nList_AO)

  ! Leave LuGamma2 open until the end. Closed by CloseP
  !close(LuGAMMA2)

  call mma_deallocate(B_t)
  call mma_deallocate(QVec)

  call mma_deallocate(AOList)

  call DaClos(Lu_Q)
  call DaClos(LUAPT2)
  close(LuGAMMA)

end do ! iSym

call CWTime(TotCPU1,TotWall1)
!write(u6,*) 'CPU/Wall Time for mult_with_q_caspt2:',totcpu1-totcpu0,totwall1-totwall0

contains

subroutine construct_AOList()

  use Index_Functions, only: iTri, nTri_Elem
  use Gateway_Info, only: CutInt
  use Gateway_global, only: force_part_c
  use iSD_data, only: iSD
  use RI_glob, only: DMLT, iBDsh, nJDens
  use SOAO_Info, only: iAOtSO
  use Constants, only: Zero
  use Definitions, only: wp, iwp

  integer(kind=iwp) :: i3, i4, iiQ, ij, ij_Shell, ijQ, ijS, iloc, indS, iS, ish, iSkal, iSO, iSym_, jjQ, jloc, jS, jsh, jSkal, &
                       kAO, kAOk, kBas, kBsInc, kCmp, klS, kS, kSO, kSOk, lAO, lAOl, lBas, lBsInc, lCmp, lMaxDens, lS, lSO, lSOl, &
                       nij_Shell, nSkal, nSkal_Auxiliary, nSkal_Valence
  real(kind=wp) :: A_int_ij, Dm_ij, ThrAO, TMax_all, XDm_ii, XDm_ij, XDm_jj, XDm_max
  logical(kind=iwp) :: DoFock, DoGrad, Indexation
  integer(kind=iwp), allocatable :: iShij(:,:)
  real(kind=wp), allocatable :: MaxDens(:), TMax_Valence(:,:), Tmp(:,:)
  integer(kind=iwp), external :: Cho_irange
  logical(kind=iwp), external :: Rsv_Tsk2

  call Set_Basis_Mode('Auxiliary')
  call Nr_Shells(nSkal_Auxiliary)
  call Set_Basis_Mode('Valence')
  call Nr_Shells(nSkal_Valence)
  call Set_Basis_Mode('WithAuxiliary')
  !call Nr_Shells(nSkal)
  !nSkal = nSkal_Valence
  call SetUp_iSD()

  Indexation = .false.
  DoFock = .false.
  DoGrad = .true.
  ThrAO = Zero
  call Setup_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)

  !nSkal_Valence = nSkal-nSkal_Auxiliary
  call mma_allocate(TMax_Valence,nSkal_Valence,nSkal_Valence,Label='TMax_Valence')

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
  call mma_deallocate(Tmp)

  ! Calculate maximum density value for each shellpair

  lMaxDens = nTri_Elem(nSkal_Valence)
  call mma_allocate(MaxDens,lMaxDens,Label='MaxDens')
  MaxDens(:) = Zero
  do iSym_=0,nSym-1
    kS = 1+nSkal_Valence*iSym_ ! note diff wrt declaration of iBDsh
    !do j=1,nBas(iSym_)
    do jloc=1,nBas(1)
      jsh = Cho_Irange(jloc,iBDsh(kS),nSkal_Valence,.true.)
      do iloc=1,jloc
        ish = Cho_Irange(iloc,iBDsh(kS),nSkal_Valence,.true.)
        ijS = nTri_Elem(jsh-1)+ish
        do iSO=1,nJDens
          if (.not. DMLT(iSO)%Active) cycle
          ij = iTri(jloc,iloc)
          Dm_ij = abs(DMLT(iSO)%SB(iSym_+1)%A1(ij))
          MaxDens(ijS) = max(MaxDens(ijS),Dm_ij)
        end do
      end do
    end do
  end do

  call qpg_iArray('LBList',Found,nSkal2_)
  call mma_allocate(iShij,2,nSkal2_,Label='iShij')
  ij_Shell = 0
  indS = 0
  do iSkal=1,nSkal_Valence
    iiQ = nTri_Elem(iSkal)
    XDm_ii = MaxDens(iiQ)
    do jSkal=1,iSkal
      indS = indS+1
      jjQ = nTri_Elem(jSkal)
      XDm_jj = MaxDens(jjQ)
      ijQ = iTri(iSkal,jSkal)
      XDm_ij = MaxDens(ijQ)
      XDm_max = max(XDm_ij,XDm_ii,XDm_jj)
      A_int_ij = TMax_Valence(iSkal,jSkal)
      if (A_int_ij*TMax_all >= CutInt) then
        if (A_int_ij*XDm_max >= CutInt) then
          ij_Shell = ij_Shell+1
          iShij(1,ij_Shell) = iSkal
          iShij(2,ij_Shell) = jSkal
#         ifdef _DEBUGPRINT_
          write(u6,*) 'ij_Shell,iSkal,jSkal=',ij_Shell,iSkal,jSkal
#         endif
        end if
      end if
    end do
  end do
  call mma_deallocate(MaxDens)
  call mma_deallocate(TMax_Valence)
  !write(u6,*) 'total pair = ',nTri_Elem(nSkal_Valence)
  !write(u6,*) 'actual pair= ',nij_shell

  !call qpg_iArray('LBList',Found,nSkal2_)
  !write(u6,*) 'nSkal2_ = ',nSkal2_
  call mma_allocate(LBList,nSkal2_,Label='LBList')
  LBList = 0
  call Get_iArray('LBList',LBList,nSkal2_)

  iOpt = 1
  call Init_Tsk2(id,ij_Shell,iOpt,LBList)
  call mma_deallocate(LBList)

  indS = 0
  do while (Rsv_Tsk2(id,klS))

    kS = iShij(1,klS)
    lS = iShij(2,klS)
    indS = indS+1

    kCmp = iSD(2,kS)
    kBas = iSD(3,kS)
    kAO = iSD(7,kS)

    lCmp = iSD(2,lS)
    lBas = iSD(3,lS)
    lAO = iSD(7,lS)

    if (force_part_c) then
      kBsInc = (kBas+1)/2
      lBsInc = (lBas+1)/2
    else
      kBsInc = kBas
      lBsInc = lBas
    end if

    do i4=1,lCmp
      lSO = iAOtSO(lAO+i4,0)
      do i3=1,kCmp
        kSO = iAOtSO(kAO+i3,0)
        do lAOl=0,min(lBsInc,lBas)-1
          lSOl = lSO+lAOl
          do kAOk=0,min(kBsInc,kBas)-1
            kSOk = kSO+kAOk
            nCalAO = nCalAO+1
            AOList(1,nCalAO) = kSOk
            AOList(2,nCalAO) = lSOl
            !write(u6,*) ncalao,ksok,lsol
          end do
        end do
      end do
    end do
  end do

  ! Construct the offset array of shell and AO
  nList_Shell = 0
  nList_Shell(myRank+1) = indS
# ifdef _MOLCAS_MPP_
  if (is_real_par()) call GAIGOP(nList_Shell,NPROCS,'+')
# endif
  kS = 1
  nij_Shell = 0
  do iloc=1,NPROCS
    lS = nList_Shell(iloc)
    nij_Shell = nij_Shell+lS
    nList_Shell(iloc) = kS
    kS = kS+lS
  end do

  nList_AO = 0
  nList_AO(myRank+1) = nCalAO
# ifdef _MOLCAS_MPP_
  if (is_real_par()) call GAIGOP(nList_AO,NPROCS,'+')
# endif
  kS = 1
  nij_Shell = 0
  do iloc=1,NPROCS
    lS = nList_AO(iloc)
    nij_Shell = nij_Shell+lS
    nList_AO(iloc) = kS
    kS = kS+lS
  end do
  nCalAO_tot = kS-1
  iAOstart = nList_AO(myRank+1)
  if (myRank+1 == NPROCS) then
    iAOlast = kS-1
  else
    iAOlast = nList_AO(myRank+2)-1
  end if

  call Free_Tsk2(id)
  call mma_deallocate(iShij)
  call Term_Ints()
  call Free_iSD()

end subroutine construct_AOList

end subroutine Mult_with_Q_CASPT2
