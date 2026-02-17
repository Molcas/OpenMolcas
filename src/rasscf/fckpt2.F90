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
! Copyright (C) 1990, Bjorn O. Roos                                    *
!***********************************************************************

subroutine FCKPT2(CMOO,CMON,FI,FP,FTR,VEC,WO,SQ,CMOX)
! Purpose: To diagonalize the inactive, active,
! and external parts of the Fock matrix (FI+FA) for CASPT2
! and order the eigenvalues and eigenvectors after energy.
! This ordering is of value in subsequent CI and CASPT2
! calculations. All diagonal elements of
! the transformed Fock matrix are collected in FDIAG for
! later printing. The new MO's are written onto JOBIPH
! in address IADR15(9). Also the inactive Fock matrix FI
! is saved on Jobiph for use in CASPT2.
! Note:these orbitals leave the CI expansion
! invariant only for CAS wave functions.
! Called from SXCTL if IFINAL=1 (after last MC iteration)
!
! Modifications: - Produce full matrix FP and write it to JOBIPH.
! B. Roos, Lund, June 1990
!
! ********** IBM-3090 Release 88 09 07 **********

use Index_Functions, only: iTri, nTri_Elem
use rasscf_global, only: NORBT, NTOT3, FDIAG, ixSym, IADR15
use PrintLevel, only: DEBUG, VERBOSE
use output_ras, only: IPRLOC
use general_data, only: JOBIPH, NASH, NBAS, NDEL, NFRO, NISH, NRS1, NRS2, NRS3, NSSH, NSYM, NTOT, NTOT2
#ifdef _HDF5_
use mh5, only: mh5_close_dset, mh5_create_dset_int, mh5_create_dset_real, mh5_create_file, mh5_put_dset
use fciqmc, only: tNonDiagStochPT2, tPrepStochCASPT2
use RASWfn, only: wfn_mocoef
use stdalloc, only: mma_allocate, mma_deallocate
#endif
#ifdef _ENABLE_CHEMPS2_DMRG_
use rasscf_global, only: NAC
#endif
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: CMOO(*), CMON(*), FI(*), FP(*), FTR(*), VEC(*), WO(*), SQ(*), CMOX(*)
integer(kind=iwp) :: i, I_F, iAd15, IB, iBas, ID, IFD, II, ioff, iPrLev, IST, ISTFCK, ISTMO, ISTMO1, iSym, j, M_IN, N_OT, NA, NA1, &
                     NAB, NABT, NAO, NAT, NB, NBF, NBT, NDNB, NDO, NEO, NEO1, NFNB, NFO, NI, NI1, NIJ, NIO, NIO1, NJ, NO1, NOC, &
                     NOO, NP, NPQ, NQ, NR1, NR11, NR2, NR21, NR3, NR31, NT, NT1, NTT, NTU, NTUT, NU, NUT
real(kind=wp) :: FMIN
#ifdef _ENABLE_CHEMPS2_DMRG_
integer(kind=iwp) :: iChMolpro(8), ifock, iiash, iOrb, jOrb, LuFck, nOrbTot
character(len=3) :: Label
integer(kind=iwp), allocatable :: OrbSym(:)
integer(kind=iwp), external :: IsFreeUnit
#endif
#ifdef _HDF5_
integer(kind=iwp) :: dset_id, file_id, idx, k, nActOrb, nOrbCount, offset
integer(kind=iwp), allocatable :: indices(:,:)
real(kind=wp), allocatable :: fockmat(:,:), vals(:), vecs(:,:)
#endif

! Local print level (if any)
IPRLEV = IPRLOC(4)
if (IPRLEV >= DEBUG) write(u6,*) ' Entering FCKPT2'

IB = 0
ISTMO1 = 1
ISTFCK = 0
ID = 0

#ifdef _HDF5_
if (tPrepStochCASPT2 .or. tNonDiagStochPT2) then
  nActOrb = 0
  do isym=1,nsym
    nActOrb = nActOrb+nAsh(isym)
  end do
  if (tPrepStochCASPT2) then
    call mma_allocate(indices,2,nActOrb)
    call mma_allocate(vals,nActOrb)
  else
    call mma_allocate(indices,2,nTri_Elem(nActOrb))
    call mma_allocate(vals,nTri_Elem(nActOrb))
    call mma_allocate(fockmat,nActOrb,nActOrb)
    call mma_allocate(vecs,nActOrb,nActOrb)
    fockmat(:,:) = Zero
    vecs(:,:) = Zero
  end if
  indices(:,:) = 0
  vals(:) = Zero
  nOrbCount = 0  ! keeps track of indices over different irreps
end if
#endif

#ifdef _ENABLE_CHEMPS2_DMRG_
ifock = 1
norbtot = 0
do iiash=1,nsym
  norbtot = norbtot+nAsh(iiash)
end do

! Get character table to convert MOLPRO symmetry format
call MOLPRO_ChTab(nSym,Label,iChMolpro)

! Convert orbital symmetry into MOLPRO format
call mma_allocate(OrbSym,NAC,Label='OrbSym')
iOrb = 1
do iSym=1,nSym
  do jOrb=1,NASH(iSym)
    OrbSym(iOrb) = iChMolpro(iSym)
    iOrb = iOrb+1
  end do
end do

LuFCK = isFreeUnit(27)
!open(unit=LuFCK,file='FOCK_CHEMPS2',action='write' status='replace')
call molcas_open(LuFCK,'FOCK_CHEMPS2')
write(LuFCK,'(1X,A12,I2,A1)') '&FOCK NACT= ',norbtot,','
write(LuFCK,'(2X,A7)',advance='NO') 'ORBSYM='
do iOrb=1,norbtot
  write(LuFCK,'(I1,A1)',advance='NO') OrbSym(iOrb),','
end do
write(LuFCK,*)
write(LuFCK,*) '/'
call mma_deallocate(OrbSym)
#endif

do ISYM=1,NSYM
  NBF = NBAS(ISYM)
  NFO = NFRO(ISYM)
  NIO = NISH(ISYM)
  NAO = NASH(ISYM)
  NR1 = NRS1(ISYM)
  NR2 = NRS2(ISYM)
  NR3 = NRS3(ISYM)
  NEO = NSSH(ISYM)
  NOO = NFO+NIO+NAO
  N_OT = NIO+NAO+NEO
  NOC = NIO+NAO
  ISTMO = ISTMO1+NFO*NBF
  !*********************************************************************
  ! Frozen orbitals (move MO's to CMON and set FDIAG to zero)
  !*********************************************************************
  if (NFO /= 0) then
    NFNB = NBF*NFO
    CMON(ISTMO1:ISTMO1+NFNB-1) = CMOO(ISTMO1:ISTMO1+NFNB-1)
    FDIAG(IB+1:IB+NFO) = Zero
  end if

  ! Clear the MO transformation matrix CMOX

  CMOX(1:N_OT**2) = Zero

  !*********************************************************************
  ! Inactive part of the Fock matrix
  !*********************************************************************

  if (NIO /= 0) then
    ! MOVE FP TO TRIANGULAR FORM
    NIJ = 0
    do NI=1,NIO
      do NJ=1,NI
        NIJ = NIJ+1
        FTR(NIJ) = FP(NIJ+ISTFCK)
        if (IXSYM(IB+NFO+NI) /= IXSYM(IB+NFO+NJ)) FTR(NIJ) = Zero
      end do
    end do
    ! DIAGONALIZE
    call unitmat(VEC,NIO)
    call Jacob(FTR,VEC,NIO,NIO)
    ! MOVE EIGENVALUES TO FDIAG.

    II = 0
    NO1 = IB+NFO
    do NI=1,NIO
      II = II+NI
      FDIAG(NO1+NI) = FTR(II)
    end do

    ! Sort eigenvalues and orbitals after energy

    if (NIO > 1) then
      NIO1 = NIO-1
      do NI=1,NIO1
        NI1 = NI+1
        M_IN = NI
        do NJ=NI1,NIO
          if (FDIAG(NO1+NJ) < FDIAG(NO1+M_IN)) M_IN = NJ
        end do
        if (M_IN == NI) GO TO 20
        FMIN = FDIAG(NO1+M_IN)
        FDIAG(NO1+M_IN) = FDIAG(NO1+NI)
        FDIAG(NO1+NI) = FMIN
        call DSWAP_(NIO,VEC(1+NIO*(NI-1)),1,VEC(1+NIO*(M_IN-1)),1)
20      continue
      end do
    end if
    call DGEACC(One,VEC,NIO,'N',CMOX,N_OT,NIO,NIO)
  end if

  !*********************************************************************
  ! Active part of the Fock matrix
  !*********************************************************************

  if (NAO /= 0) then
    !*******************************************************************
    ! RAS1 part of the Fock matrix
    !*******************************************************************
    if (NR1 /= 0) then
      ! MOVE FP TO TRIANGULAR FORM
      NTU = 0
      do NT=1,NR1
        do NU=1,NT
          NTU = NTU+1
          NTT = NT+NIO
          NUT = NU+NIO
          NTUT = ISTFCK+iTri(NTT,NUT)
          FTR(NTU) = FP(NTUT)
          if (IXSYM(IB+NFO+NTT) /= IXSYM(IB+NFO+NUT)) FTR(NTU) = Zero
        end do
      end do
      ! DIAGONALIZE
      call unitmat(VEC,NR1)
      call Jacob(FTR,VEC,NR1,NR1)

      ! Move eigenvalues to FDIAG.

      II = 0
      NO1 = IB+NFO+NIO
      do NT=1,NR1
        II = II+NT
        FDIAG(NO1+NT) = FTR(II)
      end do

      ! Sort eigenvalues and orbitals after energy

      if (NR1 > 1) then
        NR11 = NR1-1
        do NT=1,NR11
          NT1 = NT+1
          M_IN = NT
          do NU=NT1,NR1
            if (FDIAG(NO1+NU) < FDIAG(NO1+M_IN)) M_IN = NU
          end do
          if (M_IN == NT) GO TO 41
          FMIN = FDIAG(NO1+M_IN)
          FDIAG(NO1+M_IN) = FDIAG(NO1+NT)
          FDIAG(NO1+NT) = FMIN
          call DSWAP_(NR1,VEC(1+NR1*(NT-1)),1,VEC(1+NR1*(M_IN-1)),1)
41        continue
        end do
      end if
      call DGEACC(One,VEC,NR1,'N',CMOX(1+N_OT*NIO+NIO),N_OT,NR1,NR1)

#     ifdef _ENABLE_CHEMPS2_DMRG_
      II = 0
      NO1 = IB+NFO+NIO
      do NT=1,NR1
        write(LuFCK,'(1X,ES23.16E2,I4,I4)') FDIAG(NO1+NT),ifock,ifock
        ifock = ifock+1
      end do
#     endif
    end if ! NR1

    !*******************************************************************
    ! RAS2 part of the Fock matrix
    !*******************************************************************
    if (NR2 /= 0) then
      ! MOVE FP TO TRIANGULAR FORM
      NTU = 0
      do NT=1,NR2
        do NU=1,NT
          NTU = NTU+1
          NTT = NT+NIO+NR1
          NUT = NU+NIO+NR1
          NTUT = ISTFCK+iTri(NTT,NUT)
          ! decoupling test of virtual orbitals
          ! if ((NT > 12) .and. (NU < 13)) FP(NTUT) = Zero
          ! write(u6,*) "t, u, F(t,u)", NT, NU, FP(NTUT)
#         ifdef _HDF5_
          if (tNonDiagStochPT2) then
            if (iprlev >= debug) write(u6,*) 'fock(t,u)',NT+nOrbCount,NU+nOrbCount,FP(NTUT)
            fockmat(NT+nOrbCount,NU+nOrbCount) = FP(NTUT)
          end if
#         endif
          FTR(NTU) = FP(NTUT)
          if (IXSYM(IB+NFO+NTT) /= IXSYM(IB+NFO+NUT)) FTR(NTU) = Zero
        end do
      end do

      ! DIAGONALIZE
      call unitmat(VEC,NR2)
      call Jacob(FTR,VEC,NR2,NR2)

      ! Move eigenvalues to FDIAG.

      II = 0
      NO1 = IB+NFO+NIO+NR1
      do NT=1,NR2
        II = II+NT
        FDIAG(NO1+NT) = FTR(II)
      end do

      ! Sort eigenvalues and orbitals after energy

      if (NR2 > 1) then
        NR21 = NR2-1
        do NT=1,NR21
          NT1 = NT+1
          M_IN = NT
          do NU=NT1,NR2
            if (FDIAG(NO1+NU) < FDIAG(NO1+M_IN)) M_IN = NU
          end do
          if (M_IN == NT) GO TO 42
          FMIN = FDIAG(NO1+M_IN)
          FDIAG(NO1+M_IN) = FDIAG(NO1+NT)
          FDIAG(NO1+NT) = FMIN
          call DSWAP_(NR2,VEC(1+NR2*(NT-1)),1,VEC(1+NR2*(M_IN-1)),1)
42        continue
        end do
      end if
      call DGEACC(One,VEC,NR2,'N',CMOX(1+N_OT*(NIO+NR1)+NIO+NR1),N_OT,NR2,NR2)

#     ifdef _HDF5_
      if (tNonDiagStochPT2) then
        ! grab the eigenvectors of the Fock matrix as well
        do i=1,NR2**2
          j = modulo((i-1),NR2)
          k = modulo((i-1-j)/NR2,NR2)
          vecs(j+1+nOrbCount,k+1+nOrbCount) = vec(i)
        end do
        nOrbCount = nOrbCount+NR2
      end if
      if (tPrepStochCASPT2) then
        do i=1,NR2
          offset = IB+NFO+NIO+NR1  ! for frozen, inactive, RAS1
          indices(:,i+nOrbCount) = i+nOrbCount
          vals(i+nOrbCount) = fdiag(offset+i)
        end do
        ! increment by number of RAS2 orbs in this irrep
        nOrbCount = nOrbCount+NR2
      end if
#     endif

#     ifdef _ENABLE_CHEMPS2_DMRG_
      II = 0
      NO1 = IB+NFO+NIO+NR1
      do NT=1,NR2
        write(LuFCK,'(1X,ES23.16E2,I4,I4)') FDIAG(NO1+NT),ifock,ifock
        ifock = ifock+1
      end do
#     endif
    end if ! NR2

    !*******************************************************************
    ! RAS3 part of the Fock matrix
    !*******************************************************************
    if (NR3 /= 0) then
      ! MOVE FP TO TRIANGULAR FORM
      NTU = 0
      do NT=1,NR3
        do NU=1,NT
          NTU = NTU+1
          NTT = NT+NIO+NR1+NR2
          NUT = NU+NIO+NR1+NR2
          NTUT = ISTFCK+iTri(NTT,NUT)
          FTR(NTU) = FP(NTUT)
          if (IXSYM(IB+NFO+NTT) /= IXSYM(IB+NFO+NUT)) FTR(NTU) = Zero
        end do
      end do
      ! DIAGONALIZE
      call unitmat(VEC,NR3)
      call Jacob(FTR,VEC,NR3,NR3)

      ! Move eigenvalues to FDIAG.

      II = 0
      NO1 = IB+NFO+NIO+NR1+NR2
      do NT=1,NR3
        II = II+NT
        FDIAG(NO1+NT) = FTR(II)
      end do

      ! Sort eigenvalues and orbitals after energy

      if (NR3 > 1) then
        NR31 = NR3-1
        do NT=1,NR31
          NT1 = NT+1
          M_IN = NT
          do NU=NT1,NR3
            if (FDIAG(NO1+NU) < FDIAG(NO1+M_IN)) M_IN = NU
          end do
          if (M_IN == NT) GO TO 43
          FMIN = FDIAG(NO1+M_IN)
          FDIAG(NO1+M_IN) = FDIAG(NO1+NT)
          FDIAG(NO1+NT) = FMIN
          call DSWAP_(NR3,VEC(1+NR3*(NT-1)),1,VEC(1+NR3*(M_IN-1)),1)
43        continue
        end do
      end if
      call DGEACC(One,VEC,NR3,'N',CMOX(1+N_OT*(NIO+NR1+NR2)+NIO+NR1+NR2),N_OT,NR3,NR3)

#     ifdef _ENABLE_CHEMPS2_DMRG_
      II = 0
      NO1 = IB+NFO+NIO+NR1+NR2
      do NT=1,NR3
        write(LuFCK,'(1X,ES23.16E2,I4,I4)') FDIAG(NO1+NT),ifock,ifock
        ifock = ifock+1
      end do
#     endif
    end if ! NR3

  end if ! NAO

  !*********************************************************************
  ! external part of the Fock matrix
  !*********************************************************************
  if (NEO /= 0) then
    ! MOVE FP TO TRIANGULAR FORM
    NAB = 0
    do NA=1,NEO
      do NB=1,NA
        NAB = NAB+1
        NAT = NA+NIO+NAO
        NBT = NB+NIO+NAO
        NABT = ISTFCK+iTri(NAT,NBT)
        FTR(NAB) = FP(NABT)
        if (IXSYM(IB+NFO+NAT) /= IXSYM(IB+NFO+NBT)) FTR(NAB) = Zero
      end do
    end do
    ! DIAGONALIZE
    call unitmat(VEC,NEO)
    call Jacob(FTR,VEC,NEO,NEO)

    ! Move eigenvalues to FDIAG.

    II = 0
    NO1 = IB+NFO+NIO+NAO
    do NA=1,NEO
      II = II+NA
      FDIAG(NO1+NA) = FTR(II)
    end do

    ! Sort eigenvalues and orbitals after energy

    if (NEO > 1) then
      NEO1 = NEO-1
      do NA=1,NEO1
        NA1 = NA+1
        M_IN = NA
        do NB=NA1,NEO
          if (FDIAG(NO1+NB) < FDIAG(NO1+M_IN)) M_IN = NB
        end do
        if (M_IN == NA) GO TO 60
        FMIN = FDIAG(NO1+M_IN)
        FDIAG(NO1+M_IN) = FDIAG(NO1+NA)
        FDIAG(NO1+NA) = FMIN
        call DSWAP_(NEO,VEC(1+NEO*(NA-1)),1,VEC(1+NEO*(M_IN-1)),1)
60      continue
      end do
    end if
    call DGEACC(One,VEC,NEO,'N',CMOX(1+N_OT*NOC+NOC),N_OT,NEO,NEO)
  end if

  ! Transform molecular orbitals

  if (NBF*NTOT > 0) call DGEMM_('N','N',NBF,N_OT,N_OT,One,CMOO(ISTMO),NBF,CMOX,N_OT,Zero,CMON(ISTMO),NBF)

  !*********************************************************************
  ! Deleted orbitals (move MO's and set zero to FDIAG)
  !*********************************************************************
  NDO = NDEL(ISYM)
  if (NDO /= 0) then
    NDNB = NDO*NBF
    IST = ISTMO1+NBF*(NOO+NEO)
    CMON(IST:IST+NDNB-1) = CMOO(IST:IST+NDNB-1)
    FDIAG(IB+NBF-NDO+1:IB+NBF) = Zero
  end if

  ! Transform inactive Fock matrix FI and the CASPT2 matrix FP

  if (N_OT > 0) then
    call SQUARE(FI(ISTFCK+1),SQ,1,N_OT,N_OT)
    call DGEMM_('N','N',N_OT,N_OT,N_OT,One,SQ,N_OT,CMOX,N_OT,Zero,VEC,N_OT)
    call DGEMM_('T','N',N_OT,N_OT,N_OT,One,CMOX,N_OT,VEC,N_OT,Zero,SQ,N_OT)

    ! Move transformed Fock matrix back to FI

    NPQ = ISTFCK
    do NP=1,N_OT
      do NQ=1,NP
        NPQ = NPQ+1
        FI(NPQ) = SQ(N_OT*(NP-1)+NQ)
      end do
    end do

    ! The FP matrix

    call SQUARE(FP(ISTFCK+1),SQ,1,N_OT,N_OT)
    call DGEMM_('N','N',N_OT,N_OT,N_OT,One,SQ,N_OT,CMOX,N_OT,Zero,VEC,N_OT)
    call DGEMM_('T','N',N_OT,N_OT,N_OT,One,CMOX,N_OT,VEC,N_OT,Zero,SQ,N_OT)

    ! Move transformed Fock matrix back to FP

    NPQ = ISTFCK
    do NP=1,N_OT
      do NQ=1,NP
        NPQ = NPQ+1
        FP(NPQ) = SQ(N_OT*(NP-1)+NQ)
      end do
    end do
  end if

  IB = IB+NBF
  ISTFCK = ISTFCK+nTri_Elem(N_OT)
  ISTMO1 = ISTMO1+NBF**2
  ID = ID+nTri_Elem(NAO)
end do

#ifdef _HDF5_
if (tPrepStochCASPT2 .or. tNonDiagStochPT2) then
  file_id = mh5_create_file('fockdump.h5')
  if (tNonDiagStochPT2) then  ! linearise quadratic Fock matrix
    do i=1,nActOrb
      do j=1,i
        idx = iTri(i,j)
        indices(1,idx) = i
        indices(2,idx) = j
        vals(idx) = fockmat(i,j)
      end do
    end do
    call mma_deallocate(fockmat)
    dset_id = mh5_create_dset_real(file_id,'ACT_FOCK_EIGVECS',2,[nActOrb,nActOrb])
    call mh5_put_dset(dset_id,vecs)
    call mh5_close_dset(dset_id)
    call mma_deallocate(vecs)
  end if
  dset_id = mh5_create_dset_int(file_id,'ACT_FOCK_INDEX',2,[2,size(vals)])
  call mh5_put_dset(dset_id,indices)
  call mh5_close_dset(dset_id)
  dset_id = mh5_create_dset_real(file_id,'ACT_FOCK_VALUES',1,[size(vals)])
  call mh5_put_dset(dset_id,vals)
  call mh5_close_dset(dset_id)
  call mma_deallocate(indices)
  call mma_deallocate(vals)
  if (tPrepStochCASPT2) write(u6,*) 'Diagonal active Fock matrix dumped.'
  if (tNonDiagStochPT2) write(u6,*) 'Non-diagonal active Fock matrix dumped.'
end if
#endif

#ifdef _ENABLE_CHEMPS2_DMRG_
close(LuFCK)
#endif
if (IPRLEV >= VERBOSE) then
  write(u6,*) ' Diagonal elements of the Fock matrix in FCKPT2:'
  write(u6,'(1X,10F11.6)') (FDIAG(I),I=1,NTOT)
end if

!***********************************************************************
! Orthogonalise new orbitals
!***********************************************************************

call SUPSCH(WO,CMOO,CMON)
call ORTHO_RASSCF(WO,CMOX,CMON,SQ)

!***********************************************************************
! Write new orbitals to JOBIPH/rasscf.h5
!***********************************************************************

if (IPRLEV >= DEBUG) then
  write(u6,*)
  write(u6,*) ' CMO in FCKPT2 after diag and orthog'
  write(u6,*) ' ---------------------'
  write(u6,*)
  ioff = 0
  do iSym=1,nSym
    iBas = nBas(iSym)
    if (iBas /= 0) then
      write(u6,*) 'Sym =',iSym
      do i=1,iBas
        write(u6,*) (CMON(ioff+iBas*(i-1)+j),j=1,iBas)
      end do
      iOff = iOff+(iBas*iBas)
    end if
  end do
end if

IAD15 = IADR15(9)
call DDAFILE(JOBIPH,1,CMON,NTOT2,IAD15)
#ifdef _HDF5_
call mh5_put_dset(wfn_mocoef,CMON)
#endif

! Write FI, FP and FDIAG to JOBIPH
! First remove frozen and deleted part of FDIAG

I_F = 0
IFD = 0
do ISYM=1,NSYM
  NBF = NBAS(ISYM)
  do NB=1,NBF
    IFD = IFD+1
    if ((NB > NFRO(ISYM)) .and. (NB <= NBF-NDEL(ISYM))) then
      I_F = I_F+1
      SQ(I_F) = FDIAG(IFD)
    end if
  end do
end do

IAD15 = IADR15(10)
call DDAFILE(JOBIPH,1,FI,NTOT3,IAD15)
call DDAFILE(JOBIPH,1,FP,NTOT3,IAD15)
call DDAFILE(JOBIPH,1,SQ,NORBT,IAD15)

end subroutine FCKPT2
