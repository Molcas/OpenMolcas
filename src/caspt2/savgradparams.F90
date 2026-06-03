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
! Copyright (C) 2023, Yoshio Nishimoto                                 *
!***********************************************************************
!
! Save and restore many quantities that are used in CASPT2 gradient
! calculations using disk at present (hopefully)
! Save with Mode = 1, and restore with Mode = 2
! state-dependent quantities

subroutine SavGradParams(Mode,IDSAVGRD)

use Index_Functions, only: nTri_Elem
use EQSOLV, only: IDBMAT, IDSMAT, IDSTMAT, IDTMAT, IVECX
use fake_GA, only: GA_Arrays
use caspt2_global, only: do_lindep, DREF, IDBoriMat, iTasks_grad, LUGRAD, LUSBT, LUSOLV, LUSTD, NBUF1_GRAD, nTasks_grad, PREF
use caspt2_module, only: E2Tot, EASum, ERef, HZERO, jState, MxCase, nAshT, nASup, nBTri, nCases, nG1, nG2, nG3, nG3Tot, nInDep, &
                         nISup, nState, nSym, RefEne, RFPert
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, myRank
use GA_Wrapper, only: DBL_MB, GA_Destroy, GA_NodeId
use caspt2_global, only: LURHS
use caspt2_module, only: cLab10, iAdr10, IOFFRHS
#endif
use SC_NEVPT2, only: IDBMAT_NEVPT2, SC_amplitude, SC_prop
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, byte

implicit none
integer(kind=iwp), intent(in) :: Mode
integer(kind=iwp), intent(inout) :: IDSAVGRD
integer(kind=iwp) :: ICASE, ID, iLUID, IORW, ISYM, NAS, NIN, NIS, NMAX
integer(kind=iwp), allocatable :: IWRK1(:)
integer(kind=byte), allocatable :: idxG3(:,:)
real(kind=wp), allocatable :: WRK1(:)
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: I, IDISK, IEND, ISTA, JEND, JSTA, LDM, LDW, lg_S, lg_ST, lg_T, mV1, NBLOCK
logical(kind=iwp) :: bStat
#endif

!! Shift the address due to SavGradParams2
if (IDSAVGRD == 0) then
  IDSAVGRD = NSTATE+3*NSTATE**2
  if (RFpert) IDSAVGRD = IDSAVGRD+NBTRI+1
end if
#ifdef _MOLCAS_MPP_
I = 0
myRank = GA_NODEID()
#endif

!! Decide what to do
if (Mode == 1) then
  IORW = 1 !! Write
else if (Mode == 2) then
  IORW = 2 !! Read
end if

!! Save internal contractions-related quantities
!! 1. Some integers
!! - Number of independent vectors
call IDAFILE(LUGRAD,IORW,NINDEP(1,1),8*MXCASE,IDSAVGRD)
!! - Number of active indices
call IDAFILE(LUGRAD,IORW,NASUP(1,1),8*MXCASE,IDSAVGRD)
!! - Number of inactive + secondary indices
call IDAFILE(LUGRAD,IORW,NISUP(1,1),8*MXCASE,IDSAVGRD)

!! 2. RDMs
!! - NG1, NG2, NG3
call mma_allocate(IWRK1,6,Label='IWRK1')
if (IORW == 1) then
  IWRK1(1) = NG1
  IWRK1(2) = NG2
  IWRK1(3) = NG3
  IWRK1(4) = NG3TOT
  IWRK1(5) = NBUF1_GRAD
  IWRK1(6) = nTasks_grad
  call IDAFILE(LUGRAD,IORW,IWRK1,6,IDSAVGRD)
else if (IORW == 2) then
  call IDAFILE(LUGRAD,IORW,IWRK1,6,IDSAVGRD)
  NG1 = IWRK1(1)
  NG2 = IWRK1(2)
  NG3 = IWRK1(3)
  NG3TOT = IWRK1(4)
  NBUF1_GRAD = IWRK1(5)
  nTasks_grad = IWRK1(6)
  iTasks_grad(:) = 0
end if
call mma_deallocate(IWRK1)
call IDAFILE(LUGRAD,IORW,iTasks_grad,NASHT**2,IDSAVGRD)

NMAX = NG3

do ISYM=1,NSYM
  do ICASE=1,11
    NIN = NINDEP(ISYM,ICASE)
    NAS = NASUP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    NMAX = max(NMAX,nTri_Elem(NAS),NAS*NIN,NIS)
  end do
  do ICASE=12,13
    NAS = NASUP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    NMAX = max(NMAX,NAS*NIS)
  end do
end do

call mma_allocate(WRK1,NMAX,Label='WRK1')

call mma_allocate(idxG3,6,NG3,label='idxG3')
if (IORW == 1) then
  !! NG3 index
  iLUID = 0
  call I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)
  call I1DAFILE(LUGRAD,1,idxG3,6*NG3,IDSAVGRD)

  !! D and F
  call PT2_GET(NG1,' GAMMA1',WRK1)
  call DDAFILE(LUGRAD,IORW,WRK1,NG1,IDSAVGRD)
  call PT2_GET(NG2,' GAMMA2',WRK1)
  call DDAFILE(LUGRAD,IORW,WRK1,NG2,IDSAVGRD)
  call PT2_GET(NG3,' GAMMA3',WRK1)
  call DDAFILE(LUGRAD,IORW,WRK1,NG3,IDSAVGRD)

  if (HZERO /= 'DYALL') then
    call PT2_GET(NG1,' DELTA1',WRK1)
    call DDAFILE(LUGRAD,IORW,WRK1,NG1,IDSAVGRD)
    call PT2_GET(NG2,' DELTA2',WRK1)
    call DDAFILE(LUGRAD,IORW,WRK1,NG2,IDSAVGRD)
    call PT2_GET(NG3,' DELTA3',WRK1)
    call DDAFILE(LUGRAD,IORW,WRK1,NG3,IDSAVGRD)
  end if

  !! EASUM
  WRK1(1) = EASUM
  call DDAFILE(LUGRAD,IORW,WRK1,1,IDSAVGRD)
else if (IORW == 2) then
# ifdef _MOLCAS_MPP_
  if (is_real_par()) then
    !! Reset, because NG3 can be different
    !! STINI has been skipped
    IADR10(:,1) = -1
    IADR10(:,2) = 0
    CLAB10(:) = '   EMPTY'
    IADR10(1,1) = 0
  end if
# endif
  !! NG3 index
  call I1DAFILE(LUGRAD,2,idxG3,6*NG3,IDSAVGRD)
  iLUID = 0
  call I1DAFILE(LUSOLV,1,idxG3,6*NG3,iLUID)

  !! D and F
  call DDAFILE(LUGRAD,IORW,WRK1,NG1,IDSAVGRD)
  call PT2_PUT(NG1,' GAMMA1',WRK1)
  call DDAFILE(LUGRAD,IORW,WRK1,NG2,IDSAVGRD)
  call PT2_PUT(NG2,' GAMMA2',WRK1)
  call DDAFILE(LUGRAD,IORW,WRK1,NG3,IDSAVGRD)
  call PT2_PUT(NG3,' GAMMA3',WRK1)

  if (HZERO /= 'DYALL') then
    call DDAFILE(LUGRAD,IORW,WRK1,NG1,IDSAVGRD)
    call PT2_PUT(NG1,' DELTA1',WRK1)
    call DDAFILE(LUGRAD,IORW,WRK1,NG2,IDSAVGRD)
    call PT2_PUT(NG2,' DELTA2',WRK1)
    call DDAFILE(LUGRAD,IORW,WRK1,NG3,IDSAVGRD)
    call PT2_PUT(NG3,' DELTA3',WRK1)
  end if

  !! EASUM
  call DDAFILE(LUGRAD,IORW,WRK1,1,IDSAVGRD)
  EASUM = WRK1(1)
  call GETDPREF(DREF,size(DREF),PREF,size(PREF))
  EREF = REFENE(JSTATE)
end if
call mma_deallocate(idxG3)

!! 3. LUSBT matrices
do ISYM=1,NSYM
  do ICASE=1,11
    NIN = NINDEP(ISYM,ICASE)
    NAS = NASUP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    if (IORW == 1) then
      if (NIN > 0) then
        !! Active overlap
#       ifdef _MOLCAS_MPP_
        if (is_real_par() .and. ((icase == 1) .or. (icase == 4))) then
          call PSBMAT_GETMEM('S',lg_S,NAS)
          call PSBMAT_READ('S',iCase,iSym,lg_S,NAS)
          call GA_Distribution(lg_S,myRank,ISTA,IEND,JSTA,JEND)
          if ((ISTA > 0) .and. (JSTA > 0)) then
            call GA_Access(lg_S,ISTA,IEND,JSTA,JEND,mV1,LDM)
            NBLOCK = LDM*(JEND-JSTA+1)
            call DDAFILE(LUGRAD,1,DBL_MB(mV1),NBLOCK,IDSAVGRD)
            call GA_Release(lg_S,ISTA,IEND,JSTA,JEND)
          end if
          call PSBMAT_FREEMEM(lg_S)
        else
#       endif
          ID = IDSMAT(ISYM,ICASE)
          if (ID >= 0) then
            call DDAFILE(LUSBT,2,WRK1,nTri_Elem(NAS),ID)
            call DDAFILE(LUGRAD,1,WRK1,nTri_Elem(NAS),IDSAVGRD)
          end if
#       ifdef _MOLCAS_MPP_
        end if
#       endif
        if (NAS > 0) then
          !! ST matrix
#         ifdef _MOLCAS_MPP_
          if (is_real_par() .and. ((icase == 1) .or. (icase == 4))) then
            call GA_CREATE_STRIPED('H',NAS,NIN,'STMAT',lg_ST)
            call PSBMAT_READ('M',iCase,iSym,lg_ST,NAS*NIN)
            call GA_Distribution(lg_ST,myRank,ISTA,IEND,JSTA,JEND)
            if ((ISTA > 0) .and. (JSTA > 0)) then
              call GA_Access(lg_ST,ISTA,IEND,JSTA,JEND,mV1,LDM)
              NBLOCK = LDM*(JEND-JSTA+1)
              call DDAFILE(LUGRAD,1,DBL_MB(mV1),NBLOCK,IDSAVGRD)
              call GA_Release(lg_ST,ISTA,IEND,JSTA,JEND)
            end if
            bStat = GA_Destroy(lg_ST)
          else
#         endif
            ID = IDSTMAT(ISYM,ICASE)
            call DDAFILE(LUSBT,2,WRK1,NAS*NIN,ID)
            call DDAFILE(LUGRAD,1,WRK1,NAS*NIN,IDSAVGRD)
#         ifdef _MOLCAS_MPP_
          end if
#         endif
          !! Transformation matrix (eigenvector)
#         ifdef _MOLCAS_MPP_
          if (is_real_par() .and. ((icase == 1) .or. (icase == 4))) then
            call GA_CREATE_STRIPED('H',NAS,NIN,'TMAT',lg_T)
            call PSBMAT_READ('T',iCase,iSym,lg_T,NAS*NIN)
            call GA_Distribution(lg_T,myRank,ISTA,IEND,JSTA,JEND)
            if ((ISTA > 0) .and. (JSTA > 0)) then
              call GA_Access(lg_T,ISTA,IEND,JSTA,JEND,mV1,LDM)
              NBLOCK = LDM*(JEND-JSTA+1)
              call DDAFILE(LUGRAD,1,DBL_MB(mV1),NBLOCK,IDSAVGRD)
              call GA_Release(lg_T,ISTA,IEND,JSTA,JEND)
            end if
            bStat = GA_Destroy(lg_T)
          else
#         endif
            ID = IDTMAT(ISYM,ICASE)
            call DDAFILE(LUSBT,2,WRK1,NAS*NIN,ID)
            call DDAFILE(LUGRAD,1,WRK1,NAS*NIN,IDSAVGRD)
#         ifdef _MOLCAS_MPP_
          end if
#         endif
        end if
        !! Eigenvalue
        ID = IDBMAT(ISYM,ICASE)
        call DDAFILE(LUSBT,2,WRK1,NIN,ID)
        call DDAFILE(LUGRAD,1,WRK1,NIN,IDSAVGRD)
        if (SC_amplitude .and. NIS > 0) then
          call DDAFILE(LUSBT,2,WRK1,NIS,ID)
          call DDAFILE(LUGRAD,1,WRK1,NIS,IDSAVGRD)
        end if
      end if
      !! IS
      !call DDAFILE(LUSBT,2,WRK1,NIS,ID)
      !call DDAFILE(LUGRAD,1,WRK1,NIS,IDSAVGRD)
      if (do_lindep .and. (NAS > 0)) then
        ID = IDBoriMat(ISYM,ICASE)
        call DDAFILE(LUSTD,2,WRK1,nTri_Elem(NAS),ID)
        call DDAFILE(LUGRAD,1,WRK1,nTri_Elem(NAS),IDSAVGRD)
      end if
      !! Original B matrix, needed in SC-NEVPT2 gradient
      if (HZERO == 'DYALL' .and. SC_amplitude) then
        ID = IDBMAT_NEVPT2(ISYM,ICASE,1)
        call DDAFILE(LUSBT,2,WRK1,nTri_Elem(NAS),ID)
        call DDAFILE(LUGRAD,1,WRK1,nTri_Elem(NAS),IDSAVGRD)
        ID = IDBMAT_NEVPT2(ISYM,ICASE,2)
        call DDAFILE(LUSBT,2,WRK1,NIS,ID)
        call DDAFILE(LUGRAD,1,WRK1,NIS,IDSAVGRD)
      end if
    else if (IORW == 2) then
      if (NIN > 0) then
        !! Active overlap
#       ifdef _MOLCAS_MPP_
        if (is_real_par() .and. ((icase == 1) .or. (icase == 4))) then
          call PSBMAT_GETMEM('S',lg_S,NAS)
          call GA_Distribution(lg_S,myRank,ISTA,IEND,JSTA,JEND)
          if ((ISTA > 0) .and. (JSTA > 0)) then
            call GA_Access(lg_S,ISTA,IEND,JSTA,JEND,mV1,LDM)
            NBLOCK = LDM*(JEND-JSTA+1)
            call DDAFILE(LUGRAD,2,DBL_MB(mV1),NBLOCK,IDSAVGRD)
            call GA_Release(lg_S,ISTA,IEND,JSTA,JEND)
          end if
          call PSBMAT_WRITE('S',iCase,iSym,lg_S,NAS)
          call PSBMAT_FREEMEM(lg_S)
        else
#       endif
          ID = IDSMAT(ISYM,ICASE)
          if (ID >= 0) then
            call DDAFILE(LUGRAD,2,WRK1,nTri_Elem(NAS),IDSAVGRD)
            call DDAFILE(LUSBT,1,WRK1,nTri_Elem(NAS),ID)
          end if
#       ifdef _MOLCAS_MPP_
        end if
#       endif
        if (NAS > 0) then
          !! ST matrix
#         ifdef _MOLCAS_MPP_
          if (is_real_par() .and. ((icase == 1) .or. (icase == 4))) then
            call GA_CREATE_STRIPED('H',NAS,NIN,'STMAT',lg_ST)
            call GA_Distribution(lg_ST,myRank,ISTA,IEND,JSTA,JEND)
            if ((ISTA > 0) .and. (JSTA > 0)) then
              call GA_Access(lg_ST,ISTA,IEND,JSTA,JEND,mV1,LDM)
              NBLOCK = LDM*(JEND-JSTA+1)
              call DDAFILE(LUGRAD,2,DBL_MB(mV1),NBLOCK,IDSAVGRD)
              call GA_Release(lg_ST,ISTA,IEND,JSTA,JEND)
            end if
            call PSBMAT_WRITE('M',iCase,iSym,lg_ST,NAS*NIN)
            bStat = GA_Destroy(lg_ST)
          else
#         endif
            call DDAFILE(LUGRAD,2,WRK1,NAS*NIN,IDSAVGRD)
            ID = IDSTMAT(ISYM,ICASE)
            call DDAFILE(LUSBT,1,WRK1,NAS*NIN,ID)
#         ifdef _MOLCAS_MPP_
          end if
#         endif
          !! Transformation matrix (eigenvector)
#         ifdef _MOLCAS_MPP_
          if (is_real_par() .and. ((icase == 1) .or. (icase == 4))) then
            call GA_CREATE_STRIPED('H',NAS,NIN,'TMAT',lg_T)
            call GA_Distribution(lg_T,myRank,ISTA,IEND,JSTA,JEND)
            if ((ISTA > 0) .and. (JSTA > 0)) then
              call GA_Access(lg_T,ISTA,IEND,JSTA,JEND,mV1,LDM)
              NBLOCK = LDM*(JEND-JSTA+1)
              call DDAFILE(LUGRAD,2,DBL_MB(mV1),NBLOCK,IDSAVGRD)
              call GA_Release(lg_T,ISTA,IEND,JSTA,JEND)
            end if
            call PSBMAT_WRITE('T',iCase,iSym,lg_T,NAS*NIN)
            bStat = GA_Destroy(lg_T)
#           include "macros.fh"
            unused_var(bStat)
          else
#         endif
            call DDAFILE(LUGRAD,2,WRK1,NAS*NIN,IDSAVGRD)
            ID = IDTMAT(ISYM,ICASE)
            call DDAFILE(LUSBT,1,WRK1,NAS*NIN,ID)
#         ifdef _MOLCAS_MPP_
          end if
#         endif
        end if
        !! Eigenvalue
        call DDAFILE(LUGRAD,2,WRK1,NIN,IDSAVGRD)
        ID = IDBMAT(ISYM,ICASE)
        call DDAFILE(LUSBT,1,WRK1,NIN,ID)
        if (SC_amplitude .and. NIS > 0) then
          call DDAFILE(LUGRAD,2,WRK1,NIS,IDSAVGRD)
          call DDAFILE(LUSBT,1,WRK1,NIS,ID)
        end if
      end if
      !! IS
      !call DDAFILE(LUGRAD,2,WRK1,NIS,IDSAVGRD)
      !call DDAFILE(LUSBT,1,WRK1,NIS,ID)
      if (do_lindep .and. (NAS > 0)) then
        call DDAFILE(LUGRAD,2,WRK1,nTri_Elem(NAS),IDSAVGRD)
        ID = IDBoriMat(ISYM,ICASE)
        call DDAFILE(LUSTD,1,WRK1,nTri_Elem(NAS),ID)
      end if
      !! Original B matrix, needed in SC-NEVPT2 gradient
      if (HZERO == 'DYALL' .and. SC_amplitude) then
        ID = IDBMAT_NEVPT2(ISYM,ICASE,1)
        call DDAFILE(LUGRAD,2,WRK1,nTri_Elem(NAS),IDSAVGRD)
        call DDAFILE(LUSBT,1,WRK1,nTri_Elem(NAS),ID)
        ID = IDBMAT_NEVPT2(ISYM,ICASE,2)
        call DDAFILE(LUGRAD,2,WRK1,NIS,IDSAVGRD)
        call DDAFILE(LUSBT,1,WRK1,NIS,ID)
      end if
    end if
  end do
end do

if (SC_amplitude) then
  do ISYM=1,NSYM
    do ICASE=12,13
      NAS = NASUP(ISYM,ICASE)
      NIS = NISUP(ISYM,ICASE)
      if (IORW == 1) then
        !! Eigenvalue
        ID = IDBMAT(ISYM,ICASE)
        if (NAS > 0) then
          call DDAFILE(LUSBT,2,WRK1,NAS,ID)
          call DDAFILE(LUGRAD,1,WRK1,NAS,IDSAVGRD)
        end if
        !! IS
        if (NIS > 0) then
          call DDAFILE(LUSBT,2,WRK1,NIS,ID)
          call DDAFILE(LUGRAD,1,WRK1,NIS,IDSAVGRD)
        end if
      else if (IORW == 2) then
        !! Eigenvalue
        ID = IDBMAT(ISYM,ICASE)
        if (NAS > 0) then
          call DDAFILE(LUGRAD,2,WRK1,NAS,IDSAVGRD)
          call DDAFILE(LUSBT,1,WRK1,NAS,ID)
        end if
        !! IS
        if (NIS > 0) then
          call DDAFILE(LUGRAD,2,WRK1,NIS,IDSAVGRD)
          call DDAFILE(LUSBT,1,WRK1,NIS,ID)
        end if
      end if
    end do
  end do
end if

!! 4. E2TOT
if (IORW == 1) then
  WRK1(1) = E2TOT
  call DDAFILE(LUGRAD,IORW,WRK1,1,IDSAVGRD)
else if (IORW == 2) then
  call DDAFILE(LUGRAD,IORW,WRK1,1,IDSAVGRD)
  E2TOT = WRK1(1)
end if

call mma_deallocate(WRK1)

!! quasi-canonical active orbital energy
!call DDAFILE(LUGRAD,IORW,EPSA,NASHT,IDSAVGRD)

!! 5. Save T-amplitude (IRHS = LURHS(1) = 51)
!! SC-NEVPT2 does not use IVECX
if (HZERO /= 'DYALL' .or. .not. SC_prop) call SaveReadT1()
!! We at least need to create IVECX for MECI search
if (HZERO == 'DYALL' .and. SC_prop .and. IORW == 2) call RHS_ZERO(IVECX)

contains

subroutine SaveReadT1()

  integer(kind=iwp) :: ICASE_, ISYM_, lg_V1, NVEC

  !! IVECX = T (solution; not quasi-variational, before lambda-eq)
  IVECX = 2

  do ICASE_=1,NCASES
    do ISYM_=1,NSYM
      NIN = NINDEP(ISYM_,ICASE_)
      if (NIN == 0) cycle
      NIS = NISUP(ISYM_,ICASE_)
      NAS = NASUP(ISYM_,ICASE_)
      NVEC = NIN*NIS
      if ((ICASE_ == 12) .or. (ICASE_ == 13)) NVEC = NAS*NIS
      if (NVEC == 0) cycle
      if ((ICASE_ == 12) .or. (ICASE_ == 13)) then
        call RHS_ALLO(NAS,NIS,lg_V1)
        if (IORW == 1) then
#         ifdef _MOLCAS_MPP_
          if (is_real_par()) then
            call GA_Distribution(lg_V1,myRank,ISTA,IEND,JSTA,JEND)
            if ((IEND-ISTA+1 == NAS) .and. (ISTA > 0)) then
              call GA_Access(lg_V1,ISTA,IEND,JSTA,JEND,mV1,LDW)
              NVEC = (IEND-ISTA+1)*(JEND-JSTA+1)
              IDISK = IOFFRHS(ISYM_,ICASE_)
              call DDAFILE(LURHS(IVECX),2,DBL_MB(mV1),NVEC,IDISK)
              call DDAFILE(LUGRAD,IORW,DBL_MB(mV1),NVEC,IDSAVGRD)
              call GA_Release(lg_V1,ISTA,IEND,JSTA,JEND)
            end if
          else
#         endif
            if (NAS*NIS > 0) then
              call RHS_READ_SR(lg_V1,ICASE_,ISYM_,IVECX)
              call DDAFILE(LUGRAD,IORW,GA_Arrays(lg_V1)%A,NAS*NIS,IDSAVGRD)
            end if
#         ifdef _MOLCAS_MPP_
          end if
#         endif
        else if (IORW == 2) then
#         ifdef _MOLCAS_MPP_
          if (is_real_par()) then
            call GA_Distribution(lg_V1,myRank,ISTA,IEND,JSTA,JEND)
            if ((IEND-ISTA+1 == NAS) .and. (ISTA > 0)) then
              call GA_Access(lg_V1,ISTA,IEND,JSTA,JEND,mV1,LDW)
              NVEC = (IEND-ISTA+1)*(JEND-JSTA+1)
              call DDAFILE(LUGRAD,IORW,DBL_MB(mV1),NVEC,IDSAVGRD)
              IDISK = IOFFRHS(ISYM_,ICASE_)
              call DDAFILE(LURHS(IVECX),1,DBL_MB(mV1),NVEC,IDISK)
              call GA_Release(lg_V1,ISTA,IEND,JSTA,JEND)
            end if
          else
#         endif
            if (NAS*NIS > 0) then
              call DDAFILE(LUGRAD,IORW,GA_Arrays(lg_V1)%A,NAS*NIS,IDSAVGRD)
              call RHS_SAVE(NIN,NIS,lg_V1,ICASE_,ISYM_,IVECX)
            end if
#         ifdef _MOLCAS_MPP_
          end if
#         endif
        end if
        call RHS_FREE(lg_V1)
      else
        call RHS_ALLO(NIN,NIS,lg_V1)
        if (IORW == 1) then
#         ifdef _MOLCAS_MPP_
          if (is_real_par()) then
            call GA_Distribution(lg_V1,myRank,ISTA,IEND,JSTA,JEND)
            if ((IEND-ISTA+1 == NIN) .and. (ISTA > 0)) then
              call GA_Access(lg_V1,ISTA,IEND,JSTA,JEND,mV1,LDW)
              NVEC = (IEND-ISTA+1)*(JEND-JSTA+1)
              IDISK = IOFFRHS(ISYM_,ICASE_)
              call DDAFILE(LURHS(IVECX),2,DBL_MB(mV1),NVEC,IDISK)
              call DDAFILE(LUGRAD,IORW,DBL_MB(mV1),NVEC,IDSAVGRD)
              call GA_Release(lg_V1,ISTA,IEND,JSTA,JEND)
            end if
          else
#         endif
            if (NIN*NIS > 0) then
              call RHS_READ_SR(lg_V1,ICASE_,ISYM_,IVECX)
              call DDAFILE(LUGRAD,IORW,GA_Arrays(lg_V1)%A,NIN*NIS,IDSAVGRD)
            end if
#         ifdef _MOLCAS_MPP_
          end if
#         endif
        else if (IORW == 2) then
#         ifdef _MOLCAS_MPP_
          if (is_real_par()) then
            call GA_Distribution(lg_V1,myRank,ISTA,IEND,JSTA,JEND)
            if ((IEND-ISTA+1 == NIN) .and. (ISTA > 0)) then
              call GA_Access(lg_V1,ISTA,IEND,JSTA,JEND,mV1,LDW)
              NVEC = (IEND-ISTA+1)*(JEND-JSTA+1)
              call DDAFILE(LUGRAD,IORW,DBL_MB(mV1),NVEC,IDSAVGRD)
              IDISK = IOFFRHS(ISYM_,ICASE_)
              call DDAFILE(LURHS(IVECX),1,DBL_MB(mV1),NVEC,IDISK)
              call GA_Release(lg_V1,ISTA,IEND,JSTA,JEND)
            end if
          else
#         endif
            if (NIN*NIS > 0) then
              call DDAFILE(LUGRAD,IORW,GA_Arrays(lg_V1)%A,NIN*NIS,IDSAVGRD)
              call RHS_SAVE(NIN,NIS,lg_V1,ICASE_,ISYM_,IVECX)
            end if
#         ifdef _MOLCAS_MPP_
          end if
#         endif
        end if
        call RHS_FREE(lg_V1)
      end if
    end do
  end do

end subroutine SaveReadT1

end subroutine SavGradParams
