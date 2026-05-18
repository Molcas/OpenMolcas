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
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine SIGDER(IVEC,JVEC,SCAL)

use Fockof, only: FAI, FAI_Full, FAT, FAT_Full, FIA, FIA_Full, FIT, FIT_Full, FTA, FTA_Full, FTI, FTI_Full, IOFFIA, IOFFIT, IOFFTA
use EQSOLV, only: IFCOUP
#if defined(_MOLCAS_MPP_) && defined(_GA_)
use Para_Info, only: Is_Real_Par
#endif
use fake_GA, only: Allocate_GA_Array, Deallocate_GA_Array, GA_Arrays
use caspt2_global, only: FIFA, idSDMat, LISTS, LUSTD
use caspt2_module, only: CPUSGM, FockType, G1SecIn, nActEl, nAsh, nASup, nCases, nInDep, nIsh, nISup, nOrb, nSsh, nSym, TIOSGM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: iVec, jVec
real(kind=wp), intent(in) :: Scal
integer(kind=iwp) :: iCase, iCase1, idSDer, IfC, IFIFA, iLoop, IMLTOP, ISYM, ISYM1, ISYM2, lCX, lg_CX, lg_D2, lg_Sgm2, Lg_SgmX, &
                     lg_V1, lSgmX, lSgmX_Sta, Max_MESG_SIZE, NA, NAS, NAS1, NAS2, nCX, nD1, NFIA, NFIT, NFTA, NI, NIN1, NIN2, &
                     NIS1, NIS2, nLoop, NO, NS, nSgm1, nSgm2, lSgm2_Sta, nD2, iCase2, nSgm2_Blk, nSgmX, nSgmX_Blk, MaxLen
real(kind=wp) :: CPU, CPU0, CPU1, FACT, TIO, TIO0, TIO1
real(kind=wp), allocatable :: D1(:), D2(:), SDER1(:), SDER2(:), SGM1(:), SGM2(:), WRK(:)
#if defined(_MOLCAS_MPP_) && defined(_GA_)
logical(kind=iwp) :: bStat
#include "global.fh"
#endif

! Work in the MO basis
! We need both explicit and implicit overlap derivatives. The latter
! comes from the derivative of the transformation matrix.

! p,q: inactive or secondary
! y,z: active (t,u)
! a,b: internally contracted
! T1_{px} * S1_{xy} f_{yz} * T2_{pz}
! = T1pa*C1xa * S1xy * fyz * T2pb*C2zb
! Derivative of S1:
! = (T1Ct1)px * (T2Ct2*f)py * dS1xy/da
! Derivative of C1 (or, Lagrangian multiplier in MO basis):
! = T1pa*dC1xa/da * S1xy * fyz * T2pb*C2zb
!   ...
! = -1/2 (T1Ct1)pu * dS1tu/da * (T2Ct2*f*S1C1*Ct1)pt
! Derivative of C2 (or, Lagrangian multiplier in MO basis):
! = T1pa*C1xa * S1xy * fyz * T2pb*dC2zb/da
!   ...
! = -1/2 (T1Ct1St1*f*C2*Ct2)pt * (T2Ct2)pu * dS2tu/da

! About IMLTOP for the SGM subroutine
! With IMLTOP=0: the vector for the second argument has to be
! contravariant form (T*C),
! With IMLTOP=1: the vector for the first  argument has to be
! covariant form (T*SC),

! Allocate some matrices for storing overlap and transformation
! derivatives. Here constructs these derivatives in the MO basis,
! but not in the internally contracted basis.

MaxLen = 0
do iCase=1,11
  do iSym=1,nSym
    nAS = nASUP(iSym,iCase)
    MaxLen = max(MaxLen,nAS*nAS)
  end do
end do

call mma_allocate(WRK,MaxLen,Label='WRK')
call DCopy_(MaxLen,[Zero],0,WRK,1)

do iCase=1,11
  do iSym=1,nSym
    nAS = nASUP(iSym,iCase)
    idSDer = idSDMat(iSym,iCase)
    call DDAFILE(LuSTD,1,WRK,nAS*nAS,idSDer)
  end do
end do
call mma_deallocate(WRK)

! If the G1 correction to the Fock matrix is used, then the
! inactive/virtual coupling elements (which are non-zero for the
! case of average CASSCF) cannot be used in the CASPT2 equations.
if ((FOCKTYPE == 'G1') .and. (.not. G1SECIN)) then
  IFCOUP(12,5) = 0
  IFCOUP(13,5) = 0
end if

! Transform to standard representation:
call PTRTOC(0,IVEC,IVEC) !! T*C (internally contracted -> MO)
if (IVEC /= JVEC) call PTRTOC(0,JVEC,JVEC)

! Set up non-diagonal blocks of Fock matrix:
! SVC: add transposed fock matrix blocks
NFIT = 0
NFIA = 0
NFTA = 0
do ISYM=1,NSYM
  NI = NISH(ISYM)
  NA = NASH(ISYM)
  NS = NSSH(ISYM)
  IOFFIT(ISYM) = NFIT
  IOFFIA(ISYM) = NFIA
  IOFFTA(ISYM) = NFTA
  NFIT = NFIT+NA*NI
  NFIA = NFIA+NS*NI
  NFTA = NFTA+NS*NA
end do
NFIT = NFIT+1
NFIA = NFIA+1
NFTA = NFTA+1

call mma_allocate(FIT_Full,NFIT,Label='FIT_Full')
call mma_allocate(FTI_Full,NFIT,Label='FTI_Full')

call mma_allocate(FIA_Full,NFIA,Label='FIA_Full')
call mma_allocate(FAI_Full,NFIA,Label='FAI_Full')

call mma_allocate(FTA_Full,NFTA,Label='FTA_Full')
call mma_allocate(FAT_Full,NFTA,Label='FAT_Full')

IFIFA = 1
do ISYM=1,NSYM
  NI = NISH(ISYM)
  NA = NASH(ISYM)
  NS = NSSH(ISYM)
  NO = NORB(ISYM)

  FIT(ISYM)%A(1:NA*NI) => FIT_Full(IOFFIT(ISYM)+1:IOFFIT(ISYM)+NA*NI)
  FTI(ISYM)%A(1:NA*NI) => FTI_Full(IOFFIT(ISYM)+1:IOFFIT(ISYM)+NA*NI)

  FIA(ISYM)%A(1:NS*NI) => FIA_Full(IOFFIA(ISYM)+1:IOFFIA(ISYM)+NS*NI)
  FAI(ISYM)%A(1:NS*NI) => FAI_Full(IOFFIA(ISYM)+1:IOFFIA(ISYM)+NS*NI)

  FTA(ISYM)%A(1:NS*NA) => FTA_Full(IOFFTA(ISYM)+1:IOFFTA(ISYM)+NS*NA)
  FAT(ISYM)%A(1:NS*NA) => FAT_Full(IOFFTA(ISYM)+1:IOFFTA(ISYM)+NS*NA)

  call FBLOCK(FIFA(IFIFA),NO,NI,NA,NS,FIT(ISYM)%A(:),FTI(ISYM)%A(:),FIA(ISYM)%A(:),FAI(ISYM)%A(:),FTA(ISYM)%A(:),FAT(ISYM)%A(:))

  IFIFA = IFIFA+(NO*(NO+1))/2

end do

call TIMING(CPU0,CPU,TIO0,TIO)

! Is it possible to reduce to one loop? We have to compute bra and
! ket overlap and bra and ket wavefunctions are not identical, so
! it seems impossible to reduce?

NLOOP = 2
do ILOOP=1,NLOOP
  !! ILOOP1 : <T+lambda|H|T       >
  !! ILOOP2 : <T       |H|T+lambda>

  ! Loop over types and symmetry block of sigma vector:
  do ICASE1=1,11
    !do ICASE1=1,NCASES
    do ISYM1=1,NSYM
      if (NINDEP(ISYM1,ICASE1) == 0) cycle
      NIS1 = NISUP(ISYM1,ICASE1)
      NAS1 = NASUP(ISYM1,ICASE1)
      NIN1 = NINDEP(ISYM1,ICASE1)
      NSGM2 = NIS1*NAS1
      if (NSGM2 == 0) cycle

      call mma_allocate(SGM2,NSGM2,Label='SGM2')
      SGM2(:) = Zero

      NSGM1 = 0
      !LSGM1 = 1
      if (ICASE1 == 1) then
        NSGM1 = NASH(ISYM1)*NISH(ISYM1)
      else if (ICASE1 == 4) then
        NSGM1 = NASH(ISYM1)*NSSH(ISYM1)
      else if ((ICASE1 == 5) .and. (ISYM1 == 1)) then
        NSGM1 = NIS1
      end if
      call mma_allocate(SGM1,max(1,NSGM1),Label='SGM1')
      SGM1(:) = Zero

      IMLTOP = 0
      do ICASE2=ICASE1+1,NCASES
        IFC = IFCOUP(ICASE2,ICASE1)
        if (IFC == 0) cycle
        do ISYM2=1,NSYM
          if (NINDEP(ISYM2,ICASE2) == 0) cycle
          NIS2 = NISUP(ISYM2,ICASE2)
          NAS2 = NASUP(ISYM2,ICASE2)
          NCX = NIS2*NAS2
          if (NCX == 0) cycle

          call RHS_ALLO(NAS2,NIS2,lg_CX)
          call RHS_READ(NAS2,NIS2,lg_CX,ICASE2,ISYM2,IVEC)
          if ((IVEC /= JVEC) .and. (ILOOP == 2)) then
            !! T = T + \lambda
            if (SCAL /= One) call RHS_SCAL(NAS2,NIS2,lg_CX,SCAL)
            call RHS_ALLO(NAS2,NIS2,lg_V1)
            call RHS_READ(NAS2,NIS2,lg_V1,ICASE2,ISYM2,JVEC)
            call RHS_DAXPY(NAS2,NIS2,One,lg_V1,lg_CX)
            call RHS_FREE(lg_V1)
          end if
          ! SVC: for case H (12,13) we can now pass the distributed array ID to
          ! the SGM subroutines
          if ((ICASE2 == 12) .or. (ICASE2 == 13)) then
            LCX = lg_CX
            !XTST = RHS_DDOT(NAS2,NIS2,lg_CX,lg_CX)
          else
            LCX = Allocate_GA_Array(NCX,'CX')
            call RHS_GET(NAS2,NIS2,lg_CX,GA_Arrays(LCX)%A)
            call RHS_FREE(lg_CX)
            !XTST = DDOT_(NCX,GA_Arrays(LCX)%A,1,GA_Arrays(LCX)%A,1)
          end if

#         ifdef _DEBUGPRINT_
          write(u6,*) ' ISYM1,ICASE1:',ISYM1,ICASE1
          write(u6,*) ' ISYM2,ICASE2:',ISYM2,ICASE2
          write(u6,*) ' SIGMA calling SGM with IMLTOP=',IMLTOP
#         endif
          ! Compute contribution SGM2 <- CX, and SGM1 <- CX  if any
          call SGM(IMLTOP,ISYM1,ICASE1,ISYM2,ICASE2,SGM1,size(SGM1),SGM2,size(SGM2),LCX,LISTS,size(LISTS))

          if ((ICASE2 == 12) .or. (ICASE2 == 13)) then
            call RHS_FREE(lg_CX)
          else
            call Deallocate_GA_Array(LCX)
          end if
        end do
      end do

      !-SVC: sum the replicate arrays:
      MAX_MESG_SIZE = 2**27
      do LSGM2_STA=1,NSGM2,MAX_MESG_SIZE
        NSGM2_BLK = min(MAX_MESG_SIZE,NSGM2-LSGM2_STA+1)
        call GADGOP(SGM2(LSGM2_STA),NSGM2_BLK,'+')
      end do

      if (NSGM1 > 0) call GADGOP(SGM1,NSGM1,'+')

      ! If there are 1-electron contributions, add them into the 2-el
      ! part (This requires a non-empty active space.)
      if (NSGM1 > 0) then
        FACT = One/real(max(1,NACTEL),kind=wp)
        if (ICASE1 == 1) then
          call SPEC1A(IMLTOP,FACT,ISYM1,SGM2,size(SGM2),SGM1,size(SGM1))
        else if (ICASE1 == 4) then
          call SPEC1C(IMLTOP,FACT,ISYM1,SGM2,size(SGM2),SGM1,size(SGM1))
        else if ((ICASE1 == 5) .and. (ISYM1 == 1)) then
          call SPEC1D(IMLTOP,FACT,SGM2,size(SGM2),SGM1,size(SGM1))
        end if
      end if
      call mma_deallocate(SGM1)

      !-SVC: no need for the replicate arrays any more, fall back to one array
      call RHS_ALLO(NAS1,NIS1,lg_SGM2)
      call RHS_PUT(NAS1,NIS1,lg_SGM2,SGM2)
      call mma_deallocate(SGM2)

      ! Add to sigma array. Multiply by S to  lower index.
      !call RHS_ALLO(NAS1,NIS1,lg_SGMX)
      !call RHS_READ(NAS1,NIS1,lg_SGMX,ICASE1,ISYM1,JVEC)
      if (ICASE1 <= 11) then
        call RHS_ALLO(NAS1,NIS1,lg_CX)
        call RHS_READ(NAS1,NIS1,lg_CX,ICASE1,ISYM1,IVEC)
        if ((IVEC /= JVEC) .and. (ILOOP == 1)) then
          !! T = T + \lambda
          if (SCAL /= One) call RHS_SCAL(NAS1,NIS1,lg_CX,SCAL)
          call RHS_ALLO(NAS1,NIS1,lg_V1)
          call RHS_READ(NAS1,NIS1,lg_V1,ICASE1,ISYM1,JVEC)
          call RHS_DAXPY(NAS1,NIS1,One,lg_V1,lg_CX)
          call RHS_FREE(lg_V1)
        end if

        call mma_allocate(SDER1,NAS1*NAS1,Label='SDER1')
        idSDer = idSDMat(iSym1,iCase1)
        call DDAFILE(LuSTD,2,SDER1,nAS1*nAS1,idSDer)

        call C1S1DER(SDER1,NAS1*NAS1)

        idSDer = idSDMat(iSym1,iCase1)
        call DDAFILE(LuSTD,1,SDER1,nAS1*nAS1,idSDer)
        call mma_deallocate(SDER1)

        call RHS_FREE(lg_CX)
      end if

      !if ((ICASE1 /= 12) .and. (ICASE1 /= 13)) then
      !  call RHS_STRANS(NAS1,NIS1,ALPHA,lg_SGM2,lg_SGMX,ICASE1,ISYM1)
      !else
      !  call RHS_DAXPY(NAS1,NIS1,ALPHA,lg_SGM2,lg_SGMX)
      !end if
      call RHS_FREE(lg_SGM2)

      ! Write SGMX to disk.
      !call RHS_SAVE(NAS1,NIS1,lg_SGMX,ICASE1,ISYM1,JVEC)
      !call RHS_FREE(lg_SGMX)
    end do
  end do

  IMLTOP = 1
  ! Loop over types and symmetry block of CX vector:
  do ICASE1=1,11
    !do ICASE1=1,NCASES
    do ISYM1=1,NSYM
      if (NINDEP(ISYM1,ICASE1) == 0) cycle
      NIS1 = NISUP(ISYM1,ICASE1)
      NAS1 = NASUP(ISYM1,ICASE1)
      ND2 = NIS1*NAS1
      if (ND2 == 0) cycle

      call RHS_ALLO(NAS1,NIS1,lg_D2)
      call RHS_SCAL(NAS1,NIS1,lg_D2,Zero)
      ! Contract S*CX to form D2. Also form D1 from D2, if needed.

      NCX = ND2
      call RHS_ALLO(NAS1,NIS1,lg_CX)
      call RHS_READ(NAS1,NIS1,lg_CX,ICASE1,ISYM1,IVEC)

      if ((IVEC /= JVEC) .and. (ILOOP == 1)) then
        !! T = T + \lambda
        if (SCAL /= One) call RHS_SCAL(NAS1,NIS1,lg_CX,SCAL)
        call RHS_ALLO(NAS1,NIS1,lg_V1)
        call RHS_READ(NAS1,NIS1,lg_V1,ICASE1,ISYM1,JVEC)
        call RHS_DAXPY(NAS1,NIS1,One,lg_V1,lg_CX)
        call RHS_FREE(lg_V1)
      end if

      if ((ICASE1 /= 12) .and. (ICASE1 /= 13)) then
        call RHS_STRANS(NAS1,NIS1,One,lg_CX,lg_D2,ICASE1,ISYM1)
      else
        call RHS_DAXPY(NAS1,NIS1,One,lg_CX,lg_D2)
      end if
      call RHS_FREE(lg_CX)

      call mma_allocate(D2,ND2,Label='D2')
      call RHS_GET(NAS1,NIS1,lg_D2,D2)
      call RHS_FREE(lg_D2)

      ND1 = 0
      !LD1 = 1
      IMLTOP = 1
      FACT = One/real(max(1,NACTEL),kind=wp)
      if (ICASE1 == 1) then
        ND1 = NASH(ISYM1)*NISH(ISYM1)
        if (ND1 > 0) then
          call mma_allocate(D1,ND1,Label='D1')
          D1(:) = Zero
          call SPEC1A(IMLTOP,FACT,ISYM1,D2,size(D2),D1,size(D1))
        end if
      else if (ICASE1 == 4) then
        ND1 = NASH(ISYM1)*NSSH(ISYM1)
        if (ND1 > 0) then
          call mma_allocate(D1,ND1,Label='D1')
          D1(:) = Zero
          call SPEC1C(IMLTOP,FACT,ISYM1,D2,size(D2),D1,size(D1))
        end if
      else if ((ICASE1 == 5) .and. (ISYM1 == 1)) then
        ND1 = NIS1
        if (ND1 > 0) then
          call mma_allocate(D1,ND1,Label='D1')
          D1(:) = Zero
          call SPEC1D(IMLTOP,FACT,D2,size(D2),D1,size(D1))
        end if
      end if
      if (.not. allocated(D1)) call mma_allocate(D1,1,Label='D1')

      !! No need to compute for ICASE2 = 12 and 13
      do ICASE2=ICASE1+1,11 !! NCASES
        if (IFCOUP(ICASE2,ICASE1) == 0) cycle
        do ISYM2=1,NSYM
          if (NINDEP(ISYM2,ICASE2) == 0) cycle
          NIS2 = NISUP(ISYM2,ICASE2)
          NAS2 = NASUP(ISYM2,ICASE2)
          NIN2 = NINDEP(ISYM2,ICASE2)
          NSGMX = NIS2*NAS2
          if (NSGMX == 0) cycle

          if ((ICASE2 == 12) .or. (ICASE2 == 13)) then
            call RHS_ALLO(NAS2,NIS2,lg_SGMX)
            call RHS_READ(NAS2,NIS2,lg_SGMX,ICASE2,ISYM2,JVEC)
            LSGMX = lg_SGMX
          else
            LSGMX = Allocate_GA_Array(NSGMX,'SGMX')
          end if

#         ifdef _DEBUGPRINT_
          write(u6,*) ' ISYM1,ICASE1:',ISYM1,ICASE1
          write(u6,*) ' ISYM2,ICASE2:',ISYM2,ICASE2
          write(u6,*) ' SIGMA calling SGM with IMLTOP=',IMLTOP
#         endif
          ! Compute contribution SGMX <- D2, and SGMX <- D1  if any
          call SGM(IMLTOP,ISYM1,ICASE1,ISYM2,ICASE2,D1,size(D1),D2,size(D2),LSGMX,LISTS,size(LISTS))
          !if (iCase2 <= 11) then
          !end if

          if ((ICASE2 /= 12) .and. (ICASE2 /= 13)) then
            MAX_MESG_SIZE = 2**27
            do LSGMX_STA=1,NSGMX,MAX_MESG_SIZE
              NSGMX_BLK = min(MAX_MESG_SIZE,NSGMX-LSGMX_STA+1)
              call GADGOP(GA_Arrays(LSGMX)%A(LSGMX_STA),NSGMX_BLK,'+')
            end do
            !call RHS_ALLO(NAS2,NIS2,lg_SGMX)
            !call RHS_READ(NAS2,NIS2,lg_SGMX,ICASE2,ISYM2,JVEC)
            !call RHS_ADD(NAS2,NIS2,lg_SGMX,GA_Array(LSGMX)%A)
            !! do C2DER
            call mma_allocate(SDER2,NAS2*NAS2,Label='SDER2')
            idSDer = idSDMat(iSym2,iCase2)
            call DDAFILE(LuSTD,2,SDER2,nAS2*nAS2,idSDer)

            call C2DER(SDER2,NAS2*NAS2)

            idSDer = idSDMat(iSym2,iCase2)
            call DDAFILE(LuSTD,1,SDER2,nAS2*nAS2,idSDer)
            call mma_deallocate(SDER2)

            !call Deallocate_GA_Array(LSGMX)
          end if

          !-SVC: no need for the replicate arrays any more, fall back to one array
          !call RHS_SAVE(NAS2,NIS2,lg_SGMX,ICASE2,ISYM2,JVEC)
          if ((ICASE2 == 12) .or. (ICASE2 == 13)) then
            call RHS_FREE(lg_SGMX)
          else
            call Deallocate_GA_Array(LSGMX)
          end if
        end do
      end do
      call mma_deallocate(D2)
      call mma_deallocate(D1)
    end do
  end do

end do

call TIMING(CPU1,CPU,TIO1,TIO)
CPUSGM = CPUSGM+(CPU1-CPU0)
TIOSGM = TIOSGM+(TIO1-TIO0)

call mma_deallocate(FIT_Full)
call mma_deallocate(FTI_Full)
call mma_deallocate(FIA_Full)
call mma_deallocate(FAI_Full)
call mma_deallocate(FTA_Full)
call mma_deallocate(FAT_Full)
do iSym=1,nSym
  nullify(FIT(iSym)%A,FTI(iSym)%A,FIA(iSym)%A,FAI(iSym)%A,FTA(iSym)%A,FAT(iSym)%A)
end do

! Transform contrav C  to ei%Agenbasis of H0(diag):
call PTRTOSR(1,IVEC,IVEC)
if (IVEC /= JVEC) call PTRTOSR(1,JVEC,JVEC)

contains

subroutine C1S1DER(SDER,nSDER)

  integer(kind=iwp), intent(in) :: nSDER
  real(kind=wp), intent(inout) :: SDER(nSDER)
  integer(kind=iwp) :: iType
# if defined(_MOLCAS_MPP_) && defined(_GA_)
  integer(kind=iwp) :: lg_SDER
# endif

  ! (T2Ct2*f)py * (T1Ct1)pz * dS1yz/da
  ! -1/2 (T2Ct2*f*S1*C1*Ct1)pt * (T1Ct1)pu * dS1tu/da

  !! Finalize the derivative of S1
  !! 2S. (T2Ct2*f) * T1Ct1
# if defined(_MOLCAS_MPP_) && defined(_GA_)
  if (is_real_par()) then
    call GA_CREATE_STRIPED('H',NAS1,NAS1,'SDER',lg_SDER)
    call GA_PUT(lg_SDER,1,NAS1,1,NAS1,SDER,NAS1)
    call GA_DGEMM('N','T',NAS1,NAS1,NIS1,Two,lg_CX,lg_SGM2,One,lg_SDER)
  else
# endif
    call DGEMM_('N','T',NAS1,NAS1,NIS1,Two,GA_Arrays(lg_CX)%A,NAS1,GA_Arrays(lg_SGM2)%A,NAS1,One,SDER,NAS1)
# if defined(_MOLCAS_MPP_) && defined(_GA_)
  end if
# endif
  !do i=1,nas1*nis1
  !  write(u6,'(i4,2f20.10)') ,i,GA_Arrays(lg_cx)%A(i),GA_Arrays(lg_sgm2)%A(i)
  !end do

  !! Next, the derivative of C1
  !! 2C. (T2Ct2*f) * S1*C1 (MO -> IC)
  !!     lg_T * lg_V2 -> lg_V1
  call RHS_ALLO(NIN1,NIS1,lg_V1)
  ITYPE = 1
  !! SGM2 is local quantity, so put this in GA?
  call RHS_SR2C(ITYPE,1,NAS1,NIS1,NIN1,lg_V1,lg_SGM2,ICASE1,ISYM1)
  !! 3C. (T2Ct2*f) * S1*C1 * Ct1 (IC -> MO)
  !!     lg_T * lg_V1 -> lg_V2
  ITYPE = 0
  call RHS_SR2C(ITYPE,0,NAS1,NIS1,NIN1,lg_V1,lg_SGM2,ICASE1,ISYM1)
  call RHS_FREE(lg_V1)

  !! 4C. (T1Ct1*f) * (T2Ct2St2*f*C1*Ct1)
# if defined(_MOLCAS_MPP_) && defined(_GA_)
  if (is_real_par()) then
    call GA_DGEMM('N','T',NAS1,NAS1,NIS1,-One,lg_CX,lg_SGM2,One,lg_SDER)
    call GA_GET(lg_SDER,1,NAS1,1,NAS1,SDER,NAS1)
    bStat = GA_destroy(lg_SDER)
  else
# endif
    call DGEMM_('N','T',NAS1,NAS1,NIS1,-One,GA_Arrays(lg_CX)%A,NAS1,GA_Arrays(lg_SGM2)%A,NAS1,One,SDER,NAS1)
# if defined(_MOLCAS_MPP_) && defined(_GA_)
  end if
# endif

end subroutine C1S1DER

subroutine C2DER(SDER,nSDER)

  integer(kind=iwp), intent(in) :: nSDER
  real(kind=wp), intent(inout) :: SDER(nSDER)
  integer(kind=iwp) :: iType, lg_Sgm, lg_v2
# if defined(_MOLCAS_MPP_) && defined(_GA_)
  integer(kind=iwp) :: lg_SDER
# endif

  ! -1/2 (T2Ct2)pu * dS2tu/da * (T1Ct1St1*f*C2*Ct2)pt

# if defined(_MOLCAS_MPP_) && defined(_GA_)
  if (is_real_par()) then
    call GA_CREATE_STRIPED('V',NAS2,NIS2,'SDER',lg_SGMX)
    call GA_PUT(lg_SGMX,1,NAS2,1,NIS2,GA_Arrays(LSGMX)%A,NAS2)
  else
# endif
    lg_SGMX = LSGMX
# if defined(_MOLCAS_MPP_) && defined(_GA_)
  end if
# endif

  !! For icase = 12 or 13, there is no need to transform,
  !! so LTMP is always replicated, but T for A and C are
  !! distributed?, so...
  call RHS_ALLO(NIN2,NIS2,lg_V2)
  !! 2. (T1Ct1St1*f) * C2 (MO -> IC; LTMP -> LTMP2)
  ITYPE = 0
  call RHS_SR2C(ITYPE,1,NAS2,NIS2,NIN2,lg_V2,lg_SGMX,ICASE2,ISYM2)
  !! 3. (T1Ct1St1*f) * C2 * Ct2 (IC -> MO; LTMP2 -> LTMP)
  call RHS_SR2C(ITYPE,0,NAS2,NIS2,NIN2,lg_V2,lg_SGMX,ICASE2,ISYM2)
  call RHS_FREE(lg_V2)

  !! 4. (T2Ct2*f) * (T1Ct1St1*f*C2*Ct2)
  call RHS_ALLO(NAS2,NIS2,lg_SGM)
  call RHS_READ(NAS2,NIS2,lg_SGM,ICASE2,ISYM2,IVEC)
  if ((IVEC /= JVEC) .and. (ILOOP == 2)) then
    !! T = T + \lambda
    if (SCAL /= One) call RHS_SCAL(NAS2,NIS2,lg_SGM,SCAL)
    call RHS_ALLO(NAS2,NIS2,lg_V1)
    call RHS_READ(NAS2,NIS2,lg_V1,ICASE2,ISYM2,JVEC)
    call RHS_DAXPY(NAS2,NIS2,One,lg_V1,lg_SGM)
    call RHS_FREE(lg_V1)
  end if
# if defined(_MOLCAS_MPP_) && defined(_GA_)
  if (is_real_par()) then
    call GA_CREATE_STRIPED('H',NAS2,NAS2,'SDER',lg_SDER)
    call GA_PUT(lg_SDER,1,NAS2,1,NAS2,SDER,NAS2)
    call GA_DGEMM('N','T',NAS2,NAS2,NIS2,-One,lg_SGM,lg_SGMX,One,lg_SDER)
    call GA_GET(lg_SDER,1,NAS2,1,NAS2,SDER,NAS2)
    bStat = GA_destroy(lg_SGMX)
    bStat = GA_destroy(lg_SDER)
  else
# endif
    call DGEMM_('N','T',NAS2,NAS2,NIS2,-One,GA_Arrays(lg_SGM)%A,NAS2,GA_Arrays(LSGMX)%A,NAS2,One,SDER,NAS2)
# if defined(_MOLCAS_MPP_) && defined(_GA_)
  end if
# endif
  call RHS_FREE(lg_SGM)

end subroutine C2DER

end subroutine sigder
