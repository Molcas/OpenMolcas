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
! Copyright (C) 1994,1999, Per Ake Malmqvist                           *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
! 1999: GEMINAL R12 ENABLED                  *
!--------------------------------------------*

subroutine SIGMA_CASPT2(ALPHA,BETA,IVEC,JVEC)
! Compute |JVEC> := BETA* |JVEC> + ALPHA* (H0-E0)* |IVEC>
! where the vectors are represented in transformed basis and
! are  stored at positions IVEC and JVEC on the LUSOLV unit.

use Index_Functions, only: nTri_Elem
use Fockof, only: FAI, FAI_Full, FAT, FAT_Full, FIA, FIA_Full, FIT, FIT_Full, FTA, FTA_Full, FTI, FTI_Full, IOFFIA, IOFFIT, IOFFTA
use EQSOLV, only: IFCoup
use Sigma_data, only: IFTEST, NFDXP, NFMV, NFR1, NFSCA
use fake_GA, only: Allocate_GA_Array, Deallocate_GA_Array, GA_Arrays
use caspt2_global, only: FIFA, LISTS
use general_data, only: nActel, nAsh
use caspt2_module, only: CPUSGM, FockType, G1SecIn, MaxIt, nASup, nCases, nInDep, nIsh, nISup, nOrb, nSsh, nSym, &
                         ThrShn, ThrShS, TIOSGM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: ALPHA, BETA
integer(kind=iwp), intent(in) :: IVEC, JVEC
integer(kind=iwp) :: iCASE1, iCase2, IfC, IFIFA, IMLTOP, ISYM, ISYM1, ISYM2, lCX, lg_CX, lg_D2, lg_Sgm2, lg_SgmX, lSgm2_Sta, &
                     lSgmX, lSgmX_Sta, Max_MESG_Size, NA, NAS1, NAS2, nCX, ND1, ND2, NFIA, NFIT, NFTA, NI, NIS1, NIS2, NO, NS, &
                     nSgm1, nSgm2, nSgm2_Blk, nSgmX, nSgmX_blk
real(kind=wp) :: CPU, CPU0, CPU1, Fact, TIO, TIO0, TIO1, XTST
real(kind=wp), allocatable :: D1(:), D2(:), SGM1(:), SGM2(:)
real(kind=wp), external :: RHS_DDot, DDot_

#ifdef _DEBUGPRINT_
write(u6,*) ' Entering SIGMA.'
write(u6,*) ' Compute |JVEC> := Beta*|JVEC> + Alpha*(H0-E0)|IVEC>'
write(u6,'(1x,a,2f15.6)') 'Alpha,Beta:',Alpha,Beta
write(u6,'(1x,a,2i5)') 'IVEC,JVEC:',IVEC,JVEC
#endif

! If the G1 correction to the Fock matrix is used, then the
! inactive/virtual coupling elements (which are non-zero for the
! case of average CASSCF) cannot be used in the CASPT2 equations.
if ((FOCKTYPE == 'G1') .and. (.not. G1SECIN)) then
  IFCOUP(12,5) = 0
  IFCOUP(13,5) = 0
end if

IFTEST = 0
! Flop counts:
NFSCA = 0
NFDXP = 0
NFMV = 0
NFR1 = 0
! First compute diagonal block contributions:
!TEST write(u6,*) ' First, do it for (H0(diag)-E0).'
call PSGMDIA(ALPHA,BETA,IVEC,JVEC)
if (ALPHA == Zero) return
if (MAXIT == 0) return
!TEST write(u6,*) ' From now on, scaling with BETA is already done.'
!TEST write(u6,*) ' Test print after SGMDIA call in SIGMA:'
!TEST write(u6,*) ' Should be zero, in first iteration.'
!TEST call OVLPRT(JVEC,JVEC,OVLAPS)
! From now on, scaling with BETA is already done.

! Transform to standard representation:
call PTRTOC(0,IVEC,IVEC)
call PTRTOC(1,JVEC,JVEC)

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
NFIT = NFIT+1 !?
NFIA = NFIA+1 !?
NFTA = NFTA+1 !?

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

  if (NO > 0) then
    FIT(ISYM)%A(1:NA*NI) => FIT_Full(IOFFIT(ISYM)+1:IOFFIT(ISYM)+NA*NI)
    FTI(ISYM)%A(1:NA*NI) => FTI_Full(IOFFIT(ISYM)+1:IOFFIT(ISYM)+NA*NI)

    FIA(ISYM)%A(1:NS*NI) => FIA_Full(IOFFIA(ISYM)+1:IOFFIA(ISYM)+NS*NI)
    FAI(ISYM)%A(1:NS*NI) => FAI_Full(IOFFIA(ISYM)+1:IOFFIA(ISYM)+NS*NI)

    FTA(ISYM)%A(1:NS*NA) => FTA_Full(IOFFTA(ISYM)+1:IOFFTA(ISYM)+NS*NA)
    FAT(ISYM)%A(1:NS*NA) => FAT_Full(IOFFTA(ISYM)+1:IOFFTA(ISYM)+NS*NA)

    call FBLOCK(FIFA(IFIFA),NO,NI,NA,NS,FIT(ISYM)%A(:),FTI(ISYM)%A(:),FIA(ISYM)%A(:),FAI(ISYM)%A(:),FTA(ISYM)%A(:),FAT(ISYM)%A(:))

    IFIFA = IFIFA+nTri_Elem(NO)
  end if

end do

call TIMING(CPU0,CPU,TIO0,TIO)
! Loop over types and symmetry block of sigma vector:
do ICASE1=1,11
  do ISYM1=1,NSYM
    if (NINDEP(ISYM1,ICASE1) == 0) cycle
    NIS1 = NISUP(ISYM1,ICASE1)
    NAS1 = NASUP(ISYM1,ICASE1)
    NSGM2 = NIS1*NAS1
    if (NSGM2 == 0) cycle

    call mma_allocate(SGM2,NSGM2,Label='SGM2')
    SGM2(:) = Zero

    if (ICASE1 == 1) then
      NSGM1 = NASH(ISYM1)*NISH(ISYM1)
    else if (ICASE1 == 4) then
      NSGM1 = NASH(ISYM1)*NSSH(ISYM1)
    else if ((ICASE1 == 5) .and. (ISYM1 == 1)) then
      NSGM1 = NIS1
    else
      NSGM1 = 0
    end if
    call mma_allocate(SGM1,max(1,NSGM1),LABEL='SGM1')
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
        ! SVC: for case H (12,13) we can now pass the distributed array ID to
        ! the SGM subroutines
        if ((ICASE2 == 12) .or. (ICASE2 == 13)) then
          LCX = lg_CX
          XTST = RHS_DDOT(NAS2,NIS2,lg_CX,lg_CX)
        else
          LCX = Allocate_GA_Array(NCX,'CX')
          call RHS_GET(NAS2,NIS2,lg_CX,GA_Arrays(LCX)%A)
          call RHS_FREE(lg_CX)
          XTST = DDOT_(NCX,GA_Arrays(LCX)%A,1,GA_Arrays(LCX)%A,1)
        end if

        if (XTST > 1.0e12_wp) then
          write(u6,'(1x,a,6i10)') ' SIGMA A. ICASE2,ISYM2:',ICASE2,ISYM2
          call Crash()
        end if

#       ifdef _DEBUGPRINT_
        write(u6,*) ' ISYM1,ICASE1:',ISYM1,ICASE1
        write(u6,*) ' ISYM2,ICASE2:',ISYM2,ICASE2
        write(u6,*) ' SIGMA calling SGM with IMLTOP=',IMLTOP
#       endif
        ! Compute contribution SGM2 <- CX, and SGM1 <- CX  if any
        call SGM(IMLTOP,ISYM1,ICASE1,ISYM2,ICASE2,SGM1,size(SGM1),SGM2,size(SGM2),LCX,LISTS,size(LISTS))

        if ((ICASE2 == 12) .or. (ICASE2 == 13)) then
          call RHS_FREE(lg_CX)
        else
          call Deallocate_GA_Array(LCX)
        end if

        ! Check for colossal values of SGM2 and SGM1
        XTST = DDOT_(NSGM2,SGM2,1,SGM2,1)
        if (XTST > 1.0e12_wp) then
          write(u6,'(1x,a,6i10)') ' SIGMA B. ICASE1,ISYM1:',ICASE1,ISYM1
          write(u6,'(1x,a,6i10)') '          ICASE2,ISYM2:',ICASE2,ISYM2
          call Crash()
        end if

        if (NSGM1 > 0) then
          XTST = DDOT_(NSGM1,SGM1,1,SGM1,1)
          if (XTST > 1.0e12_wp) then
            write(u6,'(1x,a,6i10)') ' SIGMA B2. ICASE1,ISYM1:',ICASE1,ISYM1
            write(u6,'(1x,a,6i10)') '           ICASE2,ISYM2:',ICASE2,ISYM2
            call Crash()
          end if
        end if

      end do
    end do

    !-SVC: sum the replicate arrays:
    MAX_MESG_SIZE = 2**27
    do LSGM2_STA=1,NSGM2,MAX_MESG_SIZE
      NSGM2_BLK = min(MAX_MESG_SIZE,NSGM2-LSGM2_STA+1)
      call GADGOP(SGM2(LSGM2_STA:),NSGM2_BLK,'+')
    end do

    if (NSGM1 > 0) call GADGOP(SGM1,NSGM1,'+')

    !XTST2 = DDOT_(NSGM2,SGM2,1,SGM2,1)
    !XTST1 = Zero
    !if (NSGM1 > 0) XTST1 = DDOT_(NSGM1,SGM1,1,SGM1,1)
    !write(u6,'(1x,a,a,i2,2f16.6)') 'Contr. SGM2, SGM1, ',cases(icase1),isym1,xtst2,xtst1

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

      XTST = DDOT_(NSGM2,SGM2,1,SGM2,1)
      if (XTST > 1.0e12_wp) then
        write(u6,'(1x,a,6i10)') ' SIGMA C. ICASE1,ISYM1:',ICASE1,ISYM1
        call Crash()
      end if

    end if
    call mma_deallocate(SGM1)

    !-SVC: no need for the replicate arrays any more, fall back to one array
    call RHS_ALLO(NAS1,NIS1,lg_SGM2)
    call RHS_PUT(NAS1,NIS1,lg_SGM2,SGM2)
    call mma_deallocate(SGM2)

    ! Add to sigma array. Multiply by S to  lower index.
    NSGMX = NSGM2
    call RHS_ALLO(NAS1,NIS1,lg_SGMX)
    call RHS_READ(NAS1,NIS1,lg_SGMX,ICASE1,ISYM1,JVEC)

    XTST = RHS_DDOT(NAS1,NIS1,lg_SGMX,lg_SGMX)
    if (XTST > 1.0e12_wp) then
      write(u6,'(1x,a,6i10)') ' SIGMA D. ICASE1,ISYM1:',ICASE1,ISYM1
      write(u6,'(1x,a,6i10)') '          ICASE2,ISYM2:',ICASE2,ISYM2
      call Crash()
    end if

    !if ((ICASE1 /= 12) .and. (ICASE1 /= 13)) then
    call RHS_STRANS(NAS1,NIS1,ALPHA,lg_SGM2,lg_SGMX,ICASE1,ISYM1)
    !else
    !  call RHS_DAXPY(NAS1,NIS1,ALPHA,lg_SGM2,lg_SGMX)
    !end if
    call RHS_FREE(lg_SGM2)

    XTST = RHS_DDOT(NAS1,NIS1,lg_SGMX,lg_SGMX)
    if (XTST > 1.0e12_wp) then
      write(u6,'(1x,a,6i10)') ' SIGMA E. ICASE1,ISYM1:',ICASE1,ISYM1
      write(u6,'(1x,a,6i10)') '          ICASE2,ISYM2:',ICASE2,ISYM2
      call Crash()
    end if

    ! Write SGMX to disk.
    call RHS_SAVE(NAS1,NIS1,lg_SGMX,ICASE1,ISYM1,JVEC)
    call RHS_FREE(lg_SGMX)
  end do
end do

IMLTOP = 1
! Loop over types and symmetry block of CX vector:
do ICASE1=1,11
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

    XTST = RHS_DDOT(NAS1,NIS1,lg_CX,lg_CX)
    if (XTST > 1.0e12_wp) then
      write(u6,'(1x,a,6i10)') ' SIGMA F. ICASE1,ISYM1:',ICASE1,ISYM1
      call Crash()
    end if

    if ((ICASE1 /= 12) .and. (ICASE1 /= 13)) then
      call RHS_STRANS(NAS1,NIS1,ALPHA,lg_CX,lg_D2,ICASE1,ISYM1)
    else
      call RHS_DAXPY(NAS1,NIS1,ALPHA,lg_CX,lg_D2)
    end if
    call RHS_FREE(lg_CX)

    !PAM Sanity check:
    XTST = RHS_DDOT(NAS1,NIS1,lg_D2,lg_D2)
    if (XTST > 1.0e12_wp) then
      write(u6,'(1x,a,6i10)') ' SIGMA G1 ICASE1,ISYM1:',ICASE1,ISYM1
      write(u6,'(1x,a,6i10)') '          ICASE2,ISYM2:',ICASE2,ISYM2
      call Crash()
    end if

    call mma_allocate(D2,ND2,Label='D2')
    call RHS_GET(NAS1,NIS1,lg_D2,D2)
    call RHS_FREE(lg_D2)

    ND1 = 0
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

    if (ND1 > 0) then
      XTST = DDOT_(ND1,D1,1,D1,1)
      if (XTST > 1.0e12_wp) then
        write(u6,'(1x,a,6i10)') ' SIGMA G2 ICASE1,ISYM1:',ICASE1,ISYM1
        write(u6,'(1x,a,6i10)') '          ICASE2,ISYM2:',ICASE2,ISYM2
        call Crash()
      end if
    end if

    do ICASE2=ICASE1+1,NCASES
      if (IFCOUP(ICASE2,ICASE1) == 0) cycle
      do ISYM2=1,NSYM
        if (NINDEP(ISYM2,ICASE2) == 0) cycle
        NIS2 = NISUP(ISYM2,ICASE2)
        NAS2 = NASUP(ISYM2,ICASE2)
        NSGMX = NIS2*NAS2
        if (NSGMX == 0) cycle

        if ((ICASE2 == 12) .or. (ICASE2 == 13)) then
          call RHS_ALLO(NAS2,NIS2,lg_SGMX)
          call RHS_READ(NAS2,NIS2,lg_SGMX,ICASE2,ISYM2,JVEC)
          LSGMX = lg_SGMX
        else
          LSGMX = Allocate_GA_Array(NSGMX,'SGMX')
        end if

        ! SVC: this array is just zero....
        !XTST = DDOT_(NSGMX,GA_Array(LSGMX)%A,1,GA_Array(LSGMX)%A,1)
        !if (XTST > 1.0e12_wp) then
        !  write(u6,'(1x,a,6i10)') ' SIGMA H. ICASE1,ISYM1:',ICASE1,ISYM1
        !  write(u6,'(1x,a,6i10)') '          ICASE2,ISYM2:',ICASE2,ISYM2
        !  call Crash()
        !end if

#       ifdef _DEBUGPRINT_
        write(u6,*) ' ISYM1,ICASE1:',ISYM1,ICASE1
        write(u6,*) ' ISYM2,ICASE2:',ISYM2,ICASE2
        write(u6,*) ' SIGMA calling SGM with IMLTOP=',IMLTOP
#       endif
        ! Compute contribution SGMX <- D2, and SGMX <- D1  if any
        call SGM(IMLTOP,ISYM1,ICASE1,ISYM2,ICASE2,D1,size(D1),D2,size(D2),LSGMX,LISTS,size(LISTS))

        if ((ICASE2 == 12) .or. (ICASE2 == 13)) then
          XTST = RHS_DDOT(NAS2,NIS2,lg_SGMX,lg_SGMX)
        else
          XTST = DDOT_(NSGMX,GA_Arrays(LSGMX)%A,1,GA_Arrays(LSGMX)%A,1)
        end if

        if (XTST > 1.0e12_wp) then
          write(u6,'(1x,a,6i10)') ' SIGMA I. ICASE1,ISYM1:',ICASE1,ISYM1
          write(u6,'(1x,a,6i10)') '          ICASE2,ISYM2:',ICASE2,ISYM2
          call Crash()
        end if

        if ((ICASE2 /= 12) .and. (ICASE2 /= 13)) then
          MAX_MESG_SIZE = 2**27
          do LSGMX_STA=1,NSGMX,MAX_MESG_SIZE
            NSGMX_BLK = min(MAX_MESG_SIZE,NSGMX-LSGMX_STA+1)
            call GADGOP(GA_Arrays(LSGMX)%A(LSGMX_STA),NSGMX_BLK,'+')
          end do
          call RHS_ALLO(NAS2,NIS2,lg_SGMX)
          call RHS_READ(NAS2,NIS2,lg_SGMX,ICASE2,ISYM2,JVEC)
          call RHS_ADD(NAS2,NIS2,lg_SGMX,GA_Arrays(LSGMX)%A)
          call Deallocate_GA_Array(LSGMX)
        end if

        !-SVC: no need for the replicate arrays any more, fall back to one array
        call RHS_SAVE(NAS2,NIS2,lg_SGMX,ICASE2,ISYM2,JVEC)
        call RHS_FREE(lg_SGMX)
      end do
    end do
    call mma_deallocate(D2)
    call mma_deallocate(D1)
  end do
end do

call TIMING(CPU1,CPU,TIO1,TIO)
CPUSGM = CPUSGM+(CPU1-CPU0)
TIOSGM = TIOSGM+(TIO1-TIO0)

#ifdef _DEBUGPRINT_
write(u6,*) ' End of SIGMA. Flop counts:'
write(u6,'(a,i12)') ' In MLTSCA:',NFSCA
write(u6,'(a,i12)') ' In MLTDXP:',NFDXP
write(u6,'(a,i12)') ' In MLTMV :',NFMV
write(u6,'(a,i12)') ' In MLTR1 :',NFR1
write(u6,*)
#endif

call mma_deallocate(FIT_Full)
call mma_deallocate(FTI_Full)
call mma_deallocate(FIA_Full)
call mma_deallocate(FAI_Full)
call mma_deallocate(FTA_Full)
call mma_deallocate(FAT_Full)
do iSym=1,nSym
  nullify(FIT(iSym)%A,FTI(iSym)%A,FIA(iSym)%A,FAI(iSym)%A,FTA(iSym)%A,FAT(iSym)%A)
end do

! Transform contrav C  to eigenbasis of H0(diag):
call PTRTOSR(1,IVEC,IVEC)
! Transform covar. sigma to eigenbasis of H0(diag):
call PTRTOSR(0,JVEC,JVEC)

contains

subroutine Crash()

  write(u6,*) ' Colossal value detected in SIGMA.'
  write(u6,*) ' This implies that the thresholds used for linear'
  write(u6,*) ' dependence removal must be increased.'
  write(u6,*) ' Present values, THRSHN, THRSHS:',THRSHN,THRSHS
  write(u6,*) ' Use keyword THRESHOLD in input to increase these'
  write(u6,*) ' values and then run again.'
  call ABEND()

end subroutine Crash

end subroutine SIGMA_CASPT2
