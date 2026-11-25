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
! The SGUGA self-tests (SG_Test_EpqErs_Ref, MOLCAS_CI_VERIFY).
! (Default public, as SGUGA_Core; the interface is the module SGUGA.)
module SGUGA_Checks

use SGUGA_Core
use SGUGA_MatEl
use SGUGA_Triplet
use SGUGA_Cnf
use SGUGA_Dets
use SGUGA_Epq
use SGUGA_Tables

implicit none

contains

subroutine SG_Test_EpqErs_Ref(iState,LSM,HPsi,HPsi_RAS,Epq_Only)

  use Constants, only: Zero, One
  use Definitions, only: wp, iwp, u6
  use stdalloc, only: mma_allocate, mma_deallocate
  use, intrinsic :: iso_fortran_env, only: int64

  implicit none

! If Epq_Only is present and true, only <I|Epq|J> is verified and the
! <I|EpqErs|J> tests are skipped.
  logical(kind=iwp), intent(in), optional :: Epq_Only

  integer(kind=iwp), intent(in) :: iState
  integer(kind=iwp), intent(in) :: LSM
! the sigma routines of CAS and RAS graphs (sg_h_psi, sg_h_psi_ras of
! SGUGA_Sigma), the reference of the <I|H|J> check (Check_MatEl_H)
  procedure(SG_HPsi_Proc) :: HPsi, HPsi_RAS

!
! Operator-product topologies.
!
  integer(kind=iwp), parameter :: iUU = 1
  integer(kind=iwp), parameter :: iLL = 2
  integer(kind=iwp), parameter :: iUL = 3
  integer(kind=iwp), parameter :: iLU = 4
  integer(kind=iwp), parameter :: iUO = 5
  integer(kind=iwp), parameter :: iOU = 6
  integer(kind=iwp), parameter :: iLO = 7
  integer(kind=iwp), parameter :: iOL = 8
  integer(kind=iwp), parameter :: iOO = 9

  integer(kind=iwp), parameter :: nTopology = 9

  character(len=2), parameter :: TopologyName(nTopology) = &
                                 [character(len=2) :: 'UU','LL','UL','LU','UO','OU','LO','OL','OO']

!
! Diagonal-status classes.
!
!     CC: closed   x closed
!     DC: diagonal x closed
!     CD: closed   x diagonal
!     DD: diagonal x diagonal
!
  integer(kind=iwp), parameter :: iCC = 1
  integer(kind=iwp), parameter :: iDC = 2
  integer(kind=iwp), parameter :: iCD = 3
  integer(kind=iwp), parameter :: iDD = 4

  integer(kind=iwp), parameter :: nDiagonalCase = 4

  character(len=2), parameter :: DiagonalName(nDiagonalCase) = &
                                 [character(len=2) :: 'CC','DC','CD','DD']

!
! Orientation classes for closed x closed products.
!
!     CC: canonical x canonical
!     RC: reverse   x canonical
!     CR: canonical x reverse
!     RR: reverse   x reverse
!
  integer(kind=iwp), parameter :: iCanCan = 1
  integer(kind=iwp), parameter :: iRevCan = 2
  integer(kind=iwp), parameter :: iCanRev = 3
  integer(kind=iwp), parameter :: iRevRev = 4

  integer(kind=iwp), parameter :: nOrientationCase = 4

  character(len=2), parameter :: &
    OrientationName(nOrientationCase) = &
    [character(len=2) :: 'CC','RC','CR','RR']

  integer(kind=iwp) :: iTopology
  integer(kind=iwp) :: iDiagonalCase
  integer(kind=iwp) :: iOrientationCase

  integer(kind=iwp) :: p
  integer(kind=iwp) :: q
  integer(kind=iwp) :: r
  integer(kind=iwp) :: s

  integer(kind=iwp) :: pMin
  integer(kind=iwp) :: pMax
  integer(kind=iwp) :: qMin
  integer(kind=iwp) :: qMax
  integer(kind=iwp) :: rMin
  integer(kind=iwp) :: rMax
  integer(kind=iwp) :: sMin
  integer(kind=iwp) :: sMax

  integer(kind=iwp) :: I, J, K

  integer(kind=iwp) :: nCSF

!
! Counters indexed by operator-product topology.
!
  integer(kind=iwp) :: nOperatorPairs(nTopology)
  integer(kind=iwp) :: nMatrixElements(nTopology)
  integer(kind=iwp) :: nDiscrepancies(nTopology)

!
! Diagonal-status counters:
!
!     first index  = topology
!     second index = CC, DC, CD, or DD
!
  integer(kind=iwp) :: nPairDiagonal( &
                       nTopology,nDiagonalCase)

  integer(kind=iwp) :: nErrDiagonal( &
                       nTopology,nDiagonalCase)

!
! Orientation counters for closed x closed products:
!
!     first index  = topology
!     second index = CC, RC, CR, or RR
!
  integer(kind=iwp) :: nPairOrientation( &
                       nTopology,nOrientationCase)

  integer(kind=iwp) :: nErrOrientation( &
                       nTopology,nOrientationCase)

  real(kind=wp) :: XRef
  real(kind=wp) :: XTest
  real(kind=wp) :: Difference

  real(kind=wp), parameter :: Threshold = 1.0e-12_wp

  logical :: FirstDiagonal
  logical :: SecondDiagonal

  logical :: FirstCanonical
  logical :: SecondCanonical

  integer(kind=iwp) :: nTest, nFail
  real(kind=wp) :: MaxDiff
  integer(kind=iwp) :: iSYCI
!
! Quick zero test: orbital occupations of the CSFs of symmetry LSM,
! Occ(t,I) = occupation of level t in CSF I, and the number of matrix
! elements per topology that were screened as exactly zero.
!
  integer(kind=iwp), allocatable :: Occ(:,:)
  integer(kind=iwp) :: nScreened(nTopology)
!
! Sampling of the (I,J) pairs. For a given operator (pair), all pairs are
! tested if there are at most MaxPairs of them. Otherwise random pairs are
! drawn (at most MaxDraw, and at most twice the number of pairs) and
! accepted until QuotaNZ pairs that pass the
! quick zero test (compared with the reference) and QuotaZ pairs that are
! zero by occupation (fast routine must return exactly zero) are taken.
! The seed is taken from the clock, so each run tests a different subset,
! unless MOLCAS_RANDOM_SEED is set. The seed is printed before the first
! sampled pair, so a failed run can be repeated with
! >export MOLCAS_RANDOM_SEED=<seed>.
!
  integer(kind=iwp), parameter :: MaxPairs = 40000, QuotaNZ = 10000, &
                                  QuotaZ = 10000, MaxDraw = 2000000
  integer(kind=int64) :: iRandState, iRandSeed
  logical :: SeedFromEnv, SeedReported
  integer(kind=iwp) :: nPairsIJ, nDrawMax, iDraw, nTakeZ, nTakeNZ
  integer(kind=iwp) :: nSampledOps(nTopology)
  logical :: Sampled
  logical :: Screened
#ifdef _CHECK_QUICK_ZERO_
  real(kind=wp) :: XChk
#endif

#ifdef _DEBUGPRINT_

  do ISYCI=1,SGS(istate)%nSym
    write(u6,*) 'ISYCI=',ISYCI
    call Print_CSF_Walks(SGS(iState),CIS(iState),ISYCI)
  end do
#endif

! Nothing to test without active levels (e.g. RASSCF with an empty active space)
  if (SGS(iState)%nLev == 0) then
    write(u6,'(1X,A)') 'SG_Test_EpqErs_Ref: no active levels, nothing to test'
    return
  end if

! State and symmetry of the module-level matrix-element routines
  ME_iState = iState
  ME_LSM = LSM
  call Reset_Decode_ME()
  if (.not. allocated(EXS(iState)%IOrd)) call Mk_Coup_Index(EXS(iState))
! Verify that MatEl_Epq is correct

  call Init_Rand_Seed()
  do ISYCI=1,SGS(istate)%nSym
    call Check_MatEl_Epq(iState,ISYCI,1.0e-12_wp,nTest,nFail,MaxDiff)
    if (nFail /= 0) call Abend()
  end do
! The ASTA walk products of E_pq and of the triplet operator T_pq
  call Check_MatEl_Tpq(iState,nFail)
  if (nFail /= 0) call Abend()

  if (present(Epq_Only)) then
    if (Epq_Only) return
  end if

  nCSF = CIS(iState)%NCSF(LSM)

!
! Orbital occupations of all CSFs of symmetry LSM, for the quick zero test,
! unpacked from the step vectors (O(nCSF*nLev)).
!
  call mma_allocate(Occ,SGS(iState)%nLev,nCSF,Label='Occ')
  do I=1,nCSF
    call CSF_Occupation(LSM,I,Occ(:,I))
  end do
  nScreened = 0
  nSampledOps = 0

  nOperatorPairs = 0
  nMatrixElements = 0
  nDiscrepancies = 0

  nPairDiagonal = 0
  nErrDiagonal = 0

  nPairOrientation = 0
  nErrOrientation = 0

  write(u6,*)
  write(u6,*) '========================================='
  write(u6,*) ' GENERAL TWO-OPERATOR TEST'
! List the topologies from TopologyName, so the header stays in step with
! the topologies actually tested.
  write(u6,'(2X,*(A,:,", "))')(TopologyName(iTopology),iTopology=1,nTopology)
  write(u6,*) ' INCLUDING DIAGONAL OPERATORS'
! The products E_pq E_rs are restricted to Gamma(p)Gamma(q) = Gamma(r)Gamma(s)
! (totally symmetric, as in H); each operator may change the symmetry.
  write(u6,*) ' TOTALLY SYMMETRIC PRODUCTS (SINGLE OPERATORS OF ANY SYMMETRY)'
  write(u6,*) '========================================='
  write(u6,'(A,I0)') 'Number of graph levels: ', &
    SGS(iState)%nLev
  write(u6,'(A,I0)') 'Split-graph midlevel:   ', &
    SGS(iState)%MidLev
  write(u6,'(A,I0)') 'Number of CSFs:         ', &
    nCSF
  write(u6,*)

!
! Test all nine operator-product topologies:
!
!     UU: E_pq(upper) E_rs(upper)
!     LL: E_pq(lower) E_rs(lower)
!     UL: E_pq(upper) E_rs(lower)
!     LU: E_pq(lower) E_rs(upper)
!     UO: E_pq(upper) E_rs(open, i.e. crossing the midlevel)
!     OU: E_pq(open) E_rs(upper)
!     LO: E_pq(lower) E_rs(open)
!     OL: E_pq(open) E_rs(lower)
!     OO: E_pq(open) E_rs(open)
!
  do iTopology=1,nTopology

    ! Empty ranges by default: every valid topology sets them below and the
    ! default case aborts, but some compilers cannot see that (Abend does
    ! not return) and warn that they may be used uninitialized.
    pMin = 1; pMax = 0
    qMin = 1; qMax = 0
    rMin = 1; rMax = 0
    sMin = 1; sMax = 0

    select case (iTopology)

      case (iUU)
        pMin = SGS(iState)%MidLev+1; pMax = SGS(iState)%nLev
        qMin = pMin; qMax = pMax
        rMin = SGS(iState)%MidLev+1; rMax = SGS(iState)%nLev
        sMin = rMin; sMax = rMax
      case (iLL)
        pMin = 1; pMax = SGS(iState)%MidLev
        qMin = pMin; qMax = pMax
        rMin = 1; rMax = SGS(iState)%MidLev
        sMin = rMin; sMax = rMax
      case (iUL)
        pMin = SGS(iState)%MidLev+1; pMax = SGS(iState)%nLev
        qMin = pMin; qMax = pMax
        rMin = 1; rMax = SGS(iState)%MidLev
        sMin = rMin; sMax = rMax
      case (iLU)
        pMin = 1; pMax = SGS(iState)%MidLev
        qMin = pMin; qMax = pMax
        rMin = SGS(iState)%MidLev+1; rMax = SGS(iState)%nLev
        sMin = rMin; sMax = rMax
      case (iUO)
        pMin = SGS(iState)%MidLev+1; pMax = SGS(iState)%nLev
        qMin = pMin; qMax = pMax
        rMin = 1; rMax = SGS(iState)%nLev
        sMin = 1; sMax = SGS(iState)%nLev
      case (iOU)
        ! First operator open (filtered below), second operator upper.
        pMin = 1; pMax = SGS(iState)%nLev
        qMin = 1; qMax = SGS(iState)%nLev
        rMin = SGS(iState)%MidLev+1; rMax = SGS(iState)%nLev
        sMin = rMin; sMax = rMax
      case (iLO)
        ! First operator lower, second operator open (filtered below).
        pMin = 1; pMax = SGS(iState)%MidLev
        qMin = 1; qMax = SGS(iState)%MidLev
        rMin = 1; rMax = SGS(iState)%nLev
        sMin = 1; sMax = SGS(iState)%nLev
      case (iOL)
        ! First operator open (filtered below), second operator lower.
        pMin = 1; pMax = SGS(iState)%nLev
        qMin = 1; qMax = SGS(iState)%nLev
        rMin = 1; rMax = SGS(iState)%MidLev
        sMin = 1; sMax = SGS(iState)%MidLev
      case (iOO)
        ! Both operators open (filtered below).
        pMin = 1; pMax = SGS(iState)%nLev
        qMin = 1; qMax = SGS(iState)%nLev
        rMin = 1; rMax = SGS(iState)%nLev
        sMin = 1; sMax = SGS(iState)%nLev
!     Restrict to only canonical case r > s
!     rMin=SGS(iState)%MidLev+1; rMax=SGS(iState)%nLev
!     sMin=1; sMax=SGS(iState)%MidLev

      case default

        write(u6,'(A,I0)') &
          'Illegal operator-product topology: ',iTopology
        call Abend()

    end select

    write(u6,'(A,A)') &
      'Testing operator-product topology: ', &
      TopologyName(iTopology)

    !
    ! Enumerate every ordered operator product E_pq E_rs
    ! belonging to the current topology.
    !
    do p=pMin,pMax
      do q=qMin,qMax

        if ((iTopology == iOU) .or. (iTopology == iOL) .or. (iTopology == iOO)) then
          ! Enforce that the first operator is open
          if (p == q) cycle
          if ((p > SGS(iState)%MidLev) .eqv. &
              (q > SGS(iState)%MidLev)) cycle
        end if

        FirstDiagonal = p == q
        FirstCanonical = p > q

        !
        ! Diagonal operators are totally symmetric; the others may
        ! change the symmetry.
        !
        ! (Each operator may change the symmetry; only the product has
        ! to be totally symmetric, see the second operator.)

        do r=rMin,rMax
          do s=sMin,sMax

            if ((iTopology == iUO) .or. (iTopology == iLO) .or. (iTopology == iOO)) then
              ! Enforce that the second operator is open
              if (r == s) cycle
              if ((r > SGS(iState)%MidLev) .eqv. &
                  (s > SGS(iState)%MidLev)) cycle
            end if

            SecondDiagonal = r == s
            SecondCanonical = r > s

            ! Totally symmetric product: Gamma(p)Gamma(q) = Gamma(r)Gamma(s),
            ! as in H. The single operators may change the symmetry.
            if (Mul(SGS(iState)%ISM(r),SGS(iState)%ISM(s)) /= &
                Mul(SGS(iState)%ISM(p),SGS(iState)%ISM(q))) cycle

            nOperatorPairs(iTopology) = &
              nOperatorPairs(iTopology)+1

            !
            ! Classify according to diagonal status.
            !
            if (.not. FirstDiagonal .and. &
                .not. SecondDiagonal) then

              iDiagonalCase = iCC

              if (FirstCanonical .and. &
                  SecondCanonical) then

                iOrientationCase = iCanCan

              else if (.not. FirstCanonical .and. &
                       SecondCanonical) then

                iOrientationCase = iRevCan

              else if (FirstCanonical .and. &
                       .not. SecondCanonical) then

                iOrientationCase = iCanRev

              else

                iOrientationCase = iRevRev

              end if

              nPairOrientation( &
                iTopology,iOrientationCase) = &
                nPairOrientation( &
                iTopology,iOrientationCase)+1

            else if (FirstDiagonal .and. &
                     .not. SecondDiagonal) then

              iDiagonalCase = iDC
              iOrientationCase = 0

            else if (.not. FirstDiagonal .and. &
                     SecondDiagonal) then

              iDiagonalCase = iCD
              iOrientationCase = 0

            else

              iDiagonalCase = iDD
              iOrientationCase = 0

            end if

            nPairDiagonal( &
              iTopology,iDiagonalCase) = &
              nPairDiagonal( &
              iTopology,iDiagonalCase)+1

            !
            ! Compare the CSF matrices: all (I,J) pairs, or a random
            ! sample of them if there are more than MaxPairs.
            !
            nPairsIJ = nCSF*nCSF
            Sampled = nPairsIJ > MaxPairs
            ! Cap the draws relative to the block size, so that sampling a block
            ! just above MaxPairs never costs more than two full loops.
            nDrawMax = merge(min(MaxDraw,2*nPairsIJ),nPairsIJ,Sampled)
            nTakeZ = 0
            nTakeNZ = 0
            if (Sampled) nSampledOps(iTopology) = nSampledOps(iTopology)+1
            do iDraw=1,nDrawMax
              if (Sampled) then
                if (nTakeZ >= QuotaZ .and. nTakeNZ >= QuotaNZ) exit
                I = Rand_Index(nCSF)
                J = Rand_Index(nCSF)
              else
                I = (iDraw-1)/nCSF+1
                J = mod(iDraw-1,nCSF)+1
              end if

              ! Skip the expensive reference when the occupations
              ! already show that <I|EpqErs|J> is zero. The fast
              ! routine is still called and must then return 0.
              Screened = Quick_Zero(Occ(:,I),Occ(:,J),p,q,r,s)
              if (Sampled) then
                if (Screened) then
                  if (nTakeZ >= QuotaZ) cycle
                  nTakeZ = nTakeZ+1
                else
                  if (nTakeNZ >= QuotaNZ) cycle
                  nTakeNZ = nTakeNZ+1
                end if
              end if
              if (Screened) then
                XRef = Zero
                nScreened(iTopology) = nScreened(iTopology)+1
#ifdef _CHECK_QUICK_ZERO_
                ! Verify the quick zero test itself (as expensive as
                ! no screening at all; define it only to test it).
                XChk = MatEl_EpqErs_Ref(LSM,I,J,p,q,r,s)
                if (abs(XChk) > Threshold) then
                  write(u6,'(A,6I6,ES22.14)') &
                    'Quick_Zero failed: I,J,p,q,r,s,Ref = ', &
                    I,J,p,q,r,s,XChk
                  call Abend()
                end if
#endif
              else
                XRef = MatEl_EpqErs_Ref( &
                       LSM, &
                       I,J, &
                       p,q, &
                       r,s)
              end if

              ! Set in every valid case below; the default aborts.
              XTest = Zero
              select case (iTopology)

                case (iUU)

                  XTest = MatEl_EpqErs_UU( &
                          LSM, &
                          I,J, &
                          p,q, &
                          r,s)

                case (iLL)

                  XTest = MatEl_EpqErs_LL( &
                          LSM, &
                          I,J, &
                          p,q, &
                          r,s)

                case (iUL)

                  XTest = MatEl_EpqErs_UL( &
                          LSM, &
                          I,J, &
                          p,q, &
                          r,s)

                case (iLU)

                  XTest = MatEl_EpqErs_LU( &
                          LSM, &
                          I,J, &
                          p,q, &
                          r,s)

                case (iUO)

                  if (I == 2 .and. J == 10 .and. &
                      p == 3 .and. q == 4 .and. r == 1 .and. s == 3) then
                    write(u6,*) 'Test <I|Epq|K>'
                    do K=1,nCSF
                      if (K /= 3) cycle
                      XTest = MatEl_Epq(SGS(istate),CIS(istate),EXS(istate), &
                                        p,q,LSM,I,K,.true.)
                      write(u6,*) 'I,p,q,K,<I|Epq|K>=',I,p,q,K,XTest
                    end do
                    write(u6,*) 'Test <K|Ers|J>'
                    do K=1,nCSF
                      if (K /= 3) cycle
                      XTest = MatEl_Epq(SGS(istate),CIS(istate),EXS(istate), &
                                        r,s,LSM,K,J,.true.)
                      write(u6,*) 'K,r,s,J,<K|Ers|J>=',K,r,s,J,XTest
                    end do
                  end if

                  XTest = MatEl_EpqErs_UO(LSM,I,J,p,q,r,s)

                case (iOU)

                  XTest = MatEl_EpqErs_OU(LSM,I,J,p,q,r,s)

                case (iLO)

                  XTest = MatEl_EpqErs_LO(LSM,I,J,p,q,r,s)

                case (iOL)

                  XTest = MatEl_EpqErs_OL(LSM,I,J,p,q,r,s)

                case (iOO)

                  XTest = MatEl_EpqErs_OO(LSM,I,J,p,q,r,s)

                case default

                  write(u6,'(A,I0)') &
                    'Illegal topology in comparison: ', &
                    iTopology
                  call Abend()

              end select

              nMatrixElements(iTopology) = &
                nMatrixElements(iTopology)+1

              Difference = XRef-XTest

              ! A screened element is exactly zero, so the fast
              ! routine must return exactly zero, not just a small
              ! number.
              if (Screened) then
                if (XTest == Zero) cycle
              else
                if (abs(Difference) <= Threshold) cycle
              end if

              nDiscrepancies(iTopology) = &
                nDiscrepancies(iTopology)+1

              nErrDiagonal( &
                iTopology,iDiagonalCase) = &
                nErrDiagonal( &
                iTopology,iDiagonalCase)+1

              if (iDiagonalCase == iCC) then

                nErrOrientation( &
                  iTopology,iOrientationCase) = &
                  nErrOrientation( &
                  iTopology,iOrientationCase)+1

              end if

              !
              ! Columns:
              !
              !   topology
              !   diagonal-status class
              !   orientation class
              !   p q r s
              !   I J
              !   XRef XTest Difference
              !
              if (iDiagonalCase == iCC) then

                write(u6, &
                      '(A,1X,A,1X,A,1X,4I4,'// &
                      '1X,2I7,3(1X,ES22.14))') &
                  TopologyName(iTopology), &
                  DiagonalName(iDiagonalCase), &
                  OrientationName( &
                  iOrientationCase), &
                  p,q,r,s, &
                  I,J, &
                  XRef,XTest,Difference

              else

                write(u6, &
                      '(A,1X,A,1X,A,1X,4I4,'// &
                      '1X,2I7,3(1X,ES22.14))') &
                  TopologyName(iTopology), &
                  DiagonalName(iDiagonalCase), &
                  '  ', &
                  p,q,r,s, &
                  I,J, &
                  XRef,XTest,Difference

              end if

            end do

          end do
        end do

      end do
    end do

    write(u6,*)

  end do

!
! Print separate summaries for UU, LL, UL, and LU.
!
  write(u6,*)
  write(u6,*) '========================================='
  write(u6,*) ' TWO-OPERATOR TEST SUMMARY'
  write(u6,*) '========================================='

  do iTopology=1,nTopology

    write(u6,*)
    write(u6,'(A,A)') &
      'Topology:               ', &
      TopologyName(iTopology)

    write(u6,'(A,I0)') &
      'Operator pairs:         ', &
      nOperatorPairs(iTopology)

    write(u6,'(A,I0)') &
      'Matrix elements:        ', &
      nMatrixElements(iTopology)

    write(u6,'(A,I0)') &
      'Screened as zero:       ', &
      nScreened(iTopology)

    write(u6,'(A,I0)') &
      'Operator pairs sampled: ', &
      nSampledOps(iTopology)

    write(u6,'(A,I0)') &
      'Discrepancies:          ', &
      nDiscrepancies(iTopology)

    write(u6,*)

    write(u6,'(A,I0,A,I0)') &
      'Closed x closed:        pairs = ', &
      nPairDiagonal(iTopology,iCC), &
      ', discrepancies = ', &
      nErrDiagonal(iTopology,iCC)

    write(u6,'(A,I0,A,I0)') &
      'Diagonal x closed:      pairs = ', &
      nPairDiagonal(iTopology,iDC), &
      ', discrepancies = ', &
      nErrDiagonal(iTopology,iDC)

    write(u6,'(A,I0,A,I0)') &
      'Closed x diagonal:      pairs = ', &
      nPairDiagonal(iTopology,iCD), &
      ', discrepancies = ', &
      nErrDiagonal(iTopology,iCD)

    write(u6,'(A,I0,A,I0)') &
      'Diagonal x diagonal:    pairs = ', &
      nPairDiagonal(iTopology,iDD), &
      ', discrepancies = ', &
      nErrDiagonal(iTopology,iDD)

    write(u6,*)

    write(u6,'(A,I0,A,I0)') &
      'Canonical x canonical:  pairs = ', &
      nPairOrientation(iTopology,iCanCan), &
      ', discrepancies = ', &
      nErrOrientation(iTopology,iCanCan)

    write(u6,'(A,I0,A,I0)') &
      'Reverse x canonical:    pairs = ', &
      nPairOrientation(iTopology,iRevCan), &
      ', discrepancies = ', &
      nErrOrientation(iTopology,iRevCan)

    write(u6,'(A,I0,A,I0)') &
      'Canonical x reverse:    pairs = ', &
      nPairOrientation(iTopology,iCanRev), &
      ', discrepancies = ', &
      nErrOrientation(iTopology,iCanRev)

    write(u6,'(A,I0,A,I0)') &
      'Reverse x reverse:      pairs = ', &
      nPairOrientation(iTopology,iRevRev), &
      ', discrepancies = ', &
      nErrOrientation(iTopology,iRevRev)

  end do

  write(u6,*)
  if (SeedReported) then
    write(u6,'(A,I0,A)') 'Random sampling used, seed = ',iRandSeed, &
      '  (reproduce with >export MOLCAS_RANDOM_SEED=<seed>)'
    write(u6,*)
  end if
  write(u6,*) '========================================='
  write(u6,*) ' END OF TWO-OPERATOR TEST'
  write(u6,*) '========================================='
  write(u6,*)

! <I|H|J> from MatEl_Epq and the MatEl_EpqErs_* routines, against sg_h_psi.
  if (sum(nDiscrepancies) == 0) call Check_MatEl_H()

  call mma_deallocate(Occ)

  if (sum(nDiscrepancies) /= 0) then
    write(u6,'(A,I0)') 'SG_Test_EpqErs_Ref: total discrepancies = ',sum(nDiscrepancies)
    call Abend()
  end if

contains

#ifdef _DEBUGPRINT_
  subroutine Print_CSF_Walks(SGS,CIS,ISYCI)

    use Symmetry_Info, only: Mul
    use Definitions, only: iwp, u6

    implicit none

    type(SGStruct), intent(in) :: SGS
    type(CIStruct), intent(in) :: CIS

    integer(kind=iwp), intent(in) :: ISYCI

    integer(kind=iwp) :: MV
    integer(kind=iwp) :: ISYUP
    integer(kind=iwp) :: ISYDWN

    integer(kind=iwp) :: IUP
    integer(kind=iwp) :: IDWN

    integer(kind=iwp) :: NUP
    integer(kind=iwp) :: NDWN

    integer(kind=iwp) :: IUW0
    integer(kind=iwp) :: IDW0

    integer(kind=iwp) :: ICUPOS
    integer(kind=iwp) :: ICDPOS

    integer(kind=iwp) :: ICUP
    integer(kind=iwp) :: ICDWN

    integer(kind=iwp) :: LEV
    integer(kind=iwp) :: NNN
    integer(kind=iwp) :: IC1, I

    integer(kind=iwp) :: ICS(SGS%nLev)

    character(len=1), parameter :: CODE(0:3) = ['0','u','d','2']
    character(len=200) :: LINE

    write(u6,*)
    write(u6,*) '========================================='
    write(u6,*) ' CSF WALKS'
    write(u6,*) '========================================='
    write(u6,*)

    I = 0
    do MV=1,CIS%nMidV

      do ISYUP=1,SGS%nSym

        ISYDWN = Mul(ISYUP,ISYCI)

        NUP = CIS%NOW(1,ISYUP,MV)
        NDWN = CIS%NOW(2,ISYDWN,MV)

        if (NUP*NDWN == 0) cycle

        write(u6,*)
        write(u6,'(A,4I6)') &
          'MV,ISYUP,ISYDWN,NUP=', &
          MV,ISYUP,ISYDWN,NUP

        write(u6,'(A,I6)') 'NDWN=',NDWN

        IUW0 = 1-CIS%nIpWlk+CIS%IOW(1,ISYUP,MV)
        IDW0 = 1-CIS%nIpWlk+CIS%IOW(2,ISYDWN,MV)

        do IDWN=1,NDWN

          ICDPOS = IDW0+IDWN*CIS%nIpWlk
          ICDWN = CIS%ICase(ICDPOS)

          NNN = 0

          do LEV=1,SGS%MidLev

            NNN = NNN+1

            if (NNN == CSFG_nPack+1) then
              NNN = 1
              ICDPOS = ICDPOS+1
              ICDWN = CIS%ICase(ICDPOS)
            end if

            IC1 = ICDWN/4
            ICS(LEV) = ICDWN-4*IC1
            ICDWN = IC1

          end do

          do IUP=1,NUP

            ICUPOS = IUW0+IUP*CIS%nIpWlk
            ICUP = CIS%ICase(ICUPOS)

            NNN = 0

            do LEV=SGS%MidLev+1,SGS%nLev

              NNN = NNN+1

              if (NNN == CSFG_nPack+1) then
                NNN = 1
                ICUPOS = ICUPOS+1
                ICUP = CIS%ICase(ICUPOS)
              end if

              IC1 = ICUP/4
              ICS(LEV) = ICUP-4*IC1
              ICUP = IC1

            end do
            I = I+1

            LINE = ' '

            write(LINE(1:25),'(A,I3,A,I3,A,I3,A,I3,A,I3,A)') &
              'I=',I,'(',MV,':',ISYUP,':',IUP,'/',IDWN,')'

            do LEV=1,SGS%nLev
              LINE(30+LEV:30+LEV) = CODE(ICS(LEV))
            end do

            write(u6,'(A)') trim(LINE)

          end do

        end do

      end do
    end do

    write(u6,*)

  end subroutine Print_CSF_Walks
#endif

  function MatEl_Epq_Ref(iState,ISYI,I,J,p,q) result(value)

    use Symmetry_Info, only: Mul
    use Constants, only: Zero, One
    use Definitions, only: iwp, wp
    use stdalloc, only: mma_allocate, mma_deallocate

    integer(kind=iwp), intent(in) :: iState
    integer(kind=iwp), intent(in) :: ISYI
    integer(kind=iwp), intent(in) :: I, J
    integer(kind=iwp), intent(in) :: p, q

    real(kind=wp) :: value

    integer(kind=iwp) :: iSYJ, nCSFI, nCSFJ

    real(kind=wp), allocatable :: Ket(:)
    real(kind=wp), allocatable :: Res(:)
    logical(kind=iwp) :: Reuse_Save

    value = Zero

    Reuse_Save = EXS(iState)%ReUse_SGTmp
    EXS(iState)%ReUse_SGTmp = .false.

    ISYJ = Mul(Mul(SGS(iState)%ISm(p),SGS(iState)%iSM(q)),ISYI)

    nCSFI = CIS(iState)%NCSF(iSYI)
    nCSFJ = CIS(iState)%NCSF(iSYJ)

    call mma_allocate(Ket,nCSFI,Label='Ket')
    call mma_allocate(Res,nCSFJ,Label='Res')

    Ket(:) = Zero
    Ket(J) = One    ! C_J=1.00

    Res(:) = Zero

    call SG_Epq_Psi( &
      SGS(iState), &
      CIS(iState), &
      EXS(iState), &
      p,q,One,ISYI, &
      Ket,Res)

    value = Res(I)

    call mma_deallocate(Ket)
    call mma_deallocate(Res)

    EXS(iState)%ReUse_SGTmp = Reuse_Save

  end function MatEl_Epq_Ref

!
! <I|H|J> for I and J of symmetry LSM, in the form used by sg_h_psi:
!
!   H = Sum_{i>=j, G(ij)=1} k_ij Ehat_ij + Sum_{i>=j} Sum_{k>=l, G(kl)=G(ij)} 1/2 (kl|ij) Ehat_kl Ehat_ij
!
! with Ehat_ij = E_ij + E_ji (E_ii for i = j), k_ij = h_ij - 1/2 Sum_k (ik|kj),
! and triangularly packed TU (h) and TUVX ((kl|ij)), as in sg_h_psi.
!
  function MatEl_H(I,J,TU,nTU,TUVX,nTUVX) result(H)

    use Index_functions, only: iTri
    use Constants, only: Half

    integer(kind=iwp), intent(in) :: I, J, nTU, nTUVX
    real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX)
    real(kind=wp) :: H
    integer(kind=iwp) :: iO, jO, kO, lO, ij, kl, nL
    real(kind=wp) :: OneInt, TwoInt, G

    H = Zero
    nL = SGS(iState)%nLev
    do iO=1,nL
      do jO=1,iO
        ij = iTri(iO,jO)
        if (SGS(iState)%ISm(iO) == SGS(iState)%ISm(jO)) then
          TwoInt = Zero
          do kO=1,nL
            TwoInt = TwoInt+TUVX(iTri(iTri(iO,kO),iTri(kO,jO)))
          end do
          OneInt = TU(ij)-Half*TwoInt
          H = H+OneInt*MatEl_Epq(SGS(iState),CIS(iState),EXS(iState),iO,jO,LSM,I,J,.false.)
          if (iO /= jO) H = H+OneInt*MatEl_Epq(SGS(iState),CIS(iState),EXS(iState),jO,iO,LSM,I,J,.false.)
        end if
        do kO=1,nL
          do lO=1,kO
            kl = iTri(kO,lO)
            if (Mul(SGS(iState)%ISm(kO),SGS(iState)%ISm(lO)) /= &
                Mul(SGS(iState)%ISm(iO),SGS(iState)%ISm(jO))) cycle
            G = Half*TUVX(iTri(kl,ij))
            if (G == Zero) cycle
            ! Ehat_kl Ehat_ij
            H = H+G*MatEl_EpqErs(I,J,kO,lO,iO,jO,Occ(:,I),Occ(:,J))
            if (iO /= jO) H = H+G*MatEl_EpqErs(I,J,kO,lO,jO,iO,Occ(:,I),Occ(:,J))
            if (kO /= lO) then
              H = H+G*MatEl_EpqErs(I,J,lO,kO,iO,jO,Occ(:,I),Occ(:,J))
              if (iO /= jO) H = H+G*MatEl_EpqErs(I,J,lO,kO,jO,iO,Occ(:,I),Occ(:,J))
            end if
          end do
        end do
      end do
    end do

  end function MatEl_H

!
! <I|H|J> in a RAS graph, in the operator-ordered form of sg_h_psi_ras
! (Olsen et al. 1988), so that the products <I|E_Q E_P|J> = Sum_K
! <I|E_Q|K><K|E_P|J> over the graph CSFs K are exact:
!
!   H = Sum_ab h'_ab E_ab + Sum_{P<=Q} f_PQ (Q|P) E_Q E_P,  P = (i,j) -> (i-1)*n+j,
!   h'_ab = h_ab - Sum_{m<a} (am|mb) - [a>b] (aa|ab) - 1/2 delta_ab (aa|aa),
!
! f_PP = 1/2, else 1; TU and TUVX indexed by level, as in MatEl_H.
!
  function MatEl_H_RAS(I,J,TU,nTU,TUVX,nTUVX) result(H)

    use Index_functions, only: iTri
    use Constants, only: Half

    integer(kind=iwp), intent(in) :: I, J, nTU, nTUVX
    real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX)
    real(kind=wp) :: H
    integer(kind=iwp) :: iO, jO, kO, lO, mO, ij, kl, nL
    real(kind=wp) :: OneInt, G

    H = Zero
    nL = SGS(iState)%nLev
    do iO=1,nL
      do jO=1,nL
        if (SGS(iState)%ISm(iO) == SGS(iState)%ISm(jO)) then
          OneInt = TU(iTri(iO,jO))
          do mO=1,iO-1
            OneInt = OneInt-TUVX(iTri(iTri(iO,mO),iTri(mO,jO)))
          end do
          if (iO > jO) OneInt = OneInt-TUVX(iTri(iTri(iO,iO),iTri(iO,jO)))
          if (iO == jO) OneInt = OneInt-Half*TUVX(iTri(iTri(iO,iO),iTri(iO,iO)))
          H = H+OneInt*MatEl_Epq(SGS(iState),CIS(iState),EXS(iState),iO,jO,LSM,I,J,.false.)
        end if
        ij = (iO-1)*nL+jO
        do kO=iO,nL
          do lO=1,nL
            kl = (kO-1)*nL+lO
            if (kl < ij) cycle
            if (Mul(SGS(iState)%ISm(kO),SGS(iState)%ISm(lO)) /= &
                Mul(SGS(iState)%ISm(iO),SGS(iState)%ISm(jO))) cycle
            G = TUVX(iTri(iTri(kO,lO),iTri(iO,jO)))
            if (kl == ij) G = Half*G
            if (G == Zero) cycle
            H = H+G*MatEl_EpqErs(I,J,kO,lO,iO,jO,Occ(:,I),Occ(:,J))
          end do
        end do
      end do
    end do

  end function MatEl_H_RAS

!
! Verify MatEl_H against sg_h_psi: with reproducible pseudo-random integrals
! (any values with the permutational symmetry of the triangular packing will
! do), compare <I|H|J> with (H e_J)_I for all I and for all columns J, or a
! random sample of MaxHCol columns. Aborts on a mismatch.
!
  subroutine Check_MatEl_H()

    integer(kind=iwp), parameter :: MaxHCol = 20
    real(kind=wp), parameter :: TolH = 1.0e-10_wp
    integer(kind=iwp) :: nL, nTU, nTUVX, nC, nCol, iCol, I, J, k, nTestH, nFailH, nFailD
    integer(kind=int64) :: iInt
    real(kind=wp), allocatable :: TU(:), TUVX(:), Psi(:), Sgm(:), TU_O(:), TUVX_O(:)
    real(kind=wp) :: X, Dev, MaxDev, MaxDevD, T0, T1, W, TimH, TimD
    logical(kind=iwp) :: Reuse_Save, AllCols, RAS
    integer(kind=iwp) :: p, q, r, s, pq, rs
    integer(kind=iwp), allocatable :: L2Orb(:)

    nL = SGS(iState)%nLev
    nTU = nL*(nL+1)/2
    nTUVX = nTU*(nTU+1)/2
    nC = CIS(iState)%NCSF(LSM)

    ! Reproducible integrals in [-0.5,0.5), independent of the sampling seed.
    call mma_allocate(TU,nTU,Label='TU_H')
    call mma_allocate(TUVX,nTUVX,Label='TUVX_H')
    iInt = 271828183_int64
    do k=1,nTU
      iInt = mod(48271_int64*iInt,2147483647_int64)
      TU(k) = real(iInt,kind=wp)/2147483647.0_wp-0.5_wp
    end do
    do k=1,nTUVX
      iInt = mod(48271_int64*iInt,2147483647_int64)
      TUVX(k) = real(iInt,kind=wp)/2147483647.0_wp-0.5_wp
    end do
    ! new integrals for MatEl_H_Direct
    call Reset_MatEl_H()

    ! RAS graph: the reference is sg_h_psi_ras (sg_h_psi needs intermediates
    ! outside the RAS space); it reads the integrals in symmetry-blocked
    ! orbital order, MatEl_H(_Direct) by level.
    RAS = SGS(iState)%nRsPrt /= 1
    if (RAS) then
      call mma_allocate(L2Orb,nL,Label='L2Orb_H')
      call SG_Level_Map(SGS(iState),L2Orb)
      call mma_allocate(TU_O,nTU,Label='TU_O')
      call mma_allocate(TUVX_O,nTUVX,Label='TUVX_O')
      do p=1,nL
        do q=1,p
          pq = Tri_SG(p,q)
          TU_O(Tri_SG(L2Orb(p),L2Orb(q))) = TU(pq)
          do r=1,nL
            do s=1,r
              rs = Tri_SG(r,s)
              if (rs > pq) cycle
              TUVX_O(Tri_SG(Tri_SG(L2Orb(p),L2Orb(q)),Tri_SG(L2Orb(r),L2Orb(s)))) = TUVX(Tri_SG(pq,rs))
            end do
          end do
        end do
      end do
      call mma_deallocate(L2Orb)
    end if

    call mma_allocate(Psi,nC,Label='Psi_H')
    call mma_allocate(Sgm,nC,Label='Sgm_H')
    Reuse_Save = EXS(iState)%ReUse_SGTMP
    EXS(iState)%ReUse_SGTMP = .false.

    AllCols = nC <= MaxHCol
    nCol = merge(nC,MaxHCol,AllCols)
    nTestH = 0
    nFailH = 0
    nFailD = 0
    MaxDev = Zero
    MaxDevD = Zero
    TimH = Zero
    TimD = Zero
    do iCol=1,nCol
      J = merge(iCol,Rand_Index(nC),AllCols)
      Psi(:) = Zero
      Psi(J) = One
      if (RAS) then
        call HPsi_RAS(SGS(iState),CIS(iState),EXS(iState),Psi,nC,LSM,Sgm,TUVX_O,nTUVX,TU_O,nTU)
      else
        call HPsi(SGS(iState),CIS(iState),EXS(iState),Psi,nC,LSM,Sgm,TUVX,nTUVX,TU,nTU)
      end if
      do I=1,nC
        call CWTime(T0,W)
        if (RAS) then
          X = MatEl_H_RAS(I,J,TU,nTU,TUVX,nTUVX)
        else
          X = MatEl_H(I,J,TU,nTU,TUVX,nTUVX)
        end if
        call CWTime(T1,W)
        TimH = TimH+T1-T0
        Dev = abs(X-Sgm(I))
        MaxDev = max(MaxDev,Dev)
        nTestH = nTestH+1
        if (Dev > TolH*max(One,abs(Sgm(I)))) then
          nFailH = nFailH+1
          if (nFailH <= 10) write(u6,'(A,2I8,3ES22.13)') 'MatEl_H mismatch: I,J,<I|H|J>,sg_h_psi,diff = ', &
            I,J,X,Sgm(I),X-Sgm(I)
        end if
        call CWTime(T0,W)
        X = MatEl_H_Direct(I,J,Occ(:,I),Occ(:,J),TU,nTU,TUVX,nTUVX)
        call CWTime(T1,W)
        TimD = TimD+T1-T0
        Dev = abs(X-Sgm(I))
        MaxDevD = max(MaxDevD,Dev)
        if (Dev > TolH*max(One,abs(Sgm(I)))) then
          nFailD = nFailD+1
          if (nFailD <= 10) write(u6,'(A,2I8,3ES22.13)') 'MatEl_H_Direct mismatch: I,J,<I|H|J>,sg_h_psi,diff = ', &
            I,J,X,Sgm(I),X-Sgm(I)
        end if
      end do
    end do

    EXS(iState)%ReUse_SGTMP = Reuse_Save
    call mma_deallocate(Sgm)
    call mma_deallocate(Psi)
    call mma_deallocate(TUVX)
    call mma_deallocate(TU)
    if (RAS) then
      call mma_deallocate(TUVX_O)
      call mma_deallocate(TU_O)
    end if

    write(u6,*)
    write(u6,'(A)') merge(' MatEl_H verification against sg_h_psi_ras (pseudo-random integrals)', &
                          ' MatEl_H verification against sg_h_psi     (pseudo-random integrals)',RAS)
    write(u6,'(A,I0,A,I0,A)') 'Columns tested         : ',nCol,' of ',nC,merge(' (all)   ',' (random)',AllCols)
    write(u6,'(A,I0)') 'Matrix elements tested : ',nTestH
    write(u6,'(A,2I8)') 'Failures (H, Direct)   : ',nFailH,nFailD
    write(u6,'(A,2ES12.4)') 'Maximum abs. deviation : ',MaxDev,MaxDevD
    write(u6,'(A,2F12.3)') 'CPU time (H, Direct)   : ',TimH,TimD
    if ((nFailH /= 0) .or. (nFailD /= 0)) then
      write(u6,'(A)') 'MatEl_H consistency check FAILED.'
      call Abend()
    end if
    write(u6,'(A)') 'MatEl_H consistency check PASSED.'

  end subroutine Check_MatEl_H

!> Direct matrix element of a spin-free excitation operator between two CSFs.
!!
!! Computes EIJ = < I | E(IP,IQ) | J >, where I belongs to symmetry
!! Mul(ISm(IP),ISm(IQ),ISYCJ) and J belongs to ISYCJ.
!!
!! This is intentionally a reference implementation:
!!   * no SGTMP cache or saved state;
!!   * no sorting or mutation of ICOUP;
!!   * no packing, tiling, SIMD, GEMM, or scratch lists;
!!   * only the coupling entries needed for the requested CSFs are visited.
!!
!! I and J are 1-based CSF indices within their respective symmetry sectors.
!! CPQ is not included: the returned value is the bare <I|E_pq|J>.

!> Exhaustively compare MatEl_Epq and MatEl_Epq_Ref for one ket symmetry.
!!
!! The routine checks every orbital pair (p,q), every ket CSF J in
!! symmetry ISYCI, and every bra CSF I in the symmetry generated by E_pq.
!! It prints each inconsistency and a final summary.
  subroutine Check_MatEl_Epq(istate,ISYCI,Tol,nTest,nFail,MaxDiff)
    use Symmetry_Info, only: Mul
    use Constants, only: Zero
    use Definitions, only: wp, iwp
#include "intent.fh"

    integer(kind=iwp), intent(in) :: ISYCI, iState
    real(kind=wp), intent(in) :: Tol
    integer(kind=iwp), intent(out) :: nTest, nFail
    real(kind=wp), intent(out) :: MaxDiff

    integer(kind=iwp) :: I, J, p, q
    integer(kind=iwp) :: nCSFI, nCSFJ, ISYCJ
    integer(kind=iwp) :: iMax, jMax, pMax, qMax
    real(kind=wp) :: value, Value_Ref, Diff, Scale
    logical(kind=iwp) :: Reuse_Save, Verbose
    !
    ! Quick zero test: orbital occupations of the CSFs of every symmetry,
    ! OccAll(t,iOff(ISY)+I) = occupation of level t in CSF I of symmetry ISY.
    ! They are unpacked from the step vectors by CSF_Occupation, not taken
    ! from MatEl_Epq, which is what is being tested. (Taking them from the
    ! reference diagonal costs nLev*nCSF**2 and dominated large cases.)
    !
    integer(kind=iwp), allocatable :: OccAll(:,:)
    integer(kind=iwp) :: iOff(8), ISY, nScreen
#ifdef _CHECK_QUICK_ZERO_
    integer(kind=iwp) :: iLev
#endif
    logical(kind=iwp) :: Screened
    ! Sampling of (I,J) as in the two-operator test (parameters of the host).
    integer(kind=iwp) :: nPairsIJ, nDrawMax, iDraw, nTakeZ, nTakeNZ, nSampledPQ
    logical(kind=iwp) :: Sampled

    nTest = 0_iwp
    nFail = 0_iwp
    MaxDiff = Zero
    iMax = 0_iwp
    jMax = 0_iwp
    pMax = 0_iwp
    qMax = 0_iwp
    Verbose = .false.
! Verbose = .True.

    nCSFI = CIS(iState)%NCSF(ISYCI)

    write(u6,*) '============================'
    write(u6,*) '   MatEl_Epq Verification   '
    write(u6,*) '============================'
    write(u6,*) 'iState,iSYCI: ',iState,ISYCI
    write(u6,*) 'SGS(iState)%iSM(:): ',SGS(iState)%iSM(:)

    ! Disable reuse of SGTMP
    Reuse_Save = EXS(iState)%ReUse_SGTMP
    EXS(iState)%ReUse_SGTMP = .false.

    iOff(:) = 0
    do ISY=2,SGS(iState)%nSym
      iOff(ISY) = iOff(ISY-1)+CIS(iState)%NCSF(ISY-1)
    end do
    call mma_allocate(OccAll,SGS(iState)%nLev, &
                      iOff(SGS(iState)%nSym)+CIS(iState)%NCSF(SGS(iState)%nSym), &
                      Label='OccAll')
    do ISY=1,SGS(iState)%nSym
      do I=1,CIS(iState)%NCSF(ISY)
        call CSF_Occupation(ISY,I,OccAll(:,iOff(ISY)+I))
#ifdef _CHECK_QUICK_ZERO_
        ! Cross-check the unpacked occupations against the reference diagonal.
        do iLev=1,SGS(iState)%nLev
          if (OccAll(iLev,iOff(ISY)+I) /= nint(MatEl_Epq_Ref(istate,ISY,I,I,iLev,iLev))) then
            write(u6,'(A,3I7)') 'CSF_Occupation failed: ISY,I,level =',ISY,I,iLev
            call Abend()
          end if
        end do
#endif
      end do
    end do
    nScreen = 0_iwp
    nSampledPQ = 0_iwp

    do q=1,SGS(iState)%nLev
      do p=1,SGS(iState)%nLev

        ISYCJ = Mul(Mul(SGS(iState)%ISm(p),SGS(iState)%ISm(q)),ISYCI)
        nCSFJ = CIS(iState)%NCSF(ISYCJ)

        ! All (I,J) pairs, or a random sample of them if there are more than
        ! MaxPairs.
        nPairsIJ = nCSFI*nCSFJ
        Sampled = nPairsIJ > MaxPairs
        ! Cap the draws relative to the block size, so that sampling a block
        ! just above MaxPairs never costs more than two full loops.
        nDrawMax = merge(min(MaxDraw,2*nPairsIJ),nPairsIJ,Sampled)
        nTakeZ = 0
        nTakeNZ = 0
        if (Sampled) nSampledPQ = nSampledPQ+1_iwp
        do iDraw=1,nDrawMax
          if (Sampled) then
            if (nTakeZ >= QuotaZ .and. nTakeNZ >= QuotaNZ) exit
            I = Rand_Index(nCSFI)
            J = Rand_Index(nCSFJ)
          else
            J = (iDraw-1)/nCSFI+1
            I = mod(iDraw-1,nCSFI)+1
          end if

          ! Skip the reference when <I|Epq|J> is zero by the occupations;
          ! MatEl_Epq must then return exactly zero.
          Screened = Epq_Is_Zero(OccAll(:,iOff(ISYCI)+I),OccAll(:,iOff(ISYCJ)+J),I,J,p,q)
          if (Sampled) then
            if (Screened) then
              if (nTakeZ >= QuotaZ) cycle
              nTakeZ = nTakeZ+1
            else
              if (nTakeNZ >= QuotaNZ) cycle
              nTakeNZ = nTakeNZ+1
            end if
          end if

          ! evaluate <I| Epq |J>

          value = MatEl_Epq(SGS(iState), &
                            CIS(iState), &
                            EXS(iState), &
                            p,q,ISYCJ,I,J,Verbose)

          if (Screened) then
            Value_Ref = Zero
            nScreen = nScreen+1_iwp
#ifdef _CHECK_QUICK_ZERO_
            if (abs(MatEl_Epq_Ref(istate,ISYCJ,I,J,p,q)) > Tol) then
              write(u6,'(A,6I7)') 'Epq_Is_Zero failed: ISYCI,I,ISYJ,J,p,q =', &
                ISYCI,I,ISYCJ,J,p,q
              call Abend()
            end if
#endif
          else
            Value_Ref = MatEl_Epq_Ref(istate,ISYCJ,I,J,p,q)
          end if

          nTest = nTest+1_iwp
          Diff = abs(value-Value_Ref)
          Scale = max(1.0_wp,abs(value),abs(Value_Ref))

          if (Diff > MaxDiff) then
            MaxDiff = Diff
            iMax = I
            jMax = J
            pMax = p
            qMax = q
          end if

          if ((Screened .and. value /= Zero) .or. Diff > Tol*Scale) then
            nFail = nFail+1_iwp
            write(u6,'(A,6I7,2ES24.14)') &
              'MATEL_EPQ mismatch: ISYCI,I,ISYJ,J,p,q,value,Value_Ref =', &
              ISYCI,I,ISYCJ,J,p,q,value,Value_Ref
          end if

        end do
      end do
    end do

    EXS(iState)%ReUse_SGTMP = Reuse_Save
    call mma_deallocate(OccAll)

    write(u6,'(/,A,I0)') 'MatEl_Epq tests       : ',nTest
    write(u6,'(A,I0)') 'Screened as zero      : ',nScreen
    write(u6,'(A,I0)') 'Sampled (p,q) blocks  : ',nSampledPQ
    write(u6,'(A,I0)') 'MatEl_Epq failures    : ',nFail
    write(u6,'(A,ES12.4)') 'Relative tolerance    : ',Tol
    if (nFail > 0) then
      write(u6,'(A,ES24.14)') 'Maximum abs. diff.    : ',MaxDiff
      write(u6,'(A,4I8)') 'At I,J,p,q            : ', &
        iMax,jMax,pMax,qMax
    end if

    if (nFail == 0) then
      write(u6,'(A)') 'MatEl_Epq consistency check PASSED.'
    else
      write(u6,'(A)') 'MatEl_Epq consistency check FAILED.'
    end if
    write(u6,*) '============================'

  end subroutine Check_MatEl_Epq

!> Check of the ASTA walk products (ASTA_MatEl) for all symmetries and
!> level pairs: E_pq against MatEl_Epq and against the determinant
!> expansion of the CSFs (CSF_Det_Expand, Det_MatEl), and the triplet
!> operator T_pq = a+_pa a_qa - a+_pb a_qb against the determinant
!> expansion. All (I,J) pairs of a (p,q) block, or a random sample of
!> nSample of them if there are more than MaxPairs. Skipped if the
!> determinant expansion would be too large.
  subroutine Check_MatEl_Tpq(iState,nFail)
    integer(kind=iwp), intent(in) :: iState
    integer(kind=iwp), intent(out) :: nFail
    integer(kind=iwp), parameter :: MxDetTot = 4000000, nSample = 2000, MxPrt = 20
    real(kind=wp), parameter :: Tol = 1.0e-12_wp
    integer(kind=iwp) :: nLev, nSym, ISY, ISYI, ISYJ, I, J, k, p, q, iDraw, nDraw, nTest, iOffI, iOffJ, nDetTot, nD
    integer(kind=iwp) :: iCSF(8), nOpen, nNZ_T, nNZ_TOff
    integer(kind=iwp), allocatable :: Steps(:,:), DOff(:), DA(:), DB(:), XA(:), XB(:)
    real(kind=wp), allocatable :: DC(:), XC(:)
    real(kind=wp) :: eW, eR, eD, tW, tD, tT, dER, dED, dTD, dTT
    logical(kind=iwp) :: Sampled

    nFail = 0
    nLev = SGS(iState)%nLev
    nSym = SGS(iState)%nSym
    iCSF(:) = 0
    do ISY=2,nSym
      iCSF(ISY) = iCSF(ISY-1)+CIS(iState)%NCSF(ISY-1)
    end do
    k = iCSF(nSym)+CIS(iState)%NCSF(nSym)
    call mma_allocate(Steps,nLev,max(1,k),Label='Steps')
    call mma_allocate(DOff,k+1,Label='DOff')
    nDetTot = 0
    DOff(1) = 0
    do ISY=1,nSym
      do I=1,CIS(iState)%NCSF(ISY)
        call SG_CSF_StepVec(SGS(iState),CIS(iState),ISY,I,Steps(:,iCSF(ISY)+I))
        nOpen = count((Steps(:,iCSF(ISY)+I) == 1) .or. (Steps(:,iCSF(ISY)+I) == 2))
        nDetTot = nDetTot+2**min(nOpen,30)
      end do
    end do
    if (nDetTot > MxDetTot) then
      write(u6,'(A,I0,A)') ' MatEl_Tpq check skipped: ',nDetTot,' determinants in the expansion'
      call mma_deallocate(Steps)
      call mma_deallocate(DOff)
      return
    end if
    call mma_allocate(DA,max(1,nDetTot),Label='DA')
    call mma_allocate(DB,max(1,nDetTot),Label='DB')
    call mma_allocate(DC,max(1,nDetTot),Label='DC')
    do k=1,iCSF(nSym)+CIS(iState)%NCSF(nSym)
      call CSF_Det_Expand(nLev,Steps(:,k),nD,XA,XB,XC)
      if (abs(sum(XC(1:nD)**2)-One) > Tol) then
        write(u6,*) 'CSF_Det_Expand: CSF not normalized',k,sum(XC(1:nD)**2)
        call Abend()
      end if
      DA(DOff(k)+1:DOff(k)+nD) = XA(1:nD)
      DB(DOff(k)+1:DOff(k)+nD) = XB(1:nD)
      DC(DOff(k)+1:DOff(k)+nD) = XC(1:nD)
      DOff(k+1) = DOff(k)+nD
      call mma_deallocate(XA)
      call mma_deallocate(XB)
      call mma_deallocate(XC)
    end do

    nTest = 0
    nNZ_T = 0
    nNZ_TOff = 0
    dER = Zero
    dED = Zero
    dTD = Zero
    dTT = Zero
    do ISYJ=1,nSym
      do q=1,nLev
        do p=1,nLev
          ISYI = Mul(Mul(SGS(iState)%ISm(p),SGS(iState)%ISm(q)),ISYJ)
          if ((CIS(iState)%NCSF(ISYI) == 0) .or. (CIS(iState)%NCSF(ISYJ) == 0)) cycle
          Sampled = CIS(iState)%NCSF(ISYI)*CIS(iState)%NCSF(ISYJ) > MaxPairs
          nDraw = merge(nSample,CIS(iState)%NCSF(ISYI)*CIS(iState)%NCSF(ISYJ),Sampled)
          do iDraw=1,nDraw
            if (Sampled) then
              I = Rand_Index(CIS(iState)%NCSF(ISYI))
              J = Rand_Index(CIS(iState)%NCSF(ISYJ))
            else
              J = (iDraw-1)/CIS(iState)%NCSF(ISYI)+1
              I = mod(iDraw-1,CIS(iState)%NCSF(ISYI))+1
            end if
            iOffI = iCSF(ISYI)+I
            iOffJ = iCSF(ISYJ)+J
            eW = ASTA_MatEl(nLev,Steps(:,iOffI),Steps(:,iOffJ),p,q,.false.)
            tW = ASTA_MatEl(nLev,Steps(:,iOffI),Steps(:,iOffJ),p,q,.true.)
            eR = MatEl_Epq(SGS(iState),CIS(iState),EXS(iState),p,q,ISYJ,I,J,.false.)
            eD = Det_MatEl(DOff(iOffI+1)-DOff(iOffI),DA(DOff(iOffI)+1:DOff(iOffI+1)),DB(DOff(iOffI)+1:DOff(iOffI+1)), &
                           DC(DOff(iOffI)+1:DOff(iOffI+1)),DOff(iOffJ+1)-DOff(iOffJ),DA(DOff(iOffJ)+1:DOff(iOffJ+1)), &
                           DB(DOff(iOffJ)+1:DOff(iOffJ+1)),DC(DOff(iOffJ)+1:DOff(iOffJ+1)),p,q,One)
            tD = Det_MatEl(DOff(iOffI+1)-DOff(iOffI),DA(DOff(iOffI)+1:DOff(iOffI+1)),DB(DOff(iOffI)+1:DOff(iOffI+1)), &
                           DC(DOff(iOffI)+1:DOff(iOffI+1)),DOff(iOffJ+1)-DOff(iOffJ),DA(DOff(iOffJ)+1:DOff(iOffJ+1)), &
                           DB(DOff(iOffJ)+1:DOff(iOffJ+1)),DC(DOff(iOffJ)+1:DOff(iOffJ+1)),p,q,-One)
            tT = Table_MatEl_T(SGS(iState),CIS(iState),Steps(:,iOffI),Steps(:,iOffJ),p,q)
            dTT = max(dTT,abs(tT-tW))
            if (abs(tT-tW) > Tol) then
              nFail = nFail+1
              if (nFail <= MxPrt) write(u6,'(A,6I5,2ES20.12)') ' Table_MatEl_T mismatch: ISYI,I,ISYJ,J,p,q, table, walk =', &
                ISYI,I,ISYJ,J,p,q,tT,tW
            end if
            nTest = nTest+1
            if (abs(tD) > 1.0e-10_wp) then
              nNZ_T = nNZ_T+1
              ! bra and ket differ above the loop head (top segments)
              if (any(Steps(max(p,q)+1:nLev,iOffI) /= Steps(max(p,q)+1:nLev,iOffJ))) nNZ_TOff = nNZ_TOff+1
            end if
            dER = max(dER,abs(eW-eR))
            dED = max(dED,abs(eW-eD))
            dTD = max(dTD,abs(tW-tD))
            if ((abs(eW-eR) > Tol) .or. (abs(eW-eD) > Tol) .or. (abs(tW-tD) > Tol)) then
              nFail = nFail+1
              if (nFail <= MxPrt) then
                write(u6,'(A,6I5)') ' MatEl_Tpq mismatch: ISYI,I,ISYJ,J,p,q =',ISYI,I,ISYJ,J,p,q
                write(u6,'(A,5ES20.12)') '   E walk, MatEl_Epq, det; T walk, det:',eW,eR,eD,tW,tD
                write(u6,'(A,30I1)') '   bra steps (tail first): ',Steps(:,iOffI)
                write(u6,'(A,30I1)') '   ket steps (tail first): ',Steps(:,iOffJ)
              end if
            end if
          end do
        end do
      end do
    end do
    write(u6,'(/,A,I0)') ' MatEl_Tpq (ASTA walk products) tests: ',nTest
    write(u6,'(A,ES10.3)') '   E_pq walk vs MatEl_Epq,     max. abs. difference ',dER
    write(u6,'(A,ES10.3)') '   E_pq walk vs determinants,  max. abs. difference ',dED
    write(u6,'(A,ES10.3)') '   T_pq walk vs determinants,  max. abs. difference ',dTD
    write(u6,'(A,ES10.3)') '   T_pq segment tables vs walk, max. abs. difference ',dTT
    write(u6,'(A,I0,A,I0,A)') '   nonzero T_pq elements: ',nNZ_T,' (',nNZ_TOff,' with different walks above the loop head)'
    write(u6,'(A,I0)') '   failures: ',nFail
    call Check_Xpq_Psi(iState,nFail)
    call mma_deallocate(Steps)
    call mma_deallocate(DOff)
    call mma_deallocate(DA)
    call mma_deallocate(DB)
    call mma_deallocate(DC)
  end subroutine Check_MatEl_Tpq

!> SG_Xpq_Psi against SG_Epq_Psi (E_pq) and against the dense product with
!> ASTA_MatEl (T_pq), for all symmetries and level pairs with at most
!> MaxPairs (I,J) pairs, with a fixed test vector.
  subroutine Check_Xpq_Psi(iState,nFail)
    integer(kind=iwp), intent(in) :: iState
    integer(kind=iwp), intent(inout) :: nFail
    real(kind=wp), parameter :: Tol = 1.0e-12_wp
    integer(kind=iwp) :: nLev, ISYI, ISYJ, I, J, p, q, nI, nJ, nTest
    integer(kind=iwp), allocatable :: SI(:), SJ(:)
    real(kind=wp), allocatable :: C(:), S1(:), S2(:)
    real(kind=wp) :: dE, dT, dTR
    logical(kind=iwp) :: Reuse_Save
    nLev = SGS(iState)%nLev
    call mma_allocate(SI,nLev,Label='SI')
    call mma_allocate(SJ,nLev,Label='SJ')
    Reuse_Save = EXS(iState)%ReUse_SGTMP
    EXS(iState)%ReUse_SGTMP = .false.
    nTest = 0
    dE = Zero
    dT = Zero
    dTR = Zero
    do ISYJ=1,SGS(iState)%nSym
      nJ = CIS(iState)%NCSF(ISYJ)
      if (nJ == 0) cycle
      call mma_allocate(C,nJ,Label='C')
      do J=1,nJ
        C(J) = sin(real(7*J+3*ISYJ,kind=wp))
      end do
      do q=1,nLev
        do p=1,nLev
          ISYI = Mul(Mul(SGS(iState)%ISm(p),SGS(iState)%ISm(q)),ISYJ)
          nI = CIS(iState)%NCSF(ISYI)
          if ((nI == 0) .or. (nI*nJ > MaxPairs)) cycle
          call mma_allocate(S1,nI,Label='S1')
          call mma_allocate(S2,nI,Label='S2')
          ! E_pq
          S1(:) = Zero
          call SG_Xpq_Psi(iState,p,q,One,ISYJ,C,S1,.false.)
          S2(:) = Zero
          call SG_Epq_Psi(SGS(iState),CIS(iState),EXS(iState),p,q,One,ISYJ,C,S2)
          dE = max(dE,maxval(abs(S1-S2)))
          if (maxval(abs(S1-S2)) > Tol) then
            nFail = nFail+1
            write(u6,'(A,3I5,ES12.4)') ' SG_Xpq_Psi (E_pq) differs from SG_Epq_Psi: ISYJ,p,q =',ISYJ,p,q,maxval(abs(S1-S2))
          end if
          ! T_pq: the split-graph tables against the reference
          S1(:) = Zero
          call SG_Tpq_Psi(iState,p,q,One,ISYJ,C,S1)
          S2(:) = Zero
          call SG_Tpq_Psi_Ref(iState,p,q,One,ISYJ,C,S2)
          dTR = max(dTR,maxval(abs(S1-S2)))
          if (maxval(abs(S1-S2)) > Tol) then
            nFail = nFail+1
            write(u6,'(A,3I5,ES12.4)') ' SG_Tpq_Psi (tables) differs from SG_Tpq_Psi_Ref: ISYJ,p,q =',ISYJ,p,q,maxval(abs(S1-S2))
          end if
          S2(:) = Zero
          do I=1,nI
            call SG_CSF_StepVec(SGS(iState),CIS(iState),ISYI,I,SI)
            do J=1,nJ
              call SG_CSF_StepVec(SGS(iState),CIS(iState),ISYJ,J,SJ)
              S2(I) = S2(I)+ASTA_MatEl(nLev,SI,SJ,p,q,.true.)*C(J)
            end do
          end do
          dT = max(dT,maxval(abs(S1-S2)))
          if (maxval(abs(S1-S2)) > Tol) then
            nFail = nFail+1
            write(u6,'(A,3I5,ES12.4)') ' SG_Tpq_Psi differs from the dense product: ISYJ,p,q =',ISYJ,p,q,maxval(abs(S1-S2))
          end if
          nTest = nTest+1
          call mma_deallocate(S1)
          call mma_deallocate(S2)
        end do
      end do
      call mma_deallocate(C)
    end do
    EXS(iState)%ReUse_SGTMP = Reuse_Save
    call mma_deallocate(SI)
    call mma_deallocate(SJ)
    write(u6,'(A,I0,A)') ' SG_Xpq_Psi tests: ',nTest,' (symmetry, p, q) blocks'
    write(u6,'(A,ES10.3)') '   E_pq vs SG_Epq_Psi,           max. abs. difference ',dE
    write(u6,'(A,ES10.3)') '   T_pq vs dense ASTA_MatEl,     max. abs. difference ',dT
    write(u6,'(A,ES10.3)') '   T_pq split-graph tables vs reference, max. abs. difference ',dTR
  end subroutine Check_Xpq_Psi

!
! Random index in 1..n from a Park-Miller generator (x <- 48271 x mod
! (2**31-1)) with the state iRandState of the host. Self-contained, so the
! sequence does not depend on other users of a random number generator.
!
  function Rand_Index(n) result(iRand)
    integer(kind=iwp), intent(in) :: n
    integer(kind=iwp) :: iRand
    if (.not. SeedReported) then
      if (SeedFromEnv) then
        write(u6,'(/,A,I0,A)') 'SG_Test_EpqErs_Ref: random sampling of (I,J) pairs, seed = ', &
          iRandSeed,' (from MOLCAS_RANDOM_SEED)'
      else
        write(u6,'(/,A,I0,A)') 'SG_Test_EpqErs_Ref: random sampling of (I,J) pairs, seed = ', &
          iRandSeed,' (from the clock)'
      end if
      write(u6,'(A,I0,/)') 'To reproduce this run: >export MOLCAS_RANDOM_SEED=',iRandSeed
      SeedReported = .true.
    end if
    iRandState = mod(48271_int64*iRandState,2147483647_int64)
    iRand = 1+int(mod(iRandState,int(n,kind=int64)),kind=iwp)
  end function Rand_Index

!
! Seed of the generator: MOLCAS_RANDOM_SEED if set, otherwise the clock.
! (getSeed is not used, since inside verify it gives a deterministic seed.)
!
  subroutine Init_Rand_Seed()
    character(len=72) :: Line
    integer(kind=int64) :: iCount
    call get_environment_variable('MOLCAS_RANDOM_SEED',Line)
    SeedFromEnv = Line /= ' '
    if (SeedFromEnv) then
      read(Line,*) iRandSeed
    else
      call system_clock(iCount)
      iRandSeed = 1_int64+mod(abs(iCount),2147483646_int64)
    end if
    ! The Park-Miller state must lie in 1..2**31-2.
    iRandState = 1_int64+mod(abs(iRandSeed),2147483646_int64)
    SeedReported = .false.
  end subroutine Init_Rand_Seed

  !
  ! Quick zero test for <I|E_pq|J> from the occupations OccI, OccJ of I, J.
  ! Returns .true. only if the matrix element is exactly zero:
  !   p == q: E_pp is diagonal in the CSF basis, so zero unless I == J
  !           (the symmetries of I and J are then equal);
  !   p /= q: zero unless n_J(q) >= 1, n_J(p) <= 1 and
  !           n_I(t) = n_J(t) + d(t,p) - d(t,q) for every level t.
  !
  function Epq_Is_Zero(OccI,OccJ,I,J,p,q) result(IsZero)
    integer(kind=iwp), intent(in) :: OccI(:), OccJ(:), I, J, p, q
    logical(kind=iwp) :: IsZero
    integer(kind=iwp) :: t, nI
    IsZero = .true.
    if (p == q) then
      IsZero = I /= J
      return
    end if
    if (OccJ(q) == 0 .or. OccJ(p) == 2) return
    do t=1,size(OccJ)
      nI = OccJ(t)
      if (t == p) nI = nI+1
      if (t == q) nI = nI-1
      if (nI /= OccI(t)) return
    end do
    IsZero = .false.
  end function Epq_Is_Zero

end subroutine SG_Test_EpqErs_Ref

end module SGUGA_Checks
