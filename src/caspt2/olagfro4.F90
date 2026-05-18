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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************
      SUBROUTINE OLagFro4(NBSQT,iSym0,iSymI,iSymJ,iSymK,iSymL0,         &
     &                    DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,WRK1)

      USE CHOVEC_IO, only: NVLOC_CHOBATCH
      use Cholesky, only: InfVec, nDimRS
      use ChoCASPT2, only: NUMCHO_PT2, NCHSPC, MXNVC
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp, u6
      use Constants, only: Zero, One, Half
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
      use caspt2_module, only: NSYM, NBAS, NBTCHES

      implicit none

#include "warnings.h"
#include "intent.fh"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#endif

      integer(kind=iwp), intent(in) :: NBSQT, iSym0, iSymI, iSymJ,      &
     &                                 iSymK, iSymL0
      real(kind=wp), intent(inout) :: DPT2AO(NBSQT), DPT2CAO(NBSQT)
      real(kind=wp), intent(_OUT_) :: FPT2AO(NBSQT), FPT2CAO(NBSQT),    &
     &                                WRK1(NBSQT)

      real(kind=wp), allocatable :: CHSPC(:), WRK2(:)
      integer(kind=iwp) :: ISTLT(8), ISTSQ(8), iSkip(8), ipWRK(8),      &
     &  nnbstr(8,3), iSym, maxvec, n2, jSym, nB, nB2, nB3, nBasI, nBasJ,&
     &  iSymIJ, nBasIJ, nBasK, iSMax, iSymL, nBasL, nBasKL, IBATCH_TOT, &
     &  JRED1, JRED2, JSTART, NVECS_RED, ILOC, IRC, JRED, NBATCH, JV1,  &
     &  IBATCH, JNUM, JV2, JREDC, NUMV, MUSED, ipVecL, iVec, jVref,     &
     &  lscr, JREDL, JVEC1, iSwap, i, j
      real(kind=wp) :: tmp

      ! INFVEC(I,J,K)=IWORK(ip_INFVEC-1+MAXVEC*N2*(K-1)+MAXVEC*(J-1)+I)

      !! It shoudl be zero, but just in case
      FPT2AO(1:NBSQT) = Zero
      FPT2CAO(1:NBSQT) = Zero

#ifdef _MOLCAS_MPP_
      If (Is_Real_Par()) Then
        !! To broadcast DPT2AO and DPT2CAO
        If (.not.King()) Then
          DPT2AO(1:NBSQT) = Zero
          DPT2CAO(1:NBSQT) = Zero
        End If
        CALL GADGOP (DPT2AO,NBSQT,'+')
        CALL GADGOP (DPT2CAO,NBSQT,'+')
      End If
#endif

      iSym = iSym0
      call getritrfinfo(nnbstr,maxvec,n2)

      ISTSQ(1)=0
      ISTLT(1)=0
      Do jSym = 2, nSym
        nB  = nBas(jSym-1)
        nB2 = nB*nB
        nB3 = (nB2+nB)/2
        ISTSQ(jSym) = ISTSQ(jSym-1) + nB2
        ISTLT(jSym) = ISTLT(jSym-1) + nB3
      End Do
      Do jSym = 1, nSym
        iSkip(jSym) = 1
        ipWRK(jSym) = 1
      End Do

      nBasI  = nBas(iSymI)
      nBasJ  = nBas(iSymJ)
      iSymIJ = 1+iEor(iSymI-1,iSymJ-1)
      nBasIJ = nBasI*nBasJ
      If (iSymI == iSymJ) nBasIJ = (nBasI*(nBasI+1))/2
      If (nBasIJ == 0) Return

      nBasK  = nBas(iSymK)
      iSMax  = iSymK
      If (iSymK == iSymI) iSMax = iSymJ
      iSymL  = 1+iEor(iSymIJ-1,iSymK-1)
      IF (iSymL > iSMax) Return !! should not
      nBasL  = nBas(iSymL0)
      nBasKL = nBasK*nBasL
      IF (iSymK == iSymL0) nBasKL = (nBasK*(nBasK+1))/2
      If (nBasKL == 0) Return

      call mma_allocate(CHSPC,NCHSPC,Label='CHSPC')
      call mma_allocate(WRK2,NBSQT,Label='WRK2')

      IBATCH_TOT=NBTCHES(iSym)

      IF(NUMCHO_PT2(iSym) == 0) Return

      ! ipnt=ip_InfVec+MaxVec_PT2*(1+InfVec_N2_PT2*(iSym-1))
      ! JRED1=iWork(ipnt)
      ! JRED2=iWork(ipnt-1+NumCho_PT2(iSym))
      JRED1=InfVec(1,2,iSym)
      JRED2=InfVec(NumCho_PT2(iSym),2,iSym)
!     write(u6,*) 'jred1,jred2 = ', jred1,jred2

! Loop over JRED
      DO JRED=JRED1,JRED2

        CALL Cho_X_nVecRS(JRED,iSym,JSTART,NVECS_RED)
        IF(NVECS_RED == 0) Cycle

        ILOC=3
        CALL CHO_X_SETRED(IRC,ILOC,JRED)
! For a reduced set, the structure is known, including
! the mapping between reduced index and basis set pairs.
! The reduced set is divided into suitable batches.
! First vector is JSTART. Nr of vectors in r.s. is NVECS_RED.
        ! JEND=JSTART+NVECS_RED-1

! Determine batch length for this reduced set.
! Make sure to use the same formula as in the creation of disk
! address tables, etc, above:
        NBATCH=1+(NVECS_RED-1)/MXNVC

! Loop over IBATCH
        JV1=JSTART
        DO IBATCH=1,NBATCH
          IBATCH_TOT=IBATCH_TOT+1

          JNUM=NVLOC_CHOBATCH(IBATCH_TOT)
          JV2=JV1+JNUM-1

          JREDC=JRED
! Read a batch of reduced vectors
          CALL CHO_VECRD(CHSPC,NCHSPC,JV1,JV2,iSym,                     &
     &                            NUMV,JREDC,MUSED)

          IF(NUMV /= JNUM) THEN
            write(u6,*)' Rats! CHO_VECRD was called, assuming it to'
            write(u6,*)' read JNUM vectors. Instead it returned NUMV'
            write(u6,*)' vectors: JNUM, NUMV=',JNUM,NUMV
            write(u6,*)' Back to the drawing board?'
            CALL QUIT(_RC_INTERNAL_ERROR_)
          END IF
          IF(JREDC /= JRED) THEN
            write(u6,*)' Rats! It was assumed that the Cholesky vectors'
            write(u6,*)' in HALFTRNSF all belonged to a given reduced'
            write(u6,*)' set, but they don''t!'
            write(u6,*)' JRED, JREDC:',JRED,JREDC
            write(u6,*)' Back to the drawing board?'
            write(u6,*)' Let the program continue and see what happens.'
          END IF

          ipVecL = 1
          Do iVec = 1, NUMV
            !! (strange) reduced form -> squared AO vector (mu nu|iVec)
            jVref = 1 !! only for iSwap=1
!           lscr  = nBasI*(nBasI+1)/2
            ! If (l_NDIMRS < 1) Then
            If (size(nDimRS) < 1) Then
              lscr  = NNBSTR(iSym,3)
            Else
              JREDL = INFVEC(iVec,2,iSym)
              ! lscr  = iWork(ip_nDimRS+iSym-1+nSym*(JREDL-1)) !! JRED?
              lscr  = nDimRS(iSym,JREDL)
            End If
            JVEC1 = 1
            iSwap = 2
            WRK2(:) = Zero
            Call Cho_ReOrdr(irc,CHSPC(ipVecL),lscr,jVref,               &
     &                      JVEC1,1,1,iSym,JREDC,iSwap,ipWRK,WRK2,      &
     &                      iSkip)
            ipVecL = ipVecL + lscr
!
!           ----- Fock-like transformations -----
!
            Call FDGTRF_RI(WRK2,DPT2AO ,FPT2AO )
            Call FDGTRF_RI(WRK2,DPT2CAO,FPT2CAO)
          End Do
          JV1=JV1+JNUM
        End Do
      End Do

      call mma_deallocate(CHSPC)
      call mma_deallocate(WRK2)

      !! Have to symmetrize Fock-transformed matrices
      Do i = 1, nBasI
        Do j = 1, i-1
          tmp = (FPT2AO(i+nBasI*(j-1))+FPT2AO(j+nBasI*(i-1)))*Half
          FPT2AO(i+nBasI*(j-1)) = Tmp
          FPT2AO(j+nBasI*(i-1)) = Tmp
          tmp = (FPT2CAO(i+nBasI*(j-1))+FPT2CAO(j+nBasI*(i-1)))*Half
          FPT2CAO(i+nBasI*(j-1)) = Tmp
          FPT2CAO(j+nBasI*(i-1)) = Tmp
        End Do
      End Do

#ifdef _MOLCAS_MPP_
      If (Is_Real_Par()) Then
        CALL GADGOP (FPT2AO,NBSQT,'+')
        CALL GADGOP (FPT2CAO,NBSQT,'+')
      End If
#endif

      Return

      Contains

      Subroutine FDGTRF_RI(ChoVec,DD,FF)

      implicit none

      real(kind=wp), intent(in) :: ChoVec(nBasI**2), DD(nBasI**2)
      real(kind=wp), intent(inout) :: FF(nBasI**2)

      real(kind=wp) :: Scal
      real(kind=wp), external :: ddot_

      !! Coulomb
      Scal = DDot_(nBasI**2,DD,1,ChoVec,1)
      FF(1:nBasI**2) = FF(1:nBasI**2) + Scal*ChoVec(1:nBasI**2)

      !! Exchange
      Call DGEMM_('T','N',nBasI,nBasI,nBasI,                            &
     &            One,ChoVec,nBasI,DD,nBasI,                            &
     &            Zero,WRK1,nBasI)
      Call DGEMM_('T','T',nBasI,nBasI,nBasI,                            &
     &           -Half,ChoVec,nBasI,WRK1,nBasI,                         &
     &            One,FF,nBasI)

      End Subroutine FDGTRF_RI

      End Subroutine OLagFro4
