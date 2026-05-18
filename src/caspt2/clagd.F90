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

      SUBROUTINE CLagD(NASHT,NG3,NSTATE,G1,G2,G3,DG1,DG2,DG3,DF1,DF2,   &
     &                 DF3,DEASUM,DEPSA,VECROT)

      use caspt2_global, only: imag_shift, iVecL,                       &
     &                         sigma_p_epsilon, LUSBT,                  &
     &                         LUSOLV, ipea_shift
      use Constants, only: Zero, Half, Two, Four
      use EQSOLV, only: IDSMAT, IDBMAT, IRHS, IVECX, IVECR, IVECW
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp, byte
      use fake_GA, only: GA_Arrays
      use caspt2_module, only: IFMSCOUP, NSYM, NASH, NAES, NASUP, NISUP,&
     &                         NINDEP, EPSA, EASUM, NTUES, NTGEUES,     &
     &                         NTGTUES
#ifdef _MOLCAS_MPP_
      use caspt2_global, only: do_lindep, idSDMat, LUSTD, real_shift
      use definitions, only: u6
      USE Para_Info, ONLY: Is_Real_Par, King
#endif

      implicit none

#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

      integer(kind=iwp), intent(in) :: NASHT, NG3, NSTATE
      real(kind=wp), intent(inout) :: G1(NASHT,NASHT),                  &
     & G2(NASHT,NASHT,NASHT,NASHT),G3(NG3),DG1(NASHT,NASHT),            &
     & DG2(NASHT,NASHT,NASHT,NASHT),DG3(NG3),DF1(NASHT,NASHT),          &
     & DF2(NASHT,NASHT,NASHT,NASHT),DF3(NG3),DEASUM,DEPSA(NASHT,NASHT)
      real(kind=wp), intent(in) :: VECROT(NSTATE)

      real(kind=wp), allocatable :: LBD(:),LID(:) !!,VEC1(:),VEC2(:)
      real(kind=wp), allocatable :: SMat(:),BDER(:),SDER(:)
#ifdef _MOLCAS_MPP_
      integer(kind=byte), ALLOCATABLE :: idxG3(:,:)
      real(kind=wp), allocatable :: VEC1(:),VEC2(:),VEC3(:),VEC4(:),    &
     &                              VEC5(:)
#endif
      integer(kind=iwp) :: iCase, iSym, NIN, NIS, NVEC, NAS, ID,        &
     &                     iLUID
#ifdef _MOLCAS_MPP_
      integer(kind=iwp) :: MYRANK, lg_S, mS, LDV
      integer(kind=iwp) :: ILO, IHI, JLO, JHI
#endif
      integer(kind=iwp) :: NS, idS, idum, iTabs, iUabs, iVabs, iXabs,   &
     &  iYabs, iTU2, iTUabs, iTgeUabs, iTgtUabs, iXY2,                  &
     &  iXYabs, iXgeYabs, iXgtYabs, NSEQ
      integer(kind=iwp) :: lg_V1, lg_V2, lg_V3, lg_V4, lg_V5

      real(kind=wp) :: ScalB1, ScalB2, ScalS1, ScalS2, ET, EU, EX, EY,  &
     & ATUXY, BDERval, bsBDER, SDERval, ATYU, ATYX, ATUX, ATUY

      Do iCase = 1, 11
!       cycle
!       If (icase /= 10.and.icase /= 11) cycle ! G
!       If (icase /= 10)                 cycle ! GP
!       If (icase /=  6.and.icase /=  7) cycle ! E
!       If (icase /=  8.and.icase /=  9) cycle ! F
!       If (icase /=  8)                 cycle ! FP
!       If (icase /=  2.and.icase /=  3) cycle ! B
!       If (icase /=  5)                 cycle ! D
!       If (icase /=  4)                 cycle ! C
!       If (icase /=  1)                 cycle ! A
        Do iSym = 1, nSym
          nIN  = nINDEP(iSym,iCase)
          If (nIN == 0) Cycle
          nIS  = nISUP(iSym,iCase)
          NVEC = nIN*nIS
          nAS  = nASUP(iSym,iCase)
          If (nVec == 0) Cycle

#ifdef _MOLCAS_MPP_
          if (is_real_par()) then
            if (iCase /= 1 .and. iCase /= 4) then
              call mma_allocate(BDER,NAS**2,Label='BDER')
              call mma_allocate(SDER,NAS**2,Label='SDER')
              BDER(:) = Zero
              SDER(:) = Zero
            end if
          else
#endif
            call mma_allocate(BDER,NAS**2,Label='BDER')
            call mma_allocate(SDER,NAS**2,Label='SDER')
            BDER(:) = Zero
            SDER(:) = Zero
#ifdef _MOLCAS_MPP_
          end if
#endif

#if defined(_MOLCAS_MPP_) && defined(_GA_)
          if (is_real_par() .and. (icase == 1 .or. icase == 4)) then
            call CLagDX_MPP()
            cycle
          else
#endif

!         write(u6,*) 'for icase = ', icase
!         write(u6,*) '# of independent vecs:', nin
!         write(u6,*) '# of non-active pairs:', nis
!         write(u6,*) '# of     active pairs:', nas
!         write(u6,*) 'dimension for Vec = ', nin*nis
          !! lg_V1 = T (solution; not quasi-variational)
          Call RHS_ALLO(nIN,nIS,lg_V1)
          Call RHS_READ_SR(lg_V1,iCase,iSym,iVecX)
          !! lg_V2 = lambda (shift correction)
          Call RHS_ALLO(nIN,nIS,lg_V2)
          CALL RHS_READ_SR(lg_V2,iCase,iSym,iVecR)
          if (sigma_p_epsilon /= Zero) then
            call mma_allocate(LBD,nAS,Label='LBD')
            call mma_allocate(LID,nIS,Label='LID')
            iD = iDBMat(iSym,iCase)
            Call dDaFile(LUSBT,2,LBD,nAS,iD)
            Call dDaFile(LUSBT,2,LID,nIS,iD)
            Call CASPT2_ResD(3,nIN,nIS,lg_V2,lg_V1,LBD,LID)
            call mma_deallocate(LBD)
            call mma_deallocate(LID)
          end if

          !! lg_V3 = RHS (in IC basis)
          Call RHS_ALLO(nIN,nIS,lg_V3)
          Call RHS_READ_SR(lg_V3,iCase,iSym,iRHS)
          !! lg_V4 = RHS (in MO basis)
          Call RHS_ALLO(nAS,nIS,lg_V4)
          Call RHS_READ_C (lg_V4,iCase,iSym,iVecW)
          !! lg_V5 = RHS2 (in IC basis)
          If (IFMSCOUP) Then
            Call RHS_ALLO(nIN,nIS,lg_V5)
            Call RHS_READ_SR(lg_V5,iCase,iSym,iVecL) ! 7
          Else
            lg_V5 = lg_V3
          End If

#ifdef _MOLCAS_MPP_
          IF (Is_Real_Par()) THEN
            IF (KING()) THEN
              ! copy global array to local buffer
              call mma_allocate(VEC1,NVEC,Label='VEC1')
              CALL GA_GET(lg_V1,1,NIN,1,NIS,VEC1,NIN)
              call mma_allocate(VEC2,NVEC,Label='VEC2')
              CALL GA_GET(lg_V2,1,NIN,1,NIS,VEC2,NIN)
              call mma_allocate(VEC3,NIN*NIS,Label='VEC3')
              CALL GA_GET(lg_V3,1,NIN,1,NIS,VEC3,NIN)
              call mma_allocate(VEC4,NAS*NIS,Label='VEC4')
              CALL GA_GET(lg_V4,1,NAS,1,NIS,VEC4,NAS)
              !! is it possible to avoid allocating twice?
!             IF (IFMSCOUP) THEN
!               call mma_allocate(VEC5,NIN*NIS,Label='VEC5')
!               CALL GA_GET(lg_V5,1,NIN,1,NIS,VEC5,NIN)
!             ELSE
!               LVEC5 = LVEC3
!             END IF
              call mma_allocate(VEC5,NIN*NIS,Label='VEC5')
              CALL GA_GET(lg_V5,1,NIN,1,NIS,VEC5,NIN)

              CALL CLagDX(0,ISYM,ICASE,VEC1,VEC2,                       &
     &                    VEC3,VEC4,                                    &
     &                    nIN,nIS,nAS,nState,                           &
     &                    VECROT,VEC5,lg_V2,BDER,SDER)

              ! free local buffer
              call mma_deallocate(VEC1)
              call mma_deallocate(VEC2)
              call mma_deallocate(VEC3)
              call mma_deallocate(VEC4)
!             IF (IFMSCOUP) THEN
              call mma_deallocate(VEC5)
!             END IF
            END IF
            CALL GASYNC()
          ELSE
#endif
            CALL CLagDX(0,iSym,iCase,GA_Arrays(lg_V1)%A,                &
     &                               GA_Arrays(lg_V2)%A,                &
     &                               GA_Arrays(lg_V3)%A,                &
     &                               GA_Arrays(lg_V4)%A,                &
     &                  nIN,nIS,nAS,nState,                             &
     &                  VECROT,GA_Arrays(lg_V5)%A,lg_V2,BDER,SDER)
#ifdef _MOLCAS_MPP_
          END IF
#endif

          If (imag_shift  /= Zero .or. sigma_p_epsilon /= Zero) Then
            nAS = nASUP(iSym,iCase)
            call mma_allocate(LBD,nAS,Label='LBD')
            call mma_allocate(LID,nIS,Label='LID')
            iD = iDBMat(iSym,iCase)
            Call dDaFile(LUSBT,2,LBD,nAS,iD)
            Call dDaFile(LUSBT,2,LID,nIS,iD)

            CALL RHS_READ_SR(lg_V1,ICASE,ISYM,iVecX)
            CALL RHS_READ_SR(lg_V2,ICASE,ISYM,iVecR)
            Call CASPT2_ResD(2,nIN,nIS,lg_V1,lg_V2,LBD,LID)
            call mma_deallocate(LBD)
            call mma_deallocate(LID)

#ifdef _MOLCAS_MPP_
            IF (Is_Real_Par()) THEN
              IF (KING()) THEN
                ! copy global array to local buffer
                call mma_allocate(VEC1,NVEC,Label='VEC1')
                CALL GA_GET(lg_V1,1,NIN,1,NIS,VEC1,NIN)
                call mma_allocate(VEC2,NVEC,Label='VEC2')
                CALL GA_GET(lg_V2,1,NIN,1,NIS,VEC2,NIN)

                !! lvec3 to lvec5 are not used
                !! only lvec1 and lvec2?
                call mma_allocate(VEC3,1,Label='VEC3')
                call mma_allocate(VEC4,1,Label='VEC4')
                call mma_allocate(VEC5,1,Label='VEC5')
                CALL CLagDX(1,ISYM,ICASE,VEC1,VEC2,                     &
     &                      VEC3,VEC4,                                  &
     &                      nIN,nIS,nAS,nState,                         &
     &                      VECROT,VEC5,lg_V2,BDER,SDER)

                ! free local buffer
                call mma_deallocate(VEC1)
                call mma_deallocate(VEC2)
                call mma_deallocate(VEC3)
                call mma_deallocate(VEC4)
                call mma_deallocate(VEC5)
              END IF
              CALL GASYNC()
            ELSE
#endif
              CALL CLagDX(1,iSym,iCase,GA_Arrays(lg_V1)%A,              &
     &                                 GA_Arrays(lg_V2)%A,              &
     &                                 GA_Arrays(lg_V3)%A,              &
     &                                 GA_Arrays(lg_V4)%A,              &
     &                    nIN,nIS,nAS,nState,                           &
     &                    VECROT,GA_Arrays(lg_V5)%A,lg_V2,BDER,SDER)
#ifdef _MOLCAS_MPP_
            end if
#endif
          End If

          !! for non-separable density/derivative
          CALL RHS_READ_SR(lg_V1,ICASE,ISYM,iVecX)
          CALL RHS_READ_SR(lg_V2,ICASE,ISYM,iVecR)
#if defined(_MOLCAS_MPP_) && defined(_GA_)
          end if
#endif

          if (iCase ==  1) call CLagDXA(NAS,BDER,SDER)
          if (iCase ==  2) call CLagDXB(NAS,BDER,SDER)
          if (iCase ==  3) call CLagDXB(NAS,BDER,SDER)
          if (iCase ==  4) call CLagDXC(NAS,BDER,SDER)
          if (iCase ==  5) call CLagDXD(NAS,BDER,SDER)
          if (iCase ==  6) call CLagDXE(NAS,BDER,SDER)
          if (iCase ==  7) call CLagDXE(NAS,BDER,SDER)
          if (iCase ==  8) call CLagDXF(NAS,BDER,SDER)
          if (iCase ==  9) call CLagDXF(NAS,BDER,SDER)
          if (iCase == 10) call CLagDXG(NAS,BDER,SDER)
          if (iCase == 11) call CLagDXG(NAS,BDER,SDER)

          CALL RHS_FREE(lg_V1)
          CALL RHS_FREE(lg_V2)
          CALL RHS_FREE(lg_V3)
          CALL RHS_FREE(lg_V4)
          If (IFMSCOUP) CALL RHS_FREE(lg_V5)

#ifdef _MOLCAS_MPP_
          if (is_real_par()) then
            if (iCase /= 1 .and. iCase /= 4) then
              call mma_deallocate(BDER)
              call mma_deallocate(SDER)
            end if
          else
#endif
            call mma_deallocate(BDER)
            call mma_deallocate(SDER)
#ifdef _MOLCAS_MPP_
          end if
#endif
        End Do
      End Do

#ifdef _MOLCAS_MPP_
      if (is_real_par()) then
        iCase = 4
        MYRANK=GA_NODEID()
        Do iSym = 1, nSym
          nAS  = nASUP(iSym,iCase)
          CALL PSBMAT_GETMEM('S',lg_S,NAS)
          CALL PSBMAT_READ('S',4,iSym,lg_S,NAS)
          CALL GA_Distribution(lg_S   ,myRank,ILO,IHI,JLO,JHI)
          Call GA_Access(lg_S   ,ILO,IHI,JLO,JHI,mS   ,LDV)
          CALL mma_allocate(idxG3,6,NG3,label='idxG3')
          iLUID=0
          CALL I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)
          !! DF3 is done after icase=4
          !! construct G3 matrix in lg_S
          CALL MKSC_G3_MPP(ISYM,DBL_MB(mS),ILO,IHI,NAS,LDV,             &
     &                     NG3,G3,IDXG3)
          Call GA_Release_Update(lg_S,ILO,IHI,JLO,JHI)
          call DF3_DEPSA_MPP(NG3,NASHT,DF3,DEPSA,lg_S,idxG3)

          Call GA_Release(lg_S   ,ILO,IHI,JLO,JHI)
          CALL mma_deallocate(idxG3)
          CALL PSBMAT_FREEMEM(lg_S)
        End Do
      end if
#endif

      Return

      contains
!
!-----------------------------------------------------------------------
!
      Subroutine CLagDXA(NAS,BDER,SDER)

      implicit none

      integer(kind=iwp), intent(in) :: NAS
      real(kind=wp), intent(in) :: BDER(NAS,NAS)
      real(kind=wp), intent(inout) :: SDER(NAS,NAS)

      integer(kind=byte), allocatable :: idxG3(:,:)

      NS = NAS*(NAS+1)/2
      call mma_allocate(SMat,NS,Label='SMat')
      idS = idSMAT(iSym,1)
      CALL DDAFILE(LUSBT,2,SMat,NS,idS)

      idum=0
      Call CLagDXA_DP (iSym,nAS,nAshT,BDER,SDER,                        &
     &                 DG1,DG2,DF1,DF2,DEPSA,DEASUM,                    &
     &                 1,nAS,1,nAS,0,G1,G2,SMat,SMat,idum)

      !! G3 and F3 relevant
      CALL mma_allocate(idxG3,6,NG3,label='idxG3')
      iLUID=0
      CALL I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)
!     idS = idSMAT(iSym,4)
!     CALL DDAFILE(LUSBT,2,SMat,NS,idS)
      CALL MKSC_G3(iSym,SMat,NS,nG3,G3,idxG3)
      call CLagDXA_FG3(iSym,nAS,nAshT,NG3,NS,BDER,SDER,                 &
     &                 DF1,DF2,DF3,DG1,DG2,DG3,DEPSA,G2,                &
     &                 SMat,idxG3)
      call mma_deallocate(idxG3)

      call mma_deallocate(SMat)

      return

      End Subroutine CLagDXA
!
!-----------------------------------------------------------------------
!
      Subroutine CLagDXB(NAS,BDER,SDER)

      USE SUPERINDEX, only: MTGEU, MTGTU

      implicit none

      integer(kind=iwp), intent(in) :: NAS
      real(kind=wp), intent(in) :: BDER(NAS,NAS)
      real(kind=wp), intent(inout) :: SDER(NAS,NAS)

      real(kind=wp), allocatable :: WrkBbf(:,:,:,:),WrkSbf(:,:,:,:)
      integer(kind=iwp) :: iTU, iXY, iT, iU, iX, iY, iV

      call mma_allocate(WrkBbf,nAshT,nAshT,nAshT,nAshT,Label='WrkBbf')
      call mma_allocate(WrkSbf,nAshT,nAshT,nAshT,nAshT,Label='WrkSbf')
      WrkBbf(:,:,:,:) = Zero
      WrkSbf(:,:,:,:) = Zero

      If (ipea_shift /= Zero) Then
        NS = NAS*(NAS+1)/2
        call mma_allocate(SMat,NS,Label='SMat')
        idS = idSMAT(iSym,iCase)
        CALL DDAFILE(LUSBT,2,SMat,NS,idS)
      End If
      ScalB1 = Zero
      ScalB2 = Zero
      ScalS1 = Zero
      ScalS2 = Zero
      iTabs  = 0
      iUabs  = 0
      iXabs  = 0
      iYabs  = 0
      iTgeUabs = 0
      iTgtUabs = 0
      iXgeYabs = 0
      iXgtYabs = 0
      Do iTU = 1, nAS
        If (iCase == 2) Then
          iTgeUabs = iTU + nTgeUes(iSym)
          iTabs    = mTgeU(1,iTgeUabs)
          iUabs    = mTgeU(2,iTgeUabs)
        Else If (iCase == 3) Then
          iTgtUabs = iTU + nTgtUes(iSym)
          iTabs    = mTgtU(1,iTgtUabs)
          iUabs    = mTgtU(2,iTgtUabs)
        End If
        ET = EPSA(iTabs)
        EU = EPSA(iUabs)
        DO iXY = 1, nAS
          If (iCase == 2) Then
            iXgeYabs = iXY + nTgeUes(iSym)
            iXabs    = mTgeU(1,iXgeYabs)
            iYabs    = mTgeU(2,iXgeYabs)
          Else If (iCase == 3) Then
            iXgtYabs = iXY + nTgtUes(iSym)
            iXabs    = mTgtU(1,iXgtYabs)
            iYabs    = mTgtU(2,iXgtYabs)
          End If
          EX = EPSA(iXabs)
          EY = EPSA(iYabs)
          ATUXY = EASUM-ET-EU-EX-EY
!         iBadr = iTU + nAS*(iXY-1)
          BDERval = BDER(ITU,IXY)

          !! For IPEA shift
          If (iTU ==iXY .and. ipea_shift /= Zero) Then
!           idT=(iTabs*(iTabs+1))/2
            ! idU=(iUabs*(iUabs+1))/2
            NSEQ = iTU*(iTU+1)/2
            bsBDER = ipea_shift*Half*BDERval
!         !! ipea_shift*Half*(DREF(IDT)+DREF(IDU))*SDP(ITGEU)
            DG1(iTabs,iTabs) = DG1(iTabs,iTabs) + bsBDER*SMat(NSEQ)
            DG1(iUabs,iUabs) = DG1(iUabs,iUabs) + bsBDER*SMat(NSEQ)
            SDER(iTU,iXY) = SDER(iTU,iXY)                               &
     &        + (G1(iTabs,iTabs)+G1(iUabs,iUabs))*bsBDER
          End If
          SDERval = SDER(ITU,IXY)
          If (iTabs == iUabs) Then
            BDERval = BDERval*Two
            SDERval = SDERval*Two
          End If

          If (iCase == 2) Then
            ScalB1 = BDERval
            ScalB2 = BDERval
            ScalS1 = SDERval
            ScalS2 = SDERval
          Else If (iCase == 3) Then
            ScalB1 = BDERval
            ScalB2 =-BDERval
            ScalS1 = SDERval
            ScalS2 =-SDERval
          End If

          WRKBBF(iTabs,iUabs,iXabs,iYabs)                               &
     &      = WRKBBF(iTabs,iUabs,iXabs,iYabs) + ScalB1
          WRKBBF(iTabs,iUabs,iYabs,iXabs)                               &
     &      = WRKBBF(iTabs,iUabs,iYabs,iXabs) + ScalB2
          WRKSBF(iTabs,iUabs,iXabs,iYabs)                               &
     &      = WRKSBF(iTabs,iUabs,iXabs,iYabs) + ScalS1
          WRKSBF(iTabs,iUabs,iYabs,iXabs)                               &
     &      = WRKSBF(iTabs,iUabs,iYabs,iXabs) + ScalS2
          If (iTabs /= iUabs) Then
          WRKBBF(iUabs,iTabs,iXabs,iYabs)                               &
     &      = WRKBBF(iUabs,iTabs,iXabs,iYabs) + ScalB2
          WRKBBF(iUabs,iTabs,iYabs,iXabs)                               &
     &      = WRKBBF(iUabs,iTabs,iYabs,iXabs) + ScalB1
          WRKSBF(iUabs,iTabs,iXabs,iYabs)                               &
     &      = WRKSBF(iUabs,iTabs,iXabs,iYabs) + ScalS2
          WRKSBF(iUabs,iTabs,iYabs,iXabs)                               &
     &      = WRKSBF(iUabs,iTabs,iYabs,iXabs) + ScalS1
          End If
        End Do
      End Do

      !! it,iu,... are actually itabs,iuabs,...
      Do iT = 1, nAshT
        ET = EPSA(iT)
        DO iU = 1, nAshT
          EU = EPSA(iU)
          Do iX = 1, nAshT
            EX = EPSA(iX)
            Do iY = 1, nAshT
              EY = EPSA(iY)
              BDERval = WRKBBF(iT,iU,iX,iY)
              SDERval = WRKSBF(iT,iU,iX,iY)

              !! term 1 (w/o delta)
              ATUXY = EASUM-ET-EU-EX-EY
              !! G1 and F1 derivative
              DF2(iX,iT,iY,iU) = DF2(iX,iT,iY,iU) + BDERval
              DG2(iX,iT,iY,iU) = DG2(iX,iT,iY,iU)                       &
     &          - ATUXY*BDERval + SDERval
              !! EASUM derivative
              DEASUM = DEASUM - BDERval*G2(iX,iT,iY,iU)
              !! EPSA derivative
              Do iV = 1, NASHT
                DEPSA(iT,iV) = DEPSA(iT,iV) + BDERval*G2(iX,iV,iY,iU)
                DEPSA(iU,iV) = DEPSA(iU,iV) + BDERval*G2(iX,iT,iY,iV)
                DEPSA(iX,iV) = DEPSA(iX,iV) + BDERval*G2(iV,iT,iY,iU)
                DEPSA(iY,iV) = DEPSA(iY,iV) + BDERval*G2(iX,iT,iV,iU)
              End Do

              BDERval = BDERval*Two
              SDERval = SDERval*Two

              !! term 2 (dxt)
              If (iX == iT) Then
                ATYU = EASUM-ET-EY-EU
                !! G1 and F1 derivative
                DF1(iY,iU) = DF1(iY,iU) - BDERval
                DG1(iY,iU) = DG1(iY,iU) + ATYU*BDERval - SDERval
                !! EASUM derivative
                DEASUM = DEASUM + BDERval*G1(iY,iU)
                !! EPSA derivative
                Do iV = 1, NASHT
                  DEPSA(iY,iV) = DEPSA(iY,iV) - BDERval*G1(iV,iU)
                  DEPSA(iU,iV) = DEPSA(iU,iV) - BDERval*G1(iY,iV)
                End Do
              End If
              !! Additional EPSA derivative
              DEPSA(iX,iT) = DEPSA(iX,iT) - BDERval*G1(iY,iU)
              !! dxt*dyu term
              If (iY == iU) DEPSA(iX,iT) = DEPSA(iX,iT) + Two*BDERval
              If (iX == iT) DEPSA(iY,iU) = DEPSA(iY,iU) + Two*BDERval

              !! term 3 (dyu)
              If (iY == iU) Then
                ATYX = EASUM-ET-EY-EX
                !! G1 and F1 derivative
                DF1(iX,iT) = DF1(iX,iT) - BDERval
                DG1(iX,iT) = DG1(iX,iT) + ATYX*BDERval - SDERval
                !! EASUM derivative
                DEASUM = DEASUM + BDERval*G1(iX,iT)
                !! EPSA derivative
                Do iV = 1, NASHT
                  DEPSA(iX,iV) = DEPSA(iX,iV) - BDERval*G1(iV,iT)
                  DEPSA(iT,iV) = DEPSA(iT,iV) - BDERval*G1(iX,iV)
                End Do
              End If
              !! Additional EPSA derivative
              DEPSA(iY,iU) = DEPSA(iY,iU) - BDERval*G1(iX,iT)

              BDERval = BDERval*Half
              SDERval = SDERval*Half

              !! term 4 (dyt)
              If (iY == iT) Then
                ATUX = EASUM-ET-EU-EX
                !! G1 and F1 derivative
                DF1(iX,iU) = DF1(iX,iU) + BDERval
                DG1(iX,iU) = DG1(iX,iU) - ATUX*BDERval + SDERval
                !! EASUM derivative
                DEASUM = DEASUM - BDERval*G1(iX,iU)
                !! EPSA derivative
                Do iV = 1, NASHT
                  DEPSA(iX,iV) = DEPSA(iX,iV) + BDERval*G1(iV,iU)
                  DEPSA(iU,iV) = DEPSA(iU,iV) + BDERval*G1(iX,iV)
                End Do
              End If
              !! Additional EPSA derivative
              DEPSA(iY,iT) = DEPSA(iY,iT) + BDERval*G1(iX,iU)
              !! dxu*dyt term
              If (iY == iT) DEPSA(iX,iU) = DEPSA(iX,iU) - Two*BDERval
              If (iX == iU) DEPSA(iY,iT) = DEPSA(iY,iT) - Two*BDERval

              !! term 5 (dxu)
              If (iX == iU) Then
                ATUY = EASUM-ET-EU-EY
                !! G1 and F1 derivative
                DF1(iY,iT) = DF1(iY,iT) + BDERval
                DG1(iY,iT) = DG1(iY,iT) - ATUY*BDERval + SDERval
                !! EASUM derivative
                DEASUM = DEASUM - BDERval*G1(iY,iT)
                !! EPSA derivative
                Do iV = 1, NASHT
                  DEPSA(iY,iV) = DEPSA(iY,iV) + BDERval*G1(iV,iT)
                  DEPSA(iT,iV) = DEPSA(iT,iV) + BDERval*G1(iY,iV)
                End Do
              End If
              !! Additional EPSA derivative
              DEPSA(iX,iU) = DEPSA(iX,iU) + BDERval*G1(iY,iT)
            End Do
          End Do
        End Do
      End Do
      If (ipea_shift /= Zero) call mma_deallocate(SMat)

      call mma_deallocate(WrkBbf)
      call mma_deallocate(WrkSbf)

      return

      End Subroutine CLagDXB
!
!-----------------------------------------------------------------------
!
      Subroutine CLagDXC(NAS,BDER,SDER)

      implicit none

      integer(kind=iwp), intent(in) :: NAS
      real(kind=wp), intent(in) :: BDER(NAS,NAS)
      real(kind=wp), intent(inout) :: SDER(NAS,NAS)

      integer(kind=byte), allocatable :: idxG3(:,:)

      NS = NAS*(NAS+1)/2
      call mma_allocate(SMat,NS,Label='SMat')
      idS = idSMAT(iSym,4)
      CALL DDAFILE(LUSBT,2,SMat,NS,idS)

      idum = 0
      Call CLagDXC_DP (iSym,nAS,nAshT,BDER,SDER,                        &
     &                 DG1,DG2,DF1,DF2,DEPSA,DEASUM,                    &
     &                 1,nAS,1,nAS,0,G1,G2,SMat,SMat,idum)

      !! G3 and F3 relevant
      CALL mma_allocate(idxG3,6,NG3,label='idxG3')
      iLUID=0
      CALL I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)
!     idS = idSMAT(iSym,4)
!     CALL DDAFILE(LUSBT,2,SMat,NS,idS)
      CALL MKSC_G3(iSym,SMat,NS,nG3,G3,idxG3)
      call CLagDXC_FG3(iSym,nAS,nAshT,NG3,NS,BDER,SDER,                 &
     &                 DF1,DF2,DF3,DG1,DG2,DG3,DEPSA,G2,                &
     &                 SMat,idxG3)
      call mma_deallocate(idxG3)

      call mma_deallocate(SMat)

      return

      End Subroutine CLagDXC
!
!-----------------------------------------------------------------------
!
      Subroutine CLagDXD(NAS,BDER,SDER)

      USE SUPERINDEX, only: MTU

      implicit none

      integer(kind=iwp), intent(in) :: NAS
      real(kind=wp), intent(in) :: BDER(NAS,NAS)
      real(kind=wp), intent(inout) :: SDER(NAS,NAS)

      real(kind=wp) :: BDER1, BDER2, SDER1, SDER2, ETX
      integer(kind=iwp) :: iTU, iXY, iV

      If (ipea_shift /= Zero) Then
        NS = NAS*(NAS+1)/2
        call mma_allocate(SMat,NS,Label='SMat')
        idS = idSMAT(iSym,iCase)
        CALL DDAFILE(LUSBT,2,SMat,NS,idS)
      End If

      Do iTU = 1, nAS/2
        iTU2   = iTU + nAS/2
        iTUabs = iTU + nTUes(iSym)
        iTabs  = mTU(1,iTUabs)
        iUabs  = mTU(2,iTUabs)
        ET     = EPSA(iTabs)
        DO iXY = 1, nAS/2
          iXY2   = iXY + nAS/2
          iXYabs = iXY + nTUes(iSym)
          iXabs  = mTU(1,iXYabs)
          iYabs  = mTU(2,iXYabs)
          EX     = EPSA(iXabs)
          ETX    = ET+EX

          BDER1 = BDER(iTU ,iXY )                                       &
     &          - BDER(iTU ,iXY2)*Half                                  &
     &          - BDER(iTU2,iXY )*Half
          BDER2 = BDER(iTU2,iXY2)

          !! Derivative of B11
          DF2(iUabs,iTabs,iXabs,iYabs)                                  &
     &      = DF2(iUabs,iTabs,iXabs,iYabs) + Two*BDER1
          DG2(iUabs,iTabs,iXabs,iYabs)                                  &
     &      = DG2(iUabs,iTabs,iXabs,iYabs) + Two*(ETX-EASUM)*BDER1
          DEASUM = DEASUM - Two*G2(iUabs,iTabs,iXabs,iYabs)*BDER1
          If (iXabs == iTabs) Then
            DF1(iUabs,iYabs) = DF1(iUabs,iYabs) + Two*BDER1
            DG1(iUabs,iYabs) = DG1(iUabs,iYabs) + Two*(ET-EASUM)*BDER1
            DEASUM = DEASUM - Two*G1(iUabs,iYabs)*BDER1
          End If
          DO iV = 1, nAsh(iSym)
            IVABS=IV+NAES(ISYM)
            DEPSA(iTabs,iVabs) = DEPSA(iTabs,iVabs)                     &
     &        + Two*BDER1*G2(iUabs,iVabs,iXabs,iYabs)
            DEPSA(iXabs,iVabs) = DEPSA(iXabs,iVabs)                     &
     &        + Two*BDER1*G2(iUabs,iTabs,iVabs,iYabs)
          End Do
          DEPSA(iTabs,iXabs) = DEPSA(iTabs,iXabs)                       &
     &      + Two*G1(iUabs,iYabs)*BDER1

          !! Derivative of B22
          DF2(iXabs,iTabs,iUabs,iYabs)                                  &
     &      = DF2(iXabs,iTabs,iUabs,iYabs) - BDER2
          DG2(iXabs,iTabs,iUabs,iYabs)                                  &
     &      = DG2(iXabs,iTabs,iUabs,iYabs) - (ETX-EASUM)*BDER2
          DEASUM = DEASUM + G2(iXabs,iTabs,iUabs,iYabs)*BDER2
          If (iXabs == iTabs) Then
            DF1(iUabs,iYabs) = DF1(iUabs,iYabs) + Two*BDER2
            DG1(iUabs,iYabs) = DG1(iUabs,iYabs) + Two*(EX-EASUM)*BDER2
            DEASUM = DEASUM - Two*G1(iUabs,iYabs)*BDER2
          End If
          DO iV = 1, nAsh(iSym)
            IVABS=IV+NAES(ISYM)
            DEPSA(iTabs,iVabs) = DEPSA(iTabs,iVabs)                     &
     &        - BDER2*G2(iXabs,iVabs,iUabs,iYabs)
            DEPSA(iXabs,iVabs) = DEPSA(iXabs,iVabs)                     &
     &        - BDER2*G2(iVabs,iTabs,iUabs,iYabs)
          End Do
          DEPSA(iXabs,iTabs) = DEPSA(iXabs,iTabs)                       &
     &      + Two*G1(iUabs,iYabs)*BDER2

          If (iTU == iXY .and. ipea_shift /= Zero) Then
!      !! ipea_shift*Half*(Two-DREF(IDU)+DREF(IDT))*SD(ITU)
            bsBDER = ipea_shift*Half*BDER(iTU,iXY)
            NSEQ = iTU*(iTU+1)/2
            DG1(iTabs,iTabs) = DG1(iTabs,iTabs) + bsBDER*SMat(NSEQ)
            DG1(iUabs,iUabs) = DG1(iUabs,iUabs) - bsBDER*SMat(NSEQ)
            SDER(iTU,iXY) = SDER(iTU,iXY)                               &
     &        + bsBDER*(Two+G1(iTabs,iTabs)-G1(iUabs,iUabs))
!    !! ipea_shift*Half*(Two-DREF(IDU)+DREF(IDT))*SD(ITU+NAS)
            bsBDER = ipea_shift*Half*BDER(iTU2,iXY2)
            NSEQ = iTU2*(iTU2+1)/2
            DG1(iTabs,iTabs) = DG1(iTabs,iTabs) + bsBDER*SMat(NSEQ)
            DG1(iUabs,iUabs) = DG1(iUabs,iUabs) - bsBDER*SMat(NSEQ)
            SDER(iTU2,iXY2) = SDER(iTU2,iXY2)                           &
     &        + bsBDER*(Two+G1(iTabs,iTabs)-G1(iUabs,iUabs))
          End If

          SDER1 = SDER(iTU ,iXY )                                       &
     &          - SDER(iTU ,iXY2)*Half                                  &
     &          - SDER(iTU2,iXY )*Half
          SDER2 = SDER(iTU2,iXY2)

          !! Derivative of S11
          DG2(iUabs,iTabs,iXabs,iYabs)                                  &
     &      = DG2(iUabs,iTabs,iXabs,iYabs) + Two*SDER1
          If (iXabs == iTabs) Then
            DG1(iUabs,iYabs) = DG1(iUabs,iYabs) + Two*SDER1
          End If
          !! Derivative of S22
          DG2(iXabs,iTabs,iUabs,iYabs)                                  &
     &      = DG2(iXabs,iTabs,iUabs,iYabs) - SDER2
          If (iXabs == iTabs) Then
            DG1(iUabs,iYabs) = DG1(iUabs,iYabs) + Two*SDER2
          End If
        End Do
      End Do
      If (ipea_shift /= Zero) call mma_deallocate(SMat)

      return

      End Subroutine CLagDXD
!
!-----------------------------------------------------------------------
!
      Subroutine CLagDXE(NAS,BDER,SDER)

      implicit none

      integer(kind=iwp), intent(in) :: NAS
      real(kind=wp), intent(in) :: BDER(NAS,NAS)
      real(kind=wp), intent(inout) :: SDER(NAS,NAS)

      real(kind=wp) :: VAL
      integer(kind=iwp) :: IT, IU, iV

      If (ipea_shift /= Zero) Then
        NS = NAS*(NAS+1)/2
        call mma_allocate(SMat,NS,Label='SMat')
        idS = idSMAT(iSym,6)
        CALL DDAFILE(LUSBT,2,SMat,NS,idS)
        !! ipea_shift*Half*DREF(IDT)*SD(IT)
        DO IT = 1, NAS
          ITABS=IT+NAES(ISYM)
          VAL = ipea_shift*Half*BDER(IT,IT)
          SDER(IT,IT) = SDER(IT,IT) + G1(ITABS,ITABS)*VAL
          NSEQ = IT*(IT-1)/2+IT
          DG1(ITABS,ITABS) = DG1(ITABS,ITABS) + SMat(NSEQ)*VAL
        End Do
        call mma_deallocate(SMat)
      End If

      DO IT = 1, NAS
        ITABS=IT+NAES(ISYM)
        ET = EPSA(ITABS)
        DO IU = 1, NAS
          IUABS=IU+NAES(ISYM)
          EU = EPSA(IUABS)
          !! Derivative of the B matrix
          !! B_{tu} = -F1_{tu} + (Esum-e_t-e_u)*G1(tu)
          DG1(ITABS,IUABS) = DG1(ITABS,IUABS)                           &
     &      + (EASUM-ET-EU)*BDER(IT,IU)
          DEASUM = DEASUM + G1(ITABS,IUABS)*BDER(IT,IU)
          DF1(ITABS,IUABS) = DF1(ITABS,IUABS) - BDER(IT,IU)
          DO IV = 1, NASH(ISYM)
            IVABS=IV+NAES(ISYM)
            DEPSA(ITABS,IUABS) = DEPSA(ITABS,IUABS)                     &
     &        - G1(ITABS,IVABS)*BDER(IV,IU)                             &
     &        - G1(IUABS,IVABS)*BDER(IV,IT)
          End Do
          DEPSA(ITABS,IUABS) = DEPSA(ITABS,IUABS) + Two*BDER(IT,IU)
          !! Derivative of the S matrix
          DG1(ITABS,IUABS) = DG1(ITABS,IUABS) - SDER(IT,IU)
        END DO
      END DO

      return

      End Subroutine CLagDXE
!
!-----------------------------------------------------------------------
!
      Subroutine CLagDXF(NAS,BDER,SDER)

      USE SUPERINDEX, only: MTGTU, MTGEU

      implicit none

      integer(kind=iwp), intent(in) :: NAS
      real(kind=wp), intent(in) :: BDER(NAS,NAS)
      real(kind=wp), intent(inout) :: SDER(NAS,NAS)

      integer(kind=iwp) :: iTU, iXY

      If (ipea_shift /= Zero) Then
        NS = NAS*(NAS+1)/2
        call mma_allocate(SMat,NS,Label='SMat')
        idS = idSMAT(iSym,iCase)
        CALL DDAFILE(LUSBT,2,SMat,NS,idS)
      End If
      ScalB1 = Zero
      ScalB2 = Zero
      ScalS1 = Zero
      ScalS2 = Zero
      iXabs  = 0
      iYabs  = 0
      iTabs  = 0
      iUabs  = 0
      iTgeUabs = 0
      iTgtUabs = 0
      iXgeYabs = 0
      iXgtYabs = 0
      Do iTU = 1, nAS
        If (iCase ==  8) Then
          iTgeUabs = iTU + nTgeUes(iSym)
          iTabs    = mTgeU(1,iTgeUabs)
          iUabs    = mTgeU(2,iTgeUabs)
        Else If (iCase ==  9) Then
          iTgtUabs = iTU + nTgtUes(iSym)
          iTabs    = mTgtU(1,iTgtUabs)
          iUabs    = mTgtU(2,iTgtUabs)
        End If
        DO iXY = 1, nAS !! iTU
          If (iCase ==  8) Then
            iXgeYabs = iXY + nTgeUes(iSym)
            iXabs    = mTgeU(1,iXgeYabs)
            iYabs    = mTgeU(2,iXgeYabs)
          Else If (iCase ==  9) Then
            iXgtYabs = iXY + nTgtUes(iSym)
            iXabs    = mTgtU(1,iXgtYabs)
            iYabs    = mTgtU(2,iXgtYabs)
          End If

          BDERval = BDER(ITU,IXY)
          If (iTU == iXY .and. ipea_shift /= Zero) Then
!           idT=(iTabs*(iTabs+1))/2
            ! idU=(iUabs*(iUabs+1))/2
            NSEQ = iTU*(iTU+1)/2
            bsBDER = ipea_shift*Half*BDERval
!     !! ipea_shift*Half*(Four-DREF(IDT)-DREF(IDU))*SDP(ITGEU)
            DG1(iTabs,iTabs) = DG1(iTabs,iTabs) - SMat(NSEQ)*bsBDER
            DG1(iUabs,iUabs) = DG1(iUabs,iUabs) - SMat(NSEQ)*bsBDER
            SDER(ITU,IXY) = SDER(ITU,IXY)                               &
     &        + (Four-G1(iTabs,iTabs)-G1(iUabs,iUabs))*bsBDER
          End If
          SDERval = SDER(ITU,IXY)
          If (iTabs == iUabs) Then
            BDERval = Two*BDERval
            SDERval = Two*SDERval
          End If

          If (iCase ==  8) Then
            ScalB1 = BDERval
            ScalB2 = BDERval
            ScalS1 = SDERval
            ScalS2 = SDERval
          Else If (iCase ==  9) Then
            ScalB1 = BDERval
            ScalB2 =-BDERval
            ScalS1 = SDERval
            ScalS2 =-SDERval
          End If

          !! Derivative of the B matrix
          !! B(tuxy) -> PREF(tx,uy)
          DEASUM = DEASUM - ScalB1*G2(iTabs,iXabs,iUabs,iYabs)          &
     &                    - ScalB2*G2(iTabs,iYabs,iUabs,iXabs)
          If (iTabs /= iUabs)                                           &
     &    DEASUM = DEASUM - ScalB2*G2(iUabs,iXabs,iTabs,iYabs)          &
     &                    - ScalB1*G2(iUabs,iYabs,iTabs,iXabs)

          ! iTX = iTabs+nAshT*(iXabs-1)
          ! iUY = iUabs+nAshT*(iYabs-1)
          ! iTY = iTabs+nAshT*(iYabs-1)
          ! iUX = iUabs+nAshT*(iXabs-1)

          DF2(iTabs,iXabs,iUabs,iYabs)                                  &
     &      = DF2(iTabs,iXabs,iUabs,iYabs) + ScalB1
          DF2(iTabs,iYabs,iUabs,iXabs)                                  &
     &      = DF2(iTabs,iYabs,iUabs,iXabs) + ScalB2
          If (iTabs /= iUabs) Then
            DF2(iUabs,iXabs,iTabs,iYabs)                                &
     &        = DF2(iUabs,iXabs,iTabs,iYabs) + ScalB2
            DF2(iUabs,iYabs,iTabs,iXabs)                                &
     &        = DF2(iUabs,iYabs,iTabs,iXabs) + ScalB1
          End If
          DG2(iTabs,iXabs,iUabs,iYabs)                                  &
     &      = DG2(iTabs,iXabs,iUabs,iYabs) + ScalS1-EASUM*ScalB1
          DG2(iTabs,iYabs,iUabs,iXabs)                                  &
     &      = DG2(iTabs,iYabs,iUabs,iXabs) + ScalS2-EASUM*ScalB2
          If (iTabs /= iUabs) Then
            DG2(iUabs,iXabs,iTabs,iYabs)                                &
     &        = DG2(iUabs,iXabs,iTabs,iYabs) + ScalS2-EASUM*ScalB2
            DG2(iUabs,iYabs,iTabs,iXabs)                                &
     &        = DG2(iUabs,iYabs,iTabs,iXabs) + ScalS1-EASUM*ScalB1
          End If
        End Do
      End Do
      If (ipea_shift /= Zero) call mma_deallocate(SMat)

      return

      End Subroutine CLagDXF
!
!-----------------------------------------------------------------------
!
      Subroutine CLagDXG(NAS,BDER,SDER)

      implicit none

      integer(kind=iwp), intent(in) :: NAS
      real(kind=wp), intent(in) :: BDER(NAS,NAS)
      real(kind=wp), intent(inout) :: SDER(NAS,NAS)

      real(kind=wp) :: VAL
      integer(kind=iwp) :: IT, IU

      If (ipea_shift /= Zero) Then
        NS = NAS*(NAS+1)/2
        call mma_allocate(SMat,NS,Label='SMat')
        idS = idSMAT(iSym,10)
        CALL DDAFILE(LUSBT,2,SMat,NS,idS)
        !! ipea_shift*Half*(Two-DREF(IDT))*SD(IT)
        DO IT = 1, NAS
          ITABS=IT+NAES(ISYM)
          VAL = ipea_shift*Half*BDER(IT,IT)
          SDER(IT,IT) = SDER(IT,IT) + (Two-G1(ITABS,ITABS))*VAL
          NSEQ = IT*(IT-1)/2+IT
          DG1(ITABS,ITABS) = DG1(ITABS,ITABS) - SMat(NSEQ)*VAL
        End Do
        call mma_deallocate(SMat)
      End If

      DO IT = 1, NAS
        ITABS=IT+NAES(ISYM)
        DO IU = 1, NAS
          IUABS=IU+NAES(ISYM)
          !! Derivative of the B matrix
          DG1(ITABS,IUABS) = DG1(ITABS,IUABS) - EASUM*BDER(IT,IU)
          DEASUM = DEASUM - G1(ITABS,IUABS)*BDER(IT,IU)
          DF1(ITABS,IUABS) = DF1(ITABS,IUABS) + BDER(IT,IU)
          !! Derivative of the S matrix
          DG1(ITABS,IUABS) = DG1(ITABS,IUABS) + SDER(IT,IU)
        END DO
      END DO

      return

      End Subroutine CLagDXG
!
!-----------------------------------------------------------------------
!
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      Subroutine CLagDX_MPP()

      use caspt2_global, only: iVecL
      use caspt2_module, only: MAXIT, JSTATE
      use Constants, only: One

      implicit none

#include "global.fh"
#include "mafdecls.fh"

      integer(kind=byte), ALLOCATABLE :: idxG3(:,:)
      real(kind=wp),allocatable :: EIG(:),WRK(:,:)

      logical(kind=iwp) :: bStat, invar_act
      integer(kind=iwp) :: myrank, lg_T, lg_WRK, lg_WRK2,               &
     &                     lg_BDER, iLoV1, iHiV1, jLoV1, jHiV1, NROW,   &
     &                     NCOL, idB, mV1, LDV1, i, j, iICB, jICB,      &
     &                     lg_SDER, idSD, mBDER, mSDER
      real(kind=wp) :: SCAL, EigI, EigJ, tmp
!
!     Construct active density in NAS basis
!     Although non-GA version is also implemented, I noticed that
!     scatter operations require GA, so I should just use GA_DGEMM
!
      SCAL = One
      IF (IFMSCOUP) SCAL = VECROT(jState)
      MYRANK=GA_NODEID()

      invar_act = .true.
      if (sigma_p_epsilon /= Zero) invar_act = .false.
      !! First, distribute the transformation matrix
      CALL GA_CREATE_STRIPED ('H',NAS,NIN,'TRANS',lg_T)
      CALL PSBMAT_READ('T',iCase,iSym,lg_T,NAS*NIN)
!     CALL PSBMAT_READ('T',iCase,iSym,lg_T,NAS)
      !! only the master node has written it
!     if (ifmscoup) then
!     else
!     IF (KING()) THEN
!       call mma_allocate(TRANS,NAS*NIN,Label='TRANS')
!       IDT=IDTMAT(ISYM,ICASE)
!       CALL DDAFILE(LUSBT,2,TRANS,NAS*NIN,IDT)
!       CALL GA_PUT(lg_T,1,NAS,1,NIN,TRANS,NAS)
!       call mma_deallocate(TRANS)
!     END IF
!     endif
      CALL GA_SYNC()

      !! Allocate BDER in NIN
      CALL GA_CREATE_STRIPED ('V',NIN,NIN,'WRK',lg_WRK)
      CALL GA_ZERO(lg_WRK)
      CALL GA_SYNC()
!
!     mode = 0 operations for B derivative
!
      !! First, construct the density with NIN basis
      !! lg_V1 = T (solution; not quasi-variational)
      Call RHS_ALLO(nIN,nIS,lg_V1)
      Call RHS_READ_SR(lg_V1,iCase,iSym,iVecX)
      call GA_DGEMM ('N','T',NIN,NIN,NIS,                               &
     &               SCAL,lg_V1,lg_V1,Zero,lg_WRK)

      If ((real_shift /= Zero) .OR. (imag_shift /= Zero)                &
     &    .OR. (sigma_p_epsilon /= Zero) .OR. IFMSCOUP) Then
        if (sigma_p_epsilon /= Zero) then
          !! lg_V2 = lambda (shift correction)
          Call RHS_ALLO(nIN,nIS,lg_V2)
          Call RHS_READ_SR(lg_V2,iCase,iSym,iVecR)
          call mma_allocate(LBD,nAS,Label='LBD')
          call mma_allocate(LID,nIS,Label='LID')
          iD = iDBMat(iSym,iCase)
          Call dDaFile(LUSBT,2,LBD,nAS,iD)
          Call dDaFile(LUSBT,2,LID,nIS,iD)
          !! this scaling is needed, so GA_DGEMM cannot be used
          Call CASPT2_ResD(3,nIN,nIS,lg_V2,lg_V1,LBD,LID)
          call mma_deallocate(LBD)
          call mma_deallocate(LID)
        else
          !! lg_V2 = lambda (shift correction)
          Call RHS_ALLO(nIN,nIS,lg_V2)
          Call RHS_READ_SR(lg_V2,iCase,iSym,iVecR)
        end if
        call GA_DGEMM ('N', 'T', NIN, NIN, NIS, Half,                   &
     &                 lg_V1, lg_V2, One, lg_WRK)
        call GA_DGEMM ('N', 'T', NIN, NIN, NIS, Half,                   &
     &                 lg_V2, lg_V1, One, lg_WRK)
        Call RHS_FREE(lg_V2)
      End If

!     if (sigma_p_epsilon /= Zero) then
!       CALL RHS_READ_SR(lg_V2,ICASE,ISYM,iVecR)
!     endif

      CALL GA_SYNC()
!
!     mode = 1 operations for B derivative
!     lg_V1 is still loaded
!
      if (imag_shift /= Zero .or. sigma_p_epsilon /= Zero) then
        CALL GA_Scale (lg_WRK,-One)

        !! T*T is skipped

        !! lg_V2 = lambda (shift correction)
        Call RHS_ALLO(nIN,nIS,lg_V2)
        Call RHS_READ_SR(lg_V2,iCase,iSym,iVecR)
!       CALL GA_Distribution (lg_V1,myRank,iLoV1,iHiV1,jLoV1,jHiV1)
!       CALL GA_Distribution (lg_V2,myRank,iLoV2,iHiV2,jLoV2,jHiV2)
!       NROW1 = iHiV1-iLoV1+1
!       NCOL1 = jHiV1-jLoV1+1
!       NROW2 = iHiV2-iLoV2+1
!       NCOL2 = jHiV2-jLoV2+1
!       if (nrow1*ncol1*nrow2*ncol2 > 0) then
          call mma_allocate(LBD,nAS,Label='LBD')
          call mma_allocate(LID,nIS,Label='LID')
          iD = iDBMat(iSym,iCase)
          Call dDaFile(LUSBT,2,LBD,nAS,iD)
          Call dDaFile(LUSBT,2,LID,nIS,iD)
          !! this scaling is needed, so GA_DGEMM cannot be used
          Call CASPT2_ResD(2,nIN,nIS,lg_V1,lg_V2,LBD,LID)
          call mma_deallocate(LBD)
          call mma_deallocate(LID)
!         Call GA_Access(lg_V1,iLoV1,iHiV1,jLoV1,jHiV1,mV1,LDV1)
!         Call GA_Access(lg_V2,iLoV2,iHiV2,jLoV2,jHiV2,mV2,LDV2)
!         Call CLag_MPP_NT(Half,DBL_MB(mV1),NROW1,DBL_MB(mV2),NROW2,
!    *                     lg_WRK,NIN,NCOL1)
!         Call CLag_MPP_NT(Half,DBL_MB(mV2),NROW2,DBL_MB(mV1),NROW1,
!    *                     lg_WRK,NIN,NCOL2)
!         Call GA_Release(lg_V1,iLoV1,iHiV1,jLoV1,jHiV1)
!         Call GA_Release(lg_V2,iLoV2,iHiV2,jLoV2,jHiV2)
          CALL GA_DGEMM ('N','T',NIN,NIN,NIS,                           &
     &                   Half,lg_V1,lg_V2,One,lg_WRK)
          CALL GA_DGEMM ('N','T',NIN,NIN,NIS,                           &
     &                   Half,lg_V2,lg_V1,One,lg_WRK)
!       end if
        Call RHS_FREE(lg_V2)

        !! Restore the original T
        Call RHS_READ_SR(lg_V1,iCase,iSym,iVecX)
        CALL GA_Scale (lg_WRK,-One)
        CALL GA_SYNC()
      end if

      Call RHS_FREE(lg_V1)
!
!     B derivative in NIN completed
!
      if (invar_act) then
        CALL GA_CREATE_STRIPED ('H',NAS,NIN,'WRK2',lg_WRK2)
        !! Use the same stripe as PSBMAT_GEMEM?
        CALL GA_CREATE_STRIPED ('H',NAS,NAS,'BDER',lg_BDER)

        !! NIN -> NAS transformation of B derivative
        !! Need 4 GAs; is it possible to reduce?
        !! lg_WRK is used later, so probably not
        CALL GA_DGEMM ('N','N',NAS,NIN,NIN,                             &
     &                 One,lg_T,lg_WRK,Zero,lg_WRK2)
        CALL GA_DGEMM ('N','T',NAS,NAS,NIN,                             &
     &                 One,lg_WRK2,lg_T,Zero,lg_BDER)
!       if (king()) then
!         call mma_allocate(VEC1,NAS*NAS,Label='WRK1')
!         CALL GA_GET(lg_bder,1,NAS,1,NAS,VEC1,NAS)
!         WRITE (*,*) 'B DERIVATIVE IN NAS'
!         CALL SQPRT(VEC1,NAS)
!         call mma_deallocate(VEC1)
!       end if

        !! cannot destroy lg_WRK; it is used for overlap derivative
        bStat = GA_destroy(lg_WRK2)
      else
        !! allocate BDER temporarily with the same dimension of WRK
        CALL GA_CREATE_STRIPED ('V',NIN,NIN,'BDER',lg_BDER)
        CALL GA_COPY(lg_WRK,lg_BDER)
        CALL GA_ZERO(lg_WRK)
      end if

      !! At present, lg_T, lg_WRK, and lg_BDER are in GA.
      !! It is possible to destroy lg_T here and can reduce the max
      !! allocated GA by one, but we anyway have to restore the whole
      !! transformation matrix in memory in the master node, so we can
      !! probably assume that there is no problem in the availability
      !! of memory.

      !! It seems that it is not good to decompose the B derivative
      !! matrix to G and F derivative contributions here, because some
      !! contributions are computed as overlap derivatives, so it is
      !! best to decompose later simultaneously. Fortunately, the max
      !! number of allocated GAs remains unchanged (4 at most).

      !! at this point, lg_T, lg_WRK, and lg_BDER are distributed
!     if (king()) then
!       call mma_allocate(VEC1,NAS*NAS,Label='WRK1')
!       CALL GA_GET(lg_wrk,1,NIN,1,NIN,VEC1,NIN)
!       WRITE (*,*) 'B DERIVATIVE IN NIN'
!       CALL SQPRT(VEC1,NIN)
!       call mma_deallocate(VEC1)
!     end if
!
!     mode = 0 operations for S derivative
!
      !! Scale with the eigenvalue
      if (invar_act) then
        CALL GA_Distribution (lg_WRK,myRank,iLoV1,iHiV1,jLoV1,jHiV1)
        NROW=iHiV1-iLoV1+1
        NCOL=jHiV1-jLoV1+1
        if (NROW > 0 .and. NCOL > 0) then
          call mma_allocate(EIG,NIN,Label='EIG')
          idB  = idBMAT(iSym,iCase)
          CALL DDAFILE(LUSBT,2,EIG,NIN,IDB)
          CALL GA_Access(lg_WRK,iLoV1,iHiV1,jLoV1,jHiV1,mV1,LDV1)

          do j = 1, NCOL
            jICB = j + jLoV1 - 1
            EigJ = EIG(jICB)
            do i = 1, NROW
              iICB = i + iLoV1 - 1
              EigI = EIG(iICB)
              DBL_MB(mV1+i-1+NROW*(j-1))                                &
     &          = -DBL_MB(mV1+i-1+NROW*(j-1))*(EigI+EigJ)*Half
            end do
          end do

          CALL GA_Release_Update(lg_WRK,iLoV1,iHiV1,jLoV1,jHiV1)
          call mma_deallocate(EIG)
        end if
      end if
!     if (king()) then
!       call mma_allocate(VEC1,NIN*NIN,Label='VEC1'))
!       CALL GA_GET(lg_wrk,1,NIN,1,NIN,VEC1,NIN)
!       WRITE (*,*) 'SCALED B DERIVATIVE IN NIN'
!       CALL SQPRT(VEC1,NIN)
!       call mma_deallocate(VEC1)
!     end if

      !! Max memory allocation of 5 GAs
      !  1) Implicit overlap derivative of the 2<1|H|0> part
      !! lg_V1 = T (solution; not quasi-variational)
      Call RHS_ALLO(nIN,nIS,lg_V1)
      Call RHS_READ_SR(lg_V1,iCase,iSym,iVecX)
      !! lg_V2 = RHS2 (in IC basis)
      Call RHS_ALLO(nIN,nIS,lg_V2)
      if (ifmscoup) then
        Call RHS_READ_SR(lg_V2,iCase,iSym,iVecL)
      else
        Call RHS_READ_SR(lg_V2,iCase,iSym,iRHS)
      end if
      call GA_DGEMM ('N','T',NIN,NIN,NIS,                               &
     &              -One,lg_V2,lg_V1,One,lg_WRK)
      Call RHS_FREE(lg_V1)
      Call RHS_FREE(lg_V2)

      If ((real_shift /= Zero) .OR. (imag_shift /= Zero)                &
     &    .OR. (sigma_p_epsilon /= Zero) .OR. IFMSCOUP) Then
        !! WRK1 = -RHS*(T+lambda/2)
        !! lg_V1 = RHS (in IC basis)
        Call RHS_ALLO(nIN,nIS,lg_V1)
        Call RHS_READ_SR(lg_V1,iCase,iSym,iRHS)
        !! lg_V2 = lambda (shift correction)
        Call RHS_ALLO(nIN,nIS,lg_V2)
        Call RHS_READ_SR(lg_V2,iCase,iSym,iVecR)
        call GA_DGEMM ('N','T',NIN,NIN,NIS,                             &
     &                -Half,lg_V1,lg_V2,One,lg_WRK)
        Call RHS_FREE(lg_V1)
        Call RHS_FREE(lg_V2)
      end if
!
!     S derivative in NIN completed (some NAS operations remain)
!
      !! For sigma_p, non-canonical condition
      if (.not.invar_act) then
        call mma_allocate(WRK,NIN,NIN,Label='WRK1')
        CALL GA_GET(lg_WRK,1,NIN,1,NIN,WRK,NIN)
        CALL GA_Distribution (lg_BDER,myRank,iLoV1,iHiV1,jLoV1,jHiV1)
        NROW=iHiV1-iLoV1+1
        NCOL=jHiV1-jLoV1+1
        if (NROW > 0 .and. NCOL > 0) then
          call mma_allocate(EIG,NIN,Label='EIG')
          idB  = idBMAT(iSym,iCase)
          CALL DDAFILE(LUSBT,2,EIG,NIN,IDB)
          CALL GA_Access(lg_BDER,iLoV1,iHiV1,jLoV1,jHiV1,mBDER,LDV1)
          !! construct the off-diagonal
          do i = 1, NROW
            iICB = i + iLoV1 - 1
            EigI = EIG(iICB)
            do j = 1, NCOL
              jICB = j + jLoV1 - 1
              EigJ = EIG(jICB)
              if (iICB == jICB) cycle
              tmp = WRK(iICB,jICB) - WRK(jICB,iICB)
              tmp = tmp/(EigI-EigJ)
              DBL_MB(mBDER+i-1+NROW*(j-1)) = tmp
            end do
          end do
          !! -(e_o + e_p)*dS/da
          CALL GA_Access(lg_WRK,iLoV1,iHiV1,jLoV1,jHiV1,mV1,LDV1)
          do j = 1, NCOL
            jICB = j + jLoV1 - 1
            EigJ = EIG(jICB)
            do i = 1, NROW
              iICB = i + iLoV1 - 1
              EigI = EIG(iICB)
              DBL_MB(mV1+i-1+NROW*(j-1)) = DBL_MB(mV1+i-1+NROW*(j-1))   &
     &          - DBL_MB(mBDER+i-1+NROW*(j-1))*(EigI+EigJ)*Half
            end do
          end do
          CALL GA_Release_Update(lg_BDER,iLoV1,iHiV1,jLoV1,jHiV1)
          CALL GA_Release_Update(lg_WRK,iLoV1,iHiV1,jLoV1,jHiV1)
          call mma_deallocate(EIG)
        end if
        call mma_deallocate(WRK)

        !! IC -> MO (B matrix)
        CALL GA_CREATE_STRIPED ('H',NAS,NIN,'WRK2',lg_WRK2)
        CALL GA_DGEMM ('N','N',NAS,NIN,NIN,                             &
     &                 One,lg_T,lg_BDER,Zero,lg_WRK2)
        bStat = GA_destroy(lg_BDER)
        CALL GA_CREATE_STRIPED ('H',NAS,NAS,'BDER',lg_BDER)
        CALL GA_DGEMM ('N','T',NAS,NAS,NIN,                             &
     &                 One,lg_WRK2,lg_T,Zero,lg_BDER)
        bStat = GA_destroy(lg_WRK2)
      end if

      CALL GA_CREATE_STRIPED ('H',NAS,NIN,'WRK2',lg_WRK2)
      !! Use the same stripe as PSBMAT_GEMEM?
      CALL GA_CREATE_STRIPED ('H',NAS,NAS,'SDER',lg_SDER)

      !! NIN -> NAS transformation of S derivative
      CALL GA_DGEMM ('N','N',NAS,NIN,NIN,                               &
     &               One,lg_T,lg_WRK,Zero,lg_WRK2)
      CALL GA_DGEMM ('N','T',NAS,NAS,NIN,                               &
     &               One,lg_WRK2,lg_T,Zero,lg_SDER)
!     if (king()) then
!       call mma_allocate(VEC1,NAS*NAS,Label='VEC1')
!       CALL GA_GET(lg_sder,1,NAS,1,NAS,VEC1,NAS)
!       WRITE (*,*) 'S DERIVATIVE IN NAS (1)'
!       CALL SQPRT(VEC1,NAS)
!       call mma_deallocate(VEC1)
!     end if

      bStat = GA_destroy(lg_WRK)
      bStat = GA_destroy(lg_WRK2)

      !! Add some trivial contributions due to the dependence
      !! on the linearly independent space
      If (do_lindep .AND. nAS /= nIN) Then
        Call LinDepLag_MPP(lg_BDER,lg_SDER,nAS,nIN,iSym,iCase)
      End If

      !  2) Explicit overlap derivative of the 2<1|H|0> part
      !     Again, not for imaginary shift-specific terms
      !! E = 2<1|H|0> + <1|H0-E0|1>
      !! lg_V1 = VEC1 = T (solution; not quasi-variational)
      Call RHS_ALLO(nIN,nIS,lg_V1)
      Call RHS_READ_SR(lg_V1,iCase,iSym,iVecX)
      Call RHS_SCAL(nIN,nIS,lg_V1,SCAL)
      If ((real_shift /= Zero) .OR. (imag_shift /= Zero)                &
     &    .OR. (sigma_p_epsilon /= Zero) .OR. IFMSCOUP) Then
        !! lg_V2 = VEC2 = lambda (shift correction)
        Call RHS_ALLO(nIN,nIS,lg_V2)
        CALL RHS_READ_SR(lg_V2,iCase,iSym,iVecR)
        CALL RHS_DAXPY(NIN,NIS,Half,lg_V2,lg_V1)
        CALL RHS_FREE(lg_V2)
      end if

      CALL GA_CREATE_STRIPED ('V',NAS,NIS,'WRK',lg_WRK)
      CALL GA_DGEMM ('N','N',NAS,NIS,NIN,                               &
     &               One,lg_T,lg_V1,Zero,lg_WRK)

      Call RHS_FREE(lg_V1)
      bStat = GA_destroy(lg_T)
      !! lg_V1 = VEC4 = RHS (in MO basis)
      Call RHS_ALLO(nAS,nIS,lg_V1)
      Call RHS_READ_C (lg_V1,iCase,iSym,iVecW)

      CALL GA_DGEMM ('N','T',NAS,NAS,NIS,                               &
     &               Two,lg_WRK,lg_V1,One,lg_SDER)

      Call RHS_FREE(lg_V1)
      bStat = GA_destroy(lg_WRK)
      CALL GA_SYNC()

      !! Add the contributions from the off-diagonal coupling
      !! (i.e., CASPT2-N). Of course, this is not for imaginary shift-
      !! specific terms.
      if (MAXIT /= 0) then
        call mma_allocate(WRK,NAS,NAS,Label='WRK')
        idSD = idSDMat(iSym,iCase)
        CALL DDAFILE(LuSTD,2,WRK,nAS*nAS,idSD)
        CALL GA_Distribution(lg_SDER,myRank,ILO,IHI,JLO,JHI)
        Call GA_Access(lg_SDER,ILO,IHI,JLO,JHI,mSDER,LDV)
        NROW = IHI-ILO+1
        NCOL = JHI-JLO+1
        DO J = 1, NCOL
          DO I = 1, NROW
            DBL_MB(mSDER+I-1+NROW*(J-1))                                &
     &        = DBL_MB(mSDER+I-1+NROW*(J-1))                            &
     &        + WRK(I+ILO-1,J+JLO-1)*Half
!    *        + WORK(LWRK+I+ILO-2+NAS*(J+JLO-2))*Half
          END DO
        END DO
        Call GA_Release_Update(lg_SDER,ILO,IHI,JLO,JHI)
        call mma_deallocate(WRK)
      end if
      CALL GA_SYNC()
!     if (king()) then
!       call mma_allocate(VEC1,NAS*NAS,Label='VEC1')
!       CALL GA_GET(lg_sder,1,NAS,1,NAS,VEC1,NAS)
!       WRITE (*,*) 'S DERIVATIVE IN NAS'
!       CALL SQPRT(VEC1,NAS)
!       call mma_deallocate(VEC1)
!     end if

      !! Compute G and F derivative contributions for B and S der
      CALL PSBMAT_GETMEM('S',lg_S,NAS)
      CALL PSBMAT_READ('S',iCase,iSym,lg_S,NAS)
      CALL GA_Distribution(lg_BDER,myRank,ILO,IHI,JLO,JHI)
      CALL GA_Distribution(lg_SDER,myRank,ILO,IHI,JLO,JHI)
      CALL GA_Distribution(lg_S   ,myRank,ILO,IHI,JLO,JHI)
      Call GA_Access(lg_BDER,ILO,IHI,JLO,JHI,mBDER,LDV)
      Call GA_Access(lg_SDER,ILO,IHI,JLO,JHI,mSDER,LDV)
      Call GA_Access(lg_S   ,ILO,IHI,JLO,JHI,mS   ,LDV)
      NROW = iHi-iLo+1 + 1 !! 1 is added so that the work space is used
      NCOL = jHi-jLo+1 + 1 !! for all procs
      call mma_allocate(WRK,NROW,NCOL,Label='WRK')
      if (iCase == 1) then
        Call CLagDXA_DP(iSym,nAS,nAshT,DBL_MB(mBDER),DBL_MB(mSDER),     &
     &                  DG1,DG2,DF1,DF2,DEPSA,DEASUM,                   &
     &                  ILO,IHI,JLO,JHI,LDV,G1,G2,DBL_MB(mS),WRK,       &
     &                  lg_S)
      else if (iCase == 4) then
        Call CLagDXC_DP(iSym,nAS,nAshT,DBL_MB(mBDER),DBL_MB(mSDER),     &
     &                  DG1,DG2,DF1,DF2,DEPSA,DEASUM,                   &
     &                  ILO,IHI,JLO,JHI,LDV,G1,G2,DBL_MB(mS),WRK,       &
     &                  lg_S)
      else
        write (u6,*) 'Invalid iCase in ...'
        call abend()
      end if
      call mma_deallocate(WRK)
      Call GA_Release_Update(lg_BDER,ILO,IHI,JLO,JHI)
      Call GA_Release_Update(lg_SDER,ILO,IHI,JLO,JHI)

      call ga_sync()

      CALL mma_allocate(idxG3,6,NG3,label='idxG3')
      iLUID=0
      CALL I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)

      if (iCase == 1) then
        Call CLagDXA_FG3_MPP(iSym,NASHT,NG3,lg_BDER,lg_SDER,            &
     &                  DG1,DG2,DG3,DF1,DF2,DF3,DEPSA,G2,idxG3)
      else if (iCase == 4) then
        Call CLagDXC_FG3_MPP(iSym,NASHT,NG3,lg_BDER,lg_SDER,            &
     &                  DG1,DG2,DG3,DF1,DF2,DF3,DEPSA,G2,idxG3)

        !! DF3 is done after icase=4
        !! construct G3 matrix in lg_S
!       CALL MKSC_G3_MPP(ISYM,DBL_MB(mS),ILO,IHI,NAS,LDV,
!    &                   NG3,G3,IDXG3)
!       Call GA_Release_Update(lg_S,ILO,IHI,JLO,JHI)
!       call DF3_DEPSA_MPP(DF3,DEPSA,lg_S,idxG3)
      else
        write (u6,*) 'Invalid iCase in ...'
        call abend()
      end if
      Call GA_Release(lg_S   ,ILO,IHI,JLO,JHI)
      CALL mma_deallocate(idxG3)
      CALL PSBMAT_FREEMEM(lg_S)

      bStat = GA_destroy(lg_BDER)
      bStat = GA_destroy(lg_SDER)

#include "macros.fh"
      unused_var(bStat)

      End Subroutine CLagDX_MPP
#endif

      End Subroutine CLagD
