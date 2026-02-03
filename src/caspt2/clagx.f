************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2021, Yoshio Nishimoto                                 *
************************************************************************
      Subroutine CLagX(IFF,CLag,DEPSA,VECROT)

      use caspt2_global, only:iPrGlb
      use Constants, only: Zero
      use definitions, only: wp, iwp, u6
      use PrintLevel, only: verbose
      use stdalloc, only: mma_allocate, mma_deallocate
      use gugx, only: SGS
      use caspt2_module, only: NCONF, NASHT, NASH, ISCF, NSTATE, JSTATE,
     &                         EPSA
      use pt2_guga, only: NG1, NG2, NG3, NG3TOT
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif

      implicit none

#ifdef _MOLCAS_MPP_
#include "global.fh"
#endif

      integer(kind=iwp), intent(in) :: IFF
      real(kind=wp), intent(inout) :: CLag(nConf,nState),
     &                                DEPSA(nAshT,nAshT)
      real(kind=wp), intent(in) :: VECROT(*)

      real(kind=wp), allocatable :: G1(:), G2(:), G3(:)
      real(kind=wp), allocatable :: DG1(:), DG2(:), DG3(:), DF1(:),
     &                              DF2(:), DF3(:)

      Integer(kind=iwp) :: nLev
      real(kind=wp) :: DEASUM, iT, iU
      real(kind=wp) :: CPUT, WALLT, CPE, CPTF0, CPTF10, TIOE, TIOTF0,
     &                 TIOTF10

      nLev=SGS%nLev

      !! reduced density matrix and fock-weighted RDM
      CALL mma_allocate(G1 ,NG1, Label='G1')
      CALL mma_allocate(G2 ,NG2, Label='G2')
      CALL mma_allocate(G3 ,NG3, Label='G3')

      !! their derivative contributions
      NG3tot = NG3
      !! Use NG3tot (in pt2_guga.F90) for the moment
#ifdef _MOLCAS_MPP_
      if (is_real_par()) call gaigop_scal(ng3tot,'+')
#endif
      CALL mma_allocate(DG1,NG1,Label='DG1')
      CALL mma_allocate(DG2,NG2,Label='DG2')
      CALL mma_allocate(DG3,NG3,Label='DG3')
      CALL mma_allocate(DF1,NG1,Label='DF1')
      CALL mma_allocate(DF2,NG2,Label='DF2')
      CALL mma_allocate(DF3,NG3,Label='DF3')

      CALL PT2_GET(NG1,' GAMMA1',G1)
      CALL PT2_GET(NG2,' GAMMA2',G2)
      CALL PT2_GET(NG3,' GAMMA3',G3)
C
      !! Initialize them
      DG1(:) = Zero
      DG2(:) = Zero
      DG3(:) = Zero
      DF1(:) = Zero
      DF2(:) = Zero
      DF3(:) = Zero
      !! DEASUM is the derivative cont. of EASUM
      DEASUM = Zero

      CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
      Call CLagD(G1,G2,G3,
     *           DG1,DG2,DG3,
     *           DF1,DF2,DF3,DEASUM,
     *           DEPSA,VECROT)
      CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
      IF (IPRGLB.GE.verbose) THEN
        CPUT =CPTF10-CPTF0
        WALLT=TIOTF10-TIOTF0
        write(u6,'(a,2f10.2)')" CLagD   : CPU/WALL TIME=", cput,wallt
#ifdef _MOLCAS_MPP_
!       if (is_real_par()) CALL GADSUM ([deasum],1)
#endif
!       write(6,*) "Deasum = ", deasum
#ifdef _MOLCAS_MPP_
!       if (is_real_par()) DEASUM = DEASUM/GA_NNODES()
#endif
      END IF

      !! Some symmetrizations are likely required
      Call CLagSym(nAshT,DG1,DG2,DF1,DF2,0)

      !! Do for the derivative of EASUM
      !! EASUM=EASUM+EPSA(IT)*DREF(IT,IT)
      Do iT = 1, nAsh(1)
         DG1(iT+nAsh(1)*(iT-1)) = DG1(iT+nAsh(1)*(iT-1))
     &                          + DEASUM*EPSA(iT)
        If (ISCF.EQ.0) Then
          Do iU = 1, nAsh(1)
            DEPSA(iT,iU) = DEPSA(iT,iU)
     *        + DEASUM*G1(iT+nAsh(1)*(iU-1))
          End Do
        Else
          !! ?
        End If
      End Do
C
#ifdef _MOLCAS_MPP_
      !! the master node does the job, so distribute to slave nodes
      !! only for the G1 and G2 replicate arrays
      if (is_real_par()) then
        CALL GADSUM (DG1,NG1)
        CALL GADSUM (DG2,NG2)
        CALL GADSUM (DF1,NG1)
        CALL GADSUM (DF2,NG2)
      end if
#endif

      Call CnstCLag(IFF,CLag(1,jState),
     *              DG1,DG2,DG3,
     *              DF1,DF2,DF3,
     *              DEPSA,
     *              G1,G2,G3,nLev)

      Call mma_deallocate(G1)
      Call mma_deallocate(G2)
      Call mma_deallocate(G3)

      Call mma_deallocate(DG1)
      Call mma_deallocate(DG2)
      Call mma_deallocate(DG3)
      Call mma_deallocate(DF1)
      Call mma_deallocate(DF2)
      Call mma_deallocate(DF3)

      End Subroutine CLagX
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CLagD(G1,G2,G3,DG1,DG2,DG3,DF1,DF2,DF3,DEASUM,DEPSA,
     *                 VECROT)

      use caspt2_global, only: imag_shift, iVecL,
     *                         sigma_p_epsilon, LUSBT,
     *                         LUSOLV, ipea_shift
      use Constants, only: Zero, One, Half, Two, Four
      use EQSOLV
      use Sigma_data
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp
      use fake_GA, only: GA_Arrays
      use caspt2_module, only: IFMSCOUP, NSYM, NASH, NAES, NASHT,
     &                         JSTATE, NASUP, NISUP, NINDEP, EPSA,
     &                         NTUES, NTGEUES, NTGTUES
      use pt2_guga, only: NG3
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

      real(kind=wp), intent(inout) :: G1(NASHT,NASHT),
     & G2(NASHT,NASHT,NASHT,NASHT),G3(*),DG1(NASHT,NASHT),
     & DG2(NASHT,NASHT,NASHT,NASHT),DG3(*),DF1(NASHT,NASHT),
     & DF2(NASHT,NASHT,NASHT,NASHT),DF3(*),DEASUM,DEPSA(NASHT,NASHT)
      real(kind=wp), intent(in) :: VECROT(*)

      real(kind=wp), allocatable :: LBD(:),LID(:) !!,VEC1(:),VEC2(:)
      real(kind=wp), allocatable :: SMat(:),BDER(:),SDER(:)
#ifdef _MOLCAS_MPP_
      INTEGER*1, ALLOCATABLE :: idxG3(:,:)
      real(kind=wp), allocatable :: VEC1(:),VEC2(:),VEC3(:),VEC4(:),
     *                              VEC5(:)
#endif
      integer(kind=iwp) :: iCase, iSym, NIN, NIS, NVEC, NAS, ID,
     &                     iLUID
#ifdef _MOLCAS_MPP_
      integer(kind=iwp) :: MYRANK, mS, LDV
#endif
      integer(kind=iwp) :: NS, idS, idum, iTabs, iUabs, iVabs, iXabs,
     &  iYabs, iZabs, iTU, iTU2, iTUabs, iTgeUabs, iTgtUabs, iXY, iXY2,
     &  iXYabs, iXgeYabs, iXgtYabs, iT, iU, iV, iX, iY, NSEQ
      integer(kind=iwp) :: lg_V1, lg_V2, lg_V3, lg_V4, lg_V5, lg_S
      integer(kind=iwp) :: ILO, IHI, JLO, JHI

      real(kind=wp) :: ScalB1, ScalB2, ScalS1, ScalS2, ET, EU, EX, EY,
     & ATUXY, BDERval, bsBDER, SDERval, EASUM, ATYU, ATYX, ATUX, ATUY

      Do iCase = 1, 11
C       cycle
C       If (icase /= 10.and.icase /= 11) cycle ! G
C       If (icase /= 10)                 cycle ! GP
C       If (icase /=  6.and.icase /=  7) cycle ! E
C       If (icase /=  8.and.icase /=  9) cycle ! F
C       If (icase /=  8)                 cycle ! FP
C       If (icase /=  2.and.icase /=  3) cycle ! B
C       If (icase /=  5)                 cycle ! D
C       If (icase /=  4)                 cycle ! C
C       If (icase /=  1)                 cycle ! A
        Do iSym = 1, nSym
          nIN  = nINDEP(iSym,iCase)
          If (nIN == 0) Cycle
          nIS  = nISUP(iSym,iCase)
          NVEC = nIN*nIS
          nAS  = nASUP(iSym,iCase)
          If (nVec == 0) Cycle
C
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
            call CLagDX_MPP
            cycle
          else
#endif
C
C         write(u6,*) "for icase = ", icase
C         write(u6,*) "# of independent vecs:", nin
C         write(u6,*) "# of non-active pairs:", nis
C         write(u6,*) "# of     active pairs:", nas
C         write(u6,*) "dimension for Vec = ", nin*nis
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
C
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

              CALL CLagDX(0,ISYM,ICASE,VEC1,VEC2,
     *                    VEC3,VEC4,
     *                    nIN,nIS,nAS,
     *                    VECROT,VEC5,lg_V2,BDER,SDER)

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
            CALL CLagDX(0,iSym,iCase,GA_Arrays(lg_V1)%A,
     &                               GA_Arrays(lg_V2)%A,
     *                               GA_Arrays(lg_V3)%A,
     &                               GA_Arrays(lg_V4)%A,
     *                  nIN,nIS,nAS,
     *                  VECROT,GA_Arrays(lg_V5)%A,lg_V2,BDER,SDER)
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
                CALL CLagDX(1,ISYM,ICASE,VEC1,VEC2,
     *                      VEC3,VEC4,
     *                      nIN,nIS,nAS,
     *                      VECROT,VEC5,lg_V2,BDER,SDER)

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
              CALL CLagDX(1,iSym,iCase,GA_Arrays(lg_V1)%A,
     &                                 GA_Arrays(lg_V2)%A,
     *                                 GA_Arrays(lg_V3)%A,
     &                                 GA_Arrays(lg_V4)%A,
     *                    nIN,nIS,nAS,
     *                    VECROT,GA_Arrays(lg_V5)%A,lg_V2,BDER,SDER)
#ifdef _MOLCAS_MPP_
            end if
#endif
          End If
C
          !! for non-separable density/derivative
          CALL RHS_READ_SR(lg_V1,ICASE,ISYM,iVecX)
          CALL RHS_READ_SR(lg_V2,ICASE,ISYM,iVecR)
#if defined(_MOLCAS_MPP_) && defined(_GA_)
          end if
#endif

          if (iCase ==  1) call CLagDXA(BDER,SDER)
          if (iCase ==  2) call CLagDXB(BDER,SDER)
          if (iCase ==  3) call CLagDXB(BDER,SDER)
          if (iCase ==  4) call CLagDXC(BDER,SDER)
          if (iCase ==  5) call CLagDXD(BDER,SDER)
          if (iCase ==  6) call CLagDXE(BDER,SDER)
          if (iCase ==  7) call CLagDXE(BDER,SDER)
          if (iCase ==  8) call CLagDXF(BDER,SDER)
          if (iCase ==  9) call CLagDXF(BDER,SDER)
          if (iCase == 10) call CLagDXG(BDER,SDER)
          if (iCase == 11) call CLagDXG(BDER,SDER)

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
C
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
          CALL MKSC_G3_MPP(ISYM,DBL_MB(mS),ILO,IHI,JLO,JHI,LDV,
     &                     NG3,G3,IDXG3)
          Call GA_Release_Update(lg_S,ILO,IHI,JLO,JHI)
          call DF3_DEPSA_MPP(DF3,DEPSA,lg_S,idxG3)
C
          Call GA_Release(lg_S   ,ILO,IHI,JLO,JHI)
          CALL mma_deallocate(idxG3)
          CALL PSBMAT_FREEMEM(lg_S)
        End Do
      end if
#endif
C
      Return
C
      contains
C
C-----------------------------------------------------------------------
C
      Subroutine CLagDXA(BDER,SDER)
C
      implicit none
C
      real(kind=wp), intent(in) :: BDER(NAS,NAS), SDER(NAS,NAS)
C
      NS = NAS*(NAS+1)/2
      call mma_allocate(SMat,NS,Label='SMat')
      idS = idSMAT(iSym,1)
      CALL DDAFILE(LUSBT,2,SMat,NS,idS)
C
      idum=0
      Call CLagDXA_DP (iSym,nAS,BDER,SDER,
     *                 DG1,DG2,DF1,DF2,DEPSA,DEASUM,
     *                 1,nAS,1,nAS,0,G1,G2,SMat,SMat,idum)
C
      !! G3 and F3 relevant
      CALL mma_allocate(idxG3,6,NG3,label='idxG3')
      iLUID=0
      CALL I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)
C     idS = idSMAT(iSym,4)
C     CALL DDAFILE(LUSBT,2,SMat,NS,idS)
      CALL MKSC_G3(iSym,SMat,nG3,G3,idxG3)
      call CLagDXA_FG3(iSym,nAS,NG3,BDER,SDER,
     *                 DF1,DF2,DF3,DG1,DG2,DG3,DEPSA,G2,
     *                 SMat,idxG3)
      call mma_deallocate(idxG3)
C
      call mma_deallocate(SMat)
C
      return
C
      End Subroutine CLagDXA
C
C-----------------------------------------------------------------------
C
      Subroutine CLagDXB(BDER,SDER)
C
      USE SUPERINDEX
C
      implicit none
C
      real(kind=wp), intent(in) :: BDER(NAS,NAS)
      real(kind=wp), intent(inout) :: SDER(NAS,NAS)

      real(kind=wp), allocatable :: WrkBbf(:,:,:,:),WrkSbf(:,:,:,:)
C
      call mma_allocate(WrkBbf,nAshT,nAshT,nAshT,nAshT,Label='WrkBbf')
      call mma_allocate(WrkSbf,nAshT,nAshT,nAshT,nAshT,Label='WrkSbf')
      WrkBbf(:,:,:,:) = Zero
      WrkSbf(:,:,:,:) = Zero
C
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
C         iBadr = iTU + nAS*(iXY-1)
          BDERval = BDER(ITU,IXY)
C
          !! For IPEA shift
          If (iTU ==iXY .and. ipea_shift /= Zero) Then
C           idT=(iTabs*(iTabs+1))/2
            ! idU=(iUabs*(iUabs+1))/2
            NSEQ = iTU*(iTU+1)/2
            bsBDER = ipea_shift*Half*BDERval
!         !! ipea_shift*0.5d0*(DREF(IDT)+DREF(IDU))*SDP(ITGEU)
            DG1(iTabs,iTabs) = DG1(iTabs,iTabs) + bsBDER*SMat(NSEQ)
            DG1(iUabs,iUabs) = DG1(iUabs,iUabs) + bsBDER*SMat(NSEQ)
            SDER(iTU,iXY) = SDER(iTU,iXY)
     *        + (G1(iTabs,iTabs)+G1(iUabs,iUabs))*bsBDER
          End If
          SDERval = SDER(ITU,IXY)
          If (iTabs == iUabs) Then
            BDERval = BDERval*Two
            SDERval = SDERval*Two
          End If
C
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
C
          WRKBBF(iTabs,iUabs,iXabs,iYabs)
     *      = WRKBBF(iTabs,iUabs,iXabs,iYabs) + ScalB1
          WRKBBF(iTabs,iUabs,iYabs,iXabs)
     *      = WRKBBF(iTabs,iUabs,iYabs,iXabs) + ScalB2
          WRKSBF(iTabs,iUabs,iXabs,iYabs)
     *      = WRKSBF(iTabs,iUabs,iXabs,iYabs) + ScalS1
          WRKSBF(iTabs,iUabs,iYabs,iXabs)
     *      = WRKSBF(iTabs,iUabs,iYabs,iXabs) + ScalS2
          If (iTabs /= iUabs) Then
          WRKBBF(iUabs,iTabs,iXabs,iYabs)
     *      = WRKBBF(iUabs,iTabs,iXabs,iYabs) + ScalB2
          WRKBBF(iUabs,iTabs,iYabs,iXabs)
     *      = WRKBBF(iUabs,iTabs,iYabs,iXabs) + ScalB1
          WRKSBF(iUabs,iTabs,iXabs,iYabs)
     *      = WRKSBF(iUabs,iTabs,iXabs,iYabs) + ScalS2
          WRKSBF(iUabs,iTabs,iYabs,iXabs)
     *      = WRKSBF(iUabs,iTabs,iYabs,iXabs) + ScalS1
          End If
        End Do
      End Do
C
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
C
              !! term 1 (w/o delta)
              ATUXY = EASUM-ET-EU-EX-EY
              !! G1 and F1 derivative
              DF2(iX,iT,iY,iU) = DF2(iX,iT,iY,iU) + BDERval
              DG2(iX,iT,iY,iU) = DG2(iX,iT,iY,iU)
     *          - ATUXY*BDERval + SDERval
              !! EASUM derivative
              DEASUM = DEASUM - BDERval*G2(iX,iT,iY,iU)
              !! EPSA derivative
              Do iV = 1, NASHT
                DEPSA(iT,iV) = DEPSA(iT,iV) + BDERval*G2(iX,iV,iY,iU)
                DEPSA(iU,iV) = DEPSA(iU,iV) + BDERval*G2(iX,iT,iY,iV)
                DEPSA(iX,iV) = DEPSA(iX,iV) + BDERval*G2(iV,iT,iY,iU)
                DEPSA(iY,iV) = DEPSA(iY,iV) + BDERval*G2(iX,iT,iV,iU)
              End Do
C
              BDERval = BDERval*Two
              SDERval = SDERval*Two
C
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
C
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
C
              BDERval = BDERval*Half
              SDERval = SDERval*Half
C
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
C
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
C
      call mma_deallocate(WrkBbf)
      call mma_deallocate(WrkSbf)
C
      return
C
      End Subroutine CLagDXB
C
C-----------------------------------------------------------------------
C
      Subroutine CLagDXC(BDER,SDER)

      implicit none

      real(kind=wp), intent(in) :: BDER(NAS,NAS)
      real(kind=wp), intent(inout) :: SDER(NAS,NAS)
C
      NS = NAS*(NAS+1)/2
      call mma_allocate(SMat,NS,Label='SMat')
      idS = idSMAT(iSym,4)
      CALL DDAFILE(LUSBT,2,SMat,NS,idS)
C
      idum = 0
      Call CLagDXC_DP (iSym,nAS,BDER,SDER,
     *                 DG1,DG2,DF1,DF2,DEPSA,DEASUM,
     *                 1,nAS,1,nAS,0,G1,G2,SMat,SMat,idum)
C
      !! G3 and F3 relevant
      CALL mma_allocate(idxG3,6,NG3,label='idxG3')
      iLUID=0
      CALL I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)
C     idS = idSMAT(iSym,4)
C     CALL DDAFILE(LUSBT,2,SMat,NS,idS)
      CALL MKSC_G3(iSym,SMat,nG3,G3,idxG3)
      call CLagDXC_FG3(iSym,nAS,NG3,BDER,SDER,
     *                 DF1,DF2,DF3,DG1,DG2,DG3,DEPSA,G2,
     *                 SMat,idxG3)
      call mma_deallocate(idxG3)
C
      call mma_deallocate(SMat)
C
      return
C
      End Subroutine CLagDXC
C
C-----------------------------------------------------------------------
C
      Subroutine CLagDXD(BDER,SDER)

      USE SUPERINDEX

      implicit none

      real(kind=wp), intent(in) :: BDER(NAS,NAS)
      real(kind=wp), intent(inout) :: SDER(NAS,NAS)

      real(kind=wp) :: BDER1, BDER2, SDER1, SDER2, ETX
C
      If (ipea_shift /= Zero) Then
        NS = NAS*(NAS+1)/2
        call mma_allocate(SMat,NS,Label='SMat')
        idS = idSMAT(iSym,iCase)
        CALL DDAFILE(LUSBT,2,SMat,NS,idS)
      End If
C
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
C
          BDER1 = BDER(iTU ,iXY )
     *          - BDER(iTU ,iXY2)*Half
     *          - BDER(iTU2,iXY )*Half
          BDER2 = BDER(iTU2,iXY2)
C
          !! Derivative of B11
          DF2(iUabs,iTabs,iXabs,iYabs)
     *      = DF2(iUabs,iTabs,iXabs,iYabs) + Two*BDER1
          DG2(iUabs,iTabs,iXabs,iYabs)
     *      = DG2(iUabs,iTabs,iXabs,iYabs) + Two*(ETX-EASUM)*BDER1
          DEASUM = DEASUM - Two*G2(iUabs,iTabs,iXabs,iYabs)*BDER1
          If (iXabs == iTabs) Then
            DF1(iUabs,iYabs) = DF1(iUabs,iYabs) + Two*BDER1
            DG1(iUabs,iYabs) = DG1(iUabs,iYabs) + Two*(ET-EASUM)*BDER1
            DEASUM = DEASUM - Two*G1(iUabs,iYabs)*BDER1
          End If
          DO iV = 1, nAsh(iSym)
            IVABS=IV+NAES(ISYM)
            DEPSA(iTabs,iVabs) = DEPSA(iTabs,iVabs)
     *        + Two*BDER1*G2(iUabs,iVabs,iXabs,iYabs)
            DEPSA(iXabs,iVabs) = DEPSA(iXabs,iVabs)
     *        + Two*BDER1*G2(iUabs,iTabs,iVabs,iYabs)
          End Do
          DEPSA(iTabs,iXabs) = DEPSA(iTabs,iXabs)
     *      + Two*G1(iUabs,iYabs)*BDER1

          !! Derivative of B22
          DF2(iXabs,iTabs,iUabs,iYabs)
     *      = DF2(iXabs,iTabs,iUabs,iYabs) - BDER2
          DG2(iXabs,iTabs,iUabs,iYabs)
     *      = DG2(iXabs,iTabs,iUabs,iYabs) - (ETX-EASUM)*BDER2
          DEASUM = DEASUM + G2(iXabs,iTabs,iUabs,iYabs)*BDER2
          If (iXabs == iTabs) Then
            DF1(iUabs,iYabs) = DF1(iUabs,iYabs) + Two*BDER2
            DG1(iUabs,iYabs) = DG1(iUabs,iYabs) + Two*(EX-EASUM)*BDER2
            DEASUM = DEASUM - Two*G1(iUabs,iYabs)*BDER2
          End If
          DO iV = 1, nAsh(iSym)
            IVABS=IV+NAES(ISYM)
            DEPSA(iTabs,iVabs) = DEPSA(iTabs,iVabs)
     *        - BDER2*G2(iXabs,iVabs,iUabs,iYabs)
            DEPSA(iXabs,iVabs) = DEPSA(iXabs,iVabs)
     *        - BDER2*G2(iVabs,iTabs,iUabs,iYabs)
          End Do
          DEPSA(iXabs,iTabs) = DEPSA(iXabs,iTabs)
     *      + Two*G1(iUabs,iYabs)*BDER2
C
          If (iTU == iXY .and. ipea_shift /= Zero) Then
C      !! ipea_shift*0.5d0*(2.0d0-DREF(IDU)+DREF(IDT))*SD(ITU)
            bsBDER = ipea_shift*Half*BDER(iTU,iXY)
            NSEQ = iTU*(iTU+1)/2
            DG1(iTabs,iTabs) = DG1(iTabs,iTabs) + bsBDER*SMat(NSEQ)
            DG1(iUabs,iUabs) = DG1(iUabs,iUabs) - bsBDER*SMat(NSEQ)
            SDER(iTU,iXY) = SDER(iTU,iXY)
     *        + bsBDER*(Two+G1(iTabs,iTabs)-G1(iUabs,iUabs))
C    !! ipea_shift*0.5d0*(2.0d0-DREF(IDU)+DREF(IDT))*SD(ITU+NAS)
            bsBDER = ipea_shift*Half*BDER(iTU2,iXY2)
            NSEQ = iTU2*(iTU2+1)/2
            DG1(iTabs,iTabs) = DG1(iTabs,iTabs) + bsBDER*SMat(NSEQ)
            DG1(iUabs,iUabs) = DG1(iUabs,iUabs) - bsBDER*SMat(NSEQ)
            SDER(iTU2,iXY2) = SDER(iTU2,iXY2)
     *        + bsBDER*(Two+G1(iTabs,iTabs)-G1(iUabs,iUabs))
          End If
C
          SDER1 = SDER(iTU ,iXY )
     *          - SDER(iTU ,iXY2)*Half
     *          - SDER(iTU2,iXY )*Half
          SDER2 = SDER(iTU2,iXY2)
C
          !! Derivative of S11
          DG2(iUabs,iTabs,iXabs,iYabs)
     *      = DG2(iUabs,iTabs,iXabs,iYabs) + Two*SDER1
          If (iXabs == iTabs) Then
            DG1(iUabs,iYabs) = DG1(iUabs,iYabs) + Two*SDER1
          End If
          !! Derivative of S22
          DG2(iXabs,iTabs,iUabs,iYabs)
     *      = DG2(iXabs,iTabs,iUabs,iYabs) - SDER2
          If (iXabs == iTabs) Then
            DG1(iUabs,iYabs) = DG1(iUabs,iYabs) + Two*SDER2
          End If
        End Do
      End Do
      If (ipea_shift /= Zero) call mma_deallocate(SMat)
C
      return
C
      End Subroutine CLagDXD
C
C-----------------------------------------------------------------------
C
      Subroutine CLagDXE(BDER,SDER)

      implicit none

      real(kind=wp), intent(in) :: BDER(NAS,NAS)
      real(kind=wp), intent(inout) :: SDER(NAS,NAS)

      real(kind=wp) :: VAL
C
      If (ipea_shift /= Zero) Then
        NS = NAS*(NAS+1)/2
        call mma_allocate(SMat,NS,Label='SMat')
        idS = idSMAT(iSym,6)
        CALL DDAFILE(LUSBT,2,SMat,NS,idS)
        !! ipea_shift*0.5d0*DREF(IDT)*SD(IT)
        DO IT = 1, NAS
          ITABS=IT+NAES(ISYM)
          VAL = ipea_shift*Half*BDER(IT,IT)
          SDER(IT,IT) = SDER(IT,IT) + G1(ITABS,ITABS)*VAL
          NSEQ = IT*(IT-1)/2+IT
          DG1(ITABS,ITABS) = DG1(ITABS,ITABS) + SMat(NSEQ)*VAL
        End Do
        call mma_deallocate(SMat)
      End If
C
      DO IT = 1, NAS
        ITABS=IT+NAES(ISYM)
        ET = EPSA(ITABS)
        DO IU = 1, NAS
          IUABS=IU+NAES(ISYM)
          EU = EPSA(IUABS)
          !! Derivative of the B matrix
          !! B_{tu} = -F1_{tu} + (Esum-e_t-e_u)*G1(tu)
          DG1(ITABS,IUABS) = DG1(ITABS,IUABS)
     *      + (EASUM-ET-EU)*BDER(IT,IU)
          DEASUM = DEASUM + G1(ITABS,IUABS)*BDER(IT,IU)
          DF1(ITABS,IUABS) = DF1(ITABS,IUABS) - BDER(IT,IU)
          DO IV = 1, NASH(ISYM)
            IVABS=IV+NAES(ISYM)
            DEPSA(ITABS,IUABS) = DEPSA(ITABS,IUABS)
     *        - G1(ITABS,IVABS)*BDER(IV,IU)
     *        - G1(IUABS,IVABS)*BDER(IV,IT)
          End Do
          DEPSA(ITABS,IUABS) = DEPSA(ITABS,IUABS) + Two*BDER(IT,IU)
          !! Derivative of the S matrix
          DG1(ITABS,IUABS) = DG1(ITABS,IUABS) - SDER(IT,IU)
        END DO
      END DO
C
      return
C
      End Subroutine CLagDXE
C
C-----------------------------------------------------------------------
C
      Subroutine CLagDXF(BDER,SDER)
C
      USE SUPERINDEX

      implicit none

      real(kind=wp), intent(in) :: BDER(NAS,NAS)
      real(kind=wp), intent(inout) :: SDER(NAS,NAS)
C
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
C
          BDERval = BDER(ITU,IXY)
          If (iTU == iXY .and. ipea_shift /= Zero) Then
C           idT=(iTabs*(iTabs+1))/2
            ! idU=(iUabs*(iUabs+1))/2
            NSEQ = iTU*(iTU+1)/2
            bsBDER = ipea_shift*Half*BDERval
C     !! ipea_shift*0.5d0*(4.0d0-DREF(IDT)-DREF(IDU))*SDP(ITGEU)
            DG1(iTabs,iTabs) = DG1(iTabs,iTabs) - SMat(NSEQ)*bsBDER
            DG1(iUabs,iUabs) = DG1(iUabs,iUabs) - SMat(NSEQ)*bsBDER
            SDER(ITU,IXY) = SDER(ITU,IXY)
     *        + (Four-G1(iTabs,iTabs)-G1(iUabs,iUabs))*bsBDER
          End If
          SDERval = SDER(ITU,IXY)
          If (iTabs == iUabs) Then
            BDERval = Two*BDERval
            SDERval = Two*SDERval
          End If
C
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
C
          !! Derivative of the B matrix
          !! B(tuxy) -> PREF(tx,uy)
          DEASUM = DEASUM - ScalB1*G2(iTabs,iXabs,iUabs,iYabs)
     *                    - ScalB2*G2(iTabs,iYabs,iUabs,iXabs)
          If (iTabs /= iUabs)
     *    DEASUM = DEASUM - ScalB2*G2(iUabs,iXabs,iTabs,iYabs)
     *                    - ScalB1*G2(iUabs,iYabs,iTabs,iXabs)
C
          ! iTX = iTabs+nAshT*(iXabs-1)
          ! iUY = iUabs+nAshT*(iYabs-1)
          ! iTY = iTabs+nAshT*(iYabs-1)
          ! iUX = iUabs+nAshT*(iXabs-1)
C
          DF2(iTabs,iXabs,iUabs,iYabs)
     *      = DF2(iTabs,iXabs,iUabs,iYabs) + ScalB1
          DF2(iTabs,iYabs,iUabs,iXabs)
     *      = DF2(iTabs,iYabs,iUabs,iXabs) + ScalB2
          If (iTabs.ne.iUabs) Then
            DF2(iUabs,iXabs,iTabs,iYabs)
     *        = DF2(iUabs,iXabs,iTabs,iYabs) + ScalB2
            DF2(iUabs,iYabs,iTabs,iXabs)
     *        = DF2(iUabs,iYabs,iTabs,iXabs) + ScalB1
          End If
          DG2(iTabs,iXabs,iUabs,iYabs)
     *      = DG2(iTabs,iXabs,iUabs,iYabs) + ScalS1-EASUM*ScalB1
          DG2(iTabs,iYabs,iUabs,iXabs)
     *      = DG2(iTabs,iYabs,iUabs,iXabs) + ScalS2-EASUM*ScalB2
          If (iTabs.ne.iUabs) Then
            DG2(iUabs,iXabs,iTabs,iYabs)
     *        = DG2(iUabs,iXabs,iTabs,iYabs) + ScalS2-EASUM*ScalB2
            DG2(iUabs,iYabs,iTabs,iXabs)
     *        = DG2(iUabs,iYabs,iTabs,iXabs) + ScalS1-EASUM*ScalB1
          End If
        End Do
      End Do
      If (ipea_shift /= Zero) call mma_deallocate(SMat)
C
      return
C
      End Subroutine CLagDXF
C
C-----------------------------------------------------------------------
C
      Subroutine CLagDXG(BDER,SDER)

      implicit none

      real(kind=wp), intent(in) :: BDER(NAS,NAS)
      real(kind=wp), intent(inout) :: SDER(NAS,NAS)

      real(kind=wp) :: VAL
C
      If (ipea_shift /= Zero) Then
        NS = NAS*(NAS+1)/2
        call mma_allocate(SMat,NS,Label='SMat')
        idS = idSMAT(iSym,10)
        CALL DDAFILE(LUSBT,2,SMat,NS,idS)
        !! ipea_shift*0.5d0*(2.0d0-DREF(IDT))*SD(IT)
        DO IT = 1, NAS
          ITABS=IT+NAES(ISYM)
          VAL = ipea_shift*Half*BDER(IT,IT)
          SDER(IT,IT) = SDER(IT,IT) + (Two-G1(ITABS,ITABS))*VAL
          NSEQ = IT*(IT-1)/2+IT
          DG1(ITABS,ITABS) = DG1(ITABS,ITABS) - SMat(NSEQ)*VAL
        End Do
        call mma_deallocate(SMat)
      End If
C
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
C
      return
C
      End Subroutine CLagDXG
C
C-----------------------------------------------------------------------
C
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      Subroutine CLagDX_MPP

      use caspt2_global, only: iVecL
      use caspt2_module, only: MAXIT

      implicit none

#include "global.fh"
#include "mafdecls.fh"
C
      INTEGER*1, ALLOCATABLE :: idxG3(:,:)
      real(kind=wp),allocatable :: EIG(:),WRK(:,:)

      logical(kind=iwp) :: bStat
      integer(kind=iwp) :: myrank, nprocs, lg_T, lg_WRK, lg_WRK2,
     &                     lg_BDER, iLoV1, iHiV1, jLoV1, jHiV1, NROW,
     &                     NCOL, idB, mV1, LDV1, i, j, iICB, jICB,
     &                     lg_SDER, idSD, mBDER, mSDER
      real(kind=wp) :: SCAL, EigI, EigJ
C
C     Construct active density in NAS basis
C     Although non-GA version is also implemented, I noticed that
C     scatter operations require GA, so I should just use GA_DGEMM
C
      SCAL = One
      IF (IFMSCOUP) SCAL = VECROT(jState)
      MYRANK=GA_NODEID()
      NPROCS=GA_NNODES()
C
      !! First, distribute the transformation matrix
      CALL GA_CREATE_STRIPED ('H',NAS,NIN,'TRANS',lg_T)
      CALL PSBMAT_READ('T',iCase,iSym,lg_T,NAS*NIN)
C     CALL PSBMAT_READ('T',iCase,iSym,lg_T,NAS)
      !! only the master node has written it
C     if (ifmscoup) then
C     else
C     IF (KING()) THEN
C       call mma_allocate(TRANS,NAS*NIN,Label='TRANS')
C       IDT=IDTMAT(ISYM,ICASE)
C       CALL DDAFILE(LUSBT,2,TRANS,NAS*NIN,IDT)
C       CALL GA_PUT(lg_T,1,NAS,1,NIN,TRANS,NAS)
C       call mma_deallocate(TRANS)
C     END IF
C     endif
      CALL GA_SYNC()
C
      !! Allocate BDER in NIN
      CALL GA_CREATE_STRIPED ('V',NIN,NIN,'WRK',lg_WRK)
      CALL GA_ZERO(lg_WRK)
      CALL GA_SYNC()
C
C     mode = 0 operations for B derivative
C
      !! First, construct the density with NIN basis
      !! lg_V1 = T (solution; not quasi-variational)
      Call RHS_ALLO(nIN,nIS,lg_V1)
      Call RHS_READ_SR(lg_V1,iCase,iSym,iVecX)
      call GA_DGEMM ('N','T',NIN,NIN,NIS,
     *               SCAL,lg_V1,lg_V1,Zero,lg_WRK)
C
      If ((real_shift /= Zero) .OR. (imag_shift /= Zero)
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
        call GA_DGEMM ('N', 'T', NIN, NIN, NIS, Half,
     &                 lg_V1, lg_V2, One, lg_WRK)
        call GA_DGEMM ('N', 'T', NIN, NIN, NIS, Half,
     &                 lg_V2, lg_V1, One, lg_WRK)
        Call RHS_FREE(lg_V2)
      End If

      if (sigma_p_epsilon /= Zero) then
C       CALL RHS_READ_SR(lg_V2,ICASE,ISYM,iVecR)
      endif
C
      CALL GA_SYNC()
C
C     mode = 1 operations for B derivative
C     lg_V1 is still loaded
C
      if (imag_shift /= Zero .or. sigma_p_epsilon /= Zero) then
        CALL GA_Scale (lg_WRK,-One)
C
        !! T*T is skipped
C
        !! lg_V2 = lambda (shift correction)
        Call RHS_ALLO(nIN,nIS,lg_V2)
        Call RHS_READ_SR(lg_V2,iCase,iSym,iVecR)
C       CALL GA_Distribution (lg_V1,myRank,iLoV1,iHiV1,jLoV1,jHiV1)
C       CALL GA_Distribution (lg_V2,myRank,iLoV2,iHiV2,jLoV2,jHiV2)
C       NROW1 = iHiV1-iLoV1+1
C       NCOL1 = jHiV1-jLoV1+1
C       NROW2 = iHiV2-iLoV2+1
C       NCOL2 = jHiV2-jLoV2+1
C       if (nrow1*ncol1*nrow2*ncol2 > 0) then
          call mma_allocate(LBD,nAS,Label='LBD')
          call mma_allocate(LID,nIS,Label='LID')
          iD = iDBMat(iSym,iCase)
          Call dDaFile(LUSBT,2,LBD,nAS,iD)
          Call dDaFile(LUSBT,2,LID,nIS,iD)
          !! this scaling is needed, so GA_DGEMM cannot be used
          Call CASPT2_ResD(2,nIN,nIS,lg_V1,lg_V2,LBD,LID)
          call mma_deallocate(LBD)
          call mma_deallocate(LID)
C         Call GA_Access(lg_V1,iLoV1,iHiV1,jLoV1,jHiV1,mV1,LDV1)
C         Call GA_Access(lg_V2,iLoV2,iHiV2,jLoV2,jHiV2,mV2,LDV2)
C         Call CLag_MPP_NT(0.5D+00,DBL_MB(mV1),NROW1,DBL_MB(mV2),NROW2,
C    *                     lg_WRK,NIN,NCOL1)
C         Call CLag_MPP_NT(0.5D+00,DBL_MB(mV2),NROW2,DBL_MB(mV1),NROW1,
C    *                     lg_WRK,NIN,NCOL2)
C         Call GA_Release(lg_V1,iLoV1,iHiV1,jLoV1,jHiV1)
C         Call GA_Release(lg_V2,iLoV2,iHiV2,jLoV2,jHiV2)
          CALL GA_DGEMM ('N','T',NIN,NIN,NIS,
     *                   Half,lg_V1,lg_V2,One,lg_WRK)
          CALL GA_DGEMM ('N','T',NIN,NIN,NIS,
     *                   Half,lg_V2,lg_V1,One,lg_WRK)
C       end if
        Call RHS_FREE(lg_V2)
C
        !! Restore the original T
        Call RHS_READ_SR(lg_V1,iCase,iSym,iVecX)
        CALL GA_Scale (lg_WRK,-One)
        CALL GA_SYNC()
      end if
C
      Call RHS_FREE(lg_V1)
C
C     B derivative in NIN completed
C
      CALL GA_CREATE_STRIPED ('H',NAS,NIN,'WRK2',lg_WRK2)
      !! Use the same stripe as PSBMAT_GEMEM?
      CALL GA_CREATE_STRIPED ('H',NAS,NAS,'BDER',lg_BDER)
C
      !! NIN -> NAS transformation of B derivative
      !! Need 4 GAs; is it possible to reduce?
      !! lg_WRK is used later, so probably not
      CALL GA_DGEMM ('N','N',NAS,NIN,NIN,
     *               One,lg_T,lg_WRK,Zero,lg_WRK2)
      CALL GA_DGEMM ('N','T',NAS,NAS,NIN,
     *               One,lg_WRK2,lg_T,Zero,lg_BDER)
C     if (king()) then
C       call mma_allocate(VEC1,NAS*NAS,Label='WRK1')
C       CALL GA_GET(lg_bder,1,NAS,1,NAS,VEC1,NAS)
C       WRITE (*,*) "B DERIVATIVE IN NAS"
C       CALL SQPRT(VEC1,NAS)
C       call mma_deallocate(VEC1)
C     end if
C
      !! cannot destroy lg_WRK; it is used for overlap derivative
      bStat = GA_destroy(lg_WRK2)
C
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
C
      !! at this point, lg_T, lg_WRK, and lg_BDER are distributed
C     if (king()) then
C       call mma_allocate(VEC1,NAS*NAS,Label='WRK1')
C       CALL GA_GET(lg_wrk,1,NIN,1,NIN,VEC1,NIN)
C       WRITE (*,*) "B DERIVATIVE IN NIN"
C       CALL SQPRT(VEC1,NIN)
C       call mma_deallocate(VEC1)
C     end if
C
C     mode = 0 operations for S derivative
C
      !! Scale with the eigenvalue
      CALL GA_Distribution (lg_WRK,myRank,iLoV1,iHiV1,jLoV1,jHiV1)
      NROW=iHiV1-iLoV1+1
      NCOL=jHiV1-jLoV1+1
      if (NROW > 0 .and. NCOL > 0) then
        call mma_allocate(EIG,NIN,Label='EIG')
        idB  = idBMAT(iSym,iCase)
        CALL DDAFILE(LUSBT,2,EIG,NIN,IDB)
        CALL GA_Access(lg_WRK,iLoV1,iHiV1,jLoV1,jHiV1,mV1,LDV1)
C
        do j = 1, NCOL
          jICB = j + jLoV1 - 1
          EigJ = EIG(jICB)
          do i = 1, NROW
            iICB = i + iLoV1 - 1
            EigI = EIG(iICB)
            DBL_MB(mV1+i-1+NROW*(j-1))
     *        = -DBL_MB(mV1+i-1+NROW*(j-1))*(EigI+EigJ)*Half
          end do
        end do
C
        CALL GA_Release_Update(lg_WRK,iLoV1,iHiV1,jLoV1,jHiV1)
        call mma_deallocate(EIG)
      end if
C     if (king()) then
C       call mma_allocate(VEC1,NIN*NIN,Label='VEC1'))
C       CALL GA_GET(lg_wrk,1,NIN,1,NIN,VEC1,NIN)
C       WRITE (*,*) "SCALED B DERIVATIVE IN NIN"
C       CALL SQPRT(VEC1,NIN)
C       call mma_deallocate(VEC1)
C     end if
C
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
      call GA_DGEMM ('N','T',NIN,NIN,NIS,
     *              -One,lg_V2,lg_V1,One,lg_WRK)
      Call RHS_FREE(lg_V1)
      Call RHS_FREE(lg_V2)
C
      If ((real_shift /= Zero) .OR. (imag_shift /= Zero)
     &    .OR. (sigma_p_epsilon /= Zero) .OR. IFMSCOUP) Then
        !! WRK1 = -RHS*(T+lambda/2)
        !! lg_V1 = RHS (in IC basis)
        Call RHS_ALLO(nIN,nIS,lg_V1)
        Call RHS_READ_SR(lg_V1,iCase,iSym,iRHS)
        !! lg_V2 = lambda (shift correction)
        Call RHS_ALLO(nIN,nIS,lg_V2)
        Call RHS_READ_SR(lg_V2,iCase,iSym,iVecR)
        call GA_DGEMM ('N','T',NIN,NIN,NIS,
     *                -Half,lg_V1,lg_V2,One,lg_WRK)
        Call RHS_FREE(lg_V1)
        Call RHS_FREE(lg_V2)
      end if
C
C     S derivative in NIN completed (some NAS operations remain)
C
      CALL GA_CREATE_STRIPED ('H',NAS,NIN,'WRK2',lg_WRK2)
      !! Use the same stripe as PSBMAT_GEMEM?
      CALL GA_CREATE_STRIPED ('H',NAS,NAS,'SDER',lg_SDER)
C
      !! NIN -> NAS transformation of S derivative
      CALL GA_DGEMM ('N','N',NAS,NIN,NIN,
     *               One,lg_T,lg_WRK,Zero,lg_WRK2)
      CALL GA_DGEMM ('N','T',NAS,NAS,NIN,
     *               One,lg_WRK2,lg_T,Zero,lg_SDER)
C     if (king()) then
C       call mma_allocate(VEC1,NAS*NAS,Label='VEC1')
C       CALL GA_GET(lg_sder,1,NAS,1,NAS,VEC1,NAS)
C       WRITE (*,*) "S DERIVATIVE IN NAS (1)"
C       CALL SQPRT(VEC1,NAS)
C       call mma_deallocate(VEC1)
C     end if
C
      bStat = GA_destroy(lg_WRK)
      bStat = GA_destroy(lg_WRK2)
!
      !! Add some trivial contributions due to the dependence
      !! on the linearly independent space
      If (do_lindep .AND. nAS /= nIN) Then
        Call LinDepLag_MPP(lg_BDER,lg_SDER,nAS,nIN,iSym,iCase)
      End If
C
      !  2) Explicit overlap derivative of the 2<1|H|0> part
      !     Again, not for imaginary shift-specific terms
      !! E = 2<1|H|0> + <1|H0-E0|1>
      !! lg_V1 = VEC1 = T (solution; not quasi-variational)
      Call RHS_ALLO(nIN,nIS,lg_V1)
      Call RHS_READ_SR(lg_V1,iCase,iSym,iVecX)
      Call RHS_SCAL(nIN,nIS,lg_V1,SCAL)
      If ((real_shift /= Zero) .OR. (imag_shift /= Zero)
     &    .OR. (sigma_p_epsilon /= Zero) .OR. IFMSCOUP) Then
        !! lg_V2 = VEC2 = lambda (shift correction)
        Call RHS_ALLO(nIN,nIS,lg_V2)
        CALL RHS_READ_SR(lg_V2,iCase,iSym,iVecR)
        CALL RHS_DAXPY(NIN,NIS,Half,lg_V2,lg_V1)
        CALL RHS_FREE(lg_V2)
      end if
C
      CALL GA_CREATE_STRIPED ('V',NAS,NIS,'WRK',lg_WRK)
      CALL GA_DGEMM ('N','N',NAS,NIS,NIN,
     *               One,lg_T,lg_V1,Zero,lg_WRK)
C
      Call RHS_FREE(lg_V1)
      bStat = GA_destroy(lg_T)
      !! lg_V1 = VEC4 = RHS (in MO basis)
      Call RHS_ALLO(nAS,nIS,lg_V1)
      Call RHS_READ_C (lg_V1,iCase,iSym,iVecW)
C
      CALL GA_DGEMM ('N','T',NAS,NAS,NIS,
     *               Two,lg_WRK,lg_V1,One,lg_SDER)
C
      Call RHS_FREE(lg_V1)
      bStat = GA_destroy(lg_WRK)
      CALL GA_SYNC()
C
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
            DBL_MB(mSDER+I-1+NROW*(J-1))
     *        = DBL_MB(mSDER+I-1+NROW*(J-1))
     *        + WRK(I+ILO-1,J+JLO-1)*Half
!    *        + WORK(LWRK+I+ILO-2+NAS*(J+JLO-2))*0.5D+00
          END DO
        END DO
        Call GA_Release_Update(lg_SDER,ILO,IHI,JLO,JHI)
        call mma_deallocate(WRK)
      end if
      CALL GA_SYNC()
C     if (king()) then
C       call mma_allocate(VEC1,NAS*NAS,Label='VEC1')
C       CALL GA_GET(lg_sder,1,NAS,1,NAS,VEC1,NAS)
C       WRITE (*,*) "S DERIVATIVE IN NAS"
C       CALL SQPRT(VEC1,NAS)
C       call mma_deallocate(VEC1)
C     end if
C
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
        Call CLagDXA_DP(iSym,nAS,DBL_MB(mBDER),DBL_MB(mSDER),
     *                  DG1,DG2,DF1,DF2,DEPSA,DEASUM,
     *                  ILO,IHI,JLO,JHI,LDV,G1,G2,DBL_MB(mS),WRK,
     *                  lg_S)
      else if (iCase == 4) then
        Call CLagDXC_DP(iSym,nAS,DBL_MB(mBDER),DBL_MB(mSDER),
     *                  DG1,DG2,DF1,DF2,DEPSA,DEASUM,
     *                  ILO,IHI,JLO,JHI,LDV,G1,G2,DBL_MB(mS),WRK,
     *                  lg_S)
      else
        write (u6,*) "Invalid iCase in ..."
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
        Call CLagDXA_FG3_MPP(iSym,lg_BDER,lg_SDER,
     *                  DG1,DG2,DG3,DF1,DF2,DF3,DEPSA,G2,idxG3)
      else if (iCase == 4) then
        Call CLagDXC_FG3_MPP(iSym,lg_BDER,lg_SDER,
     *                  DG1,DG2,DG3,DF1,DF2,DF3,DEPSA,G2,idxG3)
C
        !! DF3 is done after icase=4
        !! construct G3 matrix in lg_S
C       CALL MKSC_G3_MPP(ISYM,DBL_MB(mS),ILO,IHI,JLO,JHI,LDV,
C    &                   NG3,G3,IDXG3)
C       Call GA_Release_Update(lg_S,ILO,IHI,JLO,JHI)
C       call DF3_DEPSA_MPP(DF3,DEPSA,lg_S,idxG3)
      else
        write (u6,*) "Invalid iCase in ..."
        call abend()
      end if
      Call GA_Release(lg_S   ,ILO,IHI,JLO,JHI)
      CALL mma_deallocate(idxG3)
      CALL PSBMAT_FREEMEM(lg_S)
C
      bStat = GA_destroy(lg_BDER)
      bStat = GA_destroy(lg_SDER)
C
      End Subroutine CLagDX_MPP
#endif
C
      End Subroutine CLagD
C
C-----------------------------------------------------------------------
C
      Subroutine CLagDX(Mode,iSym,iCase,VEC1,VEC2,VEC3,VEC4,nIN,nIS,nAS,
     *                  VECROT,VEC5,lg_V2,BDERmat,SDERmat)
C
      USE SUPERINDEX
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_global, only:real_shift, imag_shift,
     *                        sigma_p_epsilon
      use caspt2_global, only:do_lindep,LUSTD,idSDMat
      use caspt2_global, only: LUSBT
      use EQSOLV
      use Sigma_data
      use definitions, only: wp
      use caspt2_module
      use pt2_guga
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
C
      Implicit Real*8 (A-H,O-Z)
C
      DIMENSION VEC1(*),VEC2(*),VEC3(*),VEC4(*),VEC5(*)
      DIMENSION VECROT(*),BDERmat(*),SDERmat(*)
C
      real(kind=wp),allocatable :: WRK1(:),WRK2(:),WRK3(:),TRANS(:),
     *                             EIG(:)
C     nAshI = nAsh(iSym)
      call mma_allocate(WRK1,nAS**2,Label='WRK1')
      call mma_allocate(WRK2,MAX(nAS**2,nAS*nIS),Label='WRK2')
      call mma_allocate(WRK3,nAS**2,Label='WRK3')
      call mma_allocate(TRANS,nAS*nIN,Label='TRANS')
      call mma_allocate(EIG,nIN,Label='EIG')
C
      idT  = idTMAT(iSym,iCase)
      Call DDAFILE(LUSBT,2,TRANS,nAS*nIN,idT)
      idB  = idBMAT(iSym,iCase)
      Call DDAFILE(LUSBT,2,EIG,nIN,idB)
C
      SCAL = 1.0D+00
      IF (IFMSCOUP) SCAL = VECROT(jState)
C
      !! VEC1: solution in IC basis
      !! VEC2: lambda   in IC basis
      !! VEC3: RHS      in IC basis
      !! VEC4: RHS      in MO basis
C
      !! Form the density in internally contracted basis
      !! The G subspace is employed in the following comments
      !! as an example.
      !! i  : inactive
      !! a,b: secondary
      !! t,u: active
      !! o,p: internally contracted configuration (basis)
      !! WRK1(o,p) = \sum_{iab} T_{o,i}^{ab}*T_{p,i}^{ab}
      !! WRK1 is the effective density in the IC basis,
      !! and will be the B derivative contribution.
C
      If (Mode.eq.0) Then
        !! WRK1 = T*T
        Call DGEMM_('N','T',nIN,nIN,nIS,
     *              SCAL,VEC1,nIN,VEC1,nIN,
     *              0.0D+00,WRK1,nIN)
      Else
        Call DCopy_(nIN*nIN,[0.0D+0],0,WRK1,1)
      End If
C
      If (real_shift .NE. 0.0D+00 .OR. imag_shift .NE. 0.0D+00
     &    .OR. sigma_p_epsilon .ne. 0.0D+00 .OR. IFMSCOUP) Then
        !! WRK1 = T*T + (T*lambda+lambda*T)/2
        !! For sigma-p CASPT2, this if branch computes the pseudo-
        !! density that comes from the numerator of the shift.
        Call DGEMM_('N','T',nIN,nIN,nIS,
     *              0.5D+00,VEC2,nIN,VEC1,nIN,
     *              1.0D+00,WRK1,nIN)
        Call DGEMM_('N','T',nIN,nIN,nIS,
     *              0.5D+00,VEC1,nIN,VEC2,nIN,
     *              1.0D+00,WRK1,nIN)
      End If
      if (sigma_p_epsilon.ne.0.0d+00 .and. mode.eq.0) then
        !! the remaining is the derivative of 2<1|H|0>, so the unscaled
        !! lambda is loaded
#ifdef _MOLCAS_MPP_
        IF (Is_Real_Par()) THEN
          IF (KING()) THEN
            CALL GA_GET(lg_V2,1,NIN,1,NIS,VEC2,NIN)
          End IF
        ELSE
#endif
          CALL RHS_READ_SR(lg_V2,ICASE,ISYM,iVecR)
#ifdef _MOLCAS_MPP_
        end if
#endif
      end if
C
      !! Transform the internally contracted density to
      !! active MO basis
      !! WRK3(t,u) = ST(t,o)*WRK1(o,p)*ST(u,p)
      !! WRK3 is the derivative contribution of the B matrix
      !! in the MO basis
      Call DGEMM_('N','N',nAS,nIN,nIN,
     *            1.0D+00,TRANS,nAS,WRK1,nIN,
     *            0.0D+00,WRK2,nAS)
      Call DGEMM_('N','T',nAS,nAS,nIN,
     *            1.0D+00,WRK2,nAS,TRANS,nAS,
     *            0.0D+00,WRK3,nAS)
C     write(6,*) "B derivative in MO"
C     call sqprt(WRK3,nas)
C
      !! Implicit derivative of the IC vector. This derivative
      !! comes from the derivative of the eigenvalue only. Other
      !! contributions of the derivative of the IC vector is considered
      !! later.
      !! -(e_o + e_p)*dS/da
      Do iICB = 1, nIN
        EigI = EIG(iICB)
        Do jICB = 1, nIN
          EigJ = EIG(jICB)
          WRK1(iICB+nIN*(jICB-1))
     *      = -WRK1(iICB+nIN*(jICB-1))*(EigI+EigJ)*0.5D+00
        End Do
      End Do
C
      !! Derivative of the overlap in the IC basis.
      !! WRK1(o,p) = WRK1(o,p) - T_{o,i}^{ab}*RHS(p,i,a,b)
      !! This contribution should not be done for the imaginary
      !! shift-specific term
      !  1) Implicit overlap derivative of the 2<1|H|0> part
      If (Mode.eq.0) Then
        !! WRK1 = -RHS*T
        Call DGEMM_('N','T',nIN,nIN,nIS,
     *             -1.0D+00,VEC5,nIN,VEC1,nIN,
     *              1.0D+00,WRK1,nIN)
        If (real_shift .NE. 0.0D+00 .OR. imag_shift .NE. 0.0D+00
     &      .OR. sigma_p_epsilon .NE. 0.0D+00 .OR. IFMSCOUP) Then
          !! WRK1 = -RHS*(T+lambda/2)
          Call DGEMM_('N','T',nIN,nIN,nIS,
     *               -0.5D+00,VEC3,nIN,VEC2,nIN,
     *                1.0D+00,WRK1,nIN)
        End If
      End If
C
      !! Convert the IC basis to the MO basis
      Call DGEMM_('N','N',nAS,nIN,nIN,
     *            1.0D+00,TRANS,nAS,WRK1,nIN,
     *            0.0D+00,WRK2,nAS)
      Call DGEMM_('N','T',nAS,nAS,nIN,
     *            1.0D+00,WRK2,nAS,TRANS,nAS,
     *            0.0D+00,WRK1,nAS)
!
      !! Add some trivial contributions due to the dependence
      !! on the linearly independent space
      If (do_lindep.AND.nAS.NE.nIN) Then
        Call LinDepLag(WRK3,WRK1,nAS,nIN,iSym,iCase)
      End If
!
      !  2) Explicit overlap derivative of the 2<1|H|0> part
      !     Again, not for imaginary shift-specific terms
      If (Mode.eq.0) Then
        !! E = 2<1|H|0> + <1|H0-E0|1>
        Call DGEMM_('N','N',nAS,nIS,nIN,
     *              SCAL,TRANS,nAS,VEC1,nIN,
     *              0.0D+00,WRK2,nAS)
        If (real_shift .NE. 0.0D+00 .OR. imag_shift .NE. 0.0D+00
     &      .OR. sigma_p_epsilon .NE. 0.0D+00 .OR. IFMSCOUP) Then
          Call DGEMM_('N','N',nAS,nIS,nIN,
     *                0.5D+00,TRANS,nAS,VEC2,nIN,
     *                1.0D+00,WRK2,nAS)
        END IF
        Call DGEMM_('N','T',nAS,nAS,nIS,
     *              2.0D+00,WRK2,nAS,VEC4,nAS,
     *              1.0D+00,WRK1,nAS)
      End If
C
      !! Add the contributions from the off-diagonal coupling
      !! (i.e., CASPT2-N). Of course, this is not for imaginary shift-
      !! specific terms.
      If (MAXIT.NE.0.and.Mode.eq.0) Then
        idSD = idSDMat(iSym,iCase)
        CALL DDAFILE(LuSTD,2,WRK2,nAS*nAS,idSD)
        !! T*(T+lambda) + (T+lambda)*T is saved, so 1/2
        Call DaXpY_(nAS*nAS,0.5D+00,WRK2,1,WRK1,1)
      End If
C
      if (mode.eq.0) then
        call daxpy_(NAS**2,1.0d+00,WRK3,1,BDERmat,1)
        call daxpy_(NAS**2,1.0d+00,WRK1,1,SDERmat,1)
      else
        call daxpy_(NAS**2,-1.0d+00,WRK3,1,BDERmat,1)
        call daxpy_(NAS**2,-1.0d+00,WRK1,1,SDERmat,1)
      end if
      !! Now, convert the above contributions to derivatives of RDM,
      !! weighted Fock, etc.
      !! WRK3 is the derivative of B in the MO basis
      !! WRK1 is the derivative of S in the MO basis
      !! See and be consistent with mkbmat.f and mksmat.f
      !! Notice that F2 and G2 in mkbmat.f and mksmat.f are halved
      !! (see getdpref.f).
      !! they are moved to CLagDXA, ... CLagDXG
      IF (iCase.eq.1) Then
      Else If (iCase.eq. 2.or.iCase.eq. 3) Then !! B
      Else If (iCase.eq. 4) Then !! C
      Else If (iCase.eq. 5) Then !! D
      Else If (iCase.eq. 6.or.iCase.eq. 7) Then !! E
      Else If (iCase.eq. 8.or.iCase.eq. 9) Then !! F
      Else If (iCase.eq.10.or.iCase.eq.11) Then !! G
      End If
C
      call mma_deallocate(WRK1)
      call mma_deallocate(WRK2)
      call mma_deallocate(WRK3)
      call mma_deallocate(TRANS)
      call mma_deallocate(EIG)
C
      Return
C
      End Subroutine CLagDX
C
C-----------------------------------------------------------------------
C
      !! From poly3
      SUBROUTINE CnstCLag(IFF,CLag,DG1,DG2,DG3,DF1,DF2,DF3,DEPSA,
     *                    G1,G2,G3,nLev)

      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_global, only: iPrGlb
      use PrintLevel, only: verbose
      use gugx, only: L2ACT
      use caspt2_global, only: LUCIEX, IDTCEX, LUSOLV
      use definitions, only: wp
      use caspt2_module
      use pt2_guga
      IMPLICIT REAL*8 (A-H,O-Z)

C
      Integer, Intent(In):: nLev
      DIMENSION CLag(nConf)
      DIMENSION DG1(*),DG2(*),DG3(*),DF1(*),DF2(*),DF3(*)
      DIMENSION G1(*),G2(*),G3(*)
      DIMENSION DEPSA(*)
C
      INTEGER ILEV
      INTEGER NG3MAX
      INTEGER ILUID
      integer*1, allocatable :: idxG3(:,:)
      real(kind=wp),allocatable :: CI1(:)
C
      INTEGER IDCI
      INTEGER J
C
      INTEGER IPARDIV
C
      IF (IFF.EQ.1) THEN
C ORBITAL ENERGIES IN CI-COUPLING ORDER:
        DO ILEV=1,NLEV
          ETA(ILEV)=EPSA(L2ACT(ILEV))
        END DO
      END IF

C-SVC20100831: recompute approximate max NG3 size needed
      NG3MAX=iPARDIV(NG3TOT,NG2)

C-SVC20100831: allocate local G3 matrices
      CALL mma_allocate(idxG3,6,NG3,label='idxG3')
      iLUID=0
      CALL I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)
* NG3 will change inside subroutine MKFG3 to the actual
* number of nonzero elements, that is why here we allocate
* with NG3MAX, but we only store (PT2_PUT) the first NG3
* elements of the G3 and F3
      IF (ISCF.EQ.0) NG3=NG3MAX

      call mma_allocate(CI1,NCONF,LABEL='CI')
      If (ISCF.EQ.0) Then
        if (iff.eq.1) then
          IDCI=IDTCEX
          DO J=1,JSTATE-1
            CALL DDAFILE(LUCIEX,0,CI1,NCONF,IDCI)
          END DO
          CALL DDAFILE(LUCIEX,2,CI1,NCONF,IDCI)
        else
C         Call LoadCI_XMS('C',1,CI1,JSTATE,U0)
        end if
        IF (IPRGLB.GE.VERBOSE) THEN
          WRITE(6,*)
          IF (NSTATE.GT.1) THEN
            WRITE(6,'(A,I4)')
     &      ' With new orbitals, the CI array of state ',MSTATE(JSTATE)
          ELSE
            WRITE(6,*)' With new orbitals, the CI array is:'
          END IF
          CALL PRWF_CP2(STSYM,NCONF,CI1,CITHR)
        END IF
      Else
        CI1(1) = 1.0D+00
      End If
C
      CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
      If (ISCF.EQ.0) Then
        CALL DERFG3(CI1,CLAG,DG1,DG2,DG3,DF1,DF2,DF3,
     &              idxG3,DEPSA,G1,G2,nLev)
      Else
        CALL DERSPE(DF1,DF2,DF3,idxG3,DEPSA,G1,G2,G3)
      End If
      CALL TIMING(CPTF10,CPE,TIOTF10,TIOE)
      IF (IPRGLB.GE.verbose) THEN
        CPUT =CPTF10-CPTF0
        WALLT=TIOTF10-TIOTF0
        write(6,*)
        write(6,'(a,2f10.2)')" DERFG3  : CPU/WALL TIME=", cput,wallt
      END IF
C
      call mma_deallocate(CI1)
      call mma_deallocate(idxG3)
C
      Return
C
      End Subroutine CnstCLag
C
C-----------------------------------------------------------------------
C
      !! From poly3
      SUBROUTINE CLagEig(if_SSDMloc,force_equal,CLag,RDMEIG,nLev)
C
      use caspt2_global, only: DREF, DWGT
      use caspt2_global, only: OMGDER, Weight
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp
      use caspt2_module
      use pt2_guga
      IMPLICIT REAL*8 (A-H,O-Z)
C
      Integer, Intent(In)::nLev
      DIMENSION CLag(nConf,nState),RDMEIG(*)
      Logical   if_SSDMloc, force_equal
      real(kind=wp),allocatable :: CI1(:),WRK(:)
C
C     MODE=0: Either state-averaged or DWGT matrix
C     MODE=1: XMS-specific term, always state-averaged DM
C
      !! RDMEIG
      call mma_allocate(CI1,nConf,Label='LCI')

      Do iState = 1, nState
        If (.not.if_SSDMloc) Then
          if (force_equal .or. .not.IFSADREF) then
            WGT = 1.0D+00/nState ! force equal-weight for XMS
          else if (IFSADREF) then
            WGT = Weight(iState) ! can be unequal weight
          else
            WGT = 1.0D+00/nState ! this should not happen...
          end if
          if (abs(wgt).le.1.0d-09) cycle
          If (ISCF.EQ.0) Then
            Call LoadCI(CI1,iState)
          Else
            CI1(1) = 1.0D+00
          End If
          Call DScal_(NLEV*NLEV,WGT,RDMEIG,1)
          Call Poly1_CLag(CI1,CLag(1,iState),RDMEIG,nLev)
          Call DScal_(NLEV*NLEV,1.0D+00/WGT,RDMEIG,1)
        Else
          Wgt = DWgt(iState,jState)
          If (abs(wgt).gt.1.0d-09) Then
            If (ISCF.EQ.0) Then
              Call LoadCI(CI1,iState)
            Else
              CI1(1) = 1.0D+00
            End If
            !! how is the numerical precision?
            Call DScal_(NLEV*NLEV,WGT,RDMEIG,1)
            Call Poly1_CLag(CI1,CLag(1,iState),RDMEIG,nLev)
            Call DScal_(NLEV*NLEV,1.0D+00/WGT,RDMEIG,1)
          End If

          !! Derivative of omega for dynamically weighted density
          If (IFDW .and. zeta >= 0.0d0) Then
            If (ISCF.EQ.0) Then
              Call LoadCI(CI1,iState)
            Else
              CI1(1) = 1.0D+00
            End If
            call mma_allocate(WRK,nAshT**2,Label='WRK')
            call POLY1(CI1,nConf)
            call GETDREF(DREF,SIZE(DREF))
            Call SQUARE(DREF,WRK,1,nAshT,nAshT)
            !! probably it is doubled somewhere, so should half
            Scal = DDOT_(nAshT**2,RDMEIG,1,WRK,1)*0.5d+00
C           write (*,*) "scal = ", scal
            call mma_deallocate(WRK)
            OMGDER(iState,jState) = OMGDER(iState,jState) + Scal
          End If

        End If
      End Do
C     write(6,*) "clag before projection"
C     do istate = 1, nstate
C       write(6,*) "state = ", istate
C       do i = 1, nconf
C         write(6,'(i3,f20.10)') i,clag(i,istate)
C       end do
C     end do
C     write(6,*) "debug"
C     IF(ORBIN.EQ.'TRANSFOR') Call CLagX_TrfCI(CLAG)
C     if (proj) then
C     ovl = ddot_(nconf*nstate,ci1,1,clag,1)
C     write(6,*) "projection coeff = ",ovl
C     call daxpy_(nconf*nstate,-ovl,ci1,1,clag,1)
C     write(6,*) "clag after projection"
C     do istate = 1, nstate
C       write(6,*) "state = ", istate
C       do i = 1, nconf
C         write(6,'(i3,f20.10)') i,clag(i,istate)
C       end do
C     end do
C     end if
C
      call mma_deallocate(CI1)
C
      Return
C
      End Subroutine CLagEig
C
C-----------------------------------------------------------------------
C
      Subroutine CLagFinal(CLag,SLag)
C
      use caspt2_global, only: iPrGlb
      use PrintLevel, only: verbose
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp
      use caspt2_module
      IMPLICIT REAL*8 (A-H,O-Z)
C
      Dimension CLag(nConf,nState),SLag(*)
      real(kind=wp),allocatable :: CI1(:),CI2(:)
C

      call mma_allocate(CI1,nConf,Label='CI1')
      call mma_allocate(CI2,nConf,Label='CI2')
C
      !! Construct SLag
      ijst = 0
      do ilStat = 1, nState
        If (ISCF.EQ.0) Then
          Call LoadCI(CI1,ilStat)
        Else
          CI1(1) = 1.0D+00
        End If
        Do jlStat = 1, ilStat !! -1
          ijst = ilStat + nState*(jlStat-1)
          If (ilStat.eq.jlStat) Cycle
          If (ISCF.EQ.0) Then
            Call LoadCI(CI2,jlStat)
          Else
            CI2(1) = 1.0D+00
          End If
          Scal = DDOT_(nConf,CI1,1,CLag(1,jlStat),1)
     *         - DDOT_(nConf,CI2,1,CLag(1,ilStat),1)
          Scal = Scal/(REFENE(jlStat)-REFENE(ilStat))
          SLag(ijst) = SLag(ijst) + Scal
          IF (IPRGLB.GE.VERBOSE) THEN
            write(6,*)
            write(6,'(1x,"SLag for State ",i1,"-",i1," = ",f20.10)')
     *         ilstat,jlstat,slag(ijst)
            write(6,*)
          END IF
        end do
      end do
C
      !! This projection is required to get convergence in MCLR.
      !! also the linear equation for non-invariant CASPT2
      Do ilStat = 1, nState
        Call DCopy_(nConf,CLag(1,ilStat),1,CI1,1)
C       do i = 1, nconf
C         write(6,'(i3,f20.10)') i,clag(i,ilstat)
C       end do
        Do jlStat = 1, nState
          If (ISCF.EQ.0) Then
            Call LoadCI(CI2,jlStat)
          Else
            CI2(1) = 1.0D+00
          End If
          Ovl = DDot_(nConf,CI1,1,CI2,1)
C         write(6,*) "projection coeff = ",ovl
          Call DaXpY_(nConf,-Ovl,CI2,1,CLag(1,ilStat),1)
        End Do
C       write(6,*) "clag after projection"
C       write(6,*) "state = ", ilstat
C       do i = 1, nconf
C         write(6,'(i3,f20.10)') i,clag(i,ilstat)
C       end do
      End Do
C
      call mma_deallocate(CI1)
      call mma_deallocate(CI2)
C
      Return
C
      End Subroutine CLagFinal
C
C-----------------------------------------------------------------------
C
      SUBROUTINE POLY1_CLag(CI,CLag,RDMEIG,nLev)
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp
      use caspt2_module, only: nConf
      use pt2_guga, only: MxCI, iAdr10, cLab10
      IMPLICIT NONE
* PER-AAKE MALMQUIST, 92-12-07
* THIS PROGRAM CALCULATES THE 1-EL DENSITY
* MATRIX FOR A CASSCF WAVE FUNCTION.
      INTEGER, INTENT(IN) :: nLev
      REAL*8, INTENT(IN) :: CI(NCONF)

      INTEGER I
      REAL*8 :: CLag(*), RDMEIG(*)
      real(kind=wp),allocatable :: SGM1(:)


      IF(NLEV.GT.0) THEN
        CALL MMA_ALLOCATE(SGM1,MXCI,LABEL='SGM1')
        CALL DENS1_RPT2_CLag(CI,SGM1,CLag,RDMEIG,nLev)
      END IF
C     return !! for test purpose

* REINITIALIZE USE OF DMAT.
* The fields IADR10 and CLAB10 are kept in pt2_guga.F90
* CLAB10 replaces older field called LABEL.
      DO I=1,64
        IADR10(I,1)=-1
        IADR10(I,2)=0
        CLAB10(I)='   EMPTY'
      END DO
      IADR10(1,1)=0
* HENCEFORTH, THE CALL PUT(NSIZE,LABEL,ARRAY) WILL ENTER AN
* ARRAY ON LUDMAT AND UPDATE THE TOC.
      IF(NLEV.GT.0) THEN
        CALL MMA_DEALLOCATE(SGM1)
      END IF


      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DENS1_RPT2_CLag (CI,SGM1,CLag,RDMEIG,nLev)
      use gugx, only: SGS, L2ACT, CIS
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: iwp
      use caspt2_module, only: nConf, nState, STSym, Mul
#if defined (_MOLCAS_MPP_) && ! defined (_GA_)
      USE Para_Info, ONLY: Is_Real_Par, King, nProcs
#endif
      use pt2_guga, only: MxCI
      IMPLICIT NONE

      LOGICAL RSV_TSK
      INTEGER, INTENT(IN):: nLev
      REAL*8 CI(MXCI),SGM1(MXCI)
      REAL*8 CLag(nConf,nState),RDMEIG(NLEV,NLEV) !! Symmetry?
      integer(kind=iwp),allocatable :: TASK(:,:)

C     REAL*8 GTU

      INTEGER ID
      INTEGER IST,ISU,ISTU
      INTEGER IT,IU,LT,LU

      INTEGER ITASK,NTASKS

      INTEGER ISSG,NSGM

* Purpose: Compute the 1-electron density matrix array G1.


* For the general cases, we use actual CI routine calls, and
* have to take account of orbital order.
* We will use level inices LT,LU... in these calls, but produce
* the density matrices with usual active orbital indices.
* Translation tables L2ACT and LEVEL, in pt2_guga.F90

* SVC20100311: set up a task table with LT,LU
* SB20190319: maybe it doesn't even make sense to parallelize the 1-RDM
      nTasks=(nLev**2+nLev)/2
      nTasks = nLev**2
      CALL mma_allocate (Task,nTasks,2,Label='TASK')

      iTask=0
      ! First, IL < JL pairs.
      Do LT = 1, nLev-1
        Do LU = LT+1, nLev
          iTask = iTask + 1
          TASK(iTask,1)=LT
          TASK(iTask,2)=LU
        End Do
      End Do
      ! Then, IL = JL pairs.
      Do LT = 1, nLev
        iTask = iTask + 1
        TASK(iTask,1)=LT
        TASK(iTask,2)=LT
      End Do
      ! Last, IL > JL pairs.
      Do LT = 2, nLev
        Do LU = 1, LT-1
          iTask = iTask + 1
          TASK(iTask,1)=LT
          TASK(iTask,2)=LU
        End Do
      End Do
      IF (iTask.NE.nTasks) WRITE(6,*) "ERROR nTasks"

      Call Init_Tsk(ID, nTasks)

* SVC20100311: BEGIN SEPARATE TASK EXECUTION
 500  If (.NOT.Rsv_Tsk (ID,iTask)) GOTO 501

* Compute SGM1 = E_UT acting on CI, with T.ge.U,
* i.e., lowering operations. These are allowed in RAS.
      LT=TASK(iTask,1)
        IST=SGS%ISM(LT)
        IT=L2ACT(LT)
        LU=Task(iTask,2)
          ISU=SGS%ISM(LU)
          IU=L2ACT(LU)
          ISTU=MUL(IST,ISU)
          ISSG=MUL(ISTU,STSYM)
          NSGM=CIS%NCSF(ISSG)
          IF(NSGM.EQ.0) GOTO 500
* GETSGM2 computes E_UT acting on CI and saves it on SGM1
          CALL GETSGM2(LU,LT,STSYM,CI,SGM1)
          IF(ISTU.EQ.1) THEN
            ! Symmetry not yet
C            write(6,*) "it,iu = ", it,iu
            Call DaXpY_(NSGM,RDMEIG(IT,IU),SGM1,1,CLag,1)
C           if (IT.ne.IU)
C    *        Call DaXpY_(NSGM,2.0d+00*RDMEIG(IT,IU),SGM1,1,CLag,1)

C           GTU=DDOT_(NSGM,CI,1,SGM1,1)
C           G1(IT,IU)=GTU
C           G1(IU,IT)=GTU
          END IF

* SVC: The master node now continues to only handle task scheduling,
*     needed to achieve better load balancing. So it exits from the task
*      list. It has to do it here since each process gets at least one
*      task.
#if defined (_MOLCAS_MPP_) && ! defined (_GA_)
      IF (IS_REAL_PAR().AND.KING().AND.(NPROCS.GT.1)) GOTO 501
#endif

      GOTO 500
 501  CONTINUE

      CALL Free_Tsk(ID)

      CALL mma_deallocate(Task)

      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MLTUNF2(LST,X)
      use Sigma_data
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*)
      DIMENSION LST(4,NLST1)
      DO ILST=1,NLST1
        L1=LST(1,ILST)
        L2=LST(2,ILST)
        L3=LST(3,ILST)
        L4=LST(4,ILST)
        V=VAL1(L4)
        IY=1+INCY2*(L3-1)
        CALL DScal_(LEN1,V,X(IY),INCY1)
        write(6,'(5i4,f20.10,2i4)') ilst,l1,l2,l3,l4,v,iy,incy1
      END DO
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      !! Taken from grdctl.f
      SUBROUTINE CLagX_TrfCI(CI)
C
      use caspt2_global, only: TAT, TORB
      use caspt2_module
      IMPLICIT REAL*8 (A-H,O-Z)
C
      REAL*8 CI(*)
C
      TAT(:)=0.0D0
C
      IOFF1=0
      IOFF2=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NR1=NRAS1(ISYM)
        NR2=NRAS2(ISYM)
        NR3=NRAS3(ISYM)
        NS=NSSH(ISYM)
* Skip inactive transformation matrix:
        IOFF1=IOFF1+NI**2
* Copy RAS1 transformation matrix transposed to TAT:
        DO I=1,NR1
          DO J=1,NR1
            IJ=I+NR1*(J-1)
            JI=J+NR1*(I-1)
            TAT(IOFF2+JI)=TORB(IOFF1+IJ)
          END DO
        END DO
        IOFF1=IOFF1+NR1**2
        IOFF2=IOFF2+NR1**2
* Copy RAS2 transformation matrix transposed to TAT:
        DO I=1,NR2
          DO J=1,NR2
            IJ=I+NR2*(J-1)
            JI=J+NR2*(I-1)
            TAT(IOFF2+JI)=TORB(IOFF1+IJ)
          END DO
        END DO
        IOFF1=IOFF1+NR2**2
        IOFF2=IOFF2+NR2**2
* Copy RAS2 transformation matrix transposed to TAT:
        DO I=1,NR3
          DO J=1,NR3
            IJ=I+NR3*(J-1)
            JI=J+NR3*(I-1)
            TAT(IOFF2+JI)=TORB(IOFF1+IJ)
          END DO
        END DO
        IOFF1=IOFF1+NR3**2
        IOFF2=IOFF2+NR3**2
* Skip virtual transformation matrix:
        IOFF1=IOFF1+NS**2
      END DO
C Transform SGM to use original MO:
      ITOEND=0
      NSG=NCONF
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NR1=NRAS1(ISYM)
        NR2=NRAS2(ISYM)
        NR3=NRAS3(ISYM)
        NS=NSSH(ISYM)
        ITOSTA=ITOEND+1
        ITOEND=ITOEND+NR1**2+NR2**2+NR3**2
*        ITO=ITOSTA+NI**2
        ITO=ITOSTA
        IF(NR1.GT.0) THEN
          ISTART=NAES(ISYM)+1
          CALL TRACI_RPT2(ISTART,NR1,TAT(ITO),STSYM,
     &                                         NSG,CI)
        END IF
        ITO=ITO+NR1**2
        IF(NR2.GT.0) THEN
          ISTART=NAES(ISYM)+NR1+1
          CALL TRACI_RPT2(ISTART,NR2,TAT(ITO),STSYM,
     &                                         NSG,CI)
        END IF
        ITO=ITO+NR2**2
        IF(NR3.GT.0) THEN
          ISTART=NAES(ISYM)+NR1+NR2+1
         !! NR1 should be NR3?
          CALL TRACI_RPT2(ISTART,NR3,TAT(ITO),STSYM,
     &                                         NSG,CI)
        END IF
      END DO
C
      RETURN
C
      END SUBROUTINE CLagX_TrfCI
C
C-----------------------------------------------------------------------
C
      Subroutine CLagSym(nAshT,DG1,DG2,DF1,DF2,mode)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension DG1(nAshT,nAshT),DG2(nAshT,nAshT,nAshT,nAshT),
     *          DF1(nAshT,nAshT),DF2(nAshT,nAshT,nAshT,nAshT)
C
C     return
C     if (mode.eq.0) then
      Do iI = 1, nAshT
        Do iJ = 1, iI-1
          Val1 = DG1(iI,iJ)
          Val2 = DG1(iJ,iI)
          DG1(iI,iJ) = (Val1+Val2)*0.5D+00
          DG1(iJ,iI) = (Val1+Val2)*0.5D+00
          Val1 = DF1(iI,iJ)
          Val2 = DF1(iJ,iI)
          DF1(iI,iJ) = (Val1+Val2)*0.5D+00
          DF1(iJ,iI) = (Val1+Val2)*0.5D+00
        End Do
      End Do
C     end if
C
      If (mode.eq.0) Then
        !! Follow G2 symmetry
        Do iI = 1, nAshT
        Do iJ = 1, nAshT
        Do iK = 1, nAshT
        Do iL = 1, nAshT
            Val1 = DG2(iI,iJ,iK,iL)
            Val2 = DG2(iJ,iI,iL,iK)
            Val3 = DG2(iK,iL,iI,iJ)
            Val4 = DG2(iL,iK,iJ,iI)
            Val  = (Val1+Val2+Val3+Val4)*0.25D+00
            DG2(iI,iJ,iK,iL) = Val
            DG2(iJ,iI,iL,iK) = Val
            DG2(iK,iL,iI,iJ) = Val
            DG2(iL,iK,iJ,iI) = Val
            Val1 = DF2(iI,iJ,iK,iL)
            Val2 = DF2(iJ,iI,iL,iK)
            Val3 = DF2(iK,iL,iI,iJ)
            Val4 = DF2(iL,iK,iJ,iI)
            Val  = (Val1+Val2+Val3+Val4)*0.25D+00
            DF2(iI,iJ,iK,iL) = Val
            DF2(iJ,iI,iL,iK) = Val
            DF2(iK,iL,iI,iJ) = Val
            DF2(iL,iK,iJ,iI) = Val
        End Do
        End Do
        End Do
        End Do
      Else If (mode.eq.1) Then
        !! Follow EtuEyz symmetry
        Do iI = 1, nAshT
        Do iJ = 1, nAshT
        Do iK = 1, nAshT
        Do iL = 1, nAshT
            Val1 = DG2(iI,iJ,iK,iL)
            Val2 = DG2(iL,iK,iJ,iI)
            Val  = (Val1+Val2)*0.5D+00
          ! DG2(iI,iJ,iK,iL) = Val
          ! DG2(iL,iK,iJ,iI) = Val
C           if (ii.ne.il.and.ij.ne.ik) then
C           DG2(iI,iJ,iK,iL) = 2.0d+00*val
C           DG2(iL,iK,iJ,iI) = 0.0d+00
C           end if
            Val1 = DF2(iI,iJ,iK,iL)
            Val2 = DF2(iL,iK,iJ,iI)
            Val  = (Val1+Val2)*0.5D+00
          ! DF2(iI,iJ,iK,iL) = Val
          ! DF2(iL,iK,iJ,iI) = Val
C           if (ii.ne.il.and.ij.ne.ik) then
C           DF2(iI,iJ,iK,iL) = 2.0d+00*val
C           DF2(iL,iK,iJ,iI) = 0.0d+00
C           end if
        End Do
        End Do
        End Do
        End Do
C       write(6,*) "asdf"
C       Do iI = 1, nAshT
C       Do iJ = 1, iI
C       Do iK = 1, iJ
C       Do iL = 1, iK
C       write(6,'(4i3,f20.10)') ii,ij,ik,il,dg2(ii,ij,ik,il)
C       DG2(iI,iJ,iK,iL) = DG2(iI,iJ,iK,iL) + DG2(iL,iK,iJ,iI)
C       if (ii.ne.il.and.ij.ne.ik) DG2(iL,iK,iJ,iI) = 0.0D+00
C       DF2(iI,iJ,iK,iL) = DF2(iI,iJ,iK,iL) + DF2(iL,iK,iJ,iI)
C       if (ii.ne.il.and.ij.ne.ik) DF2(iL,iK,iJ,iI) = 0.0D+00
C       End Do
C       End Do
C       End Do
C       End Do
C       write(6,*) "asdf end"
      end if
C
      Do iI = 1, nAshT
        Do iJ = 1, nAshT ! iI
          Do iK = 1, nAshT ! iJ
            Do iL = 1, nAshT ! iK
C             Val1 = DG2(iI,iJ,iK,iL)
C             Val2 = DG2(iJ,iI,iL,iK)
C             Val3 = DG2(iK,iL,iI,iJ)
C             Val4 = DG2(iL,iK,iJ,iI)
C             Val  = (Val1+Val2+Val3+Val4)*0.25D+00
C             DG2(iI,iJ,iK,iL) = Val
C             DG2(iJ,iI,iL,iK) = Val
C             DG2(iK,iL,iI,iJ) = Val
C             DG2(iL,iK,iJ,iI) = Val
C             Val1 = DF2(iI,iJ,iK,iL)
C             Val2 = DF2(iJ,iI,iL,iK)
C             Val3 = DF2(iK,iL,iI,iJ)
C             Val4 = DF2(iL,iK,iJ,iI)
C             Val  = (Val1+Val2+Val3+Val4)*0.25D+00
C             DF2(iI,iJ,iK,iL) = Val
C             DF2(iJ,iI,iL,iK) = Val
C             DF2(iK,iL,iI,iJ) = Val
C             DF2(iL,iK,iJ,iI) = Val
C             Val1 = DF2(iI,iJ,iK,iL)
C             Val2 = DF2(iL,iK,iJ,iI)
C             Val  = (Val1+Val2)*0.5D+00
C             DF2(iI,iJ,iK,iL) = Val
C             DF2(iL,iK,iJ,iI) = Val
            End DO
          End Do
        End Do
      End Do
C
      Return
C
      End Subroutine CLagSym
C
C-----------------------------------------------------------------------
C
      Subroutine CLagDXA_FG3(iSym,nAS,NG3,BDER,SDER,
     *                       DF1,DF2,DF3,DG1,DG2,DG3,DEPSA,
     *                       G2,SC,idxG3)
C
      USE SUPERINDEX
      use EQSOLV
      use caspt2_module
C
      Implicit Real*8 (A-H,O-Z)
C
C
      Dimension BDER(nAS,nAS),SDER(nAS,nAS),DF3(*),DG3(*)
      Dimension DF1(nAshT,nAshT),DF2(nAshT,nAshT,nAshT,nAshT),
     *          DG1(nAshT,nAshT),DG2(nAshT,nAshT,nAshT,nAshT),
     *          DEPSA(nAshT,nAshT)
      Dimension G2(nAshT,nAshT,nAshT,nAshT)
      DIMENSION SC(*)
      INTEGER*1 idxG3(6,NG3)
C
      DO iG3=1,NG3
        iT=idxG3(1,iG3)
        iU=idxG3(2,iG3)
        iV=idxG3(3,iG3)
        iX=idxG3(4,iG3)
        iY=idxG3(5,iG3)
        iZ=idxG3(6,iG3)
        iST=IASYM(iT)
        iSU=IASYM(iU)
        iSV=IASYM(iV)
        iSX=IASYM(iX)
        iSY=IASYM(iY)
        iSZ=IASYM(iZ)
        ituvs=MUL(IST,MUL(ISU,ISV))
        ixyzs=MUL(ISX,MUL(ISY,ISZ))
        F3VAL=0.0D+00
        G3VAL=0.0D+00
        if(ituvs.ne.ixyzs) goto 500
        iTU=iT+NASHT*(iU-1)
        iVX=iV+NASHT*(iX-1)
        iYZ=iY+NASHT*(iZ-1)
C-SVC20100829: 12 equivalent cases, of which the second
C  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
C  - G(tuvxyz) -> SA(xut,vyz)
        jSYM=MUL(IASYM(iX),MUL(IASYM(iU),IASYM(iT)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iX,iU,iT)-nTUVES(jSYM)
          JSUP=KTUV(iV,iY,iZ)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
        if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 300
        if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 200
C  - G(vxtuyz) -> SA(uxv,tyz)
        jSYM=MUL(IASYM(iU),MUL(IASYM(iX),IASYM(iV)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iU,iX,iV)-nTUVES(jSYM)
          JSUP=KTUV(iT,iY,iZ)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
C  - G(yzvxtu) -> SA(xzy,vtu)
        jSYM=MUL(IASYM(iX),MUL(IASYM(iZ),IASYM(iY)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iX,iZ,iY)-nTUVES(jSYM)
          JSUP=KTUV(iV,iT,iU)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
C  - G(tuyzvx) -> SA(zut,yvx)
        jSYM=MUL(IASYM(iZ),MUL(IASYM(iU),IASYM(iT)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iZ,iU,iT)-nTUVES(jSYM)
          JSUP=KTUV(iY,iV,iX)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
 200   CONTINUE
C  - G(yztuvx) -> SA(uzy,tvx)
        jSYM=MUL(IASYM(iU),MUL(IASYM(iZ),IASYM(iY)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iU,iZ,iY)-nTUVES(jSYM)
          JSUP=KTUV(iT,iV,iX)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
C  - G(vxyztu) -> SA(zxv,ytu)
        jSYM=MUL(IASYM(iZ),MUL(IASYM(iX),IASYM(iV)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iZ,iX,iV)-nTUVES(jSYM)
          JSUP=KTUV(iY,iT,iU)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
 300   CONTINUE
        if (iT.eq.iU.and.iV.eq.iX.and.iY.eq.iZ) go to 500
        if (iT.eq.iU.and.iV.eq.iZ.and.iX.eq.iY) go to 500
        if (iX.eq.iV.and.iT.eq.iZ.and.iU.eq.iY) go to 500
        if (iZ.eq.iY.and.iV.eq.iU.and.iX.eq.iT) go to 500
C  - G(utxvzy) -> SA(vtu,xzy)
        jSYM=MUL(IASYM(iV),MUL(IASYM(iT),IASYM(iU)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iV,iT,iU)-nTUVES(jSYM)
          JSUP=KTUV(iX,iZ,iY)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
        if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 500
        if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 400
C  - G(xvutzy) -> SA(tvx,uzy)
        jSYM=MUL(IASYM(iT),MUL(IASYM(iV),IASYM(iX)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iT,iV,iX)-nTUVES(jSYM)
          JSUP=KTUV(iU,iZ,iY)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
C  - G(zyxvut) -> SA(vyz,xut)
        jSYM=MUL(IASYM(iV),MUL(IASYM(iY),IASYM(iZ)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iV,iY,iZ)-nTUVES(jSYM)
          JSUP=KTUV(iX,iU,iT)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
C  - G(utzyxv) -> SA(ytu,zxv)
        jSYM=MUL(IASYM(iY),MUL(IASYM(iT),IASYM(iU)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iY,iT,iU)-nTUVES(jSYM)
          JSUP=KTUV(iZ,iX,iV)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
 400   CONTINUE
C  - G(zyutxv) -> SA(tyz,uxv)
        jSYM=MUL(IASYM(iT),MUL(IASYM(iY),IASYM(iZ)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iT,iY,iZ)-nTUVES(jSYM)
          JSUP=KTUV(iU,iX,iV)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
C  - G(xvzyut) -> SA(yvx,zut)
        jSYM=MUL(IASYM(iY),MUL(IASYM(iV),IASYM(iX)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iY,iV,iX)-nTUVES(jSYM)
          JSUP=KTUV(iZ,iU,iT)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
 500   CONTINUE
C
        F3VAL = -F3VAL
        G3VAL = -G3VAL
C
        !! last line of F3 transformation in mkfg3.f
        G3VAL = G3VAL - (EPSA(iU)+EPSA(iY))*F3VAL
        Do iW = 1, nAshT
          ISUP=KTUV(iV,iW,iT)-nTUVES(iSYM)
          JSUP=KTUV(iX,iY,iZ)-nTUVES(iSYM)
          NSEQ=MAX(iSup,jSup)*(MAX(iSup,jSup)-1)/2 + MIN(iSup,jSup)
          DEPSA(iW,iU) = DEPSA(iW,iU) - F3VAL*SC(NSEQ)
C
          ISUP=KTUV(iV,iU,iT)-nTUVES(iSYM)
          JSUP=KTUV(iX,iW,iZ)-nTUVES(iSYM)
          NSEQ=MAX(iSup,jSup)*(MAX(iSup,jSup)-1)/2 + MIN(iSup,jSup)
          DEPSA(iW,iY) = DEPSA(iW,iY) - F3VAL*SC(NSEQ)
        End Do
C
        !! derivative of <0|EtuEwv,xwEyz|0>*fww
        DF3(iG3) = DF3(iG3) + F3VAL
        !! derivative of <0|EtuEvxEyz|0>
        DG3(iG3) = DG3(iG3) + G3VAL
C
        !! remaining F3 and G3 transformation in mkfg3.f
        If (iY.eq.iX) Then
          DF2(iT,iU,iV,iZ) = DF2(iT,iU,iV,iZ) - F3VAL
          DG2(iT,iU,iV,iZ) = DG2(iT,iU,iV,iZ) - EPSA(iU)*F3VAL
          Do iW = 1, nAshT
            DEPSA(iU,iW) = DEPSA(iU,iW) - F3VAL*G2(iT,iW,iV,iZ)
          End Do
          DG2(iT,iU,iV,iZ) = DG2(iT,iU,iV,iZ) - G3VAL
        End If
        If (iV.eq.iU) Then
          DF2(iT,iX,iY,iZ) = DF2(iT,iX,iY,iZ) - F3VAL
          DG2(iT,iX,iY,iZ) = DG2(iT,iX,iY,iZ) - EPSA(iY)*F3VAL
          Do iW = 1, nAshT
            DEPSA(iW,iY) = DEPSA(iW,iY) - F3VAL*G2(iT,iX,iW,iZ)
          End Do
          DG2(iT,iX,iY,iZ) = DG2(iT,iX,iY,iZ) - G3VAL
        End If
        If (iY.eq.iU) Then
          DF2(iV,iX,iT,iZ) = DF2(iV,iX,iT,iZ) - F3VAL
          DG2(iV,iX,iT,iZ) = DG2(iV,iX,iT,iZ) - EPSA(iU)*F3VAL
          DG2(iV,iX,iT,iZ) = DG2(iV,iX,iT,iZ) - G3VAL
        End If
        DEPSA(iY,iU) = DEPSA(iY,iU) - F3VAL*G2(iV,iX,iT,iZ)
        If (iY.eq.iX.and.iV.eq.iU) Then
          DF1(iT,iZ) = DF1(iT,iZ) - F3VAL
          DG1(iT,iZ) = DG1(iT,iZ) - G3VAL
        End If
      END DO
C
      Return
C
      End Subroutine CLagDXA_FG3
C
C-----------------------------------------------------------------------
C
#ifdef _MOLCAS_MPP_
      SUBROUTINE CLagDXA_FG3_MPP(ISYM,lg_BDER,lg_SDER,DG1,DG2,DG3,
     *                           DF1,DF2,DF3,DEPSA,G2,idxG3)
C
      USE MPI
      USE SUPERINDEX
      USE EQSOLV
      use definitions, only: iwp,RtoB,wp
      use stdalloc, only: mma_allocate, mma_deallocate, mma_MaxDBLE
      USE Para_Info, ONLY: Is_Real_Par, nProcs
      use caspt2_module
      use pt2_guga
      IMPLICIT REAL*8 (A-H,O-Z)

#include "global.fh"
#include "mafdecls.fh"

      DIMENSION DG1(NASHT,NASHT),DG2(NASHT,NASHT,NASHT,NASHT),DG3(*),
     *          DF1(NASHT,NASHT),DF2(NASHT,NASHT,NASHT,NASHT),DF3(*),
     *          DEPSA(NASHT,NASHT)
      DIMENSION G2(NASHT,NASHT,NASHT,NASHT)
      INTEGER*1 idxG3(6,NG3)

      integer(kind=iwp), ALLOCATABLE :: INDI(:), INDJ(:), NELBsav(:),
     *                                  NELSsav(:)
      real(kind=wp),  ALLOCATABLE :: BUFFB(:), BUFFS(:)

      ! Since we are stuck with collective calls to MPI_Alltoallv in
      ! order to gather the elements, each process needs to loop over
      ! the same number of blocks.
      NG3MAX=NG3
      CALL GAIGOP_SCAL(NG3MAX,'max')
      IF (NG3MAX.EQ.0) RETURN

      ! The global SC matrix has already been allocated, so we need to
      ! find out how much memory is left for buffering (4 equally sized
      ! buffers for sending and receiving values and indices)
      CALL mma_MaxDBLE(MAXMEM)
      ! we need two real and two integer values per element
      iscal = (iwp*2 + wp*2)/RtoB
      !MAXBUF=MIN(NINT(0.95D0*MAXMEM)/4,2000000000/8)
      MAXBUF=MIN(NINT(0.95D0*MAXMEM)/iscal,2000000000/8)
      MAXBUF=MAXBUF-2*NG3 !! for NELBsav and NELSsav

      ! Loop over blocks NG3B of NG3, so that 12*NG3B < MAXBUF/NPROCS.
      ! This guarantees that e.g. if all processes send all their data
      ! to one other, that process receives NPROCS*NG3B*12 elements
      ! in the receive buffer.
      NG3B=MAXBUF/(NPROCS*16) !! originally 12, but not sure why
      NG3B=MIN(NG3B,NG3MAX)
      CALL GAIGOP_SCAL(NG3B,'min')
      NBUF=16*NG3B

      call mma_allocate(BUFFB,NBUF,Label='BUFFB')
      call mma_allocate(BUFFS,NBUF,Label='BUFFS')
      call mma_allocate(INDI,NBUF,Label='INDI')
      call mma_allocate(INDJ,NBUF,Label='INDJ')

      NBLOCKS=(NG3MAX-1)/NG3B+1
      DO IBLOCK=1,NBLOCKS
        IG3STA=1+(IBLOCK-1)*NG3B
        IG3END=MIN(IG3STA+NG3B-1,NG3)
        call mma_allocate(NELBsav,IG3END-IG3STA+1,Label='NELBsav')
        call mma_allocate(NELSsav,IG3END-IG3STA+1,Label='NELSsav')

        ! Second pass fills the buffers with values and indices
        NtotELB = 0 ! Number of total elements of BDER
        NtotELS = 0 ! Number of total elements of SDER
        DO iG3=IG3STA,IG3END
          NELB = 0 ! Number of elements of BDER for iG3
          NELS = 0 ! Number of elements of SDER for iG3
          iT=idxG3(1,iG3)
          iU=idxG3(2,iG3)
          iV=idxG3(3,iG3)
          iX=idxG3(4,iG3)
          iY=idxG3(5,iG3)
          iZ=idxG3(6,iG3)
          iST=IASYM(iT)
          iSU=IASYM(iU)
          iSV=IASYM(iV)
          iSX=IASYM(iX)
          iSY=IASYM(iY)
          iSZ=IASYM(iZ)
          ituvs=MUL(IST,MUL(ISU,ISV))
          ixyzs=MUL(ISX,MUL(ISY,ISZ))
          if(ituvs.ne.ixyzs) goto 501
          iTU=iT+NASHT*(iU-1)
          iVX=iV+NASHT*(iX-1)
          iYZ=iY+NASHT*(iZ-1)
         !F3VAL=F3(iG3)
C-SVC20100829: 12 equivalent cases, of which the second
C  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
C  - G(tuvxyz) -> SA(xut,vyz)
        jSYM=MUL(IASYM(iX),MUL(IASYM(iU),IASYM(iT)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iU,iT)-nTUVES(jSYM)
            ICOL=KTUV(iV,iY,iZ)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
          if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 301
          if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 201
C  - G(vxtuyz) -> SA(uxv,tyz)
        jSYM=MUL(IASYM(iU),MUL(IASYM(iX),IASYM(iV)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iX,iV)-nTUVES(jSYM)
            ICOL=KTUV(iT,iY,iZ)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
C  - G(yzvxtu) -> SA(xzy,vtu)
        jSYM=MUL(IASYM(iX),MUL(IASYM(iZ),IASYM(iY)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iZ,iY)-nTUVES(jSYM)
            ICOL=KTUV(iV,iT,iU)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
C  - G(tuyzvx) -> SA(zut,yvx)
        jSYM=MUL(IASYM(iZ),MUL(IASYM(iU),IASYM(iT)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iU,iT)-nTUVES(jSYM)
            ICOL=KTUV(iY,iV,iX)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
 201   CONTINUE
C  - G(yztuvx) -> SA(uzy,tvx)
        jSYM=MUL(IASYM(iU),MUL(IASYM(iZ),IASYM(iY)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iZ,iY)-nTUVES(jSYM)
            ICOL=KTUV(iT,iV,iX)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
C  - G(vxyztu) -> SA(zxv,ytu)
        jSYM=MUL(IASYM(iZ),MUL(IASYM(iX),IASYM(iV)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iX,iV)-nTUVES(jSYM)
            ICOL=KTUV(iY,iT,iU)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
 301   CONTINUE
          if (iT.eq.iU.and.iV.eq.iX.and.iY.eq.iZ) go to 501
          if (iT.eq.iU.and.iV.eq.iZ.and.iX.eq.iY) go to 501
          if (iX.eq.iV.and.iT.eq.iZ.and.iU.eq.iY) go to 501
          if (iZ.eq.iY.and.iV.eq.iU.and.iX.eq.iT) go to 501
C  - G(utxvzy) -> SA(vtu,xzy)
        jSYM=MUL(IASYM(iV),MUL(IASYM(iT),IASYM(iU)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iT,iU)-nTUVES(jSYM)
            ICOL=KTUV(iX,iZ,iY)-nTUVES(jSYM)
            NELB  = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
          if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 501
          if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 401
C  - G(xvutzy) -> SA(tvx,uzy)
        jSYM=MUL(IASYM(iT),MUL(IASYM(iV),IASYM(iX)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iV,iX)-nTUVES(jSYM)
            ICOL=KTUV(iU,iZ,iY)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
C  - G(zyxvut) -> SA(vyz,xut)
        jSYM=MUL(IASYM(iV),MUL(IASYM(iY),IASYM(iZ)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iY,iZ)-nTUVES(jSYM)
            ICOL=KTUV(iX,iU,iT)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
C  - G(utzyxv) -> SA(ytu,zxv)
        jSYM=MUL(IASYM(iY),MUL(IASYM(iT),IASYM(iU)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iY,iT,iU)-nTUVES(jSYM)
            ICOL=KTUV(iZ,iX,iV)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
 401   CONTINUE
C  - G(zyutxv) -> SA(tyz,uxv)
        jSYM=MUL(IASYM(iT),MUL(IASYM(iY),IASYM(iZ)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iY,iZ)-nTUVES(jSYM)
            ICOL=KTUV(iU,iX,iV)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
C  - G(xvzyut) -> SA(yvx,zut)
        jSYM=MUL(IASYM(iY),MUL(IASYM(iV),IASYM(iX)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iY,iV,iX)-nTUVES(jSYM)
            ICOL=KTUV(iZ,iU,iT)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
 501   CONTINUE
       nels = nelb
          NELBsav(iG3-IG3STA+1) = NELB
          NELSsav(iG3-IG3STA+1) = NELS
        END DO
        ntotels = ntotelb

        CALL GA_GATHER(lg_BDER,BUFFB,INDI,INDJ,NtotELB)
        CALL GA_GATHER(lg_SDER,BUFFS,INDI,INDJ,NtotELS)

        ! Finally, fill the local chunk of the SC matrix (block of rows)
        ! with the received values at their appropriate place.
        NtotELB = 0
        NtotELS = 0
        DO iG3=IG3STA,IG3END
          iT=idxG3(1,iG3)
          iU=idxG3(2,iG3)
          iV=idxG3(3,iG3)
          iX=idxG3(4,iG3)
          iY=idxG3(5,iG3)
          iZ=idxG3(6,iG3)
          NELB = NELBsav(iG3-IG3STA+1)
          NELS = NELSsav(iG3-IG3STA+1)

          F3VAL = 0.0d+00
          DO i = 1, NELB
            NtotELB = NtotELB + 1
            F3VAL = F3VAL - BUFFB(NtotELB)
          END DO

          G3VAL = 0.0d+00
          DO i = 1, NELS
            NtotELS = NtotELS + 1
            G3VAL = G3VAL - BUFFS(NtotELS)
          END DO
          G3VAL = G3VAL - (EPSA(iU)+EPSA(iY))*F3VAL

          DF3(iG3) = DF3(iG3) + F3VAL
          DG3(iG3) = DG3(iG3) + G3VAL

          !! DEPSA is done in DF3_DEPSA_MPP

          !! remaining F3 and G3 transformation in mkfg3.f
          If (iY.eq.iX) Then
            DF2(iT,iU,iV,iZ) = DF2(iT,iU,iV,iZ) - F3VAL
            DG2(iT,iU,iV,iZ) = DG2(iT,iU,iV,iZ) - EPSA(iU)*F3VAL
            Do iW = 1, nAshT
              DEPSA(iU,iW) = DEPSA(iU,iW) - F3VAL*G2(iT,iW,iV,iZ)
            End Do
            DG2(iT,iU,iV,iZ) = DG2(iT,iU,iV,iZ) - G3VAL
          End If
          If (iV.eq.iU) Then
            DF2(iT,iX,iY,iZ) = DF2(iT,iX,iY,iZ) - F3VAL
            DG2(iT,iX,iY,iZ) = DG2(iT,iX,iY,iZ) - EPSA(iY)*F3VAL
            Do iW = 1, nAshT
              DEPSA(iW,iY) = DEPSA(iW,iY) - F3VAL*G2(iT,iX,iW,iZ)
            End Do
            DG2(iT,iX,iY,iZ) = DG2(iT,iX,iY,iZ) - G3VAL
          End If
          If (iY.eq.iU) Then
            DF2(iV,iX,iT,iZ) = DF2(iV,iX,iT,iZ) - F3VAL
            DG2(iV,iX,iT,iZ) = DG2(iV,iX,iT,iZ) - EPSA(iU)*F3VAL
            DG2(iV,iX,iT,iZ) = DG2(iV,iX,iT,iZ) - G3VAL
          End If
          DEPSA(iY,iU) = DEPSA(iY,iU) - F3VAL*G2(iV,iX,iT,iZ)
          If (iY.eq.iX.and.iV.eq.iU) Then
            DF1(iT,iZ) = DF1(iT,iZ) - F3VAL
            DG1(iT,iZ) = DG1(iT,iZ) - G3VAL
          End If
        END DO

        call mma_deallocate(NELBsav)
        call mma_deallocate(NELSsav)
      END DO ! end loop over blocks of G3 values

      call mma_deallocate(BUFFB)
      call mma_deallocate(BUFFS)
      call mma_deallocate(INDI)
      call mma_deallocate(INDJ)

      RETURN
C
      End Subroutine CLagDXA_FG3_MPP
#endif
C
C-----------------------------------------------------------------------
C
      Subroutine CLagDXA_DP(iSym,nAS,BDER,SDER,DG1,DG2,DF1,DF2,
     *                      DEPSA,DEASUM,iLo,iHi,jLo,jHi,LDA,
     *                      G1,G2,SA,SA2,lg_S)
C
      USE SUPERINDEX
      use caspt2_global, only:ipea_shift
      use EQSOLV
      use caspt2_module
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, nProcs
#endif
C
      Implicit Real*8 (A-H,O-Z)
C
#ifdef _MOLCAS_MPP_
#include "global.fh"
#endif
C
      Dimension BDER(*),SDER(*),
     *          DF2(nAshT,nAshT,nAshT,nAshT),
     *          DG2(nAshT,nAshT,nAshT,nAshT),
     *          DF1(nAshT,nAshT),DG1(nAshT,nAshT),DEPSA(nAshT,nAshT)
      dimension g1(nAshT,nAshT),G2(nAshT,nAshT,nasht,nasht),SA(*),SA2(*)
C     INTEGER*1 idxG3(6,NG3)
C
C     LDA == 0, if not parallel; SA is triangular
C     LDA /= 0, if parallel    ; SA is square
C     In both cases, BDER and SDER are square
C
      ISADR=0
      if (isadr.ne.0) write (6,*) lda !! just for avoid compiling error
      if (isadr.ne.0) write (6,*) lg_S !! just for avoid compiling error
      NROW = 0
      iLoS = 0
      jLoS = 0
#ifdef _MOLCAS_MPP_
      if (is_real_par()) then
        irank = 0
        CALL GA_Distribution (lg_S,iRank,iLoS,iHiS,jLoS,jHiS)
        NROW = jHiS-jLoS+1 !! = NAS
        CALL GA_GET(lg_S,jLoS,jHiS,iLoS,iHiS,SA2,NROW)
      end if
#endif
      DO 100 IXYZ=jLo,jHi
        IXYZABS=IXYZ+NTUVES(ISYM)
        IXABS=MTUV(1,IXYZABS)
        IYABS=MTUV(2,IXYZABS)
        IZABS=MTUV(3,IXYZABS)
        EX=EPSA(IXABS)
        EY=EPSA(IYABS)
        DO 101 ITUV=iLo,iHi
          ITUVABS=ITUV+NTUVES(ISYM)
          ITABS=MTUV(1,ITUVABS)
          IUABS=MTUV(2,ITUVABS)
          IVABS=MTUV(3,ITUVABS)
          ET=EPSA(ITABS)
          EU=EPSA(IUABS)
          ETU=ET+EU
          FACT=EY+EU+EX+ET-EASUM
          IF (LDA.EQ.0) THEN
            ISADR=1+iTUV-iLo+NAS*(iXYZ-jLo)
          ELSE
            ISADR=1+iTUV-iLo+LDA*(iXYZ-jLo)
          END IF
          ValB=BDER(ISADR)
C
          If (iTUV.eq.iXYZ.and.ipea_shift.ne.0.0D+00) Then
C           !! BA in the next equation refers to the active overlap
C       ipea_shift*0.5d0*BA(ISADR)*(2.0d0-DREF(IDV)+DREF(IDT)+DREF(IDU))
            bsBDER = ipea_shift*0.5D+00*ValB
            SDER(iSAdr) = SDER(iSAdr) + bsBDER*(2.0D+00
     *        +G1(iTabs,iTabs)+G1(iUabs,iUabs)-G1(iVabs,iVabs))
            IF (LDA.EQ.0) THEN
              iSAdr2 = iTUV*(iTUV+1)/2
            ELSE
              ISADR2 = 1+iTUV-iLo+LDA*(iTUV-jLo)
            END IF
            DG1(iTabs,iTabs) = DG1(iTabs,iTabs) + bsBDER*SA(iSAdr2)
            DG1(iUabs,iUabs) = DG1(iUabs,iUabs) + bsBDER*SA(iSAdr2)
            DG1(iVabs,iVabs) = DG1(iVabs,iVabs) - bsBDER*SA(iSAdr2)
          End If
C
          !! First VALUE contribution in MKBC_DP (FACT)
          SDER(ISADR) = SDER(ISADR) + FACT*ValB
          ValS=SDER(ISADR)
C
          if (lda.eq.0) then
            Do iWabs = 1, nAshT
              !! EU derivative
              iTWV = iTabs+nAshT*(iWabs-1)+nAshT**2*(iVabs-1)
              iSAdr2 = Max(iTWV,iXYZ)*(Max(iTWV,iXYZ)-1)/2
     *               + Min(iTWV,iXYZ)
              DEPSA(iWabs,iUabs) = DEPSA(iWabs,iUabs)
     *          + ValB*SA(iSAdr2)
C
              !! EY derivative
              iXWZ = iXabs+nAshT*(iWabs-1)+nAshT**2*(iZabs-1)
              iSAdr2 = Max(iTUV,iXWZ)*(Max(iTUV,iXWZ)-1)/2
     *               + Min(iTUV,iXWZ)
              DEPSA(iWabs,iYabs) = DEPSA(iWabs,iYabs)
     *          + ValB*SA(iSAdr2)
C
              !! EX derivative
              iWYZ = iWabs+nAshT*(iYabs-1)+nAshT**2*(iZabs-1)
              iSAdr2 = Max(iTUV,iWYZ)*(Max(iTUV,iWYZ)-1)/2
     *               + Min(iTUV,iWYZ)
              DEPSA(iWabs,iXabs) = DEPSA(iWabs,iXabs)
     *          + ValB*SA(iSAdr2)
C
              !! ET derivative
              iWUV = iWabs+nAshT*(iUabs-1)+nAshT**2*(iVabs-1)
              iSAdr2 = Max(iWUV,iXYZ)*(Max(iWUV,iXYZ)-1)/2
     *               + Min(iWUV,iXYZ)
              DEPSA(iWabs,iTabs) = DEPSA(iWabs,iTabs)
     *          + ValB*SA(iSAdr2)
            End Do
          else
            Do iWabs = 1, nAshT
              !! EU derivative
              iTWV = iTabs+nAshT*(iWabs-1)+nAshT**2*(iVabs-1)
              ISADR2 = 1+iTWV-jLoS+NROW*(iXYZ-iLoS)
              DEPSA(iWabs,iUabs) = DEPSA(iWabs,iUabs)
     *          + ValB*SA2(iSAdr2)
C
              !! EY derivative
              iXWZ = iXabs+nAshT*(iWabs-1)+nAshT**2*(iZabs-1)
              ISADR2 = 1+iTUV-iLo+LDA*(iXWZ-jLo)
              DEPSA(iWabs,iYabs) = DEPSA(iWabs,iYabs)
     *          + ValB*SA(iSAdr2)
C
              !! EX derivative
              iWYZ = iWabs+nAshT*(iYabs-1)+nAshT**2*(iZabs-1)
              ISADR2 = 1+iTUV-iLo+LDA*(iWYZ-jLo)
              DEPSA(iWabs,iXabs) = DEPSA(iWabs,iXabs)
     *          + ValB*SA(iSAdr2)
C
              !! ET derivative
              iWUV = iWabs+nAshT*(iUabs-1)+nAshT**2*(iVabs-1)
              ISADR2 = 1+iWUV-jLoS+NROW*(iXYZ-iLoS)
              DEPSA(iWabs,iTabs) = DEPSA(iWabs,iTabs)
     *          + ValB*SA2(iSAdr2)
            End Do
          end if

          IF (LDA.EQ.0) THEN
            iSAdr = Max(iTUV,iXYZ)*(Max(iTUV,iXYZ)-1)/2
     *            + Min(iTUV,iXYZ)
          END IF
          DEASUM = DEASUM - ValB*SA(iSAdr)
C
C         2dtx ( Fvuyz-Et*Gvuyz )
C         2 dtx Gvuyz + 2 dtx dyu Gvz
          If (iTabs.eq.iXabs) Then
            !! VALUE=VALUE+4.0D0*(FP(IP)-ET*PREF(IP))
            DF2(iVabs,iUabs,iYabs,iZabs)
     *        = DF2(iVabs,iUabs,iYabs,iZabs) + 2.0D+00*ValB
            DG2(iVabs,iUabs,iYabs,iZabs)
     *        = DG2(iVabs,iUabs,iYabs,iZabs) - 2.0D+00*ET*ValB
C
            !! VALUE=VALUE+4.0D0*PREF(IP)
            DG2(iVabs,iUabs,iYabs,iZabs)
     *        = DG2(iVabs,iUabs,iYabs,iZabs) + 2.0D+00*ValS
            If (iYabs.eq.iUabs) Then
              !! VALUE=VALUE+2.0D0*DREF(ID)
              DG1(iVabs,iZabs) = DG1(iVabs,iZabs) + 2.0D+00*ValS
            End If
          End If
          DEPSA(iTabs,iXabs) = DEPSA(iTabs,iXabs)
     *      -2.0D+00*ValB*G2(iVabs,iUabs,iYabs,iZabs)
C
C         dxu ( -Fvtyz + Eu*Gvtyz )
C         -dxu Gvtyz -dxu dyt Gvz
          If (iXabs.eq.iUabs) Then
            !! VALUE=VALUE-2.0D0*(FP(IP)-EU*PREF(IP))
            DF2(iVabs,iTabs,iYabs,iZabs)
     *        = DF2(iVabs,iTabs,iYabs,iZabs) - ValB
            DG2(iVabs,iTabs,iYabs,iZabs)
     *        = DG2(iVabs,iTabs,iYabs,iZabs) + EU*ValB
            !! VALUE=VALUE - 2.0D0*PREF(IP)
            DG2(iVabs,iTabs,iYabs,iZabs)
     *        = DG2(iVabs,iTabs,iYabs,iZabs) - ValS
            If (iYabs.eq.iTabs) Then
              !! VALUE=VALUE - DREF(ID)
              DG1(iVabs,iZabs) = DG1(iVabs,iZabs) - ValS
            End If
          End If
          DEPSA(iXabs,iUabs) = DEPSA(iXabs,iUabs)
     *      + ValB*G2(iVabs,iTabs,iYabs,iZabs)
C
C         dyt ( -Fvuxz + Et*Gvuxz +dxu (-Fvz+(Et+Eu)*Gvz))
C         -dyt Gvuxz
          If (iYabs.eq.iTabs) Then
            !! VALUE=VALUE-2.0D0*(FP(IP)-ET*PREF(IP))
            DF2(iVabs,iUabs,iXabs,iZabs)
     *        = DF2(iVabs,iUabs,iXabs,iZabs) - ValB
            DG2(iVabs,iUabs,iXabs,iZabs)
     *        = DG2(iVabs,iUabs,iXabs,iZabs) + ET*ValB
            If (iXabs.eq.iUabs) Then
              !! VALUE=VALUE - (FD(ID)-ETU*DREF(ID))
              DF1(iVabs,iZabs) = DF1(iVabs,iZabs) - ValB
              DG1(iVabs,iZabs) = DG1(iVabs,iZabs) + ETU*ValB
            End If
C
            !! VALUE=VALUE - 2.0D0*PREF(IP)
            DG2(iVabs,iUabs,iXabs,iZabs)
     *        = DG2(iVabs,iUabs,iXabs,iZabs) - ValS
          End If
          DEPSA(iYabs,iTabs) = DEPSA(iYabs,iTabs)
     *      + ValB*G2(iVabs,iUabs,iXabs,iZabs)
          If (iYabs.eq.iTabs)
     *    DEPSA(iXabs,iUabs) = DEPSA(iXabs,iUabs) + ValB*G1(iVabs,iZabs)
          If (iXabs.eq.iUabs)
     *    DEPSA(iYabs,iTabs) = DEPSA(iYabs,iTabs) + ValB*G1(iVabs,iZabs)
C
C         -dyu Gvzxt
          If (iYabs.eq.iUabs) Then
            !! VALUE=VALUE-2.0D0*(FP(IP)-EU*PREF(IP))
            DF2(iVabs,iZabs,iXabs,iTabs)
     *        = DF2(iVabs,iZabs,iXabs,iTabs) - ValB
            DG2(iVabs,iZabs,iXabs,iTabs)
     *        = DG2(iVabs,iZabs,iXabs,iTabs) + EU*ValB
            If (iXabs.eq.iTabs) Then
              !! VALUE=VALUE+2.0D0*(FD(ID)-ETU*DREF(ID))
              DF1(iVabs,iZabs) = DF1(iVabs,iZabs) + 2.0D+00*ValB
              DG1(iVabs,iZabs) = DG1(iVabs,iZabs) - 2.0D+00*ETU*ValB
            End If
C
            !! VALUE=VALUE - 2.0D0*PREF(IP)
            DG2(iVabs,iZabs,iXabs,iTabs)
     *        = DG2(iVabs,iZabs,iXabs,iTabs) - ValS
          End If
          DEPSA(iYabs,iUabs) = DEPSA(iYabs,iUabs)
     *      + ValB*G2(iVabs,iZabs,iXabs,iTabs)
          If (iYabs.eq.iUabs) DEPSA(iXabs,iTabs) = DEPSA(iXabs,iTabs)
     *      - 2.0D+00*ValB*G1(iVabs,iZabs)
          If (iXabs.eq.iTabs) DEPSA(iYabs,iUabs) = DEPSA(iYabs,iUabs)
     *      - 2.0D+00*ValB*G1(iVabs,iZabs)
 101    CONTINUE
#ifdef _MOLCAS_MPP_
        if (is_real_par() .and. IXYZ==iHiS .and. iRank/=NPROCS-1) then
          irank = irank + 1
          CALL GA_Distribution (lg_S,iRank,iLoS,iHiS,jLoS,jHiS)
          CALL GA_GET(lg_S,jLoS,jHiS,iLoS,iHiS,SA2,NROW)
        end if
#endif
 100  CONTINUE
C
      Return
C
      End Subroutine CLagDXA_DP
C
C-----------------------------------------------------------------------
C
      Subroutine CLagDXC_FG3(iSym,nAS,NG3,BDER,SDER,
     *                       DF1,DF2,DF3,DG1,DG2,DG3,DEPSA,
     *                       G2,SC,idxG3)
C
      USE SUPERINDEX
      use EQSOLV
      use caspt2_module
C
      Implicit Real*8 (A-H,O-Z)
C
C
      Dimension BDER(nAS,nAS),SDER(nAS,nAS),DF3(*),DG3(*)
      Dimension DF1(nAshT,nAshT),DF2(nAshT,nAshT,nAshT,nAshT),
     *          DG1(nAshT,nAshT),DG2(nAshT,nAshT,nAshT,nAshT),
     *          DEPSA(nAshT,nAshT)
      Dimension G2(nAshT,nAshT,nAshT,nAshT)
      DIMENSION SC(*)
      INTEGER*1 idxG3(6,NG3)
C
      DO iG3=1,NG3
        iT=idxG3(1,iG3)
        iU=idxG3(2,iG3)
        iV=idxG3(3,iG3)
        iX=idxG3(4,iG3)
        iY=idxG3(5,iG3)
        iZ=idxG3(6,iG3)
        iST=IASYM(iT)
        iSU=IASYM(iU)
        iSV=IASYM(iV)
        iSX=IASYM(iX)
        iSY=IASYM(iY)
        iSZ=IASYM(iZ)
        ituvs=MUL(IST,MUL(ISU,ISV))
        ixyzs=MUL(ISX,MUL(ISY,ISZ))
        F3VAL=0.0D+00
        G3VAL=0.0D+00
        if(ituvs.ne.ixyzs) goto 500
        iTU=iT+NASHT*(iU-1)
        iVX=iV+NASHT*(iX-1)
        iYZ=iY+NASHT*(iZ-1)
C-SVC20100829: 12 equivalent cases, of which the second
C  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
C  - G(tuvxyz) -> SC(vut,xyz)
        jSYM=MUL(IASYM(iV),MUL(IASYM(iU),IASYM(iT)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iV,iU,iT)-nTUVES(jSYM)
          JSUP=KTUV(iX,iY,iZ)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
        if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 300
        if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 200
C  - G(vxtuyz) -> SC(txv,uyz)
        jSYM=MUL(IASYM(iT),MUL(IASYM(iX),IASYM(iV)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iT,iX,iV)-nTUVES(jSYM)
          JSUP=KTUV(iU,iY,iZ)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
C  - G(yzvxtu) -> SC(vzy,xtu)
        jSYM=MUL(IASYM(iV),MUL(IASYM(iZ),IASYM(iY)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iV,iZ,iY)-nTUVES(jSYM)
          JSUP=KTUV(iX,iT,iU)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
C  - G(tuyzvx) -> SC(yut,zvx)
        jSYM=MUL(IASYM(iY),MUL(IASYM(iU),IASYM(iT)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iY,iU,iT)-nTUVES(jSYM)
          JSUP=KTUV(iZ,iV,iX)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
 200   CONTINUE
C  - G(yztuvx) -> SC(tzy,uvx)
        jSYM=MUL(IASYM(iT),MUL(IASYM(iZ),IASYM(iY)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iT,iZ,iY)-nTUVES(jSYM)
          JSUP=KTUV(iU,iV,iX)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
C  - G(vxyztu) -> SC(yxv,ztu)
        jSYM=MUL(IASYM(iY),MUL(IASYM(iX),IASYM(iV)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iY,iX,iV)-nTUVES(jSYM)
          JSUP=KTUV(iZ,iT,iU)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
 300   CONTINUE
        if (iT.eq.iU.and.iV.eq.iX.and.iY.eq.iZ) go to 500
        if (iT.eq.iU.and.iV.eq.iZ.and.iX.eq.iY) go to 500
        if (iX.eq.iV.and.iT.eq.iZ.and.iU.eq.iY) go to 500
        if (iZ.eq.iY.and.iV.eq.iU.and.iX.eq.iT) go to 500
C  - G(utxvzy) -> SC(xtu,vzy)
        jSYM=MUL(IASYM(iX),MUL(IASYM(iT),IASYM(iU)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iX,iT,iU)-nTUVES(jSYM)
          JSUP=KTUV(iV,iZ,iY)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
        if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 500
        if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 400
C  - G(xvutzy) -> SC(uvx,tzy)
        jSYM=MUL(IASYM(iU),MUL(IASYM(iV),IASYM(iX)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iU,iV,iX)-nTUVES(jSYM)
          JSUP=KTUV(iT,iZ,iY)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
C  - G(zyxvut) -> SC(xyz,vut)
        jSYM=MUL(IASYM(iX),MUL(IASYM(iY),IASYM(iZ)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iX,iY,iZ)-nTUVES(jSYM)
          JSUP=KTUV(iV,iU,iT)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
C  - G(utzyxv) -> SC(ztu,yxv)
        jSYM=MUL(IASYM(iZ),MUL(IASYM(iT),IASYM(iU)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iZ,iT,iU)-nTUVES(jSYM)
          JSUP=KTUV(iY,iX,iV)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
 400   CONTINUE
C  - G(zyutxv) -> SC(uyz,txv)
        jSYM=MUL(IASYM(iU),MUL(IASYM(iY),IASYM(iZ)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iU,iY,iZ)-nTUVES(jSYM)
          JSUP=KTUV(iT,iX,iV)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
C  - G(xvzyut) -> SC(zvx,yut)
        jSYM=MUL(IASYM(iZ),MUL(IASYM(iV),IASYM(iX)))
        IF (jSYM.EQ.iSYM) THEN
          ISUP=KTUV(iZ,iV,iX)-nTUVES(jSYM)
          JSUP=KTUV(iY,iU,iT)-nTUVES(jSYM)
          F3VAL = F3VAL + BDER(iSup,jSup)
          G3VAL = G3VAL + SDER(iSup,jSup)
        ENDIF
 500   CONTINUE
C
        !! last line of F3 transformation in mkfg3.f
C     g3val=0.0d+00 ! asdf
        G3VAL = G3VAL - (EPSA(iU)+EPSA(iY))*F3VAL
        Do iW = 1, nAshT
          ISUP=KTUV(iV,iW,iT)-nTUVES(iSYM)
          JSUP=KTUV(iX,iY,iZ)-nTUVES(iSYM)
          NSEQ=MAX(iSup,jSup)*(MAX(iSup,jSup)-1)/2 + MIN(iSup,jSup)
          DEPSA(iW,iU) = DEPSA(iW,iU) - F3VAL*SC(NSEQ)
C
          ISUP=KTUV(iV,iU,iT)-nTUVES(iSYM)
          JSUP=KTUV(iX,iW,iZ)-nTUVES(iSYM)
          NSEQ=MAX(iSup,jSup)*(MAX(iSup,jSup)-1)/2 + MIN(iSup,jSup)
          DEPSA(iW,iY) = DEPSA(iW,iY) - F3VAL*SC(NSEQ)
        End Do
C
        !! derivative of <0|EtuEwv,xwEyz|0>*fww
        DF3(iG3) = DF3(iG3) + F3VAL
        !! derivative of <0|EtuEvxEyz|0>
        DG3(iG3) = DG3(iG3) + G3VAL
C
        !! remaining F3 and G3 transformation in mkfg3.f
        If (iY.eq.iX) Then
          DF2(iT,iU,iV,iZ) = DF2(iT,iU,iV,iZ) - F3VAL
          DG2(iT,iU,iV,iZ) = DG2(iT,iU,iV,iZ) - EPSA(iU)*F3VAL
          Do iW = 1, nAshT
            DEPSA(iU,iW) = DEPSA(iU,iW) - F3VAL*G2(iT,iW,iV,iZ)
          End Do
          DG2(iT,iU,iV,iZ) = DG2(iT,iU,iV,iZ) - G3VAL
        End If
        If (iV.eq.iU) Then
          DF2(iT,iX,iY,iZ) = DF2(iT,iX,iY,iZ) - F3VAL
          DG2(iT,iX,iY,iZ) = DG2(iT,iX,iY,iZ) - EPSA(iY)*F3VAL
          Do iW = 1, nAshT
            DEPSA(iW,iY) = DEPSA(iW,iY) - F3VAL*G2(iT,iX,iW,iZ)
          End Do
          DG2(iT,iX,iY,iZ) = DG2(iT,iX,iY,iZ) - G3VAL
        End If
        If (iY.eq.iU) Then
          DF2(iV,iX,iT,iZ) = DF2(iV,iX,iT,iZ) - F3VAL
          DG2(iV,iX,iT,iZ) = DG2(iV,iX,iT,iZ) - EPSA(iU)*F3VAL
          DG2(iV,iX,iT,iZ) = DG2(iV,iX,iT,iZ) - G3VAL
        End If
        DEPSA(iY,iU) = DEPSA(iY,iU) - F3VAL*G2(iV,iX,iT,iZ)
        If (iY.eq.iX.and.iV.eq.iU) Then
          DF1(iT,iZ) = DF1(iT,iZ) - F3VAL
          DG1(iT,iZ) = DG1(iT,iZ) - G3VAL
        End If
      END DO
C
      Return
C
      End Subroutine CLagDXC_FG3
C
C-----------------------------------------------------------------------
C
#ifdef _MOLCAS_MPP_
      SUBROUTINE CLagDXC_FG3_MPP(ISYM,lg_BDER,lg_SDER,DG1,DG2,DG3,
     *                           DF1,DF2,DF3,DEPSA,G2,idxG3)
      USE MPI
      USE SUPERINDEX
      USE EQSOLV
      USE Para_Info, only: nProcs
      use definitions, only: iwp,RtoB,wp
      use stdalloc, only: mma_allocate, mma_deallocate, mma_MaxDBLE
      use caspt2_module
      use pt2_guga

      IMPLICIT REAL*8 (A-H,O-Z)


#include "global.fh"
#include "mafdecls.fh"

      DIMENSION DG1(NASHT,NASHT),DG2(NASHT,NASHT,NASHT,NASHT),DG3(*),
     *          DF1(NASHT,NASHT),DF2(NASHT,NASHT,NASHT,NASHT),DF3(*),
     *          DEPSA(NASHT,NASHT)
      DIMENSION G2(NASHT,NASHT,NASHT,NASHT)
      INTEGER*1 idxG3(6,NG3)

      integer(kind=iwp), allocatable :: INDI(:), INDJ(:), NELBsav(:),
     *                                  NELSsav(:)
      real(kind=wp),  ALLOCATABLE :: BUFFB(:), BUFFS(:)

      ! Since we are stuck with collective calls to MPI_Alltoallv in
      ! order to gather the elements, each process needs to loop over
      ! the same number of blocks.
      NG3MAX=NG3
      CALL GAIGOP_SCAL(NG3MAX,'max')
      IF (NG3MAX.EQ.0) RETURN

      ! The global SC matrix has already been allocated, so we need to
      ! find out how much memory is left for buffering (4 equally sized
      ! buffers for sending and receiving values and indices)
      CALL mma_MaxDBLE(MAXMEM)
      ! we need two real and two integer values per element
      iscal = (iwp*2 + wp*2)/RtoB
      !MAXBUF=MIN(NINT(0.95D0*MAXMEM)/4,2000000000/8)
      MAXBUF=MIN(NINT(0.95D0*MAXMEM)/iscal,2000000000/8)
      MAXBUF=MAXBUF-2*NG3 !! for NELBsav and NELSsav

      ! Loop over blocks NG3B of NG3, so that 12*NG3B < MAXBUF/NPROCS.
      ! This guarantees that e.g. if all processes send all their data
      ! to one other, that process receives NPROCS*NG3B*12 elements
      ! in the receive buffer.
      NG3B=MAXBUF/(NPROCS*16) !! originally 12, but not sure why
      NG3B=MIN(NG3B,NG3MAX)
      CALL GAIGOP_SCAL(NG3B,'min')
      NBUF=16*NG3B

      call mma_allocate(BUFFB,NBUF,Label='BUFFB')
      call mma_allocate(BUFFS,NBUF,Label='BUFFS')
      call mma_allocate(INDI,NBUF,Label='INDI')
      call mma_allocate(INDJ,NBUF,Label='INDJ')

      NBLOCKS=(NG3MAX-1)/NG3B+1
      DO IBLOCK=1,NBLOCKS
        IG3STA=1+(IBLOCK-1)*NG3B
        IG3END=MIN(IG3STA+NG3B-1,NG3)
        call mma_allocate(NELBsav,IG3END-IG3STA+1,Label='NELBsav')
        call mma_allocate(NELSsav,IG3END-IG3STA+1,Label='NELSsav')

        ! Second pass fills the buffers with values and indices
        NtotELB = 0 ! Number of total elements of BDER
        NtotELS = 0 ! Number of total elements of SDER
        DO iG3=IG3STA,IG3END
          NELB = 0 ! Number of elements of BDER for iG3
          NELS = 0 ! Number of elements of SDER for iG3
          iT=idxG3(1,iG3)
          iU=idxG3(2,iG3)
          iV=idxG3(3,iG3)
          iX=idxG3(4,iG3)
          iY=idxG3(5,iG3)
          iZ=idxG3(6,iG3)
          iST=IASYM(iT)
          iSU=IASYM(iU)
          iSV=IASYM(iV)
          iSX=IASYM(iX)
          iSY=IASYM(iY)
          iSZ=IASYM(iZ)
          ituvs=MUL(IST,MUL(ISU,ISV))
          ixyzs=MUL(ISX,MUL(ISY,ISZ))
          if(ituvs.ne.ixyzs) goto 501
          iTU=iT+NASHT*(iU-1)
          iVX=iV+NASHT*(iX-1)
          iYZ=iY+NASHT*(iZ-1)
         !F3VAL=F3(iG3)
C-SVC20100829: 12 equivalent cases, of which the second
C  half reflects the S(tuv,xyz)=S(xyz,tuv) symmetry:
C  - F(tuvxyz) -> BC(vut,xyz)
          jSYM=MUL(IASYM(iV),MUL(IASYM(iU),IASYM(iT)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iU,iT)-nTUVES(jSYM)
            ICOL=KTUV(iX,iY,iZ)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
          if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 301
          if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 201
C  - F(vxtuyz) -> BC(txv,uyz)
          jSYM=MUL(IASYM(iT),MUL(IASYM(iX),IASYM(iV)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iX,iV)-nTUVES(jSYM)
            ICOL=KTUV(iU,iY,iZ)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
C  - F(yzvxtu) -> BC(vzy,xtu)
          jSYM=MUL(IASYM(iV),MUL(IASYM(iZ),IASYM(iY)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iV,iZ,iY)-nTUVES(jSYM)
            ICOL=KTUV(iX,iT,iU)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
C  - F(tuyzvx) -> BC(yut,zvx)
          jSYM=MUL(IASYM(iY),MUL(IASYM(iU),IASYM(iT)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iY,iU,iT)-nTUVES(jSYM)
            ICOL=KTUV(iZ,iV,iX)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
 201   CONTINUE
C  - F(yztuvx) -> BC(tzy,uvx)
          jSYM=MUL(IASYM(iT),MUL(IASYM(iZ),IASYM(iY)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iT,iZ,iY)-nTUVES(jSYM)
            ICOL=KTUV(iU,iV,iX)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
C  - F(vxyztu) -> BC(yxv,ztu)
          jSYM=MUL(IASYM(iY),MUL(IASYM(iX),IASYM(iV)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iY,iX,iV)-nTUVES(jSYM)
            ICOL=KTUV(iZ,iT,iU)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
 301   CONTINUE
          if (iT.eq.iU.and.iV.eq.iX.and.iY.eq.iZ) go to 501
          if (iT.eq.iU.and.iV.eq.iZ.and.iX.eq.iY) go to 501
          if (iX.eq.iV.and.iT.eq.iZ.and.iU.eq.iY) go to 501
          if (iZ.eq.iY.and.iV.eq.iU.and.iX.eq.iT) go to 501
C  - F(utxvzy) -> BC(xtu,vzy)
          jSYM=MUL(IASYM(iX),MUL(IASYM(iT),IASYM(iU)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iT,iU)-nTUVES(jSYM)
            ICOL=KTUV(iV,iZ,iY)-nTUVES(jSYM)
            NELB  = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
          if (iTU.eq.iVX.and.iVX.eq.iYZ) go to 501
          if (iTU.eq.iVX.or.iTU.eq.iYZ.or.iVX.eq.iYZ) go to 401
C  - F(xvutzy) -> BC(uvx,tzy)
          jSYM=MUL(IASYM(iU),MUL(IASYM(iV),IASYM(iX)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iV,iX)-nTUVES(jSYM)
            ICOL=KTUV(iT,iZ,iY)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
C  - F(zyxvut) -> BC(xyz,vut)
          jSYM=MUL(IASYM(iX),MUL(IASYM(iY),IASYM(iZ)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iX,iY,iZ)-nTUVES(jSYM)
            ICOL=KTUV(iV,iU,iT)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
C  - F(utzyxv) -> BC(ztu,yxv)
          jSYM=MUL(IASYM(iZ),MUL(IASYM(iT),IASYM(iU)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iT,iU)-nTUVES(jSYM)
            ICOL=KTUV(iY,iX,iV)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
 401   CONTINUE
C  - F(zyutxv) -> BC(uyz,txv)
          jSYM=MUL(IASYM(iU),MUL(IASYM(iY),IASYM(iZ)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iU,iY,iZ)-nTUVES(jSYM)
            ICOL=KTUV(iT,iX,iV)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
C  - F(xvzyut) -> BC(zvx,yut)
          jSYM=MUL(IASYM(iZ),MUL(IASYM(iV),IASYM(iX)))
          IF (jSYM.EQ.iSYM) THEN
            IROW=KTUV(iZ,iV,iX)-nTUVES(jSYM)
            ICOL=KTUV(iY,iU,iT)-nTUVES(jSYM)
            NELB = NELB + 1
            NtotELB = NtotELB + 1
            INDI(NtotELB) = IROW
            INDJ(NtotELB) = ICOL
          ENDIF
 501   CONTINUE
       nels = nelb
          NELBsav(iG3-IG3STA+1) = NELB
          NELSsav(iG3-IG3STA+1) = NELS
        END DO
        ntotels = ntotelb

        CALL GA_GATHER(lg_BDER,BUFFB,INDI,INDJ,NtotELB)
        CALL GA_GATHER(lg_SDER,BUFFS,INDI,INDJ,NtotELS)

        ! Finally, fill the local chunk of the SC matrix (block of rows)
        ! with the received values at their appropriate place.
        NtotELB = 0
        NtotELS = 0
        DO iG3=IG3STA,IG3END
          iT=idxG3(1,iG3)
          iU=idxG3(2,iG3)
          iV=idxG3(3,iG3)
          iX=idxG3(4,iG3)
          iY=idxG3(5,iG3)
          iZ=idxG3(6,iG3)
          NELB = NELBsav(iG3-IG3STA+1)
          NELS = NELSsav(iG3-IG3STA+1)

          F3VAL = 0.0d+00
          DO i = 1, NELB
            NtotELB = NtotELB + 1
            F3VAL = F3VAL + BUFFB(NtotELB)
          END DO
C         BUFFB(iG3-IG3STA+1) = F3VAL

          G3VAL = 0.0d+00
          DO i = 1, NELS
            NtotELS = NtotELS + 1
            G3VAL = G3VAL + BUFFS(NtotELS)
          END DO
          G3VAL = G3VAL - (EPSA(iU)+EPSA(iY))*F3VAL

          DF3(iG3) = DF3(iG3) + F3VAL
          DG3(iG3) = DG3(iG3) + G3VAL

          !! do depsa
C         ISUP=KTUV(iV,iU,iT)-nTUVES(iSYM)
C         If (ISUP.GE.iLo .and. ISUP.LE.iHi) THen
C           Do iW = 1, nAshT
C             JSUP=KTUV(iX,iW,iZ)-nTUVES(iSYM)
C             NSEQ=1+ISUP-iLo+LDA*(JSUP-jLo)
C             DEPSA(iW,iY) = DEPSA(iW,iY) - F3VAL*SC(ISUP-iLo,JSUP)
C           End Do
C         End If

          !! remaining F3 and G3 transformation in mkfg3.f
          If (iY.eq.iX) Then
            DF2(iT,iU,iV,iZ) = DF2(iT,iU,iV,iZ) - F3VAL
            DG2(iT,iU,iV,iZ) = DG2(iT,iU,iV,iZ) - EPSA(iU)*F3VAL
            Do iW = 1, nAshT
              DEPSA(iU,iW) = DEPSA(iU,iW) - F3VAL*G2(iT,iW,iV,iZ)
            End Do
            DG2(iT,iU,iV,iZ) = DG2(iT,iU,iV,iZ) - G3VAL
          End If
          If (iV.eq.iU) Then
            DF2(iT,iX,iY,iZ) = DF2(iT,iX,iY,iZ) - F3VAL
            DG2(iT,iX,iY,iZ) = DG2(iT,iX,iY,iZ) - EPSA(iY)*F3VAL
            Do iW = 1, nAshT
              DEPSA(iW,iY) = DEPSA(iW,iY) - F3VAL*G2(iT,iX,iW,iZ)
            End Do
            DG2(iT,iX,iY,iZ) = DG2(iT,iX,iY,iZ) - G3VAL
          End If
          If (iY.eq.iU) Then
            DF2(iV,iX,iT,iZ) = DF2(iV,iX,iT,iZ) - F3VAL
            DG2(iV,iX,iT,iZ) = DG2(iV,iX,iT,iZ) - EPSA(iU)*F3VAL
            DG2(iV,iX,iT,iZ) = DG2(iV,iX,iT,iZ) - G3VAL
          End If
          DEPSA(iY,iU) = DEPSA(iY,iU) - F3VAL*G2(iV,iX,iT,iZ)
          If (iY.eq.iX.and.iV.eq.iU) Then
            DF1(iT,iZ) = DF1(iT,iZ) - F3VAL
            DG1(iT,iZ) = DG1(iT,iZ) - G3VAL
          End If
        END DO

        call mma_deallocate(NELBsav)
        call mma_deallocate(NELSsav)
      END DO ! end loop over blocks of G3 values

      call mma_deallocate(BUFFB)
      call mma_deallocate(BUFFS)
      call mma_deallocate(INDI)
      call mma_deallocate(INDJ)

      RETURN
C
      End Subroutine CLagDXC_FG3_MPP
#endif
C
C-----------------------------------------------------------------------
C
      Subroutine CLagDXC_DP(iSym,nAS,BDER,SDER,DG1,DG2,DF1,DF2,
     *                      DEPSA,DEASUM,iLo,iHi,jLo,jHi,LDC,
     *                      G1,G2,SC,SC2,lg_S)
C
      use caspt2_global, only:ipea_shift
      USE SUPERINDEX
      use EQSOLV
      use caspt2_module
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, nProcs
#endif
C
      Implicit Real*8 (A-H,O-Z)
C
#ifdef _MOLCAS_MPP_
#include "global.fh"
#endif
C
      Dimension BDER(*),SDER(*),
     *          DF2(nAshT,nAshT,nAshT,nAshT),
     *          DG2(nAshT,nAshT,nAshT,nAshT),
     *          DF1(nAshT,nAshT),DG1(nAshT,nAshT),DEPSA(nAshT,nAshT)
      dimension g1(nAshT,nAshT),G2(nAshT,nAshT,nasht,nasht),SC(*),SC2(*)
C     INTEGER*1 idxG3(6,NG3)
C
C     LDC == 0, if not parallel; SC is triangular
C     LDC /= 0, if parallel    ; SC is square
C     In both cases, BDER and SDER are square
C
      ISADR=0
      if (isadr.ne.0) write (6,*) ldc !! just for avoid compiling error
      if (isadr.ne.0) write (6,*) lg_S !! just for avoid compiling error
      NROW = 0
      iLoS = 0
      jLoS = 0
#ifdef _MOLCAS_MPP_
      if (is_real_par()) then
        irank = 0
        CALL GA_Distribution (lg_S,iRank,iLoS,iHiS,jLoS,jHiS)
        NROW = jHiS-jLoS+1 !! = NAS
        CALL GA_GET(lg_S,jLoS,jHiS,iLoS,iHiS,SC2,NROW)
      end if
#endif
      DO 100 IXYZ=jLo,jHi
        IXYZABS=IXYZ+NTUVES(ISYM)
        IXABS=MTUV(1,IXYZABS)
        IYABS=MTUV(2,IXYZABS)
        IZABS=MTUV(3,IXYZABS)
        EY=EPSA(IYABS)
        DO 101 ITUV=iLo,iHi
          ITUVABS=ITUV+NTUVES(ISYM)
          ITABS=MTUV(1,ITUVABS)
          IUABS=MTUV(2,ITUVABS)
          IVABS=MTUV(3,ITUVABS)
          EU=EPSA(IUABS)
          EYU=EY + EU
          FACT=EYU-EASUM
          IF (LDC.EQ.0) THEN
            ISADR=1+iTUV-iLo+NAS*(iXYZ-jLo)
          ELSE
            ISADR=1+iTUV-iLo+LDC*(iXYZ-jLo)
          END IF
          ValB=BDER(ISADR)
C
          If (iTUV.eq.iXYZ.and.ipea_shift.ne.0.0D+00) Then
C           !! BC in the next equation refers to the active overlap
C    !! ipea_shift*0.5d0*BC(ISADR)*(4.0d0-DREF(IDT)-DREF(IDV)+DREF(IDU))
            bsBDER = ipea_shift*0.5D+00*ValB
            SDER(iSAdr) = SDER(iSAdr) + bsBDER*(4.0D+00
     *       -G1(iTabs,iTabs)+G1(iUabs,iUabs)-G1(iVabs,iVabs))
            IF (LDC.EQ.0) THEN
              iSAdr2 = iTUV*(iTUV+1)/2
            ELSE
              ISADR2 = 1+iTUV-iLo+LDC*(iTUV-jLo)
            END IF
            DG1(iTabs,iTabs) = DG1(iTabs,iTabs) - bsBDER*SC(iSAdr2)
            DG1(iUabs,iUabs) = DG1(iUabs,iUabs) + bsBDER*SC(iSAdr2)
            DG1(iVabs,iVabs) = DG1(iVabs,iVabs) - bsBDER*SC(iSAdr2)
          End If
C
          !! First VALUE contribution in MKBC_DP (FACT)
C         IF (LDC.EQ.0) ISADR=ITUV*(ITUV-1)/2+IXYZ
          SDER(ISADR) = SDER(ISADR) + FACT*ValB
          ValS=SDER(ISADR)
C
          !! DEPSA can be computed simultaneously with the F3
          !! contributions, but remember that SC here is the overlap,
          !! whereas SC in FG3 subroutines are just G3 matrix
          if (ldc.eq.0) then
            Do iWabs = 1, nAshT
              iTWV = iTabs+nAshT*(iWabs-1)+nAshT**2*(iVabs-1)
              iSAdr2 = Max(iTWV,iXYZ)*(Max(iTWV,iXYZ)-1)/2
     *               + Min(iTWV,iXYZ)
              DEPSA(iWabs,iUabs) = DEPSA(iWabs,iUabs)
     *          + ValB*SC(iSAdr2)
C
              iXWZ = iXabs+nAshT*(iWabs-1)+nAshT**2*(iZabs-1)
              iSAdr2 = Max(iTUV,iXWZ)*(Max(iTUV,iXWZ)-1)/2
     *               + Min(iTUV,iXWZ)
              DEPSA(iWabs,iYabs) = DEPSA(iWabs,iYabs)
     *          + ValB*SC(iSAdr2)
            End Do
          else
            !! assume SC has all columns (which is reasonable)
C           iTWV = iTabs+nAshT**2*(iVabs-1)
C           ISADR2 = 1+iTWV-jLoS+NROW*(iXYZ-iLoS)
C           Call DaXpY_(nAshT,ValB,SC2(iSAdr2),nAshT,DEPSA(1,iUabs),1)
C           iXWZ = iXabs+nAshT**2*(iZabs-1)
C           ISADR2 = 1+iTUV-iLo+LDC*(iXWZ-jLo)
C           Call DaXpY_(nAshT,ValB,SC(iSAdr2),LDC*nAshT,
C    *                             DEPSA(1,iYabs),1)
            Do iWabs = 1, nAshT
              !! we do not have all elements, so use distributed memory
              iTWV = iTabs+nAshT*(iWabs-1)+nAshT**2*(iVabs-1)
              ISADR2 = 1+iTWV-jLoS+NROW*(iXYZ-iLoS)
              DEPSA(iWabs,iUabs) = DEPSA(iWabs,iUabs)
     *          + ValB*SC2(iSAdr2)

              !! we have all elements, so local memory
              iXWZ = iXabs+nAshT*(iWabs-1)+nAshT**2*(iZabs-1)
              ISADR2 = 1+iTUV-iLo+LDC*(iXWZ-jLo)
              DEPSA(iWabs,iYabs) = DEPSA(iWabs,iYabs)
     *          + ValB*SC(iSAdr2)
            End Do
          end if

          !! If non-parallel, overlap is triangular
          !! If parallel, the index is the same to that of SDER
          IF (LDC.EQ.0) THEN
            iSAdr = Max(iTUV,iXYZ)*(Max(iTUV,iXYZ)-1)/2
     *            + Min(iTUV,iXYZ)
          END IF
          DEASUM = DEASUM - ValB*SC(iSAdr)
C
C         dyu ( Fvztx - EPSA(u)*Gvztx )
C         dyu Gvztx
          IF(IYABS.EQ.IUABS) THEN
            !! VALUE=VALUE+2.0D0*(FP(IP)-EU*PREF(IP))
            DF2(iVabs,iZabs,iTabs,iXabs)
     *        = DF2(iVabs,iZabs,iTabs,iXabs) + ValB
            DG2(iVabs,iZabs,iTabs,iXabs)
     *        = DG2(iVabs,iZabs,iTabs,iXabs) - EU*ValB
C
            !! VALUE=VALUE+2.0D0*PREF(IP)
            DG2(iVabs,iZabs,iTabs,iXabs)
     *        = DG2(iVabs,iZabs,iTabs,iXabs) + ValS
          END IF
          DEPSA(iYabs,iUabs) = DEPSA(iYabs,iUabs)
     *      - ValB*G2(iVabs,iZabs,iTabs,iXabs)
C
C         dyx ( Fvutz - EPSA(y)*Gvutz )
C         dyx Gvutz -> dut Gzyxv
          IF(IYABS.EQ.IXABS) THEN
            !! VALUE=VALUE+2.0D0*(FP(IP)-EY*PREF(IP))
            DF2(iVabs,iUabs,iTabs,iZabs)
     *        = DF2(iVabs,iUabs,iTabs,iZabs) + ValB
            DG2(iVabs,iUabs,iTabs,iZabs)
     *        = DG2(iVabs,iUabs,iTabs,iZabs) - EY*ValB
C
            !! VALUE=VALUE+2.0D0*PREF(IP)
            DG2(iVabs,iUabs,iTabs,iZabs)
     *        = DG2(iVabs,iUabs,iTabs,iZabs) + ValS
          END IF
          DEPSA(iYabs,iXabs) = DEPSA(iYabs,iXabs)
     *      - ValB*G2(iVabs,iUabs,iTabs,iZabs)

C         dtu ( Fvxyz - EPSA(u)*Gvxyz + dyx Fvz -
C                (EPSA(u)+EPSA(y)*dyz Gvz)
C         dtu Gvxyz + dtu dyx Gvz
          IF(ITABS.EQ.IUABS) THEN
            !! VALUE=VALUE+2.0D0*(FP(IP)-EU*PREF(IP))
            DF2(iVabs,iXabs,iYabs,iZabs)
     *        = DF2(iVabs,iXabs,iYabs,iZabs) + ValB
            DG2(iVabs,iXabs,iYabs,iZabs)
     *        = DG2(iVabs,iXabs,iYabs,iZabs) - EU*ValB
C
            !! VALUE=VALUE+2.0D0*PREF(IP)
            DG2(iVabs,iXabs,iYabs,iZabs)
     *        = DG2(iVabs,iXabs,iYabs,iZabs) + ValS
            IF(IYABS.EQ.IXABS) THEN
              !! VALUE=VALUE+FD(ID)-EYU*DREF(ID)
              DF1(iVabs,iZabs) = DF1(iVabs,iZabs) + ValB
              DG1(iVabs,iZabs) = DG1(iVabs,iZabs) - EYU*ValB
C
              !! VALUE=VALUE+DREF((ID1*(ID1-1))/2+ID2)
              DG1(iVabs,iZabs) = DG1(iVabs,iZabs) + ValS
            END IF
          END IF
          DEPSA(iTabs,iUabs) = DEPSA(iTabs,iUabs)
     *      - ValB*G2(iVabs,iXabs,iYabs,iZabs)
          If (iYabs.eq.iXabs)
     *    DEPSA(iTabs,iUabs) = DEPSA(iTabs,iUabs) - ValB*G1(iVabs,iZabs)
          If (iTabs.eq.iUabs)
     *    DEPSA(iYabs,iXabs) = DEPSA(iYabs,iXabs) - ValB*G1(iVabs,iZabs)
 101    CONTINUE
#ifdef _MOLCAS_MPP_
        if (is_real_par() .and. IXYZ==iHiS .and. iRank/=NPROCS-1) then
          irank = irank + 1
          CALL GA_Distribution (lg_S,iRank,iLoS,iHiS,jLoS,jHiS)
          CALL GA_GET(lg_S,jLoS,jHiS,iLoS,iHiS,SC2,NROW)
        end if
#endif
 100  CONTINUE
C
      Return
C
      End Subroutine CLagDXC_DP
C
C-----------------------------------------------------------------------
C
#ifdef _MOLCAS_MPP_
      Subroutine DF3_DEPSA_MPP(DF3,DEPSA,lg_S,idxG3)
      USE MPI
      USE SUPERINDEX
      USE EQSOLV
      USE Para_Info, only: nProcs
      use definitions, only: wp
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module
      use pt2_guga
C
      Implicit Real*8 (A-H,O-Z)

#include "global.fh"
#include "mafdecls.fh"

      DIMENSION DF3(*),DEPSA(NASHT,NASHT)
      INTEGER*1 idxG3(6,NG3)
      real(kind=wp), allocatable :: WRK(:)

      !! do depsa
      isym=1
      Do iRank = 0, NPROCS-1
        CALL GA_Distribution(lg_S,iRank,ILO,IHI,JLO,JHI)
        NROW=iHi-iLo+1
        NCOL=jHi-jLo+1
        Call mma_allocate(WRK,NROW*NCOL,Label='WRK')
        Call GA_Get(lg_S,iLo,iHi,jLo,jHi,WRK,NROW)
        DO iG3=1,NG3
          iT=idxG3(1,iG3)
          iU=idxG3(2,iG3)
          iV=idxG3(3,iG3)
          iX=idxG3(4,iG3)
          iY=idxG3(5,iG3)
          iZ=idxG3(6,iG3)
          iST=IASYM(iT)
          iSU=IASYM(iU)
          iSV=IASYM(iV)
          iSX=IASYM(iX)
          iSY=IASYM(iY)
          iSZ=IASYM(iZ)
          ituvs=MUL(IST,MUL(ISU,ISV))
          ixyzs=MUL(ISX,MUL(ISY,ISZ))
          if(ituvs.ne.ixyzs) cycle

          F3VAL = DF3(iG3)

          ISUP1=KTUV(iV,iU,iT)-nTUVES(iSYM)
C         write (*,'("fval=",6i3,f20.10)') it,iu,iv,ix,iy,iz,f3val
          If (ISUP1.GE.iLo .and. ISUP1.LE.iHi) Then
            Do iW = 1, nAshT
              JSUP1=KTUV(iX,iW,iZ)-nTUVES(iSYM)
              NSEQ=1+ISUP1-iLo+NROW*(JSUP1-jLo)
              DEPSA(iW,iY) = DEPSA(iW,iY) - F3VAL*WRK(NSEQ)
C       write (*,'(8i3,1f20.10)')
C    *   iv,iu,it,ix,iw,iz,isup1,jsup1,wrk(nseq)
            End Do
          End If
          JSUP2=KTUV(iX,iY,iZ)-nTUVES(iSYM)
          If (JSUP2.GE.iLo .and. JSUP2.LE.iHi) Then
            Do iW = 1, nAshT
              ISUP2=KTUV(iV,iW,iT)-nTUVES(iSYM)
              NSEQ=1+JSUP2-iLo+NROW*(ISUP2-jLo)
              DEPSA(iW,iU) = DEPSA(iW,iU) - F3VAL*WRK(NSEQ)
C       write (*,'(8i3,1f20.10)')
C    *   iv,iu,it,ix,iw,iz,isup1,jsup1,wrk(nseq)
            End Do
          End If
        END DO
        call mma_deallocate(WRK)
      End Do
C
      Return
C
      End Subroutine DF3_DEPSA_MPP
#endif
C
C-----------------------------------------------------------------------
C
      Subroutine DEPSAOffC(CLag,DEPSA,FIFA,FIMO,WRK1,WRK2,U0)
C
      use caspt2_global, only:IPrGlb
      use PrintLevel, only: verbose
      use caspt2_global, only: ConvInvar,SLag
      use gugx, only: SGS, CIS
      use caspt2_global, only: LUCIEX, IDCIEX, IDTCEX
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: iwp,wp
      use caspt2_module
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      Implicit Real*8 (A-H,O-Z)
C
C
      Dimension CLag(nConf,nState),DEPSA(nAshT,nAshT),FIFA(*),FIMO(*),
     *          WRK1(nBasT,nBasT),WRK2(*),U0(nState,nState)
      Dimension Eact(nState)
      Integer :: nLev, nMidV
      real(kind=wp),allocatable :: VecST(:,:),VecS1(:,:),VecS2(:,:),
     *                             VecCID(:,:),VecPre(:),VecFancy(:),
     *                             VecCIT(:,:),INT1(:),INT2(:),G2(:)

      nLev = SGS%nLev
      nMidV= CIS%nMidV
C
      Thres = ConvInvar !! 1.0d-07
C
      If (IPRGLB.GE.verbose) Then
        Write (6,*)
        Write (6,'(3X,"Linear Equation for Non-Invariant CASPT2",
     *                " (threshold =",ES9.2,")")') Thres
        Write (6,*)
        CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
      End If
C
C     If CASPT2 energy is not invariant with respect to rotations within
C     active space (with IPEA shift and/or with RAS reference), the
C     active density obtained in constructing CI derivative is no longer
C     correct... well, it may be correct, but orbital rotations in the
C     active space cannot be parametrized in Z-vector, so analytic
C     derivatives cannot be computed with the existing module. So,
C     the active density is computed in a differnt way.
C
C     See J. Chem. Phys. 2023, 158, 174112. for details, in particular,
C     Section II C 4 "Non-invariance with respect to active MOs"
C     To be more specific, this subroutine solves the linear equation
C     (Eq. (71)) and computes the second term in Eq. (70) later.
C     CLag corresponds to the RHS in Eq. (71).
C
      !! Some post-processing of CI derivative
      !! Somehow, this has to be done in the XMS basis
      Call CLagFinalOffC(SLag)
C
C     ----- Solve the linear equation -----
C     A_{IS,JR}*X_{JR} = CLag_{IS}, where A_{IS,JR} is the CI-CI Hessian
C     which may be seen in Z-vector
C
      call mma_allocate(VecST,nConf,nState,Label='VecST')
      call mma_allocate(VecS1,nConf,nState,Label='VecS1')
      call mma_allocate(VecS2,nConf,nState,Label='VecS2')
      call mma_allocate(VecCID,nConf,nState,Label='VecCID')
C     call mma_allocate(VecS,nState*(nState-1)/2,Label='VecS')
      call mma_allocate(VecPre,nConf,Label='VecPre')
      call mma_allocate(VecFancy,nState**3,Label='VecFancy')
C
      call mma_allocate(VecCIT,nConf,nState,Label='VecCIT')
      call mma_allocate(INT1,nAshT**2,Label='INT1')
      call mma_allocate(INT2,nAshT**4,Label='INT2')
C
      !! We do not have Cholesky vectors for frozen orbitals,
      !! so may be it is not possible to get inactive energies?
      !! It can be computed with TimesE2
      iSym = 1
      Call CnstInt(0,INT1,INT2)
      ID = IDTCEX !! IDCIEX !! this parameter is hacked
      Do iState = 1, nState
        If (ISCF.EQ.0) Then
          !! quasi-canonical, XMS
          Call DDaFile(LUCIEX,2,VecCIT(1,iState),nConf,ID)
        Else
          VecCIT(1,iState) = 1.0d+00
        End If
        !! The second term should be removed
        Eact(iState)=0.0d+00
      End Do
      if (ifxms.or.ifrms) then
        !! Transform the CLag and CI vector from XMS to SCF basis
        !! Maybe, in order to define Eact
        Call DGEMM_('N','T',nConf,nState,nState,
     &              1.0D+00,CLag,nConf,U0,nState,
     &              0.0D+00,VecST,nConf)
        Call DCopy_(nConf*nState,VecST,1,CLag,1)
        Call DGEMM_('N','T',nConf,nState,nState,
     &              1.0D+00,VecCIT,nConf,U0,nState,
     &              0.0D+00,VecST,nConf)
        Call DCopy_(nConf*nState,VecST,1,VecCIT,1)
      end if
      Call TimesE2(0,VecCIT,VecS1,INT1,INT2)
      Do iState = 1, nState
        !! scaling with nState is due to the division in TimesE2
        Eact(iState) = -0.5D+00*nState*
     *    DDot_(nConf,VecS1(1,iState),1,VecCIT(1,iState),1)
      End Do
      isyci = 1
C
      !! Precondition
      Call CnstInt(2,INT1,INT2)
      Call CnstPrec(CIS%NOCSF,CIS%IOCSF,CIS%NOW,
     *              CIS%IOW,ISYCI,VecPre,VecCIT,
     *              INT1,INT2,VecFancy,nLev,
     *              nMidV)
      Call CnstInt(0,INT1,INT2)
C
      !! Begin!
      Call DCopy_(nConf*nState,CLag,1,VecST,1)
C
      !! z0 = M^{-1}*r0
      Call DCopy_(nConf*nState,VecST,1,VecS2,1)
      Call DoPrec(VecST,VecS2,VecS1,VecPre,VecFancy)
      !! p0 = z0
      Call DCopy_(nConf*nState,VecS2,1,VecCId,1)
      MaxIter = 100
      Iter    = 1
      iSym    = 1
      ! jspin   = 0
      ! r^T dot z
      ! r (residue) = ipST
      ! z (prec. r) = ipS2
      ! p (...)     = ipCId
      ! x (solution)= ipCIT
      ! Ap          = ipS1
      ! r_{k}z_{k}  = ipST*ipS2 = deltaC
      DeltaC = DDot_(nConf*nState,VecST,1,VecS2,1)
      Delta  = DeltaC
      Delta0 = Delta
C
      If (IPRGLB.GE.verbose) Write(6,*)
     &      ' Iteration       Delta           Res(CI)        '//
     &      '  DeltaC'
      Call DCopy_(nConf*nState,[0.0D+00],0,VecCIT,1)
      If (Delta0.le.Abs(Thres)) Go To 100
      Do Iter = 1, MaxIter
        If (nConf.EQ.1) Then
          Do iState = 1, nState
            VecCIT(1,iState) = 1.0d+00
          End Do
          Exit
        End If
        !! Compute Ap
        !! ipS2 is used as a workind array
        Call TimesE2(1,VecCId,VecS1,INT1,INT2)
C
        !! AlphaC = p^T*A*p
        AlphaC= DDot_(nConf*nState,VecS1,1,VecCId,1)
        !! Alpha = r^T*z / AlphaC
        Alpha = Delta/(AlphaC)
        ! new x of CI
        Call DaXpY_(nConf*nState,Alpha,VecCId,1,VecCIT,1)
        ! new r of CI
        Call DaXpY_(nConf*nState,-Alpha,VecS1,1,VecST,1)
        ResCI=sqrt(DDot_(nConf*nState,VecST,1,VecST,1))
        !! z = M^{-1}*r
        Call DCopy_(nConf*nState,VecST,1,VecS2,1)
        Call DoPrec(VecST,VecS2,VecS1,VecPre,VecFancy)
C
        !! Append new vectors
        DeltaC= Ddot_(nConf*nState,VecST,1,VecS2,1)
        Beta  = DeltaC/Delta
        Delta = DeltaC
        Call DScal_(nConf*nState,   Beta,VecCId,1)
        Call DaXpY_(nConf*nState,1.0D+00,VecS2,1,VecCId,1)
C
        If (IPRGLB.GE.verbose)
     *  Write(6,'(I7,4X,ES17.9,ES17.9,ES17.9)')
     &         iter,delta/delta0,resci,deltac
C
        Res = ResCI
        If (Res.le.Abs(Thres)) Exit
      End Do
C
      If (Iter.eq.MaxIter+1) Then
        write(6,*)
     *  "CI iteration for non-invariant CASPT2 did not converge..."
        call abend
      End If
C
  100 CONTINUE
C
      If (IPRGLB.GE.verbose) Then
        CALL TIMING(CPTF1,CPE,TIOTF1,TIOE)
        CPUT =CPTF1-CPTF0
        WALLT=TIOTF1-TIOTF0
        Write (6,*)
        Write (6,'(3X,"Linear equation converged in ",I3," steps")')
     *         iter-1
        Write (6,'(3X,"CPU and wall time (in s) = ",2F8.2)') CPUT,WALLT
        Write (6,*)
      End If
C
      If (IFXMS.OR.IFRMS) Then
        !! Transform back the CLag from CAS to XMS
        Call DGEMM_('N','N',NConf,nState,nState,
     &              1.0D+00,CLag,nConf,U0,nState,
     &              0.0D+00,VecST,nConf)
        Call DCopy_(nConf*nState,VecST,1,CLag,1)
      End If
C
      call mma_deallocate(VecS1)
      call mma_deallocate(VecS2)
      call mma_deallocate(VecCId)
C     call mma_deallocate(VecS)
      call mma_deallocate(VecPre)
      call mma_deallocate(VecFancy)
C
C     ----- Construct (a part of) the true active density -----
C     Compute the second term in Eq. (70) = Eq. (72)
C     The SCF, not XMS, basis is used
C
      ID = IDCIEX !! idtcex?
      Do iState = 1, nState
        If (ISCF.EQ.0) Then
          If (IFXMS.OR.IFRMS) THen
            !! Use unrotated (SCF) CI vector
            Call LoadCI_XMS('C',1,VecST(1,iState),iState,U0)
          Else
            Call DDaFile(LUCIEX,2,VecST(1,iState),nConf,ID)
          End If
        Else
          VecST(1,iState) = 1.0d+00
        End If
      End Do
      call mma_allocate(G2,nAshT**4,Label='G2')
      Call CnstInt(1,INT1,INT2)
      Call CnstDEPSA(VecST,VecCIT,INT1,G2,INT2)
      call mma_deallocate(G2)
C
      If (IPRGLB.GE.verbose) Then
        CALL TIMING(CPTF2,CPE,TIOTF2,TIOE)
        CPUT =CPTF2-CPTF1
        WALLT=TIOTF2-TIOTF1
        Write (6,'(3X,"Off-diagonal density is constructed")')
        Write (6,'(3X,"CPU and wall time (in s) = ",2F8.2)') CPUT,WALLT
        Write (6,*)
      End If
C
      call mma_deallocate(VecST)
      call mma_deallocate(VecCIT)
      call mma_deallocate(INT1)
      call mma_deallocate(INT2)
C
      Contains
C
C-----------------------------------------------------------------------
C
      Subroutine CLagFinalOffC(SLag)
C
C     Orthogonalize the partial derivative with respect to CI coeff
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      Dimension SLag(*)
      real(kind=wp),allocatable :: CI1(:),CI2(:)
C
      call mma_allocate(VecST,nConf,nState,Label='VecST')
      Call DGEMM_('N','T',nConf,nState,nState,
     &            1.0D+00,CLag,nConf,U0,nState,
     &            0.0D+00,VecST,nConf)
      Call DCopy_(nConf*nState,VecST,1,CLag,1)

      call mma_allocate(CI1,nConf,Label='CI1')
      call mma_allocate(CI2,nConf,Label='CI2')
C
      !! Construct SLag
      ijst = 0
      do ilStat = 1, nState
        If (ISCF.EQ.0) Then
          Call LoadCI_XMS('C',1,CI1,ilStat,U0)
        Else
          CI1(1) = 1.0D+00
        End If
        Do jlStat = 1, ilStat !! -1
          ijst = ilStat + nState*(jlStat-1)
          If (ilStat.eq.jlStat) Cycle
          If (ISCF.EQ.0) Then
            Call LoadCI_XMS('C',1,CI2,jlStat,U0)
          Else
            CI2(1) = 1.0D+00
          End If
          Scal = DDOT_(nConf,CI1,1,CLag(1,jlStat),1)
     *         - DDOT_(nConf,CI2,1,CLag(1,ilStat),1)
          Scal = Scal/(REFENE(jlStat)-REFENE(ilStat))
          SLag(ijst) = SLag(ijst) + Scal
          IF (IPRGLB.GE.verbose) THEN
            write(6,'(1x,"SLag for State ",i1,"-",i1," = ",f20.10)')
     *         ilstat,jlstat,slag(ijst)
          END IF
        end do
      end do
C
      !! Projection
      Do ilStat = 1, nState
        Call DCopy_(nConf,CLag(1,ilStat),1,CI1,1)
        Do jlStat = 1, nState
          If (ISCF.EQ.0) Then
            Call LoadCI_XMS('C',1,CI2,jlStat,U0)
          Else
            CI2(1) = 1.0D+00
          End If
          Ovl = DDot_(nConf,CI1,1,CI2,1)
          Call DaXpY_(nConf,-Ovl,CI2,1,CLag(1,ilStat),1)
        End Do
      End Do
C
      Call DGEMM_('N','N',nConf,nState,nState,
     &            1.0D+00,CLag,nConf,U0,nState,
     &            0.0D+00,VecST,nConf)
      Call DCopy_(nConf*nState,VecST,1,CLag,1)
C
      call mma_deallocate(CI1)
      call mma_deallocate(CI2)
      call mma_deallocate(VecST)
C
      Return
C
      End Subroutine CLagFinalOffC
C
C-----------------------------------------------------------------------
C
      Subroutine CnstInt(Mode,INT1,INT2)
C
      Use CHOVEC_IO
      use ChoCASPT2
C
      Implicit Real*8 (A-H,O-Z)
C
      Real*8 INT1(nAshT,nAshT),INT2(nAshT,nAshT,nAshT,nAshT)
      integer(kind=iwp),allocatable :: BGRP(:,:)
      real(kind=wp),allocatable :: KET(:)
C
      Integer Active, Inactive, Virtual
      Parameter (Inactive=1, Active=2, Virtual=3)
      Integer nSh(8,3)
C
      Call DCopy_(nAshT**2,[0.0D+00],0,INT1,1)
      Call DCopy_(nAshT**4,[0.0D+00],0,INT2,1)
C
      nFroI = nFro(iSym)
      nIshI = nIsh(iSym)
      nCorI = nFroI+nIshI
      nBasI = nBas(iSym)
C
C     --- One-Electron Integral
C
      !! Read H_{\mu \nu}
C     IRC=-1
C     IOPT=6
C     ICOMP=1
C     ISYLBL=1
C     CALL RDONE(IRC,IOPT,'OneHam  ',ICOMP,WRK2,ISYLBL)
C     !! triangular -> square transformation
C     Call Square(WRK2,WRK1,1,nBasT,nBasT)
C     !! AO -> MO transformation
C     Call DGemm_('T','N',nBasT,nBasT,nBasT,
C    *            1.0D+00,CMOPT2,nBasT,WRK1,nBasT,
C    *            0.0D+00,WRK2,nBasT)
C     Call DGemm_('N','N',nBasT,nBasT,nBasT,
C    *            1.0D+00,WRK2,nBasT,CMOPT2,nBasT,
C    *            0.0D+00,WRK1,nBasT)
      !! Inactive energy
C     Do iCorI = 1, nFro(iSym)+nIsh(iSym)
C       RIn_Ene = RIn_Ene + 2.0d+00*WRK1(iCorI,iCorI)
C     End Do
      !! Put in INT1
C     Do iAshI = 1, nAsh(iSym)
C       Do jAshI = 1, nAsh(iSym)
C         Val = WRK1(nCorI+iAshI,nCorI+jAshI)
C         INT1(iAshI,jAshI) = INT1(iAshI,jAshI) + Val
C       End Do
C     End Do
      Do iAshI = 1, nAsh(iSym)
        Do jAshI = 1, nAsh(iSym)
          Val = FIMO(nCorI+iAshI+nBasI*(nCorI+jAshI-1))
          INT1(iAshI,jAshI) = INT1(iAshI,jAshI) + Val
        End Do
      End Do
C
C     --- Two-Electron Integral
C
      iSymA = 1
      iSymI = 1
      iSymB = 1
      iSymJ = 1
C     If (.not.IfChol) Then
C       Do iCorI = 1, nFro(iSym)+nIsh(iSym)
C         iOrb = iCorI
C         jOrb = iCorI
C         Call Coul(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
C         Do jCorI = 1, nFro(iSym)+nIsh(iSym)
C           RIn_Ene = RIn_Ene + 2.0d+00*WRK1(jCorI,jCorI)
C         End Do
C         Call Exch(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
C         Do jCorI = 1, nFro(iSym)+nIsh(iSym)
C           RIn_Ene = RIn_Ene - WRK1(jCorI,jCorI)
C         End Do
C       End Do
C     End If
C
      If (IfChol) Then
        Call ICopy(NSYM,NISH,1,nSh(1,Inactive),1)
        Call ICopy(NSYM,NASH,1,nSh(1,Active  ),1)
        Call ICopy(NSYM,NSSH,1,nSh(1,Virtual ),1)
        DO JSYM=1,NSYM
          IB1=NBTCHES(JSYM)+1
          IB2=NBTCHES(JSYM)+NBTCH(JSYM)
C
          MXBGRP=IB2-IB1+1
          IF (MXBGRP.LE.0) CYCLE
          call mma_allocate(BGRP,2,MXBGRP,Label='BGRP')
          IBGRP=1
          DO IB=IB1,IB2
           BGRP(1,IBGRP)=IB
           BGRP(2,IBGRP)=IB
           IBGRP=IBGRP+1
          END DO
          NBGRP=MXBGRP
C
          CALL MEMORY_ESTIMATE(JSYM,BGRP,NBGRP,
     &                         NCHOBUF,MXPIQK,NADDBUF)
          call mma_allocate(KET,NCHOBUF,Label='KETBUF')
C         write(6,*) "nchobuf= ", nchobuf
C         write(6,*) "nbgrp= ", nbgrp
C         write(6,*) "nbtch= ", nbtch(jsym)
          Do IBGRP=1,NBGRP
C
            IBSTA=BGRP(1,IBGRP)
            IBEND=BGRP(2,IBGRP)
C           write(6,*) ibsta,ibend
C
            NV=0
            DO IB=IBSTA,IBEND
              NV=NV+NVLOC_CHOBATCH(IB)
            END DO
C
            !! int2(tuvx) = (tu|vx)/2
            !! This can be computed without frozen orbitals
            Call Get_Cholesky_Vectors(Active,Active,JSYM,
     &                                KET,nKet,
     &                                IBSTA,IBEND)
C
            Call DGEMM_('N','T',NASH(JSYM)**2,NASH(JSYM)**2,NV,
     *                  0.5D+00,KET,NASH(JSYM)**2,
     *                          KET,NASH(JSYM)**2,
     *                  0.0D+00,INT2,NASH(JSYM)**2)
          End Do
          call mma_deallocate(KET)
          call mma_deallocate(BGRP)
        End Do
      Else
        Do iAshI = 1, nAsh(iSym)
          iOrb = nCorI+iAshI
          Do jAshI = 1, nAsh(iSym)
            jOrb = nCorI+jAshI
C
            Call Coul(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
            !! Put in INT1
C           Do iCorI = 1, nFro(iSym)+nIsh(iSym)
C             INT1(iAshI,jAshI) = INT1(iAshI,jAshI)
C    *          + 2.0d+00*WRK1(iCorI,iCorI)
C           End Do
            !! Put in INT2
            Do kAshI = 1, nAsh(iSym)
              Do lAshI = 1, nAsh(iSym)
                INT2(iAshI,jAshI,kAshI,lAshI)
     *        = INT2(iAshI,jAshI,kAshI,lAshI)
     *        + WRK1(nCorI+kAshI,nCorI+lAshI)*0.5d+00
              End Do
            End Do
C
C           Call Exch(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
            !! Put in INT1
C           Do iCorI = 1, nFro(iSym)+nIsh(iSym)
C             INT1(iAshI,jAshI) = INT1(iAshI,jAshI) - WRK1(iCorI,iCorI)
C           End Do
          End Do
        End Do
      End If
#ifdef _MOLCAS_MPP_
      call GADSUM(INT2,nAshT**4)
#endif
C     write(6,*) "int2"
C     call sqprt(int2,25)
C     call sqprt(int1,5)
C     call sqprt(fimo,12)
      If (Mode.eq.0) Then
      Do IT = 1, nAshT
        Do iU = 1, nAshT
          iTU = iT + nAshT*(iU-1)
          Do iV = 1, nAshT
            Do iX = 1, nAshT
              iVX = iV + nAshT*(iX-1)
              If (iVX.gt.iTU) Then
               INT2(iT,iU,iV,IX) = INT2(iT,iU,iV,iX) + INT2(iV,iX,iT,iU)
               INT2(iV,iX,iT,iU) = 0.0D+00
              End If
            End Do
          End Do
        End Do
      End Do
      End If
C
      if (mode.eq.0.or.mode.eq.1) then
      Do IT = 1, nAshT
        Do iU = 1, nAshT
          Do iX = 1, nAshT
            INT1(IT,IU) = INT1(IT,IU) - INT2(IT,IX,IX,IU)
          End Do
        End Do
      End Do
      endif
C
      Return
C
      End Subroutine CnstInt
C
C-----------------------------------------------------------------------
C
      !! dens2_rpt2.f
      Subroutine TimesE2(Mode,CIin,CIout,INT1,INT2)
      use gugx, only: SGS, L2ACT, CIS
      Implicit Real*8 (A-H,O-Z)

      Dimension CIin(nConf,nState),CIout(nConf,nState)
      Real*8    INT1(nAshT,nAshT),INT2(nAshT,nAshT,nAshT,nAshT)
      LOGICAL   RSV_TSK
      real(kind=wp),allocatable :: SGM1(:),SGM2(:)
      integer(kind=iwp),allocatable :: TASK(:,:)
      Integer :: nLev
      nLev=SGS%nLev
      ! logical tras,uras,vras,xras
C
C     --- H_{IJ}*P_J
C    <CI1|EtuEvx|CI2>=<CI1|Evx
C
      nTasks= nLev**2
      CALL mma_allocate (Task,nTasks,2,Label='TASK')
C
      iTask=0
      DO LT=1,nLev
        DO LU=1,nLev
          iTask=iTask+1
          TASK(iTask,1)=LT
          TASK(iTask,2)=LU
        ENDDO
      ENDDO
      IF (iTask.NE.nTasks) WRITE(6,*) "ERROR nTasks"
C
      call mma_allocate(SGM1,nConf,Label='SGM1')
      call mma_allocate(SGM2,nConf,Label='SGM2')
C
      Call DCopy_(nConf*nState,[0.0D+00],0,CIout,1)
      Do kState = 1, nState
        !! Start the actual part of dens2_rpt2
        Call Init_Tsk(ID, nTasks)
C
 500    If (.NOT.Rsv_Tsk (ID,iTask)) GOTO 501
C
        LT=TASK(iTask,1)
        ! tras=.false.
        ! if (lt.le.nras1(1)) tras=.true.
          IST=SGS%ISM(LT)
          IT=L2ACT(LT)
          LU=Task(iTask,2)
          ! uras=.false.
          ! if (lu.gt.nras1(1)+nras2(1)) uras=.true.
C         if (tras.and.uras) go to 500
            ! LTU=iTask
            ISU=SGS%ISM(LU)
            IU=L2ACT(LU)
            ISTU=MUL(IST,ISU)
            ISSG=MUL(ISTU,STSYM)
            NSGM=CIS%NCSF(ISSG)
            IF(NSGM.EQ.0) GOTO 500
            !! <CIin|Etu
            CALL GETSGM2(LU,LT,STSYM,CIin(1,kState),SGM1)
            IF(ISTU.EQ.1) THEN
              !! <CIin|Etu|CIout>*I1tu
              Call DaXpY_(NSGM,INT1(IT,IU),SGM1,1,CIout(1,kState),1)
            END IF
            LVX=0
            DO LV=1,NLEV
              ISV=SGS%ISM(LV)
              IV=L2ACT(LV)
              ! vras=.false.
              ! if (lv.le.nras1(1)) vras=.true.
              DO LX=1,NLEV
                LVX=LVX+1
                ISX=SGS%ISM(LX)
                ISVX=MUL(ISV,ISX)
                ! xras=.false.
                ! if (lx.gt.nras1(1)+nras2(1)) xras=.true.
C               if (vras.and.xras) go to 110
                IF(ISVX.NE.ISTU) GOTO 110
                IX=L2ACT(LX)
                CALL GETSGM2(LX,LV,ISSG,SGM1,SGM2)
                Call DaXpY_(NSGM,INT2(IT,IU,IV,IX),SGM2,1,
     *                      CIout(1,kState),1)
 110          CONTINUE
              END DO
            END DO

C
        GOTO 500
 501    CONTINUE
        CALL Free_Tsk(ID)
        !! End the actual part of dens2_rpt2
      End Do
C
      call mma_deallocate(Task)
C
#ifdef _MOLCAS_MPP_
      CALL GAdSUM(CIout,nConf*nState)
#endif
C
C     --- -E_{S}*CJ + zL_{KL}
C
      Do kState = 1, nState
        EC=Eact(kState)
        Call Daxpy_(nConf,EC,CIin(1,kState),1,CIout(1,kState),1)
      End Do
C
      !! Project out the reference vector, just in case
      If (Mode.eq.1) Then
        Do ilStat = 1, nState
          Call DCopy_(nConf,CIout(1,ilStat),1,SGM1,1)
          Do jlStat = 1, nState
            Call LoadCI_XMS('C',1,SGM2,jlStat,U0)
            Ovl = DDot_(nConf,SGM1,1,SGM2,1)
            Call DaXpY_(nConf,-Ovl,SGM2,1,CIout(1,ilStat),1)
          End Do
        End Do
      End If
C
      call mma_deallocate(SGM1)
      call mma_deallocate(SGM2)
C
      Call DScal_(nConf*nState,2.0D+00/nState,CIout,1)
C
      Return
C
      End Subroutine TimesE2
C
C-----------------------------------------------------------------------
C
      Subroutine CnstDEPSA(CI,CIT,G1,G2,INT2)
C
      use gugx, only: SGS
      use pt2_guga
      Implicit Real*8 (A-H,O-Z)
C
C
      Dimension CI(nConf,nState),CIT(nConf,nState),G1(nAshT,nAshT),
     *          G2(nAshT,nAshT,nAshT,nAshT)
      Real*8    INT2(nAshT,nAshT,nAshT,nAshT)
      real(kind=wp),allocatable :: SGM1(:),SGM2(:),G1T(:),G2T(:),
     *                             Fock(:),FockOut(:)
      Integer :: nLev
      nLev=SGS%nLev
C
C     LOGICAL   RSV_TSK
C
C     This subroutine computes the second term in Eq. (70) or the RHS of
C     Eq. (72) in the CASPT2-IPEA gradient paper
C     CIT corresponds to \overline{Q}, if I remember correctly
C
      Call DCopy_(nAshT**2,[0.0D+00],0,G1,1)
      Call DCopy_(nAshT**4,[0.0D+00],0,G2,1)
C
      !! Construct transition(?) density matrix
      !! (<CI|Etu|CIT>+<CIT|Etu|CI>)/2, where CIT is the solution
      call mma_allocate(SGM1,nConf,Label='SGM1')
      call mma_allocate(SGM2,nConf,Label='SGM2')
      call mma_allocate(G1T,NG1,Label='GT1')
      call mma_allocate(G2T,NG2,Label='GT2')
C
C  !! This is for CASSCF orbital Lagrangian, but this may not contribute
C     Call Dens2T_RPT2(CI(1,jState),CI(1,jState),
C    *                 SGM1,SGM2,G1T,G2T,nLev)
C     Call DaXpY_(NG1,-0.5D+00,G1T,1,G1,1)
C     Call DaXpY_(NG2,-0.5D+00,G2T,1,G2,1)
C
      Do kState = 1, nState
C       Wgt = DWgt(iState,iState)
        Wgt = 1.0D+00/nState
C
        !! <CI|Etu|CIT>+<CIT|Etu|CI> and the t+ u+ x v variant
        Call Dens2T_RPT2(CI(1,kState),CIT(1,kState),
     *                   SGM1,SGM2,G1T,G2T,nLev)
        Call DaXpY_(NG1,WGT,G1T,1,G1,1)
        Call DaXpY_(NG2,WGT,G2T,1,G2,1)
C
        !! For the orbital contribution of CASSCF Lagrangian
        !! Just add the SLag rotation contributions
        ilState = kState
        Do jlState = 1, ilState-1
C         If (ilState.eq.jlState) Cycle
          vSLag = -0.5D+00*SLag(ilState,jlState)
          If (abs(vSLag).le.1.0D-08) Cycle
          Call Dens2T_RPT2(CI(1,ilState),CI(1,jlState),
     *                     SGM1,SGM2,G1T,G2T,nLev)
          Call DaXpY_(NG1,vSLag,G1T,1,G1,1)
          Call DaXpY_(NG2,vSLag,G2T,1,G2,1)
        End Do
      End Do
C
      call mma_deallocate(SGM1)
      call mma_deallocate(SGM2)
      call mma_deallocate(G1T)
      call mma_deallocate(G2T)
C
      !! Finally, construct the Fock matrix only for active-active
      !! Should be equivalent to FockGen in MCLR
      call mma_allocate(Fock,nAshT**2,Label='Fock')
      Fock(:) = 0.0d+00
C
      !! 1) FIMO term
      Do iS=1,nSym
        If (nBas(iS).gt.0) Then
          jS=iEOr(is-1,iSym-1)+1
          Do iA=1,nAsh(is)
            Do jA=1,nAsh(js)
C             rd=rDens1(iA+nA(iS),jA+nA(js))
C             ip1=nBas(iS)*(nIsh(is)+iA-1)+ipCM(is)-1
C             ip2=nBas(iS)*(nIsh(js)+jA-1) +ipmat(is,js)
              rd=G1(iA,jA)
              ip1= 1+nFro(jS)+nIsh(jS)
     *           + nBas(iS)*(nFro(iS)+nIsh(iS)+iA-1)
              ip2=1+nAsh(iS)*(jA-1)
              Call DaXpY_(nAsh(iS),Rd,FIMO(ip1),1,Fock(ip2),1)
            End Do
          End Do
        End If
      End Do
C     write(6,*) "after 1"
C     call sqprt(fock,nasht)
C
      !! 2) two-electron term (only CreQADD part)
      Do iS=1,nSym
        ipS=iEOr(is-1,isym-1)+1
        if (norb(ips).ne.0) Then
          Do jS=1,nsym
            ijS=iEOR(is-1,js-1)+1
            Do kS=1,nSym
              ls=iEOr(ijs-1,iEor(ks-1,isym-1))+1
*                                                                      *
************************************************************************
*                                                                      *
               Do kAsh=1,nAsh(kS)
                kAA=kAsh+nFro(kS)+nIsh(kS)
                Do lAsh=1,nAsh(lS)
                  lAA=lAsh+nFro(lS)+nIsh(lS)
*
*                 Pick up (pj|kl)
*
                  Call Coul(ipS,jS,kS,lS,kAA,lAA,WRK1,WRK2)
*
                  Do iAsh=1,nAsh(iS)
                    ipQ=nAsh(ipS)*(iAsh-1)
                    Do jAsh=1,nAsh(jS)
                      ipM=nFro(ipS)+nIsh(ipS)
     *                   +(nFro(jS)+nIsh(jS)+jAsh-1)*nBas(ipS)
                      Call DaXpY_(nAsh(ipS),G2(iAsh,jAsh,kAsh,lAsh)*2,
     &                            INT2(1,jAsh,kAsh,lAsh),1,
     *                            Fock(1+ipQ),1)
                      ipM=ipM+nOrb(ipS)
*
                    End Do
                  End Do
*
                End Do
              End Do
*                                                                      *
************************************************************************
*                                                                      *
            End Do  ! kS
          End Do     ! jS
        End If
      End Do           ! iS
C
      !! 3) anti-symmetrize
      !! 4) Divide by the difference of orbital energies
      call mma_allocate(FockOut,nAshT**2,Label='FockOut')
      Do iS=1,nSym
        jS=iEOR(iS-1,iSym-1)+1
        If (nAsh(is)*nAsh(jS).ne.0) Then
          !! Anti-symmetrize
          Call DGeSub(Fock,nAsh(iS),'N',
     &                Fock,nAsh(jS),'T',
     &                FockOut,nAsh(iS),
     &                nAsh(iS),nAsh(jS))


          !! Divide
          imo=1
          Do iAsh = 1, nAsh(iSym)
            iOrb = iAsh + nFro(iSym) + nIsh(iSym)
            EigI = FIFA(iMO+iOrb-1+nBas(iSym)*(iOrb-1))
            Do jAsh = 1, iAsh-1
              jOrb = jAsh + nFro(iSym) + nIsh(iSym)
              EigJ = FIFA(iMO+jOrb-1+nBas(iSym)*(jOrb-1))
              OLagIJ = FockOut(iAsh+nAsh(iSym)*(jAsh-1))
              Tmp = OLagIJ/(EigI-EigJ)
              DEPSA(iAsh,jAsh) = DEPSA(iAsh,jAsh) + Tmp
              DEPSA(jAsh,iAsh) = DEPSA(jAsh,iAsh) + Tmp
            End Do
          End Do
        End If
      End Do
C
      call mma_deallocate(FockOut)
      call mma_deallocate(Fock)
C
      Return
C
      End Subroutine CnstDEPSA
C
C-----------------------------------------------------------------------
C
      !! PRWF1_CP2
      SUBROUTINE CnstPrec(NOCSF,IOCSF,NOW,IOW,ISYCI,PRE,ci,
     *                    INT1,INT2,Fancy,nLev,nMidV)
      use gugx, only: SGS, CIS
      use caspt2_module
      use pt2_guga
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER, INTENT(IN) :: nLev
      DIMENSION NOCSF(NSYM,NMIDV,NSYM),IOCSF(NSYM,NMIDV,NSYM)
      DIMENSION NOW(2,NSYM,NMIDV),IOW(2,NSYM,NMIDV)
      DIMENSION PRE(*)
      dimension ci(*)
      REAL*8 INT1(NLEV,NLEV),INT2(NLEV,NLEV,NLEV,NLEV)
      REAL*8 Fancy(nRoots,nRoots,nRoots)

      DIMENSION ICS(MXLEV)
      Integer :: nIpWlk
      nIpWlk = CIS%nIpWlk
C
C     Construct (approximate?) preconditioner for the active linear
C     equation that should be solved for non-invariant CASPT2 methods
C     (with IPEA shift)
C
C -- NOTE: THIS PRWF ROUTINE USES THE CONVENTION THAT CI BLOCKS
C -- ARE MATRICES CI(I,J), WHERE THE   F I R S T   INDEX I REFERS TO
C -- THE   U P P E R   PART OF THE WALK.

C SVC: set up a CSF string length as LENCSF
C     LINE=' '
      LENCSF=0
      ISY=0
      DO LEV=1,NLEV
        IF(ISY.NE.SGS%ISM(LEV)) THEN
          ISY=SGS%ISM(LEV)
          LENCSF=LENCSF+1
        END IF
        LENCSF=LENCSF+1
      END DO
      LENCSF=MIN(LENCSF,256)
      LENCSF=MAX(LENCSF,10)


C     SVC2010:
C     allocate scratch memory for determinant expansion
      IF (PRSD) THEN
C       CALL MMA_ALLOCATE(LEX,NLEV,LABEL='LEX')
      END IF

C     LINE=' '

C -- THE MAIN LOOP IS OVER BLOCKS OF THE ARRAY CI
C    WITH SPECIFIED MIDVERTEX MV, AND UPPERWALK SYMMETRY ISYUP.
      DO 40 MV=1,NMIDV
        DO 41 ISYUP=1,NSYM
          NCI=NOCSF(ISYUP,MV,ISYCI)
          IF(NCI.EQ.0) GOTO 41
          NUP=NOW(1,ISYUP,MV)
          ISYDWN=MUL(ISYUP,ISYCI)
          NDWN=NOW(2,ISYDWN,MV)
          ICONF=IOCSF(ISYUP,MV,ISYCI)
          IUW0=1-NIPWLK+IOW(1,ISYUP,MV)
          IDW0=1-NIPWLK+IOW(2,ISYDWN,MV)
          IDWNSV=0
          DO 30 IDWN=1,NDWN
            DO 31 IUP=1,NUP
              ICONF=ICONF+1
C             COEF=CI(ICONF)
C -- SKIP OR PRINT IT OUT?
C             IF(ABS(COEF).LT.THR) GOTO  31
              IF(IDWNSV.NE.IDWN) THEN
                ICDPOS=IDW0+IDWN*NIPWLK
                ICDWN=CIS%ICASE(ICDPOS)
C -- UNPACK LOWER WALK.
                NNN=0
                DO 10 LEV=1,SGS%MIDLEV
                  NNN=NNN+1
                  IF(NNN.EQ.16) THEN
                    NNN=1
                    ICDPOS=ICDPOS+1
                    ICDWN=CIS%ICASE(ICDPOS)
                  END IF
                  IC1=ICDWN/4
                  ICS(LEV)=ICDWN-4*IC1
                  ICDWN=IC1
  10            CONTINUE
                IDWNSV=IDWN
              END IF
              ICUPOS=IUW0+NIPWLK*IUP
              ICUP=CIS%ICASE(ICUPOS)
C -- UNPACK UPPER WALK:
              NNN=0
              DO LEV=SGS%MIDLEV+1,NLEV
                NNN=NNN+1
                IF(NNN.EQ.16) THEN
                  NNN=1
                  ICUPOS=ICUPOS+1
                  ICUP=CIS%ICASE(ICUPOS)
                END IF
                IC1=ICUP/4
                ICS(LEV)=ICUP-4*IC1
                ICUP=IC1
              END DO
C -- PRINT IT!
              K=0
              ISY=0
              PRE(ICONF) = 0.0D+00
              DO LEV=1,NLEV
                IF(ISY.NE.SGS%ISM(LEV)) THEN
                  ISY=SGS%ISM(LEV)
                  K=K+1
C                 LINE(K:K)=' '
                END IF
                K=K+1
C               LINE(K:K)=CODE(ICS(LEV))
                IF (ICS(LEV).EQ.0) THEN
                  VAL = 0.0D+00
                ELSE IF (ICS(LEV).EQ.3) THEN
                  VAL = 2.0D+00*INT1(LEV,LEV)
C                 L=0
C                 JSY=0
                  DO LEV2=1,NLEV
                    IF (ICS(LEV2).EQ.0) THEN
                    val = val + 0.0d+00
                    ELSE IF ((LEV.EQ.LEV2.AND.ICS(LEV2).EQ.3).OR.
     *                       (LEV.NE.LEV2.AND.ICS(LEV2).EQ.1).OR.
     *                       (LEV.NE.LEV2.AND.ICS(LEV2).EQ.2)) THEN
                      val2 =  4.0d+00*int2(lev,lev ,lev2,lev2)
     *                      - 2.0d+00*int2(lev,lev2,lev ,lev2)
                      val = val + val2
                    ELSE IF (LEV.NE.LEV2.AND.ICS(LEV2).EQ.3) THEN
                      val2 =  4.0d+00*int2(lev,lev ,lev2,lev2)
     *                      - 2.0d+00*int2(lev,lev2,lev ,lev2)
                      val = val + val2
                    END IF
                  END DO
                ELSE
                  VAL = INT1(LEV,LEV)
C                 L=0
C                 JSY=0
                  DO LEV2=1,NLEV
                    IF (ICS(LEV2).EQ.0.OR.LEV.EQ.LEV2) THEN
                    ELSE IF (ICS(LEV2).EQ.3) THEN
C      val2 =  2.0d+00*int2(lev,lev ,lev2,lev2)
C    *       - 1.0d+00*int2(lev,lev2,lev ,lev2)
C      val2 = 0.0d+00
C      val = val + val2*0.5d+00
                    ELSE
                      val2 = int2(lev,lev ,lev2,lev2)
     *                     + int2(lev,lev2,lev ,lev2)
                      if (ics(lev).eq.ics(lev2)) then
                      val2 = int2(lev,lev ,lev2,lev2)
     *                     - int2(lev,lev2,lev ,lev2)
                      end if
                      val = val + val2
                    END IF
                  END DO
                END IF
                PRE(ICONF) = PRE(ICONF) + VAL
              END DO
C     SVC2010 experimental: add determinant expansion
C             IF (PRSD) THEN
c     Specify projected spin in half integer units
C     Default: use maximum spin projection
C              IMS = ISPIN-1
C              WRITE(6,*)
C              CALL EXPCSF (ICS, NLEV, IMS, LEX)
C              WRITE(6,*)
C             ENDIF
  31        CONTINUE
  30      CONTINUE
  41    CONTINUE
  40  CONTINUE
C
      !! mclr/sa_prec.f
      !! Prepare so-called fancy preconditioner
      Do iSt = 1, nRoots
        Ene = Eact(iSt)
        Do jSt = 1, nRoots
          Do kSt = 1, nRoots
            Fancy(jSt,kSt,iSt) = 0.0D+00
            Do iConf = 1, nConf
              dnum=PRE(iConf)+Ene
              dnum=Sign(Max(Abs(dnum),1.0d-16),dnum)
              Fancy(jSt,kSt,iSt) = Fancy(jSt,kSt,iSt)
     *          + CI(iConf+nConf*(jSt-1))*CI(iConf+nConf*(kSt-1))/dnum
            End Do
          End Do
        End Do
        Call MatInvert(Fancy(1,1,iSt),nRoots)
      End Do
C
C     SVC2010: free scratch for determinant expansion
      IF (PRSD) THEN
C       CALL MMA_DEALLOCATE(LEX)
      END IF
C     WRITE(6,*)
      RETURN
      END SUBROUTINE CnstPrec
C
C-----------------------------------------------------------------------
C
      !! mclr/dminvci_sa.f
      Subroutine DoPrec(VecIN,VecOUT,CI,Pre,Fancy)
C
C     Apply precondition to CI vectors, taken from the MCLR module
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension VecIN(*),VecOUT(*),CI(*),Pre(*),
     *          Fancy(nRoots,nRoots,nRoots)
      Dimension rcoeff(nRoots),alpha(nRoots)
C
      !! Standard inverse of the diagonal elements
      Do iRoots = 1, nRoots
        Do iConf = 1, nConf
          VecOUT(iConf+nConf*(iRoots-1))
     *      = VecIN(iConf+nConf*(iRoots-1))/(Pre(iConf)+Eact(iRoots))
        End Do
      End Do
C
      !! Construct reference CI vectors
      Do iRoots = 1, nRoots
        Call LoadCI_XMS('C',1,CI(1+nConf*(iRoots-1)),iRoots,U0)
      End Do
C
      !! The so-called fancy precondioner
      Do iRoots = 1, nRoots
        Do jRoots = 1, nRoots
          rcoeff(jRoots) = DDot_(nconf,VecOUT(1+nConf*(iRoots-1)),1,
     *                                 CI(1+nConf*(jRoots-1)),1)
        End Do
C
        Do jRoots = 1, nRoots
          alpha(jRoots) = 0.0D+00
          Do kRoots = 1, nRoots
            alpha(jRoots) = alpha(jRoots)
     *        + Fancy(jRoots,kRoots,iRoots)*rcoeff(kRoots)
          End Do
        End Do
C
        Do jRoots = 1, nRoots
          Do iConf = 1, nConf
            VecOUT(iConf+nConf*(iRoots-1))
     *        = VecOUT(iConf+nConf*(iRoots-1))
     *        - CI(iConf+nConf*(jRoots-1))*alpha(jRoots)
     *          /(Pre(iConf)+Eact(iRoots))
          End Do
        End Do
      End Do
C
      End Subroutine DoPrec
C
      End Subroutine DEPSAOffC
C
C-----------------------------------------------------------------------
C
      Subroutine DEPSAOffO(OLag,DEPSA,FIFA)
      use caspt2_module
C
      Implicit Real*8 (A-H,O-Z)
C
C
      Dimension OLag(*),DEPSA(nAshT,nAshT),FIFA(*)
C
C     This is much easier; similar to the frozen core approximation.
C     Corresponds to the first term in Eq. (70)
C
      iMO  = 1
      DO iSym = 1, nSym
        nAshI = nAsh(iSym)
        nOrbI = 0
        If (nAshI.ne.0) Then
          nOrbI = nBas(iSym)-nDel(iSym)
          nFroI = nFro(iSym)
          nIshI = nIsh(iSym)
          nBasI = nBas(iSym)
C
          Do iAsh = 1, nAshI
            iOrb = iAsh + nFroI + nIshI
            EigI = FIFA(iMO+iOrb-1+nBasI*(iOrb-1))
            Do jAsh = 1, iAsh-1
              jOrb = jAsh + nFroI + nIshI
              EigJ = FIFA(iMO+jOrb-1+nBasI*(jOrb-1))
              OLagIJ = OLag(iMO+iOrb-1+nOrbI*(jOrb-1))
              OLagJI = OLag(iMO+jOrb-1+nOrbI*(iOrb-1))
              Tmp = -(OLagIJ-OLagJI)/(EigI-EigJ)*0.5D+00
C
              DEPSA(iAsh,jAsh) = DEPSA(iAsh,jAsh) + Tmp
              DEPSA(jAsh,iAsh) = DEPSA(jAsh,iAsh) + Tmp
            End Do
          End Do
        End If
        iMO  = iMO  + nOrbI*nOrbI
      End Do
C
      End Subroutine DEPSAOffO
C
C-----------------------------------------------------------------------
C
      Subroutine LinDepLag(BDer,SDer,nAS,nIN,iSym,iCase)
C
      use caspt2_global, only: LUSTD, idBoriMat
      use caspt2_global, only: LUSBT
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp
      use caspt2_module
C
      Implicit Real*8 (A-H,O-Z)
C
C
      DIMENSION WGRONK(2)
      Dimension BDer(*),SDer(*)
      real(kind=wp),allocatable :: S(:),SS(:,:),VEC(:,:),EIG(:),SCA(:),
     *                             SCRATCH(:),LAG(:,:),B(:),F(:,:)
C
C     Compute contributions that arise from the non-invariance effect
C     in non-orthogonal -> orthogonal ICB rotations
C     See J. Chem. Phys. 2023, 158, 174112. for details, in particular,
C     Section II C 3 "Non-invariance with respect to orthogonal..."
C
      !! Obtain the X matrix
      !! First, read S
      NS = NAS*(NAS+1)/2
      call mma_allocate(S,NS,Label='S')
      call mma_allocate(SS,NAS,NAS,Label='SS')
      idS = idSMAT(iSym,iCase)
      CALL DDAFILE(LUSBT,2,S,NS,idS)
      IJ=0
      DO J=1,NAS
        DO I=1,J
          IJ=IJ+1
          SS(I,J)=S(IJ)
          SS(J,I)=S(IJ)
        END DO
      END DO
C
      call mma_allocate(VEC,NAS,NAS,Label='VEC')
      call mma_allocate(EIG,NAS,Label='EIG')
      call mma_allocate(SCA,NAS,Label='SCA')
      IDIAG=0
      DO I=1,NAS
        IDIAG=IDIAG+I
        SD=S(IDIAG)
        If (IFDORTHO) then
          SCA(I)=1.0D+00
        Else
          IF(SD.GT.THRSHN) THEN
* Small variations of the scale factor were beneficial
              SCA(I)=(1.0D00+DBLE(I)*3.0D-6)/SQRT(SD)
          ELSE
            SCA(I)=0.0D0
          END IF
        End If
      END DO
      IJ=0
      DO J=1,NAS
        DO I=1,J
          IJ=IJ+1
          S(IJ)=S(IJ)*SCA(I)*SCA(J)
        END DO
      END DO
C
      IJ=0
      DO J=1,NAS
        DO I=1,J
          IJ=IJ+1
          VEC(I,J)=S(IJ)
        END DO
      END DO
      INFO=0
      call dsyev_('V','L',NAS,VEC,NAS,EIG,WGRONK,-1,INFO)
      NSCRATCH=INT(WGRONK(1))
      call mma_allocate(SCRATCH,NSCRATCH,Label='SCRATCH')
      call dsyev_('V','U',NAS,VEC,NAS,EIG,SCRATCH,NSCRATCH,INFO)
      call mma_deallocate(SCRATCH)
C
      DO I=1,NAS
        SCAL=SCA(I)
C       CALL DSCAL_(NIN,SCA,VEC(1,I),NAS)
        CALL DSCAL_(NAS,SCAL,VEC(I,1),NAS)
      END DO
      call mma_deallocate(SCA)
      call mma_deallocate(S)
C
      !! Scale only the independent vectors to avoid
      !! any numerically unstable computation
      DO I=1,NAS
        EVAL=EIG(I)
        IF(EVAL.LT.THRSHS) CYCLE
        FACT=1.0D00/SQRT(EVAL)
        Call DScal_(nAS,FACT,VEC(1,I),1)
      END DO
C
      call mma_allocate(LAG,NAS,NAS,Label='LAG')
      IDB=IDBoriMat(ISYM,ICASE)
      NB=NS
      call mma_allocate(B,NB,Label='B')
      CALL DDAFILE(LUSTD,2,B,NB,IDB)
      call mma_allocate(F,NAS,NAS,Label='F')
      IJ=0
      DO J=1,NAS
        DO I=1,J
          IJ=IJ+1
          F(I,J)=B(IJ)
          F(J,I)=B(IJ)
        END DO
      END DO
C
      !! Compute the partial derivative
      !! F   : B
      !! BDER: D
      !! VEC : X^0 and X
      Call DGEMM_('N','T',NAS,NAS,NAS,
     *            2.0D+00,F,NAS,BDER,NAS,
     *            0.0D+00,LAG,NAS)
      Call DGEMM_('N','N',NAS,NAS,NAS,
     *            1.0D+00,LAG,NAS,VEC,NAS,
     *            0.0D+00,F,NAS)
      Call DCopy_(NAS*NAS,F,1,LAG,1)
C
      CALL DGEMM_('T','N',NAS,NAS,NAS,
     *            1.0D+00,VEC,NAS,LAG,NAS,
     *            0.0D+00,F,NAS)
      !! At this point,
      !! F = 2 \mathcal{X}^0 * B * D * \mathcal{X}
C
      !! remove dependent part
      !! (linearly indep-indep and dep-dep)
      Do I = 1, nAS-nIN
        Do J = 1, nAS-nIN
          F(I,J) = 0.0D+00
        End Do
      End Do
      Do I = nAS-nIN+1, nAS
        Do J = nAS-nIN+1, nAS
          F(I,J) = 0.0D+00
        End Do
      End Do
C
      !! orthogonal -> non-orthogonal
      !! Finalize Eq. (62)
      CALL DGEMM_('N','N',NAS,NAS,NAS,
     *            1.0D+00,VEC,NAS,F,NAS,
     *            0.0D+00,LAG,NAS)
      CALL DGEMM_('N','T',NAS,NAS,NAS,
     *            1.0D+00,LAG,NAS,VEC,NAS,
     *            0.0D+00,F,NAS)
C
      Call DaXpY_(nAS*nAS,1.0D+00,F,1,SDER,1)
C
      call mma_deallocate(LAG)
      call mma_deallocate(B)
      call mma_deallocate(F)
C
      call mma_deallocate(SS)
      call mma_deallocate(EIG)
      call mma_deallocate(VEC)
C
      Return
C
      End Subroutine LinDepLag
C
C-----------------------------------------------------------------------
C
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      Subroutine LinDepLag_MPP(lg_BDER,lg_SDER,nAS,nIN,iSym,iCase)
C
      use caspt2_global, only: LUSTD, idBoriMat
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp
      use caspt2_module
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
C
      Implicit Real*8 (A-H,O-Z)
C

#include "global.fh"
#include "mafdecls.fh"
C
#if ! defined (_SCALAPACK_)
      DIMENSION WGRONK(2)
#endif
      logical :: bStat
      real(kind=wp),allocatable :: EIG(:),VEC(:),SCRATCH(:)
C
C     Parallel LinDepLag
C     We always use the canonical orthonormalization.
C
      !! Obtain the X matrix
      !! First, read S
      CALL PSBMAT_GETMEM ('S',lg_S,NAS)
      CALL PSBMAT_READ ('S',iCase,iSym,lg_S,NAS)

      call mma_allocate(EIG,NAS,Label='EIG')
      EIG(:) = 0.0d+00

      myRank = GA_NodeID()
#ifdef _SCALAPACK_
      CALL PSBMAT_GETMEM('VMAT',lg_Vec,NAS)
      CALL GA_PDSYEVX_ (lg_S, lg_Vec, EIG, 0)
      bSTAT = GA_Destroy (lg_S)
#else
C here for the non-ScaLAPACK version: copy matrix to master process,
C diagonalize using the serial DSYEV routine, and copy the resulting
C eigenvectors back to a global array.  Then distribute the eigenvalues.
      IF (myRank.EQ.0) THEN
        call mma_allocate(VEC,NAS*NAS,Label='VEC')
        CALL GA_Get (lg_S, 1, NAS, 1, NAS, VEC, NAS)
      END IF
      bSTAT = GA_Destroy (lg_S)
      IF (myRank.EQ.0) THEN
        CALL DSYEV_('V','L',NAS,VEC,NAS,EIG,WGRONK,-1,INFO)
        NSCRATCH=INT(WGRONK(1))
        call mma_allocate(SCRATCH,NSCRATCH,Label='SCRATCH')
        CALL DSYEV_('V','L',NAS,VEC,NAS,EIG,SCRATCH,NSCRATCH,INFO)
        call mma_deallocate(SCRATCH)
      END IF
      CALL PSBMAT_GETMEM('VMAT',lg_Vec,NAS)
      IF (myRank.EQ.0) THEN
        CALL GA_Put (lg_Vec, 1, NAS, 1, NAS, VEC, NAS)
        call mma_deallocate(VEC)
      END IF
      CALL GADSUM(EIG,NAS)
#endif
C
      !! Scale only the independent vectors to avoid
      !! any numerically unstable computation
C Form orthonormal transformation vectors by scaling the eigenvectors.
      call GA_Sync()
      call GA_Distribution (lg_Vec, myRank, iLo, iHi, jLo, jHi)
      IF (iLo.NE.0) THEN
        call GA_Access (lg_Vec, iLo, iHi, jLo, jHi, mV, LDV)
        IF ((jHi-jLo+1).NE.NAS) THEN
          WRITE(6,*) 'SBDIAG_MPP: error in striping of lg_Vec, ABORT'
          CALL ABEND()
        END IF
        DO I=1,jHi-jLo+1 ! NAS
          EVAL=EIG(I)
          IF(EVAL.LT.THRSHS) CYCLE
          FACT=1.0D00/SQRT(EVAL)
          Call DScal_(iHi-iLo+1,FACT,DBL_MB(mV+LDV*(I-1)),1)
        END DO
        call GA_Release_Update (lg_Vec, iLo, iHi, jLo, jHi)
      END IF
      call GA_Sync()

      call mma_deallocate(EIG)
C
C     X matrix has been prpeared
C
      CALL GA_CREATE_STRIPED ('H',NAS,NAS,'Lag',lg_Lag)
      CALL PSBMAT_GETMEM ('B',lg_B,NAS)

      call GA_Distribution (lg_B, myRank, iLo, iHi, jLo, jHi)
      call GA_Access (lg_B, iLo, iHi, jLo, jHi, mB, LDB)
      IDB=IDBoriMat(ISYM,ICASE)
      CALL DDAFILE(LUSTD,2,DBL_MB(mB),(iHi-iLo+1)*(jHi-jLo+1),IDB)
      call GA_Release_Update (lg_B, iLo, iHi, jLo, jHi)
C
      !! Compute the partial derivative
      !! Work(LF)  : B --> lg_B
      !! BDER      : D --> lg_BDER
      !! Work(LVEC): X^0 and X --> lg_Vec
      Call GA_DGEMM ('N','T',NAS,NAS,NAS,
     *               2.0d+00,lg_B,lg_BDER,
     *               0.0D+00,lg_Lag)
      Call GA_DGEMM ('N','N',NAS,NAS,NAS,
     *               1.0d+00,lg_Lag,lg_VEC,
     *               0.0D+00,lg_B)
      Call GADupl(lg_B,lg_Lag)
C
      Call GA_DGEMM ('T','N',NAS,NAS,NAS,
     *               1.0d+00,lg_Vec,lg_Lag,
     *               0.0D+00,lg_B)
      !! At this point,
      !! Work(LF) = 2 \mathcal{X}^0 * B * D * \mathcal{X}
C
      !! remove dependent part
      !! (linearly indep-indep and dep-dep)
      call GA_Distribution (lg_B, myRank, iLo, iHi, jLo, jHi)
      if (iLo.ne.0) then
        call GA_Access (lg_B, iLo, iHi, jLo, jHi, mV, LDV)
C       if (ilo.le.nas-nin) then
          Do J = 1, nAS-nIN
            Do I = 1, min(iHi-iLo+1,nAS-nIN-iLo+1)
              DBL_MB(mV+i-1+LDV*(j-1)) = 0.0d+00
            End Do
          End Do
C       end if
C       if (ilo.ge.nas-nin+1) then
          Do J = nAS-nIN+1, nAS
            Do I = max(1,nAS-nIN+1-iLo+1), LDV
              DBL_MB(mV+i-1+LDV*(j-1)) = 0.0d+00
            End Do
          End Do
C       end if
        call GA_Release_Update (lg_B, iLo, iHi, jLo, jHi)
      end if
C
      !! orthogonal -> non-orthogonal
      !! Finalize Eq. (62)
      Call GA_DGEMM ('N','N',NAS,NAS,NAS,
     *               1.0d+00,lg_Vec,lg_B,
     *               0.0D+00,lg_Lag)
      Call GA_DGEMM ('N','T',NAS,NAS,NAS,
     *               1.0d+00,lg_Lag,lg_Vec,
     *               1.0D+00,lg_SDER)

      CALL PSBMAT_FREEMEM(lg_Vec)
      bSTAT = GA_Destroy (lg_Lag)
      CALL PSBMAT_FREEMEM (lg_B)
C
      Return
C
      End Subroutine LinDepLag_MPP
#endif
