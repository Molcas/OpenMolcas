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
! SavGradParams : state-dependent quantities
! SavGradParams2: state-independent quantities
!
Subroutine SavGradParams(Mode,IDSAVGRD)

  use caspt2_global, only: LUGRAD, LUSTD, do_lindep, IDBoriMat, &
                             NBUF1_GRAD, iTasks_grad, nTasks_grad
#ifdef _MOLCAS_MPP_
  use caspt2_global, only: LURHS
  use caspt2_module, only: IOFFRHS
#endif
  use caspt2_global, only: DREF, PREF
  use caspt2_global, only: LUSOLV, LUSBT
  use definitions, only: iwp,wp,byte
  use stdalloc, only: mma_allocate, mma_deallocate
  use EQSOLV, only: IDSMAT, IDBMAT, IDSTMAT, IVECX, IDTMAT
  use fake_GA, only: GA_Arrays
#ifdef _MOLCAS_MPP_
  USE Para_Info, ONLY: Is_Real_Par, King, myRank
      use pt2_guga, only: iAdr10, cLab10
#endif

  use caspt2_module, only: E2Tot, EASum, ERef, jState, MxCase, nAshT, nBTri, nState, nSym, RFPert, nCases, &
                           nInDep, nISup, nASup, RefEne

      use pt2_guga, only: nG1, nG2, nG3, nG3Tot
  Implicit None

#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
logical(kind=iwp) bStat
#endif

  integer(kind=iwp), intent(in) :: Mode
  integer(kind=iwp), intent(inout) :: IDSAVGRD

  integer(kind=iwp) :: IORW,ID,NIN,NAS,NIS,NNN,NMAX,ISYM,ICASE,iLUID
#ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: I,lg_ST,lg_S,lg_T,ISTA,IEND,JSTA,JEND,mV1,LDW,IDISK,LDM,NBLOCK
#endif

  real(kind=wp), allocatable :: WRK1(:)
  integer(kind=iwp), allocatable :: IWRK1(:)
  integer(kind=byte), allocatable :: idxG3(:,:)

! character(len=80) :: Label

  !! Shift the address due to SavGradParams2
  If (IDSAVGRD == 0) Then
    IDSAVGRD = NSTATE + 3*NSTATE**2
    if (RFpert) IDSAVGRD = IDSAVGRD + NBTRI + 1
  End If
#ifdef _MOLCAS_MPP_
  I = 0
  myRank = GA_NODEID()
#endif

  !! Decide what to do
  If (Mode == 1) Then
    IORW = 1 !! Write
  Else If (Mode == 2) Then
    IORW = 2 !! Read
  End If

  !! Save internal contractions-related quantities
  !! 1. Some integers
  !! - Number of independent vectors
  CALL IDAFILE(LUGRAD,IORW,NINDEP(1,1),8*MXCASE,IDSAVGRD)
  !! - Number of active indices
  CALL IDAFILE(LUGRAD,IORW,NASUP(1,1),8*MXCASE,IDSAVGRD)
  !! - Number of inactive + secondary indices
  CALL IDAFILE(LUGRAD,IORW,NISUP(1,1),8*MXCASE,IDSAVGRD)

  !! 2. RDMs
  !! - NG1, NG2, NG3
  Call mma_allocate(IWRK1,6,Label='IWRK1')
  If (IORW == 1) Then
    IWRK1(1) = NG1
    IWRK1(2) = NG2
    IWRK1(3) = NG3
    IWRK1(4) = NG3TOT
    IWRK1(5) = NBUF1_GRAD
    IWRK1(6) = nTasks_grad
    CALL IDAFILE(LUGRAD,IORW,IWRK1,6,IDSAVGRD)
  Else If (IORW == 2) Then
    CALL IDAFILE(LUGRAD,IORW,IWRK1,6,IDSAVGRD)
    NG1    = IWRK1(1)
    NG2    = IWRK1(2)
    NG3    = IWRK1(3)
    NG3TOT = IWRK1(4)
    NBUF1_GRAD = IWRK1(5)
    nTasks_grad= IWRK1(6)
    iTasks_grad(:) = 0
  End If
  Call mma_deallocate(IWRK1)
  CALL IDAFILE(LUGRAD,IORW,iTasks_grad,NASHT**2,IDSAVGRD)

  NMAX = 0

  Do ISYM = 1, NSYM
    Do ICASE = 1, 11
      NIN = NINDEP(ISYM,ICASE)
      NAS = NASUP(ISYM,ICASE)
      NIS = NISUP(ISYM,ICASE)
      NNN = NAS*(NAS+1)/2
      NNN = MAX(NNN,NAS*NIN)
      NNN = MAX(NNN,NIS)
      NMAX= MAX(NNN,NMAX)
    End Do
    Do ICASE = 12, 13
      NAS = NASUP(ISYM,ICASE)
      NIS = NISUP(ISYM,ICASE)
      NMAX= MAX(NAS*NIS,NMAX)
    End Do
  End Do
  NMAX = MAX(NMAX,NG3)

  Call mma_allocate(WRK1,NMAX,Label='WRK1')

  CALL mma_allocate(idxG3,6,NG3,label='idxG3')
  If (IORW == 1) Then
    !! NG3 index
    iLUID = 0
    CALL I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)
    CALL I1DAFILE(LUGRAD,1,idxG3,6*NG3,IDSAVGRD)

    !! D and F
    CALL PT2_GET(NG1,' GAMMA1',WRK1)
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG1,IDSAVGRD)
    CALL PT2_GET(NG2,' GAMMA2',WRK1)
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG2,IDSAVGRD)
    CALL PT2_GET(NG3,' GAMMA3',WRK1)
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG3,IDSAVGRD)

    CALL PT2_GET(NG1,' DELTA1',WRK1)
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG1,IDSAVGRD)
    CALL PT2_GET(NG2,' DELTA2',WRK1)
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG2,IDSAVGRD)
    CALL PT2_GET(NG3,' DELTA3',WRK1)
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG3,IDSAVGRD)

    !! EASUM
    WRK1(1) = EASUM
    CALL DDAFILE(LUGRAD,IORW,WRK1,1,IDSAVGRD)
  Else If (IORW == 2) Then
#ifdef _MOLCAS_MPP_
    if (is_real_par()) then
      !! Reset, because NG3 can be different
      !! STINI has been skipped
      DO I=1,64
        IADR10(I,1)=-1
        IADR10(I,2)=0
        CLAB10(I)='   EMPTY'
      END DO
      IADR10(1,1)=0
    end if
#endif
    !! NG3 index
    CALL I1DAFILE(LUGRAD,2,idxG3,6*NG3,IDSAVGRD)
    iLUID = 0
    CALL I1DAFILE(LUSOLV,1,idxG3,6*NG3,iLUID)

    !! D and F
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG1,IDSAVGRD)
    CALL PT2_PUT(NG1,' GAMMA1',WRK1)
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG2,IDSAVGRD)
    CALL PT2_PUT(NG2,' GAMMA2',WRK1)
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG3,IDSAVGRD)
    CALL PT2_PUT(NG3,' GAMMA3',WRK1)

    CALL DDAFILE(LUGRAD,IORW,WRK1,NG1,IDSAVGRD)
    CALL PT2_PUT(NG1,' DELTA1',WRK1)
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG2,IDSAVGRD)
    CALL PT2_PUT(NG2,' DELTA2',WRK1)
    CALL DDAFILE(LUGRAD,IORW,WRK1,NG3,IDSAVGRD)
    CALL PT2_PUT(NG3,' DELTA3',WRK1)

    !! EASUM
    CALL DDAFILE(LUGRAD,IORW,WRK1,1,IDSAVGRD)
    EASUM = WRK1(1)
    CALL GETDPREF(DREF,SIZE(DREF),PREF,SIZE(PREF))
    EREF=REFENE(JSTATE)
  End If
  Call mma_deallocate(idxG3)

  !! 3. LUSBT matrices
  Do ISYM = 1, NSYM
    Do ICASE = 1, 11
      NIN = NINDEP(ISYM,ICASE)
      NAS = NASUP(ISYM,ICASE)
      NIS = NISUP(ISYM,ICASE)
      If (IORW == 1) Then
        !! Active overlap
#ifdef _MOLCAS_MPP_
        if (is_real_par() .and. (icase==1 .or. icase==4)) then
          CALL PSBMAT_GETMEM('S',lg_S,NAS)
          CALL PSBMAT_READ('S',iCase,iSym,lg_S,NAS)
          CALL GA_Distribution (lg_S,myRank,ISTA,IEND,JSTA,JEND)
          IF (ISTA.GT.0 .AND. JSTA.GT.0) THEN
            CALL GA_Access (lg_S,ISTA,IEND,JSTA,JEND,mV1,LDM)
            NBLOCK=LDM*(JEND-JSTA+1)
            CALL DDAFILE(LUGRAD,1,DBL_MB(mV1),NBLOCK,IDSAVGRD)
            CALL GA_Release (lg_S,ISTA,IEND,JSTA,JEND)
          END IF
          CALL PSBMAT_FREEMEM(lg_S)
        else
#endif
          if (NAS > 0) then
            ID = IDSMAT(ISYM,ICASE)
            CALL DDAFILE(LUSBT,2,WRK1,NAS*(NAS+1)/2,ID)
            CALL DDAFILE(LUGRAD,1,WRK1,NAS*(NAS+1)/2,IDSAVGRD)
          end if
#ifdef _MOLCAS_MPP_
        end if
#endif
        !! ST matrix
#ifdef _MOLCAS_MPP_
        if (is_real_par() .and. (icase==1 .or. icase==4)) then
          CALL GA_CREATE_STRIPED ('H',NAS,NIN,'STMAT',lg_ST)
          CALL PSBMAT_READ ('M',iCase,iSym,lg_ST,NAS*NIN)
          CALL GA_Distribution (lg_ST,myRank,ISTA,IEND,JSTA,JEND)
          IF (ISTA.GT.0 .AND. JSTA.GT.0) THEN
            CALL GA_Access (lg_ST,ISTA,IEND,JSTA,JEND,mV1,LDM)
            NBLOCK=LDM*(JEND-JSTA+1)
            CALL DDAFILE(LUGRAD,1,DBL_MB(mV1),NBLOCK,IDSAVGRD)
            CALL GA_Release (lg_ST,ISTA,IEND,JSTA,JEND)
          END IF
          bStat = GA_Destroy (lg_ST)
        else
#endif
          if (NAS*NIN > 0) then
            ID = IDSTMAT(ISYM,ICASE)
            CALL DDAFILE(LUSBT,2,WRK1,NAS*NIN,ID)
            CALL DDAFILE(LUGRAD,1,WRK1,NAS*NIN,IDSAVGRD)
          end if
#ifdef _MOLCAS_MPP_
        end if
#endif
        !! Transformation matrix (eigenvector)
#ifdef _MOLCAS_MPP_
        if (is_real_par() .and. (icase==1 .or. icase==4)) then
          CALL GA_CREATE_STRIPED ('H',NAS,NIN,'TMAT',lg_T)
          CALL PSBMAT_READ ('T',iCase,iSym,lg_T,NAS*NIN)
          CALL GA_Distribution (lg_T,myRank,ISTA,IEND,JSTA,JEND)
          IF (ISTA.GT.0 .AND. JSTA.GT.0) THEN
            CALL GA_Access (lg_T,ISTA,IEND,JSTA,JEND,mV1,LDM)
            NBLOCK=LDM*(JEND-JSTA+1)
            CALL DDAFILE(LUGRAD,1,DBL_MB(mV1),NBLOCK,IDSAVGRD)
            CALL GA_Release (lg_T,ISTA,IEND,JSTA,JEND)
          END IF
          bStat = GA_Destroy (lg_T)
        else
#endif
          if (NAS*NIN > 0) then
            ID = IDTMAT(ISYM,ICASE)
            CALL DDAFILE(LUSBT,2,WRK1,NAS*NIN,ID)
            CALL DDAFILE(LUGRAD,1,WRK1,NAS*NIN,IDSAVGRD)
          end if
#ifdef _MOLCAS_MPP_
        end if
#endif
        !! Eigenvalue
        if (NIN > 0) then
          ID  = IDBMAT(ISYM,ICASE)
          CALL DDAFILE(LUSBT,2,WRK1,NIN,ID)
          CALL DDAFILE(LUGRAD,1,WRK1,NIN,IDSAVGRD)
        end if
        !! IS
!       CALL DDAFILE(LUSBT,2,WRK1,NIS,ID)
!       CALL DDAFILE(LUGRAD,1,WRK1,NIS,IDSAVGRD)
        if (do_lindep .and. NAS > 0) then
          ID  = IDBoriMat(ISYM,ICASE)
          CALL DDAFILE(LUSTD,2,WRK1,NAS*(NAS+1)/2,ID)
          CALL DDAFILE(LUGRAD,1,WRK1,NAS*(NAS+1)/2,IDSAVGRD)
        end if
      Else If (IORW == 2) Then
        !! Active overlap
#ifdef _MOLCAS_MPP_
        if (is_real_par() .and. (icase==1 .or. icase==4)) then
          CALL PSBMAT_GETMEM('S',lg_S,NAS)
          CALL GA_Distribution (lg_S,myRank,ISTA,IEND,JSTA,JEND)
          IF (ISTA.GT.0 .AND. JSTA.GT.0) THEN
            CALL GA_Access (lg_S,ISTA,IEND,JSTA,JEND,mV1,LDM)
            NBLOCK=LDM*(JEND-JSTA+1)
            CALL DDAFILE(LUGRAD,2,DBL_MB(mV1),NBLOCK,IDSAVGRD)
            CALL GA_Release (lg_S,ISTA,IEND,JSTA,JEND)
          END IF
          CALL PSBMAT_WRITE('S',iCase,iSym,lg_S,NAS)
          CALL PSBMAT_FREEMEM(lg_S)
        else
#endif
          if (NAS > 0) then
            CALL DDAFILE(LUGRAD,2,WRK1,NAS*(NAS+1)/2,IDSAVGRD)
            ID = IDSMAT(ISYM,ICASE)
            CALL DDAFILE(LUSBT,1,WRK1,NAS*(NAS+1)/2,ID)
          end if
#ifdef _MOLCAS_MPP_
        end if
#endif
        !! ST matrix
#ifdef _MOLCAS_MPP_
        if (is_real_par() .and. (icase==1 .or. icase==4)) then
          CALL GA_CREATE_STRIPED ('H',NAS,NIN,'STMAT',lg_ST)
          CALL GA_Distribution (lg_ST,myRank,ISTA,IEND,JSTA,JEND)
          IF (ISTA.GT.0 .AND. JSTA.GT.0) THEN
            CALL GA_Access (lg_ST,ISTA,IEND,JSTA,JEND,mV1,LDM)
            NBLOCK=LDM*(JEND-JSTA+1)
            CALL DDAFILE(LUGRAD,2,DBL_MB(mV1),NBLOCK,IDSAVGRD)
            CALL GA_Release (lg_ST,ISTA,IEND,JSTA,JEND)
          END IF
          CALL PSBMAT_WRITE ('M',iCase,iSym,lg_ST,NAS*NIN)
          bStat = GA_Destroy (lg_ST)
        else
#endif
          if (NAS*NIN > 0) then
            CALL DDAFILE(LUGRAD,2,WRK1,NAS*NIN,IDSAVGRD)
            ID = IDSTMAT(ISYM,ICASE)
            CALL DDAFILE(LUSBT,1,WRK1,NAS*NIN,ID)
          end if
#ifdef _MOLCAS_MPP_
        end if
#endif
        !! Transformation matrix (eigenvector)
#ifdef _MOLCAS_MPP_
        if (is_real_par() .and. (icase==1 .or. icase==4)) then
          CALL GA_CREATE_STRIPED ('H',NAS,NIN,'TMAT',lg_T)
          CALL GA_Distribution (lg_T,myRank,ISTA,IEND,JSTA,JEND)
          IF (ISTA.GT.0 .AND. JSTA.GT.0) THEN
            CALL GA_Access (lg_T,ISTA,IEND,JSTA,JEND,mV1,LDM)
            NBLOCK=LDM*(JEND-JSTA+1)
            CALL DDAFILE(LUGRAD,2,DBL_MB(mV1),NBLOCK,IDSAVGRD)
            CALL GA_Release (lg_T,ISTA,IEND,JSTA,JEND)
          END IF
          CALL PSBMAT_WRITE ('T',iCase,iSym,lg_T,NAS*NIN)
          bStat = GA_Destroy (lg_T)
        else
#endif
          if (NAS*NIN > 0) then
            CALL DDAFILE(LUGRAD,2,WRK1,NAS*NIN,IDSAVGRD)
            ID = IDTMAT(ISYM,ICASE)
            CALL DDAFILE(LUSBT,1,WRK1,NAS*NIN,ID)
          end if
#ifdef _MOLCAS_MPP_
        end if
#endif
        !! Eigenvalue
        if (NIN > 0) then
          CALL DDAFILE(LUGRAD,2,WRK1,NIN,IDSAVGRD)
          ID  = IDBMAT(ISYM,ICASE)
          CALL DDAFILE(LUSBT,1,WRK1,NIN,ID)
        end if
        !! IS
!       CALL DDAFILE(LUGRAD,2,WRK1,NIS,IDSAVGRD)
!       CALL DDAFILE(LUSBT,1,WRK1,NIS,ID)
        if (do_lindep .and. NAS > 0) then
          CALL DDAFILE(LUGRAD,2,WRK1,NAS*(NAS+1)/2,IDSAVGRD)
          ID  = IDBoriMat(ISYM,ICASE)
          CALL DDAFILE(LUSTD,1,WRK1,NAS*(NAS+1)/2,ID)
        end if
      End If
    End Do
  End Do

  !! 4. E2TOT
  If (IORW == 1) Then
    WRK1(1) = E2TOT
    CALL DDAFILE(LUGRAD,IORW,WRK1,1,IDSAVGRD)
  Else If (IORW == 2) Then
    CALL DDAFILE(LUGRAD,IORW,WRK1,1,IDSAVGRD)
    E2TOT = WRK1(1)
  End If

  call mma_deallocate(WRK1)

  !! quasi-canonical active orbital energy
! CALL DDAFILE(LUGRAD,IORW,EPSA,NASHT,IDSAVGRD)

  !! 5. Save T-amplitude (IRHS = LURHS(1) = 51)
  Call SaveReadT1()

Contains

  Subroutine SaveReadT1()

    Implicit None

    integer(kind=iwp) :: lg_V1,NVEC,ICASE_,ISYM_

    !! IVECX = T (solution; not quasi-variational, before lambda-eq)
    IVECX = 2

    Do ICASE_ = 1, NCASES
      Do ISYM_ = 1, NSYM
        NIN = NINDEP(ISYM_,ICASE_)
        If (NIN == 0) Cycle
        NIS = NISUP(ISYM_,ICASE_)
        NAS = NASUP(ISYM_,ICASE_)
        NVEC = NIN*NIS
        If ((ICASE_ == 12) .OR. (ICASE_ == 13)) NVEC = NAS*NIS
        If (NVEC == 0) Cycle
        If ((ICASE_ == 12) .OR. (ICASE_ == 13)) Then
          Call RHS_ALLO(NAS,NIS,lg_V1)
          If (IORW == 1) Then
#ifdef _MOLCAS_MPP_
            if (is_real_par()) then
              CALL GA_Distribution (lg_V1,myRank,ISTA,IEND,JSTA,JEND)
              IF (IEND-ISTA+1.EQ.NAS .AND. ISTA.GT.0) THEN
                CALL GA_Access (lg_V1,ISTA,IEND,JSTA,JEND,mV1,LDW)
                NVEC=(IEND-ISTA+1)*(JEND-JSTA+1)
                IDISK=IOFFRHS(ISYM_,ICASE_)
                CALL DDAFILE(LURHS(IVECX),2,DBL_MB(mV1),NVEC,IDISK)
                CALL DDAFILE(LUGRAD,IORW,DBL_MB(mV1),NVEC,IDSAVGRD)
                CALL GA_Release (lg_V1,ISTA,IEND,JSTA,JEND)
              END IF
            else
#endif
              if (NAS*NIS > 0) then
                Call RHS_READ_SR(lg_V1,ICASE_,ISYM_,IVECX)
                CALL DDAFILE(LUGRAD,IORW,GA_Arrays(lg_V1)%A,NAS*NIS,IDSAVGRD)
              end if
#ifdef _MOLCAS_MPP_
            end if
#endif
          Else If (IORW == 2) Then
#ifdef _MOLCAS_MPP_
            if (is_real_par()) then
              CALL GA_Distribution (lg_V1,myRank,ISTA,IEND,JSTA,JEND)
              IF (IEND-ISTA+1.EQ.NAS .AND. ISTA.GT.0) THEN
                CALL GA_Access (lg_V1,ISTA,IEND,JSTA,JEND,mV1,LDW)
                NVEC=(IEND-ISTA+1)*(JEND-JSTA+1)
                CALL DDAFILE(LUGRAD,IORW,DBL_MB(mV1),NVEC,IDSAVGRD)
                IDISK=IOFFRHS(ISYM_,ICASE_)
                CALL DDAFILE(LURHS(IVECX),1,DBL_MB(mV1),NVEC,IDISK)
                CALL GA_Release (lg_V1,ISTA,IEND,JSTA,JEND)
              END IF
            else
#endif
              if (NAS*NIS > 0) then
                CALL DDAFILE(LUGRAD,IORW,GA_Arrays(lg_V1)%A,NAS*NIS,IDSAVGRD)
                CALL RHS_SAVE_SR(lg_V1,ICASE_,ISYM_,IVECX)
              end if
#ifdef _MOLCAS_MPP_
            end if
#endif
          End If
          CALL RHS_FREE(lg_V1)
        Else
          Call RHS_ALLO(NIN,NIS,lg_V1)
          If (IORW == 1) Then
#ifdef _MOLCAS_MPP_
            if (is_real_par()) then
              CALL GA_Distribution (lg_V1,myRank,ISTA,IEND,JSTA,JEND)
              IF (IEND-ISTA+1.EQ.NIN .AND. ISTA.GT.0) THEN
                CALL GA_Access (lg_V1,ISTA,IEND,JSTA,JEND,mV1,LDW)
                NVEC=(IEND-ISTA+1)*(JEND-JSTA+1)
                IDISK=IOFFRHS(ISYM_,ICASE_)
                CALL DDAFILE(LURHS(IVECX),2,DBL_MB(mV1),NVEC,IDISK)
                CALL DDAFILE(LUGRAD,IORW,DBL_MB(mV1),NVEC,IDSAVGRD)
                CALL GA_Release (lg_V1,ISTA,IEND,JSTA,JEND)
              END IF
            else
#endif
              if (NIN*NIS > 0) then
                Call RHS_READ_SR(lg_V1,ICASE_,ISYM_,IVECX)
                CALL DDAFILE(LUGRAD,IORW,GA_Arrays(lg_V1)%A,NIN*NIS,IDSAVGRD)
              end if
#ifdef _MOLCAS_MPP_
            end if
#endif
          Else If (IORW == 2) Then
#ifdef _MOLCAS_MPP_
            if (is_real_par()) then
              CALL GA_Distribution (lg_V1,myRank,ISTA,IEND,JSTA,JEND)
              IF (IEND-ISTA+1.EQ.NIN .AND. ISTA.GT.0) THEN
                CALL GA_Access (lg_V1,ISTA,IEND,JSTA,JEND,mV1,LDW)
                NVEC=(IEND-ISTA+1)*(JEND-JSTA+1)
                CALL DDAFILE(LUGRAD,IORW,DBL_MB(mV1),NVEC,IDSAVGRD)
                IDISK=IOFFRHS(ISYM_,ICASE_)
                CALL DDAFILE(LURHS(IVECX),1,DBL_MB(mV1),NVEC,IDISK)
                CALL GA_Release (lg_V1,ISTA,IEND,JSTA,JEND)
              END IF
            else
#endif
              if (NIN*NIS > 0) then
                CALL DDAFILE(LUGRAD,IORW,GA_Arrays(lg_V1)%A,NIN*NIS,IDSAVGRD)
                CALL RHS_SAVE_SR(lg_V1,ICASE_,ISYM_,IVECX)
              end if
#ifdef _MOLCAS_MPP_
            end if
#endif
          End If
          CALL RHS_FREE(lg_V1)
        End If
      End Do
    End Do

  End Subroutine SaveReadT1

End Subroutine SavGradParams

Subroutine SavGradParams2(Mode,UEFF,U0,H0)
!
! It seems that values that are unchanged during the gradient loop
! have to be separately saved and restored
! If this subroutine is updated, the shift at the beginning of
! the SavGradParams subroutine should also be updated
!
  use caspt2_global, only: LUGRAD
  use definitions, only: iwp,wp
  use stdalloc, only: mma_allocate, mma_deallocate
  use caspt2_module, only: Energy, ERFSelf, nBTri, nState, RFPert

  Implicit None

  integer(kind=iwp), intent(in) :: Mode
  real(kind=wp)    , intent(inout) :: UEFF(*),U0(*),H0(*)

  integer(kind=iwp) :: IORW,ID
  real(kind=wp), allocatable :: lTemp(:)

  !! Decide what to do
  If (Mode == 1) Then
    IORW = 1 !! Write
  Else If (Mode == 2) Then
    IORW = 2 !! Read
  End If

  ID = 0
  CALL DDAFILE(LUGRAD,IORW,ENERGY,NSTATE   ,ID)
  CALL DDAFILE(LUGRAD,IORW,UEFF  ,NSTATE**2,ID)
  CALL DDAFILE(LUGRAD,IORW,U0    ,NSTATE**2,ID)
  CALL DDAFILE(LUGRAD,IORW,H0    ,NSTATE**2,ID)

  if (RFpert) then
    Call mma_allocate(lTemp,NBTRI+1,Label='lTemp')
    if (Mode == 1) then
      Call Get_dScalar('RF Self Energy',lTemp(1+NBTRI))
      Call Get_dArray('Reaction field',lTemp,NBTRI)
      CALL DDAFILE(LUGRAD,IORW,lTemp,NBTRI+1,ID)
    else if (Mode == 2) then
      CALL DDAFILE(LUGRAD,IORW,lTemp,NBTRI+1,ID)
      Call Put_dScalar('RF Self Energy',lTemp(1+NBTRI))
      ERFSelf = lTemp(1+NBTRI)
      Call Put_dArray('Reaction field',lTemp,NBTRI)
    end if
    Call mma_deallocate(lTemp)
  end if

End Subroutine SavGradParams2
