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
! Copyright (C) 1994-1996, Jeppe Olsen                                 *
!***********************************************************************
! Density matrices between L and R
!
! I12 = 1 => only one-body density
! I12 = 2 => one- and two-body density matrices
!
! Jeppe Olsen,      Oct 1994
! GAS modifications Aug 1995
! Two body density added, 1996
!
! Two-body density is stored as rho2(ijkl)=<l!e(ij)e(kl)-delta(jk)e(il)!r>
! ijkl = ij*(ij-1)/2+kl, ij.ge.kl
!
      SUBROUTINE DENSI2_MCLR(I12,RHO1,RHO2,L,R,LUL,LUR,ieaw,n1,n2)
      use Str_Info, only: STR,MXNSTR,IATPM1,IATPM2,IBTPM1,IBTPM2,
     &                    ITYP_DUMMY,NELEC,NOCTYP
      use stdalloc, only: mma_allocate, mma_deallocate
      use MCLR_Data, only: IPRCIX,IPRDIA
      use Constants, only: Zero
      use MCLR_Data, only: IDC,PSSIGN
      use MCLR_Data, only: MXSB,MXSOOB,IASTFI,IBSTFI,ISMOST,MNR1IC,
     &                     MXR3IC
      use MCLR_Data, only: MAXI,MAXK,ICISTR
      use MCLR_Data, only: NACOB,IBTSOB,NOBPTS,NTSOB
      use DetDim, only: MXPOBS,MXINKA,MXPNGAS
      use CandS, only: ICSM,ISSM,ISSPC,ICSPC
      use input_mclr, only: nIrrep,nsMOB
      IMPLICIT None
      Integer I12
*.Output
      REAL*8 RHO1(*),RHO2(*)
*. Specific input
      REAL*8 L(*),R(*)

      INTEGER LUL,LUR,ieaw,n1,n2
*
* =====
*.Input
* =====
*
*.Definition of L and R is picked up from CANDS
* with L being S and  R being C
*. Before I forget it :
      INTEGER iSXSTSM(1),IDUMMY(1)
      Integer, Allocatable:: SIOIO(:), CIOIO(:), SBLTP(:), CBLTP(:)
      Integer, Allocatable:: IX(:,:), OOS(:,:)
      Real*8, Allocatable:: CB(:), SB(:), INSCR(:), C2(:), XIXS(:,:)
      Real*8, Allocatable:: RHO1S(:), RHO1P(:), XNATO(:)
      Integer idum(1)
      Integer IPRDEN,NGAS,IATP,IBTP,JATP,JBTP,NOCTPA,NOCTPB,NAEL,NBEL,
     &        IOCTPA,IOCTPB,MXSTBL0,MAXA,MAXA1,MAXB,MAXB1,MXSTBL,MXTSOB,
     &        IOBTP,IOBSM,LSCR1,INTSCR,IATP2,IBTP2,LSCR2,LSCR12,MAXIK,
     &        LSCR3,NOOS,IMNMX,MXCIJA,MXCIJAB,MXCIJB,MXCJ,MXIJST,
     &        MXIJSTF,MXSXBL

      IDUM = 0
      IPRDEN=1
      NGAS=3

      RHO1(1:NACOB**2) = Zero
      RHO2(1:NACOB**2*(NACOB**2+1)/2) = Zero
*
* Info for this internal space
*
      IATP = IASTFI(ISSPC)
      IBTP = IBSTFI(ISSPC)
      JATP = IASTFI(ICSPC)
      JBTP = IBSTFI(ICSPC)
      IF(IATP.NE.JATP.OR.IBTP.NE.JBTP) THEN
        WRITE(6,*) ' My world is falling apart'
        WRITE(6,*) ' C and sigma belongs to different types of strings'
        WRITE(6,*) ' IATP IBTP JATP JBTP ',IATP,IBTP,JATP,JBTP
        Call Abend( )
      END IF
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*. Offsets for supergroups
      IOCTPA = 1
      IOCTPB = 1
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)

*. Largest block of strings in zero order space
      MXSTBL0 = MXNSTR
*. Largest number of strings of given symmetry and type
      MAXA = 0
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(Str(IATPM1)%NSTSO,nIrrep*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(Str(IATPM2)%NSTSO,nIrrep*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(Str(IBTPM1)%NSTSO,nIrrep*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(Str(IBTPM2)%NSTSO,nIrrep*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MXSTBL = MAX(MAXA,MAXB)
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
*. replace with MXINKA !!!
      MAXI = MIN(MXINKA,MXSTBL)
      MAXK = MIN(MXINKA,MXSTBL)
*Largest active orbital block belonging to given type and symmetry
      MXTSOB = 0
      DO IOBTP = 1, NGAS
      DO IOBSM = 1, NSMOB
       MXTSOB = MAX(MXTSOB,NOBPTS(IOBTP,IOBSM))
      END DO
      END DO
*.Local scratch arrays for blocks of C and sigma
      LSCR1 = 0
      IF(ICISTR.LE.2) THEN
        LSCR1 = MXSB
      ELSE IF(ICISTR.EQ.3) THEN
        LSCR1 = MXSOOB
      END IF
      IF(ICISTR.EQ.1) THEN
        CALL mma_allocate(CB,LSCR1,Label='CB')
        CALL mma_allocate(SB,LSCR1,Label='SB')
      END IF

*.SCRATCH space for block of two-electron density matrix
* A 4 index block with four indices belonging OS class

      INTSCR = MXTSOB ** 4

      CALL mma_allocate(INSCR,INTSCR,Label='INSCR')

*
*. Arrays giving allowed type combinations
      CALL mma_allocate(SIOIO,NOCTPA*NOCTPB,Label='SIOIO')
      CALL mma_allocate(CIOIO,NOCTPA*NOCTPB,Label='CIOIO')

      CALL IAIBCM_MCLR(MNR1IC(ISSPC),MXR3IC(ISSPC),NOCTPA,NOCTPB,
     &            Str(IATP)%EL1,Str(IATP)%EL3,
     &            Str(IBTP)%EL1,Str(IBTP)%EL3,
     &            SIOIO,IPRDEN)

*
      CALL IAIBCM_MCLR(MNR1IC(ICSPC),MXR3IC(ICSPC),NOCTPA,NOCTPB,
     &            Str(IATP)%EL1,Str(IATP)%EL3,
     &            Str(IBTP)%EL1,Str(IBTP)%EL3,
     &            CIOIO,IPRDEN)

*
* Get memory requirements
*
*
      IATP2 = MIN(IATP+2,ITYP_Dummy)
      IBTP2 = MIN(IBTP+2,ITYP_Dummy)
      CALL MXRESC(CIOIO,IATP,IBTP,NOCTPA,NOCTPB,nIrrep,
     &            Str(IATP)%NSTSO,Str(IBTP)%NSTSO,
     &            IATP+1,Str(IATP+1)%NSTSO,NOCTYP(IATP+1),
     &            Str(IBTP+1)%NSTSO,NOCTYP(IBTP+1),
     &            NSMOB,3,3,NTSOB,IPRCIX,MAXK,
     &            Str(IATP2)%NSTSO,NOCTYP(IATP2),
     &            Str(IBTP2)%NSTSO,NOCTYP(IBTP2),
     &            Str(IATP)%EL123,Str(IBTP)%EL123,
     &            MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXIJST,
     &            MXIJSTF)


      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB)
      LSCR12 = MAX(LSCR1,2*LSCR2)
      CALL mma_allocate(C2,LSCR12,Label='C2')
*
*. Space for annihilation/creation mappings
      MAXIK = MAX(MAXI,MAXK)
      LSCR3 = MAX(MXSTBL*MXTSOB,MXIJST,MAXIK*MXTSOB*MXTSOB,MXSTBL0)
      CALL mma_allocate(IX,LSCR3,4,Label='IX')
      CALL mma_allocate(XIXS,LSCR3,4,Label='XIXS')
*. Arrays giving block type
      CALL mma_allocate(SBLTP,nIrrep,Label='SBLTP')
      CALL mma_allocate(CBLTP,nIrrep,Label='CBLTP')
*. Arrays for additional symmetry operation
      CALL ZBLTP(ISMOST(1,ISSM),nIrrep,IDC,SBLTP,idum)
      CALL ZBLTP(ISMOST(1,ICSM),nIrrep,IDC,CBLTP,idum)
*.10 OOS arrayy
      NOOS = NOCTPA*NOCTPB*nIrrep
      CALL mma_allocate(OOS,NOOS,10,Label='OSS')
* scratch space containing active one body
      CALL mma_allocate(RHO1S,NACOB ** 2,Label='RHO1S')
*. For natural orbitals
      CALL mma_allocate(RHO1P,NACOB*(NACOB+1)/2,Label='RHO1P')
      CALL mma_allocate(XNATO,NACOB **2,Label='XNATO')
*. Natural orbitals in symmetry blocks
*     CALL mma_allocate(RHO1SM,NACOB ** 2,Label='RHO1SM')
*     CALL mma_allocate(XNATSM,NACOB ** 2,Label='XNATSM')
*     CALL mma_allocate(OCCSM,NACOB,Label='OCCSM')
*
*
*. Transform from combination scaling to determinant scaling
*
      IF(IDC.NE.1.AND.ICISTR.EQ.1) THEN
*. Left CI vector
        CALL SCDTC2_MCLR(L,ISMOST(1,ISSM),SBLTP,nIrrep,
     &              NOCTPA,NOCTPB,Str(IATP)%NSTSO,
     &              Str(IBTP)%NSTSO,SIOIO,IDC,
     &              2,IDUMMY,IPRDIA)
*. Right CI vector
        CALL SCDTC2_MCLR(R,ISMOST(1,ICSM),CBLTP,nIrrep,
     &              NOCTPA,NOCTPB,Str(IATP)%NSTSO,
     &              Str(IBTP)%NSTSO,CIOIO,IDC,
     &              2,IDUMMY,IPRDIA)
      END IF

      IF(ICISTR.EQ.1) THEN
        CALL GASDN2_MCLR(I12,RHO1,RHO2,R,L,CB,SB,C2,
     &              CIOIO,SIOIO,ISMOST(1,ICSM),
     &       ISMOST(1,ISSM),CBLTP,SBLTP,
     &       NACOB,
     &       Str(IATP)%NSTSO,Str(IATP)%ISTSO,
     &       Str(IBTP)%NSTSO,Str(IBTP)%ISTSO,
     &       NAEL,IATP,NBEL,IBTP,
     &       IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &       nIrrep,NSMOB,nIrrep,nIrrep,
     &       MXPNGAS,NTSOB,IBTSOB,
     &       MAXK,MAXI,LSCR1,LSCR1,
     &       C2(LSCR2+1:LSCR12),C2(1:LSCR2),
     &       iSXSTSM,
     &       NGAS,
     &       Str(IATP)%EL123,Str(IBTP)%EL123,IDC,
     &       OOS(:,1), OOS(:,2), OOS(:,3), OOS(:,4), OOS(:,5),
     &       OOS(:,6), OOS(:,7), OOS(:,8), OOS(:,9), OOS(:,10),
     &       IX(:,1),XIXS(:,1),IX(:,2),XIXS(:,2),
     &       IX(:,3),XIXS(:,3),IX(:,4),XIXS(:,4),INSCR,
     &       MXPOBS,IPRDEN,RHO1S,LUL,LUR,
     &       PSSIGN,PSSIGN,RHO1P,XNATO,ieaw,n1,n2)
      ELSE IF(ICISTR.GE.2) THEN
        CALL GASDN2_MCLR(I12,RHO1,RHO2,R,L,R,L,C2,
     &              CIOIO,SIOIO,ISMOST(1,ICSM),
     &       ISMOST(1,ISSM),CBLTP,SBLTP,
     &       NACOB,
     &       Str(IATP)%NSTSO,Str(IATP)%ISTSO,
     &       Str(IBTP)%NSTSO,Str(IBTP)%ISTSO,
     &       NAEL,IATP,NBEL,IBTP,
     &       IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &       nIrrep,NSMOB,nIrrep,nIrrep,
     &       MXPNGAS,NTSOB,IBTSOB,
     &       MAXK,MAXI,LSCR1,LSCR1,
     &       C2(LSCR2+1:LSCR12),C2(1:LSCR2),
     &       iSXSTSM,
     &       NGAS,
     &       Str(IATP)%EL123,Str(IBTP)%EL123,IDC,
     &       OOS(:,1), OOS(:,2), OOS(:,3), OOS(:,4), OOS(:,5),
     &       OOS(:,6), OOS(:,7), OOS(:,8), OOS(:,9), OOS(:,10),
     &       IX(:,1),XIXS(:,1),IX(:,2),XIXS(:,2),
     &       IX(:,3),XIXS(:,3),IX(:,4),XIXS(:,4),INSCR,
     &       MXPOBS,IPRDEN,RHO1S,LUL,LUR,
     &       PSSIGN,PSSIGN,RHO1P,XNATO,ieaw,n1,n2)
      END IF

      IF(IDC.NE.1.AND.ICISTR.EQ.1) THEN
*. Transform from combination scaling to determinant scaling
*
        CALL SCDTC2_MCLR(L,ISMOST(1,ISSM),SBLTP,nIrrep,
     &              NOCTPA,NOCTPB,Str(IATP)%NSTSO,
     &              Str(IBTP)%NSTSO,SIOIO,IDC,
     &              1,IDUMMY,IPRDIA)
        CALL SCDTC2_MCLR(R,ISMOST(1,ICSM),CBLTP,nIrrep,
     &              NOCTPA,NOCTPB,Str(IATP)%NSTSO,
     &              Str(IBTP)%NSTSO,CIOIO,IDC,
     &              1,IDUMMY,IPRDIA)
      END IF
*
*     Free memory
*
      If (ICISTR.eq.1) Then
         Call mma_deallocate(CB)
         Call mma_deallocate(SB)
      End If
      Call mma_deallocate(INSCR)
      Call mma_deallocate(SIOIO)
      Call mma_deallocate(CIOIO)
      Call mma_deallocate(C2)
      Call mma_deallocate(IX)
      Call mma_deallocate(XIXS)
      Call mma_deallocate(SBLTP)
      Call mma_deallocate(CBLTP)
      Call mma_deallocate(OOS)
      Call mma_deallocate(RHO1S)
      Call mma_deallocate(RHO1P)
      Call mma_deallocate(XNATO)

      END SUBROUTINE DENSI2_MCLR
