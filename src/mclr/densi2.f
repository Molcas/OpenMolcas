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
* Copyright (C) 1994-1996, Jeppe Olsen                                 *
************************************************************************
       SUBROUTINE DENSI2(I12,RHO1,RHO2,L,R,LUL,LUR,ieaw,n1,n2)
       use Str_Info
*
* Density matrices between L and R
*
* I12 = 1 => only one-body density
* I12 = 2 => one- and two-body density matrices
*
* Jeppe Olsen,      Oct 1994
* GAS modifications Aug 1995
* Two body density added, 1996
*
* Two-body density is stored as rho2(ijkl)=<l!e(ij)e(kl)-delta(jk)e(il)!r>
* ijkl = ij*(ij-1)/2+kl, ij.ge.kl
*
      IMPLICIT REAL*8(A-H,O-Z)

*
* =====
*.Input
* =====
*
*.Definition of L and R is picked up from CANDS
* with L being S and  R being C
#include "cands.fh"
#include "detdim.fh"
#include "orbinp_mclr.fh"
#include "cicisp_mclr.fh"
#include "cstate_mclr.fh"
#include "csm.fh"
#include "stdalloc.fh"
#include "crun_mclr.fh"
#include "cprnt_mclr.fh"
#include "spinfo_mclr.fh"

#include "Input.fh"
#include "csmprd.fh"
*. Specific input
      REAL*8 L
      DIMENSION L(*),R(*)
*.Output
      DIMENSION RHO1(*),RHO2(*)
*. Before I forget it :
      DIMENSION iSXSTSM(1),IDUMMY(1)
      Integer, Allocatable:: SIOIO(:), CIOIO(:), SBLTP(:), CBLTP(:)
      Integer, Allocatable:: STSTS(:), STSTD(:), IX(:,:), OOS(:,:)
      Real*8, Allocatable:: CB(:), SB(:), INSCR(:), C2(:), XIXS(:,:)
      Real*8, Allocatable:: RHO1S(:), RHO1P(:), XNATO(:)
*     Real*8, Allocatable:: RHO1SM(:), XNATSM(:), OCCSM(:)
      Integer idum(1)

      IDUM = 0
CFUE  IPRDEN=0
      IPRDEN=1
      ZERO = 0.0D0
      NGAS=3

      CALL SETVEC(RHO1,ZERO ,NACOB ** 2 )
      CALL SETVEC(RHO2,ZERO ,NACOB ** 2 *(NACOB**2+1)/2)
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

*
* string sym, string sym => sx sym
* string sym, string sym => dx sym
      CALL mma_allocate(STSTS,NSMST ** 2,Label='STSTS')
      CALL mma_allocate(STSTD,NSMST ** 2,Label='STSTD')
      CALL STSTSM_MCLR(STSTS,STSTD,NSMST)
*. Largest block of strings in zero order space
      MXSTBL0 = MXNSTR
*. Largest number of strings of given symmetry and type
      MAXA = 0
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(Str(IATPM1)%NSTSO,NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(Str(IATPM2)%NSTSO,NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(Str(IBTPM1)%NSTSO,NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(Str(IBTPM2)%NSTSO,NSMST*NOCTYP(IBTPM2),2)
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
* A 4 index block with four indeces belonging OS class

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
      CALL MXRESC(CIOIO,IATP,IBTP,NOCTPA,NOCTPB,NSMST,
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
      CALL mma_allocate(SBLTP,NSMST,Label='SBLTP')
      CALL mma_allocate(CBLTP,NSMST,Label='CBLTP')
*. Arrays for additional symmetry operation
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,SBLTP,idum)
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,CBLTP,idum)
*.10 OOS arrayy
      NOOS = NOCTPA*NOCTPB*NSMST
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
        CALL SCDTC2_MCLR(L,ISMOST(1,ISSM),SBLTP,NSMST,
     &              NOCTPA,NOCTPB,Str(IATP)%NSTSO,
     &              Str(IBTP)%NSTSO,SIOIO,IDC,
     &              2,IDUMMY,IPRDIA)
*. Right CI vector
        CALL SCDTC2_MCLR(R,ISMOST(1,ICSM),CBLTP,NSMST,
     &              NOCTPA,NOCTPB,Str(IATP)%NSTSO,
     &              Str(IBTP)%NSTSO,CIOIO,IDC,
     &              2,IDUMMY,IPRDIA)
      END IF

      IF(ICISTR.EQ.1) THEN
        CALL GASDN2(I12,RHO1,RHO2,R,L,CB,SB,C2,
     &              CIOIO,SIOIO,ISMOST(1,ICSM),
     &       ISMOST(1,ISSM),CBLTP,SBLTP,
     &       NACOB,
     &       Str(IATP)%NSTSO,Str(IATP)%ISTSO,
     &       Str(IBTP)%NSTSO,Str(IBTP)%ISTSO,
     &       NAEL,IATP,NBEL,IBTP,
     &       IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &       NSMST,NSMOB,NSMSX,NSMDX,
     &       MXPNGAS,NTSOB,IBTSOB,
     &       MAXK,MAXI,LSCR1,LSCR1,
     &       C2(LSCR2+1:LSCR12),C2(1:LSCR2),
     &       iSXSTSM,STSTS,STSTD,SXDXSX,
     &       ADSXA,ASXAD,NGAS,
     &       Str(IATP)%EL123,Str(IBTP)%EL123,IDC,
     &       OOS(:,1), OOS(:,2), OOS(:,3), OOS(:,4), OOS(:,5),
     &       OOS(:,6), OOS(:,7), OOS(:,8), OOS(:,9), OOS(:,10),
     &       IX(:,1),XIXS(:,1),IX(:,2),XIXS(:,2),
     &       IX(:,3),XIXS(:,3),IX(:,3),XIXS(:,4),INSCR,
     &       MXPOBS,IPRDEN,RHO1S,LUL,LUR,
     &       PSSIGN,PSSIGN,RHO1P,XNATO,ieaw,n1,n2)
      ELSE IF(ICISTR.GE.2) THEN
        CALL GASDN2(I12,RHO1,RHO2,R,L,R,L,C2,
     &              CIOIO,SIOIO,ISMOST(1,ICSM),
     &       ISMOST(1,ISSM),CBLTP,SBLTP,
     &       NACOB,
     &       Str(IATP)%NSTSO,Str(IATP)%ISTSO,
     &       Str(IBTP)%NSTSO,Str(IBTP)%ISTSO,
     &       NAEL,IATP,NBEL,IBTP,
     &       IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &       NSMST,NSMOB,NSMSX,NSMDX,
     &       MXPNGAS,NTSOB,IBTSOB,
     &       MAXK,MAXI,LSCR1,LSCR1,
     &       C2(LSCR2+1:LSCR12),C2(1:LSCR2),
     &       iSXSTSM,STSTS,STSTD,SXDXSX,
     &       ADSXA,ASXAD,NGAS,
     &       Str(IATP)%EL123,Str(IBTP)%EL123,IDC,
     &       OOS(:,1), OOS(:,2), OOS(:,3), OOS(:,4), OOS(:,5),
     &       OOS(:,6), OOS(:,7), OOS(:,8), OOS(:,9), OOS(:,10),
     &       IX(:,1),XIXS(:,1),IX(:,2),XIXS(:,2),
     &       IX(:,3),XIXS(:,3),IX(:,3),XIXS(:,4),INSCR,
     &       MXPOBS,IPRDEN,RHO1S,LUL,LUR,
     &       PSSIGN,PSSIGN,RHO1P,XNATO,ieaw,n1,n2)
      END IF

      IF(IDC.NE.1.AND.ICISTR.EQ.1) THEN
*. Transform from combination scaling to determinant scaling
*
        CALL SCDTC2_MCLR(L,ISMOST(1,ISSM),SBLTP,NSMST,
     &              NOCTPA,NOCTPB,Str(IATP)%NSTSO,
     &              Str(IBTP)%NSTSO,SIOIO,IDC,
     &              1,IDUMMY,IPRDIA)
        CALL SCDTC2_MCLR(R,ISMOST(1,ICSM),CBLTP,NSMST,
     &              NOCTPA,NOCTPB,Str(IATP)%NSTSO,
     &              Str(IBTP)%NSTSO,CIOIO,IDC,
     &              1,IDUMMY,IPRDIA)
      END IF
*
*     Free memory
*
      Call mma_deallocate(STSTS)
      Call mma_deallocate(STSTD)
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
*     Call mma_deallocate(RHO1SM)
*     Call mma_deallocate(XNATSM)
*     Call mma_deallocate(OCCSM)

      RETURN
      END
