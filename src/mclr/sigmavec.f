************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE SigmaVec(C,HC,kic)
      Use Str_Info
*
* Outer routine for sigma vector generation
* RAS space
* This is a driver routine change in some rows just so it s
* should work in my mclr program. Of course I have to change
* the name on it but if you look in Jeppes mv7old(relaci) you
* will find a subroutine that looks very much like this
* This routine is just a setup routine for memory etc
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* =====
*.Input
* =====
*
#include "cands.fh"
#include "detdim.fh"
#include "orbinp_mclr.fh"
#include "cicisp_mclr.fh"
#include "cstate_mclr.fh"
#include "csm.fh"
#include "csfbas_mclr.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "crun_mclr.fh"

#include "Input.fh"
#include "cprnt_mclr.fh"
#include "glbbas_mclr.fh"
#include "genop.fh"
#include "csmprd.fh"
      Logical allokc,allokc2
      Dimension C(*),HC(*),kic(2)
      Integer sxstsm(1)
      Dimension idummy(1)
      Integer, Allocatable:: STSTS(:), STSTD(:), SVST(:),
     &                       SIOIO(:), CIOIO(:),
     &                       SBLTP(:), CBLTP(:),
     &                       I1(:), I2(:), I3(:), I4(:),
     &                       OOS(:,:)
      Real*8, Allocatable:: SB(:), CB(:), INSCR(:),
     &                      XI1S(:), XI2S(:), XI3S(:), XI4S(:),
     &                      C2(:), CJRES(:), SIRES(:)

      LUC=0
      LUHC=0
      iprnt=200
      ZERO = 0.0D0
      IDUM = 0
      NSDET = NINT(XISPSM(ISSM,ISSPC))
*
*. The story of MV7 : I started out from nothing, absolutely zero,
*
      call dcopy_(NSDET,[Zero],0,HC,1)
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
        Call Abend
      END IF
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
* string sym, string sym => sx sym
* string sym, string sym => dx sym
      Call mma_allocate(STSTS,NSMST**2,Label='STSTS')
      Call mma_allocate(STSTD,NSMST**2,Label='STSTD')

      CALL STSTSM_MCLR(STSTS,STSTD,NSMST)

*. Largest block of strings in zero order space
      MAXA0 = IMNMX(Str(IATP)%NSTSO,NSMST*NOCTYP(IATP),2)
      MAXB0 = IMNMX(Str(IBTP)%NSTSO,NSMST*NOCTYP(IBTP),2)
      MXSTBL0 = MAX(MAXA0,MAXB0)
*. Largest number of strings of given symmetry and type
      MAXA = 0
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(Str(IATP+1)%NSTSO,NSMST*NOCTYP(IATP+1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(Str(IATP+2)%NSTSO,NSMST*NOCTYP(IATP+2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(Str(IBTP+1)%NSTSO,NSMST*NOCTYP(IBTP+1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(Str(IBTP+2)%NSTSO,NSMST*NOCTYP(IBTP+2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MXSTBL = MAX(MAXA,MAXB)
      IF(IPRCIX.GE.2 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
      MAXI = MIN( MXINKA,MXSTBL)
      MAXK = MIN( MXINKA,MXSTBL)
*Largest active orbital block belonging to given type and symmetry
      MXTSOB = IMNMX(NTSOB,3*NSMOB,2)
      MAXIJ = MXTSOB ** 2
*.Local scratch arrays for blocks of C and sigma
      LSCR1 = 0
      IF(ICISTR.LE.2) THEN
        LSCR1 = MXSB
      ELSE IF(ICISTR.EQ.3) THEN
        LSCR1 = MXSOOB
      END IF
      IF(ICISTR.EQ.1) THEN
        Call mma_allocate(CB,LSCR1,Label='CB')
        Call mma_allocate(SB,LSCR1,Label='SB')
      END IF
*.SCRATCH space for integrals
* A 4 index integral block with four indeces belonging OS class
      INTSCR = MXTSOB ** 4

      Call mma_allocate(INSCR,INTSCR,Label='INSCR')
*. Arrays giving allowed type combinations

      Call mma_allocate(SIOIO,NOCTPA*NOCTPB,Label='SIOIO')
      Call mma_allocate(CIOIO,NOCTPA*NOCTPB,Label='CIOIO')
*
*      SIOIO,CIOIO [NOCTPA x NOCTPB]
*      Allowed combinations of alpha and beta strings
*      Sigma and CI vector respectively
*
*      (RAS1 & RAS3 constrains)
*
      CALL IAIBCM_MCLR(MNR1IC(ISSPC),MXR3IC(ISSPC),NOCTPA,NOCTPB,
     &            Str(IATP)%EL1,Str(IATP)%EL3,
     &            Str(IBTP)%EL1,Str(IBTP)%EL3,
     &            SIOIO,IPRCIX)
*
      CALL IAIBCM_MCLR(MNR1IC(ICSPC),MXR3IC(ICSPC),NOCTPA,NOCTPB,
     &            Str(IATP)%EL1,Str(IATP)%EL3,
     &            Str(IBTP)%EL1,Str(IBTP)%EL3,
     &            CIOIO,IPRCIX)
*
*
*. Arrays giving block type
      Call mma_allocate(SBLTP,NSMST,Label='SBLTP')
      Call mma_allocate(CBLTP,NSMST,Label='CBLTP')
*. Arrays for additional symmetry operation
*
      Call mma_allocate(SVST,NSMST,Label='SVST')
*
*.scratch space for projected matrices and a CI block
*
*. Scratch space for CJKAIB resolution matrices
*. Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
      IF(NOPART.EQ.0) THEN
         MAXPK = MAXK
       ELSE
         MAXPK = 0
      END IF
*
* Get memory requirements
*
*     MXCJ:MXCIJA:MXCIJB:MXCIJAB:MXSXBL:MXIJST:MXIJSTF
*
#ifdef _WARNING_WORKAROUND_
      IATP2=MIN(IATP+2,NSTTYP)
      IBTP2=MIN(IbTP+2,NSTTYP)
#else
      IATP2=IATP+2
      IBTP2=IBTP+2
#endif
      CALL MXRESC(CIOIO,IATP,IBTP,NOCTPA,NOCTPB,NSMST,
     &            Str(IATP)%NSTSO,Str(IBTP)%NSTSO,
     &            IATP+1,Str(IATP+1)%NSTSO,NOCTYP(IATP+1),
     &            Str(IBTP+1)%NSTSO,NOCTYP(IBTP+1),
     &            NSMOB,3,3,NTSOB,IPRCIX,MAXpK,
     &            Str(IATP2)%NSTSO,NOCTYP(IATP2),
     &            Str(IBTP2)%NSTSO,NOCTYP(IBTP2),
     &            Str(IATP)%EL123,Str(IBTP)%EL123,
     &            MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXIJST,
     &            MXIJSTF)
*
*
*.vectors able to hold strings of given sym and type
      MAXIK = MAX(MAXI,MAXK)
*. I1 and Xi1s must also be able to hold largest st block
* and (i,j,kstring) array

      LSCR3 = MAX(MXIJST,MAXIK*MXTSOB,MXSTBL0)

      Call mma_allocate(I1,LSCR3,Label='I1')
      Call mma_allocate(I2,LSCR3,Label='I2')
      Call mma_allocate(I3,MAXIK*MXTSOB,Label='I3')
      Call mma_allocate(I4,MAXIK*MXTSOB,Label='I4')
      Call mma_allocate(XI1S,LSCR3,Label='XI1S')
      Call mma_allocate(XI2S,LSCR3,Label='XI2S')
      Call mma_allocate(XI3S,MAXIK*MXTSOB,Label='XI3S')
      Call mma_allocate(XI4S,MAXIK*MXTSOB,Label='XI4S')
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB)
*
*. Merethe Interface
      MOCAA = 0
      IF(MOCAA.NE.0) THEN
*. These blocks will also be used for single excitations so,
*. To be removed   Largest block of SX excitations
        MAXE = MAX(NAEL,NBEL)
        MAXEl3 = MIN(MAXE,MXTSOB)
        MXSXST = (MXTSOB+1)*MAXEL3
        MXSXBL = MXSXST*MXSTBL0
        IF(IPRCIX.GE.2)
     &  WRITE(6,*) ' MXSXST,MXSXBL = ', MXSXST,MXSXBL
        LSCR2 = MAX(4*MXSXBL,LSCR2)
      END IF
*
      LSCR12 = MAX(LSCR1,2*LSCR2)
*
      AlloKC=.false.
      AlloKc2=.false.
      IF(IDIAG .EQ. 2 ) THEN
*. PICO diagonalizer uses block KVEC3, use this as scratch block
        KC2 = KVEC3
        IF ( 2 * LSCR2 .GT. LSCR1 ) THEN
           AlloKC=.true.
           Call mma_allocate(CJRES,LSCR2,Label='CJRES')
           KCJRES = ip_of_Work(CJRES(1))
           Call mma_allocate(SIRES,LSCR2,Label='SIRES')
           KSIRES = ip_of_Work(SIRES(1))
        ELSE
           KCJRES = KC2
           KSIRES = KC2 + LSCR2
        END IF
      ELSE
        AlloKC2=.true.
        Call mma_allocate(C2,LSCR12,Label='C2')
        KC2 = ip_of_Work(C2(1))
        KCJRES = KC2
        KSIRES = KC2 + LSCR2
      END IF
      KSSCR = KSIRES
      KCSCR = KCJRES
*
*     Symmetry handling symmetry allowed/forbidden
*
*     Out KSBLTP [NSMST]
*
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,SBLTP,SVST)
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,CBLTP,SVST)
*.10 OOS arrays
      nOOS = NOCTPA*NOCTPB*NSMST
      Call mma_allocate(OOS,nOOS,10,Label='OOS')
*
      iiCOPY=1
      IF(NOCSF.EQ.0) THEN
* Transform C vector from CSF to SD basis
        CALL CSDTVC_MCLR(C,HC,1,wORK(KDTOC),iWORK(KICTS(kic(1))),
     &                   icsm,iiCOPY,IPRDIA)
      END IF
*
      IF(IDC.NE.1.AND.ICISTR.EQ.1) THEN

*. Transform from combination scaling to determinant scaling
        CALL SCDTC2_MCLR(C,ISMOST(1,ICSM),CBLTP,NSMST,
     &                   NOCTPA,NOCTPB,Str(IATP)%NSTSO,
     &                   Str(IBTP)%NSTSO,CIOIO,IDC,
     &                   2,IDUMMY,IPRDIA)
      END IF

*
*     Goto 987
      IF(I12.EQ.2) THEN
        IDOH2 = 1
      ELSE
        IDOH2 = 0
      END IF
*
      IF(ICISTR.EQ.1) THEN
       LLUC = 0
       LLUHC = 0
      ELSE
       LLUC = LUC
       LLUHC = LUHC
      END IF
      call dcopy_(NSDET,[ZERO],0,HC,1)
*
      IF(ICISTR.EQ.1) THEN
      CALL RASSG4(C,HC,CB,SB,wORK(KC2),
     &            CIOIO,SIOIO,ISMOST(1,ICSM),
     &            ISMOST(1,ISSM),CBLTP,SBLTP,
     &            NORB1,NORB2,NORB3,NACOB,
     &            Str(IATP)%NSTSO,Str(IATP)%ISTSO,
     &            Str(IBTP)%NSTSO,Str(IBTP)%ISTSO,
     &            NAEL,IATP,NBEL,IBTP,NOCTPA,NOCTPB,
     &            NSMST,NSMOB,NSMSX,NSMDX,NTSOB,IBTSOB,ITSOB,
     &            MAXIJ,MAXK,MAXI,ICISTR,IINSTR,INTSCR,LSCR1,
     &            LSCR1,
     &            INSCR,wORK(KCSCR),wORK(KSSCR),
     &            SXSTSM,STSTS,STSTD,SXDXSX,
     &            ADSXA,ASXAD,
     &            Str(IATP)%EL1,Str(IATP)%EL3,
     &            Str(IBTP)%EL1,Str(IBTP)%EL3,IDC,
     &            OOS(:,1), OOS(:,2), OOS(:,3), OOS(:,4),
     &            OOS(:,5), OOS(:,6), OOS(:,7), OOS(:,8),
     &            OOS(:,9), OOS(:,10),
     &            I1,XI1S,I2,XI2S,
     &            I3,XI3S,I4,XI4S,
     &            IDOH2,
     &            SVST,PSSIGN,IPRDIA,LLUC,LLUHC,IST,
     &            WORK(KCJRES),WORK(KSIRES),NOPARt,TimeDep)

      Else
       Call SysHalt('sigmavec')
*
*      IF WE USE DISK REPLACE THE FIRST VARIBLES LIKE THIS
*      CALL RASSG4(C,HC,C,HC!,
*
      End If
*
      IF(IDC.NE.1.AND.ICISTR.EQ.1) THEN
*. Transform from combination scaling to determinant scaling

        CALL SCDTC2_MCLR(HC,ISMOST(1,ISSM),SBLTP,NSMST,
     &              NOCTPA,NOCTPB,Str(IATP)%NSTSO,
     &              Str(IBTP)%NSTSO,CIOIO,IDC,
     &              1,IDUMMY,IPRDIA)
      END IF

*
      IF(NOCSF.EQ.0) THEN
* Transform HC vector from SD to CSF basis
        CALL CSDTVC_MCLR(C,HC,2,WORK(KDTOC),iWORK(KICTS(kic(2))),
     &                   ISSM,1,IPRDIA)
      END IF

*. Eliminate local memory

      Call mma_deallocate(STSTS)
      Call mma_deallocate(STSTD)
      Call mma_deallocate(INSCR)
      Call mma_deallocate(SIOIO)
      Call mma_deallocate(CIOIO)
      Call mma_deallocate(SBLTP)
      Call mma_deallocate(CBLTP)

      Call mma_deallocate(I4)
      Call mma_deallocate(I3)
      Call mma_deallocate(I2)
      Call mma_deallocate(I1)
      Call mma_deallocate(XI4S)
      Call mma_deallocate(XI3S)
      Call mma_deallocate(XI2S)
      Call mma_deallocate(XI1S)
      Call mma_deallocate(OOS)
      IF (ICISTR.EQ.1) THEN
        Call mma_deallocate(CB)
        Call mma_deallocate(SB)
      End If
      Call mma_deallocate(SVST)
      IF (AlloKc ) THEN
        Call mma_deallocate(CJRES)
        Call mma_deallocate(SIRES)
      End If
      IF (AlloKc2) Call mma_deallocate(C2)
*
      RETURN
      END
