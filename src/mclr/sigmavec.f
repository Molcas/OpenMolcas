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
#include "strbas_mclr.fh"
#include "cstate_mclr.fh"
#include "strinp_mclr.fh"
#include "stinf_mclr.fh"
#include "csm.fh"
#include "csfbas_mclr.fh"
#include "WrkSpc.fh"
#include "crun_mclr.fh"

#include "Input.fh"
#include "cprnt_mclr.fh"
#include "glbbas_mclr.fh"
#include "genop.fh"
#include "csmprd.fh"
      Logical allokc,allokc2
      Dimension C(*),HC(*),kic(2)
      integer sxstsm(1)
      dimension idummy(1)
      LUC=0
      LUHC=0
      iprnt=200
      ZERO = 0.0D0
      IDUM = 0
      NSDET = NINT(XISPSM(ISSM,ISSPC))
*
*. The story of MV7 : I started out from nothing, absolutely zero,
*
      call dcopy_(NSDET,Zero,0,HC,1)

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
      Call Getmem('KSTSTS','ALLO','INTE',KSTSTS,NSMST ** 2)
      Call Getmem('KSTSTD','ALLO','INTE',KSTSTD,NSMST**2)
      CALL STSTSM_MCLR(iwORK(KSTSTS),iwORK(KSTSTD),NSMST)
*. Largest block of strings in zero order space
      MAXA0 = IMNMX(iWORK(KNSTSO(IATP)),NSMST*NOCTYP(IATP),2)
      MAXB0 = IMNMX(iWORK(KNSTSO(IBTP)),NSMST*NOCTYP(IBTP),2)
      MXSTBL0 = MAX(MAXA0,MAXB0)
*. Largest number of strings of given symmetry and type
      MAXA = 0
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(iWORK(KNSTSO(IATP+1)),NSMST*NOCTYP(IATP+1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(iWORK(KNSTSO(IATP+2)),NSMST*NOCTYP(IATP+2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(iWORK(KNSTSO(IBTP+1)),NSMST*NOCTYP(IBTP+1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(iWORK(KNSTSO(IBTP+2)),NSMST*NOCTYP(IBTP+2),2)
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
      IF(ICISTR.LE.2) THEN
        LSCR1 = MXSB
      ELSE IF(ICISTR.EQ.3) THEN
        LSCR1 = MXSOOB
      END IF
      IF(ICISTR.EQ.1) THEN
        Call Getmem('KCB   ','ALLO','REAL',KCB,LSCR1)
        Call Getmem('KSB   ','ALLO','REAL' ,KSB,LSCR1)
      END IF
*.SCRATCH space for integrals
* A 4 index integral block with four indeces belonging OS class
      INTSCR = MXTSOB ** 4

      Call Getmem('INSCR ','ALLO','REAL' ,KINSCR,INTSCR)
*. Arrays giving allowed type combinations '

      Call Getmem('SIOIO ','ALLO','INTE' ,KSIOIO,NOCTPA*NOCTPB)
      Call Getmem('CIOIO ','ALLO','INTE' ,KCIOIO,NOCTPA*NOCTPB)
*
*      SIOIO,CIOIO [NOCTPA x NOCTPB]
*      Allowed combinations of alpha and beta strings
*      Sigma and CI vector respectively
*
*      (RAS1 & RAS3 constrains)
*
      CALL IAIBCM_MCLR(MNR1IC(ISSPC),MXR3IC(ISSPC),NOCTPA,NOCTPB,
     &            iWORK(KEL1(IATP)),iWORK(KEL3(IATP)),
     &            iWORK(KEL1(IBTP)),iWORK(KEL3(IBTP)),
     &            iwORK(KSIOIO),IPRCIX)
*
      CALL IAIBCM_MCLR(MNR1IC(ICSPC),MXR3IC(ICSPC),NOCTPA,NOCTPB,
     &            iWORK(KEL1(IATP)),iWORK(KEL3(IATP)),
     &            iWORK(KEL1(IBTP)),iWORK(KEL3(IBTP)),
     &            iwORK(KCIOIO),IPRCIX)
*
*
*. Arrays giving block type
      Call Getmem('SBLTP ','ALLO','INTE' ,KSBLTP,NSMST)
      Call Getmem('CBLTP ','ALLO','INTEGER' ,KCBLTP,NSMST)
*. Arrays for additional symmetry operation
*
      KSVST = 1

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
      CALL MXRESC(iwORK(KCIOIO),IATP,IBTP,NOCTPA,NOCTPB,NSMST,
     &            iWORK(KNSTSO(IATP)),iWORK(KNSTSO(IBTP)),
     &            IATP+1,iWORK(KNSTSO(IATP+1)),NOCTYP(IATP+1),
     &            iWORK(KNSTSO(IBTP+1)),NOCTYP(IBTP+1),
     &            NSMOB,3,3,NTSOB,IPRCIX,MAXpK,
     &            iWORK(KNSTSO(IATP+2)),NOCTYP(IATP+2),
     &            iWORK(KNSTSO(IBTP+2)),NOCTYP(IBTP+2),
     &            iWORK(KEL123(IATP)),iWORK(KEL123(IBTP)),
     &            MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXIJST,
     &            MXIJSTF)
*
*
*.vectors able to hold strings of given sym and type
      MAXIK = MAX(MAXI,MAXK)
*. I1 and Xi1s must also be able to hold largest st block
* and (i,j,kstring) array

      LSCR3 = MAX(MXIJST,MAXIK*MXTSOB,MXSTBL0)

      Call Getmem('I1    ','ALLO','INTE',KI1,LSCR3)
      Call Getmem('I2    ','ALLO','INTE',KI2,LSCR3)
      Call Getmem('I3    ','ALLO','INTE',KI3,MAXIK*MXTSOB)
      Call Getmem('I4    ','ALLO','INTE',KI4,MAXIK*MXTSOB)
      Call Getmem('XI1S  ','ALLO','REAL',KXI1S,LSCR3)
      Call Getmem('XI2S    ','ALLO','REAL',KXI2S,LSCR3)
      Call Getmem('XI3S    ','ALLO','REAL',KXI3S,MAXIK*MXTSOB)
      Call Getmem('XI4S  ','ALLO','REAL',KXI4S,MAXIK*MXTSOB)
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB)
*
*. Merethe Interface
      MOCAA = 0
      IF(MOCAA.NE.0) THEN
*. These blocks will also be used for single excitations so,
*. TO be removed   Largest block of SX excitations
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
        IF(2*LSCR2 .GT. LSCR1 ) THEN
           AlloKC=.true.
           Call Getmem('KCRJES','ALLO','REAL',KCJRES,LSCR2)
           Call Getmem('KSIRES','ALLO','REAL',KSIRES,LSCR2)
        ELSE
           KCJRES = KC2
           KSIRES = KC2 + LSCR2
        END IF
      ELSE
        AlloKc2=.true.
        Call Getmem('KC2','ALLO','REAL',KC2,LSCR12)
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
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,iwORK(KSBLTP),
     &           iwORK(KSVST))
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,iwORK(KCBLTP),
     &           iwORK(KSVST))
*.10 OOS arrayy
      NOOS = NOCTPA*NOCTPB*NSMST
      Call Getmem('OOS1','ALLO','INTE',KOOS1,NOOS)
      Call Getmem('OOS2','ALLO','INTE',KOOS2,NOOS)
      Call Getmem('OOS3','ALLO','INTE',KOOS3,NOOS)
      Call Getmem('OOS4','ALLO','INTE',KOOS4,NOOS)
      Call Getmem('OOS5','ALLO','INTE',KOOS5,NOOS)
      Call Getmem('OOS6','ALLO','INTE',KOOS6,NOOS)
      Call Getmem('OOS7','ALLO','INTE',KOOS7,NOOS)
      Call Getmem('OOS8','ALLO','INTE',KOOS8,NOOS)
      Call Getmem('OOS9','ALLO','INTE',KOOS9,NOOS)
      Call Getmem('OOS10','ALLO','INTE',KOOS10,NOOS)
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
        CALL SCDTC2_MCLR(C,ISMOST(1,ICSM),iwORK(KCBLTP),NSMST,
     &              NOCTPA,NOCTPB,iWORK(KNSTSO(IATP)),
     &              iWORK(KNSTSO(IBTP)),iwORK(KCIOIO),IDC,
     &              2,IDUMMY,IPRDIA)
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
      call dcopy_(NSDET,ZERO,0,HC,1)
*
      IF(ICISTR.EQ.1) THEN
      CALL RASSG4(C,HC,wORK(KCB),wORK(KSB),wORK(KC2),
     &            iwORK(KCIOIO),iwORK(KSIOIO),ISMOST(1,ICSM),
     &            ISMOST(1,ISSM),iwORK(KCBLTP),iwORK(KSBLTP),
     &            NORB1,NORB2,NORB3,NACOB,
     &            iWORK(KNSTSO(IATP)),iWORK(KISTSO(IATP)),
     &            iWORK(KNSTSO(IBTP)),iWORK(KISTSO(IBTP)),
     &            NAEL,IATP,NBEL,IBTP,NOCTPA,NOCTPB,
     &            NSMST,NSMOB,NSMSX,NSMDX,NTSOB,IBTSOB,ITSOB,
     &            MAXIJ,MAXK,MAXI,ICISTR,IINSTR,INSCR,LSCR1,
     &            LSCR1,
     &            wORK(KINSCR),wORK(KCSCR),wORK(KSSCR),
     &            SXSTSM,iwORK(KSTSTS),iwORK(KSTSTD),SXDXSX,
     &            ADSXA,ASXAD,
     &            iWORK(KEL1(IATP)),iWORK(KEL3(IATP)),
     &            iWORK(KEL1(IBTP)),iWORK(KEL3(IBTP)),IDC,
     &            iwORK(KOOS1),iwORK(KOOS2),iwORK(KOOS3),
     &            iwORK(KOOS4),iwORK(KOOS5),iwORK(KOOS6),
     &            iwORK(KOOS7),iwORK(KOOS8),
     &            iwORK(KOOS9),iwORK(KOOS10),
     &            iWORK(KI1),wORK(KXI1S),iWORK(KI2),wORK(KXI2S),
     &            iWORK(KI3),wORK(KXI3S),iWORK(KI4),wORK(KXI4S),
     &            IDOH2,
     &            iwORK(KSVST),PSSIGN,IPRDIA,LLUC,LLUHC,IST,
     &            wORK(KCJRES),wORK(KSIRES),NOPARt,TimeDep)

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

        CALL SCDTC2_MCLR(HC,ISMOST(1,ISSM),iwORK(KSBLTP),NSMST,
     &              NOCTPA,NOCTPB,iWORK(KNSTSO(IATP)),
     &              iWORK(KNSTSO(IBTP)),iwORK(KCIOIO),IDC,
     &              1,IDUMMY,IPRDIA)
      END IF

*
      IF(NOCSF.EQ.0) THEN
* Transform HC vector from SD to CSF basis
        CALL CSDTVC_MCLR(C,HC,2,WORK(KDTOC),iWORK(KICTS(kic(2))),
     &                   ISSM,1,IPRDIA)
      END IF

*. Eliminate local memory

      Call Getmem('KSTSTS','FREE','INTE',KSTSTS,NSMST ** 2)
      Call Getmem('KSTSTD','FREE','INTE',KSTSTD,NSMST)
      Call Getmem('INSCR ','FREE','REAL' ,KINSCR,INTSCR)
      Call Getmem('SIOIO ','FREE','INTE' ,KSIOIO,NOCTPA*NOCTPB)
      Call Getmem('CIOIO ','FREE','INTE' ,KCIOIO,NOCTPA*NOCTPB)
      Call Getmem('SBLTP ','FREE','REAL' ,KSBLTP,NSMST) !origonal
      Call Getmem('CBLTP ','FREE','REAL' ,KCBLTP,NSMS)

! =====================================================================
!      YMA testing
!      Call Getmem('KSTSTS','FREE','INTE',KSTSTS,1000 ** 2) YMA testing
!      Call Getmem('KSTSTD','FREE','INTE',KSTSTD,10000)
!      Call Getmem('INSCR ','FREE','REAL' ,KINSCR,10000)
!      Call Getmem('SIOIO ','FREE','INTE' ,KSIOIO,10000)
!      Call Getmem('CIOIO ','FREE','INTE' ,KCIOIO,10000)
!      Call Getmem('SBLTP ','FREE','REAL' ,KSBLTP,10000) !origonal
!      Call Getmem('CBLTP ','FREE','REAL' ,KCBLTP,10000)
! =====================================================================
      Call Getmem('I1    ','FREE','INTE',KI1,LSCR3)
      Call Getmem('I2    ','FREE','INTE',KI2,LSCR3)
      Call Getmem('I3    ','FREE','INTE',KI3,MAXIK*MXTSOB)
      Call Getmem('I4    ','FREE','INTE',KI4,MAXIK*MXTSOB)
      Call Getmem('XI1S  ','FREE','REAL',KXI1S,LSCR3)
      Call Getmem('XI2S  ','FREE','REAL',KXI2S,LSCR3)
      Call Getmem('XI3S  ','FREE','REAL',KXI3S,MAXIK*MXTSOB)
      Call Getmem('XI4S  ','FREE','REAL',KXI4S,MAXIK*MXTSOB)
      Call Getmem('OOS1  ','FREE','REAL',KOOS1,NOOS)
      Call Getmem('OOS2  ','FREE','REAL',KOOS2,NOOS)
      Call Getmem('OOS3  ','FREE','REAL',KOOS3,NOOS)
      Call Getmem('OOS4  ','FREE','REAL',KOOS4,NOOS)
      Call Getmem('OOS5  ','FREE','REAL',KOOS5,NOOS)
      Call Getmem('OOS6  ','FREE','REAL',KOOS6,NOOS)
      Call Getmem('OOS7  ','FREE','REAL',KOOS7,NOOS)
      Call Getmem('OOS8  ','FREE','REAL',KOOS8,NOOS)
      Call Getmem('OOS9  ','FREE','REAL',KOOS9,NOOS)
      Call Getmem('OOS10 ','FREE','REAL',KOOS10,NOOS)
      IF(ICISTR.EQ.1) THEN
        Call Getmem('KCB   ','FREE','REAL',KCB,LSCR1)
        Call Getmem('KSB   ','FREE','REAL' ,KSB,LSCR1)
      End If
      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
        Call Getmem('SVST  ','FREE','INTEGER' ,KSVST,NSMST)
      End If
      IF(AlloKc ) THEN
        Call Getmem('KCRJES','FREE','REAL',KCJRES,LSCR2)
        Call Getmem('KSIRES','FREE','REAL',KSIRES,LSCR2)
      End If
      IF(AlloKc2) THEN
        Call Getmem('KC2','FREE','REAL',KC2,LSCR12)
      End If
*
      RETURN
      END
