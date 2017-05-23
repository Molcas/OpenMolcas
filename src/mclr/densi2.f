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
*
* Density matrices between L and R
*
* I12 = 1 => only one-body density
* I12 = 2 => one- and two-body density matrices
*
* Jeppe Olsen,      Oct 94
* GAS modifications Aug 95
* Two body density added, '96
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
#include "strbas_mclr.fh"
#include "cstate_mclr.fh"
#include "strinp_mclr.fh"
#include "stinf_mclr.fh"
#include "csm.fh"
#include "WrkSpc.fh"
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
      CALL GETMEM('KSTSTS','ALLO','REAL',KSTSTS,NSMST ** 2)
      CALL GETMEM('KSTSTD','ALLO','REAL',KSTSTD,NSMST ** 2)
      CALL STSTSM_MCLR(WORK(KSTSTS),WORK(KSTSTD),NSMST)
*. Largest block of strings in zero order space
      MAXA0 = IMNMX(iWORK(KNSTSO(IATP)),NSMST*NOCTYP(IATP),2)
      MAXB0 = IMNMX(iWORK(KNSTSO(IBTP)),NSMST*NOCTYP(IBTP),2)
      MXSTBL0 = MXNSTR
*. Largest number of strings of given symmetry and type
      MAXA = 0
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(iWORK(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(iWORK(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(iWORK(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(iWORK(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
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
      MAXIJ = MXTSOB ** 2
*.Local scratch arrays for blocks of C and sigma
      IF(ICISTR.LE.2) THEN
        LSCR1 = MXSB
      ELSE IF(ICISTR.EQ.3) THEN
        LSCR1 = MXSOOB
      END IF
      IF(ICISTR.EQ.1) THEN
        CALL GETMEM('KCB','ALLO','REAL',KCB,LSCR1)
        CALL GETMEM('KSB','ALLO','REAL',KSB,LSCR1)
      END IF

*.SCRATCH space for block of two-electron density matrix
* A 4 index block with four indeces belonging OS class

      INTSCR = MXTSOB ** 4

      CALL GETMEM('INSCR','ALLO','REAL',KINSCR,INTSCR)

*
*. Arrays giving allowed type combinations '
      CALL GETMEM('SIOIO','ALLO','REAL',KSIOIO,NOCTPA*NOCTPB)
      CALL GETMEM('CIOIO','ALLO','REAL',KCIOIO,NOCTPA*NOCTPB)

      CALL IAIBCM_MCLR(MNR1IC(ISSPC),MXR3IC(ISSPC),NOCTPA,NOCTPB,
     &            iWORK(KEL1(IATP)),iWORK(KEL3(IATP)),
     &            iWORK(KEL1(IBTP)),iWORK(KEL3(IBTP)),
     &            WORK(KSIOIO),IPRDEN)

*
      CALL IAIBCM_MCLR(MNR1IC(ICSPC),MXR3IC(ICSPC),NOCTPA,NOCTPB,
     &            iWORK(KEL1(IATP)),iWORK(KEL3(IATP)),
     &            iWORK(KEL1(IBTP)),iWORK(KEL3(IBTP)),
     &            WORK(KCIOIO),IPRDEN)

*
* Get memory requirements
*
*
      CALL MXRESC(WORK(KcIOIO),IATP,IBTP,NOCTPA,NOCTPB,NSMST,
     &            iWORK(KNSTSO(IATP)),iWORK(KNSTSO(IBTP)),
     &            IATP+1,iWORK(KNSTSO(IATP+1)),NOCTYP(IATP+1),
     &            iWORK(KNSTSO(IBTP+1)),NOCTYP(IBTP+1),
     &            NSMOB,3,3,NTSOB,IPRCIX,MAXK,
     &            iWORK(KNSTSO(IATP+2)),NOCTYP(IATP+2),
     &            iWORK(KNSTSO(IBTP+2)),NOCTYP(IBTP+2),
     &            iWORK(KEL123(IATP)),iWORK(KEL123(IBTP)),
     &            MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXIJST,
     &            MXIJSTF)


      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB)
      LSCR12 = MAX(LSCR1,2*LSCR2)
      CALL GETMEM('KC2','ALLO','REAL',KC2,LSCR12)
      KSSCR = KC2
      KCSCR = KC2 + LSCR2
*
*. Space for annihilation/creation mappings
      MAXIK = MAX(MAXI,MAXK)
      LSCR3 = MAX(MXSTBL*MXTSOB,MXIJST,MAXIK*MXTSOB*MXTSOB,MXSTBL0)
      CALL GETMEM('I1','ALLO','REAL',KI1,LSCR3)
      CALL GETMEM('I2','ALLO','REAL',KI2,LSCR3)
      CALL GETMEM('I3','ALLO','REAL',KI3,LSCR3)
      CALL GETMEM('I4','ALLO','REAL',KI4,LSCR3)
      CALL GETMEM('XI1S','ALLO','REAL',KXI1S,LSCR3)
      CALL GETMEM('XI2S','ALLO','REAL',KXI2S,LSCR3)
      CALL GETMEM('XI3S','ALLO','REAL',KXI3S,LSCR3)
      CALL GETMEM('XI4S','ALLO','REAL',KXI4S,LSCR3)
*. Arrays giving block type
      CALL GETMEM('SBLTP','ALLO','REAL',KSBLTP,NSMST)
      CALL GETMEM('CBLTP','ALLO','REAL',KCBLTP,NSMST)
*. Arrays for additional symmetry operation
      KSVST = 1
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,WORK(KSBLTP),WORK(KSVST))
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,WORK(KCBLTP),WORK(KSVST))
*.10 OOS arrayy
      NOOS = NOCTPA*NOCTPB*NSMST
      CALL GETMEM('OOS1','ALLO','REAL',KOOS1,NOOS)
      CALL GETMEM('OOS2','ALLO','REAL',KOOS2,NOOS)
      CALL GETMEM('OOS3','ALLO','REAL',KOOS3,NOOS)
      CALL GETMEM('OOS4','ALLO','REAL',KOOS4,NOOS)
      CALL GETMEM('OOS5','ALLO','REAL',KOOS5,NOOS)
      CALL GETMEM('OOS6','ALLO','REAL',KOOS6,NOOS)
      CALL GETMEM('OOS7','ALLO','REAL',KOOS7,NOOS)
      CALL GETMEM('OOS8','ALLO','REAL',KOOS8,NOOS)
      CALL GETMEM('OOS9','ALLO','REAL',KOOS9,NOOS)
      CALL GETMEM('OOS10','ALLO','REAL',KOOS10,NOOS)
* scratch space containing active one body
      CALL GETMEM('RHO1S','ALLO','REAL',KRHO1S,NACOB ** 2)
*. For natural orbitals
      CALL GETMEM('RHO1P','ALLO','REAL',KRHO1P,NACOB*(NACOB+1)/2)
      CALL GETMEM('XNATO','ALLO','REAL',KXNATO,NACOB **2)
*. Natural orbitals in symmetry blocks
      CALL GETMEM('RHO1SM','ALLO','REAL',KRHO1SM,NACOB ** 2)
      CALL GETMEM('XNATSM','ALLO','REAL',KXNATSM,NACOB ** 2)
      CALL GETMEM('OCCSM','ALLO','REAL',KOCCSM,NACOB)
*
*
*. Transform from combination scaling to determinant scaling
*
      IF(IDC.NE.1.AND.ICISTR.EQ.1) THEN
*. Left CI vector
        CALL SCDTC2_MCLR(L,ISMOST(1,ISSM),WORK(KSBLTP),NSMST,
     &              NOCTPA,NOCTPB,iWORK(KNSTSO(IATP)),
     &              iWORK(KNSTSO(IBTP)),WORK(KSIOIO),IDC,
     &              2,IDUMMY,IPRDIA)
*. Right CI vector
        CALL SCDTC2_MCLR(R,ISMOST(1,ICSM),WORK(KCBLTP),NSMST,
     &              NOCTPA,NOCTPB,iWORK(KNSTSO(IATP)),
     &              iWORK(KNSTSO(IBTP)),WORK(KCIOIO),IDC,
     &              2,IDUMMY,IPRDIA)
      END IF

      IF(ICISTR.EQ.1) THEN
        CALL GASDN2(I12,RHO1,RHO2,
     &       R,L,WORK(KCB),WORK(KSB),WORK(KC2),
     &       WORK(KCIOIO),WORK(KSIOIO),ISMOST(1,ICSM),
     &       ISMOST(1,ISSM),WORK(KCBLTP),WORK(KSBLTP),
     &       NACOB,
     &       iWORK(KNSTSO(IATP)),iWORK(KISTSO(IATP)),
     &       iWORK(KNSTSO(IBTP)),iWORK(KISTSO(IBTP)),
     &       NAEL,IATP,NBEL,IBTP,
     &       IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &       NSMST,NSMOB,NSMSX,NSMDX,
     &       MXPNGAS,NTSOB,IBTSOB,
     &       MAXK,MAXI,LSCR1,LSCR1,
     &       WORK(KCSCR),WORK(KSSCR),
     &       SXSTSM,WORK(KSTSTS),WORK(KSTSTD),SXDXSX,
     &       ADSXA,ASXAD,NGAS,
     &       iWORK(KEL123(iatp)),iWORK(KEL123(ibtp)),IDC,
     &       WORK(KOOS1),WORK(KOOS2),WORK(KOOS3),WORK(KOOS4),
     &       WORK(KOOS5),WORK(KOOS6),WORK(KOOS7),WORK(KOOS8),
     &       WORK(KOOS9),WORK(KOOS10),
     &       WORK(KI1),WORK(KXI1S),WORK(KI2),WORK(KXI2S),
     &       WORK(KI3),WORK(KXI3S),WORK(KI4),WORK(KXI4S),WORK(KINSCR),
     &       MXPOBS,IPRDEN,WORK(KRHO1S),LUL,LUR,
     &       PSSIGN,PSSIGN,WORK(KRHO1P),WORK(KXNATO),ieaw,n1,n2)
      ELSE IF(ICISTR.GE.2) THEN
        CALL GASDN2(I12,RHO1,RHO2,
     &       R,L,R,L,WORK(KC2),
     &       WORK(KCIOIO),WORK(KSIOIO),ISMOST(1,ICSM),
     &       ISMOST(1,ISSM),WORK(KCBLTP),WORK(KSBLTP),
     &       NACOB,
     &       iWORK(KNSTSO(IATP)),iWORK(KISTSO(IATP)),
     &       iWORK(KNSTSO(IBTP)),iWORK(KISTSO(IBTP)),
     &       NAEL,IATP,NBEL,IBTP,
     &       IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &       NSMST,NSMOB,NSMSX,NSMDX,
     &       MXPNGAS,NTSOB,IBTSOB,
     &       MAXK,MAXI,LSCR1,LSCR1,
     &       WORK(KCSCR),WORK(KSSCR),
     &       SXSTSM,WORK(KSTSTS),WORK(KSTSTD),SXDXSX,
     &       ADSXA,ASXAD,NGAS,
     &       iWORK(KEL123(iatp)),iWORK(KEL123(ibtp)),IDC,
     &       WORK(KOOS1),WORK(KOOS2),WORK(KOOS3),WORK(KOOS4),
     &       WORK(KOOS5),WORK(KOOS6),WORK(KOOS7),WORK(KOOS8),
     &       WORK(KOOS9),WORK(KOOS10),
     &       WORK(KI1),WORK(KXI1S),WORK(KI2),WORK(KXI2S),
     &       WORK(KI3),WORK(KXI3S),WORK(KI4),WORK(KXI4S),WORK(KINSCR),
     &       MXPOBS,IPRDEN,WORK(KRHO1S),LUL,LUR,
     &       PSSIGN,PSSIGN,WORK(KRHO1P),WORK(KXNATO),ieaw,n1,n2)
      END IF

      IF(IDC.NE.1.AND.ICISTR.EQ.1) THEN
*. Transform from combination scaling to determinant scaling
*
        CALL SCDTC2_MCLR(L,ISMOST(1,ISSM),WORK(KSBLTP),NSMST,
     &              NOCTPA,NOCTPB,iWORK(KNSTSO(IATP)),
     &              iWORK(KNSTSO(IBTP)),WORK(KSIOIO),IDC,
     &              1,IDUMMY,IPRDIA)
        CALL SCDTC2_MCLR(R,ISMOST(1,ICSM),WORK(KCBLTP),NSMST,
     &              NOCTPA,NOCTPB,iWORK(KNSTSO(IATP)),
     &              iWORK(KNSTSO(IBTP)),WORK(KCIOIO),IDC,
     &              1,IDUMMY,IPRDIA)
      END IF
*
*     Free memory
*
      Call GetMem('KSTSTS','FREE','Real',KSTSTS,IDUM)
      Call GetMem('KSTSTD','FREE','Real',KSTSTD,IDUM)
      If (ICISTR.eq.1) Then
        Call GetMem('KCB','FREE','Real',KCB,IDUM)
        Call GetMem('KSB','FREE','Real',KSB,IDUM)
      End If
      Call GetMem('INTSCR','FREE','Real',KINSCR,IDUM)
      Call GetMem('SIOIO','FREE','Real',KSIOIO,IDUM)
      Call GetMem('CIOIO','FREE','Real',KCIOIO,IDUM)
      Call GetMem('KC2','FREE','Real',KC2,IDUM)
      Call GetMem('I1','FREE','Real',KI1,IDUM)
      Call GetMem('I2','FREE','Real',KI2,IDUM)
      Call GetMem('I3','FREE','Real',KI3,IDUM)
      Call GetMem('I4','FREE','Real',KI4,IDUM)
      Call GetMem('XI1S','FREE','Real',KXI1S,IDUM)
      Call GetMem('XI2S','FREE','Real',KXI2S,IDUM)
      Call GetMem('XI3S','FREE','Real',KXI3S,IDUM)
      Call GetMem('XI4S','FREE','Real',KXI4S,IDUM)
      Call GetMem('SBLTP','FREE','Real',KSBLTP,IDUM)
      Call GetMem('CBLTP','FREE','Real',KCBLTP,IDUM)
      Call GetMem('OOS1','FREE','Real',KOOS1,IDUM)
      Call GetMem('OOS2','FREE','Real',KOOS2,IDUM)
      Call GetMem('OOS3','FREE','Real',KOOS3,IDUM)
      Call GetMem('OOS4','FREE','Real',KOOS4,IDUM)
      Call GetMem('OOS5','FREE','Real',KOOS5,IDUM)
      Call GetMem('OOS6','FREE','Real',KOOS6,IDUM)
      Call GetMem('OOS7','FREE','Real',KOOS7,IDUM)
      Call GetMem('OOS8','FREE','Real',KOOS8,IDUM)
      Call GetMem('OOS9','FREE','Real',KOOS9,IDUM)
      Call GetMem('OOS10','FREE','Real',KOOS10,IDUM)
      Call GetMem('RHO1S','FREE','Real',KRHO1S,IDUM)
      Call GetMem('RHO1P','FREE','Real',KRHO1P,IDUM)
      Call GetMem('XNATO','FREE','Real',KXNATO,IDUM)
      Call GetMem('RHO1SM','FREE','Real',KRHO1SM,IDUM)
      Call GetMem('XNATSM','FREE','Real',KXNATSM,IDUM)
      Call GetMem('OCCSM','FREE','Real',KOCCSM,IDUM)

      RETURN
      END
