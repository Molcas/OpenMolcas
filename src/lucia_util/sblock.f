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
      SUBROUTINE SBLOCK(  NBLOCK,  IBLOCK,   IBOFF,      CB,     HCB,
     &                       LUC,IRESTRICT, LUCBLK,ICBAT_RES,ICBAT_INI,
     &                  ICBAT_END)
*
* Generate a set of sigma blocks,
* The NBLOCK specified in IBLOCK starting from IBOFF,
* be more specific.
*
* The blocks are delivered in HCB
*
* The blocks are scaled and reformed to combination order
* If LUCBLK.GT.0, the blocks of C corresponding to IBLOCK
* are stored on LUCBLK
*
* ICONSPA,ICONSPB  added October 1996
* ICBAT_RES, ICBAT_INI, IBBAT_END added august 1997
*
* If ICBAT_RES .eq.1 then it as assumed that only
* Cbatches ICBAT_INI to ICBAT_END are stored on  LUC
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.fh"
*
* =====
*.Input
* =====
*
*.Definition of c and sigma spaces
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*. Sigma blocks require
      INTEGER IBLOCK(8,*)
*
*./ORBINP/ : NACOB used
#include "orbinp.fh"
#include "cicisp.fh"
#include "strbas.fh"
#include "cstate.fh"
#include "strinp.fh"
#include "stinf.fh"
#include "csm.fh"
#include "WrkSpc.fh"
#include "crun.fh"
#include "gasstr.fh"
#include "cgas.fh"
#include "lucinp.fh"
#include "cprnt.fh"
#include "glbbas.fh"
#include "oper.fh"
#include "io_util.fh"
*
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
#include "cintfo.fh"
      DIMENSION CB(*),HCB(*)
*
*     IDUM = 0
*     CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SBLOCK')
*
C?    WRITE(6,*) ' IPERTOP in SBLOCK = ', IPERTOP
      NTEST = 00
      IF(NTEST.GE.5)
     &WRITE(6,*) ' SBLOCK : ISSPC,ICSPC ', ISSPC,ICSPC
C?    WRITE(6,*) ' LUC in SBLOCK ', LUC
C?    WRITE(6,*) ' I12 in SBLOCK = ', I12
      IF(LUCBLK.GT.0) THEN
        IDISK(LUCBLK)=0
      END IF
*
* Info for this internal space
*. type of alpha and beta strings
      IATP = 1
      IBTP = 2
*. alpha and beta strings with an electron removed
      IATPM1 = 3
      IBTPM1 = 4
*. alpha and beta strings with two electrons removed
      IATPM2 = 5
      IBTPM2 = 6
*
*. Number of supergroups
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Offset for supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
*. connection matrices for supergroups
*
      CALL GETMEM('CONSPA','ALLO','INTE',KCONSPA,NOCTPA**2)
      CALL GETMEM('CONSPB','ALLO','INTE',KCONSPB,NOCTPB**2)
C     SPGRPCON(IOFSPGRP,NSPGRP,NGAS,MXPNGAS,IELFSPGRP,ISPGRPCON,IPRNT)
      CALL SPGRPCON(   IOCTPA,   NOCTPA,     NGAS,  MXPNGAS, NELFSPGP,
     &              iWORK(KCONSPA),IPRCIX)
      CALL SPGRPCON(   IOCTPB,   NOCTPB,     NGAS,  MXPNGAS, NELFSPGP,
     &              iWORK(KCONSPB),IPRCIX)
*
* string sym, string sym => sx sym
* string sym, string sym => dx sym
      CALL GETMEM('KSTSTS','ALLO','INTE',KSTSTS,NSMST ** 2)
      CALL GETMEM('KSTSTD','ALLO','INTE',KSTSTD,NSMST ** 2)
      CALL STSTSM(iWORK(KSTSTS),iWORK(KSTSTD),NSMST)
*. Largest block of strings in zero order space
      MXSTBL0 = MXNSTR
*. Largest number of strings of given symmetry and type
      MAXA = 0
      MAXA0 = IMNMX(iWORK(KNSTSO(IATP)),NSMST*NOCTYP(IATP),2)
      MAXA = MAX(MAXA,MAXA0)
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(iWORK(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(iWORK(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
*
      MAXB = 0
      MAXB0 = IMNMX(iWORK(KNSTSO(IBTP)),NSMST*NOCTYP(IBTP),2)
      MAXB = MAX(MAXB,MAXB0)
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(iWORK(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(iWORK(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MXSTBL = MAX(MAXA,MAXB)
      IF(IPRCIX.GE.3 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
      MAXI = MIN( MXINKA,MXSTBL)
      MAXK = MIN( MXINKA,MXSTBL)
*Largest active orbital block belonging to given type and symmetry
      MXTSOB = 0
      DO IOBTP = 1, NGAS
      DO IOBSM = 1, NSMOB
       MXTSOB = MAX(MXTSOB,NOBPTS(IOBTP,IOBSM))
      END DO
      END DO
C?    WRITE(6,*) ' MXTSOB = ', MXTSOB
*.Local scratch arrays for blocks of C and sigma
c      IF(ISIMSYM.EQ.0) THEN
        LSCR1 = MXSOOB
c      ELSE
c        LSCR1 = MXSOOB_AS
c      END IF
      LSCR1 = MAX(LSCR1,LCSBLK)
      IF(IPRCIX.GE.3)
     &WRITE(6,*) ' ICISTR,LSCR1 ',ICISTR,LSCR1
*.SCRATCH space for integrals
* A 4 index integral block with four indeces belonging OS class
      INTSCR = MAX(MXTSOB ** 4, NTOOB**2)
      IF(IPRCIX.GE.3)
     &WRITE(6,*) ' Integral scratch space ',INTSCR
      CALL GETMEM('INSCR ','ALLO','REAL',KINSCR,INTSCR)
      CALL GETMEM('INSCR2','ALLO','REAL',KINSCR2,INTSCR)
*. Arrays giving allowed type combinations '
      CALL GETMEM('CIOIO ','ALLO','INTE',KCIOIO,NOCTPA*NOCTPB)
      CALL GETMEM('SIOIO ','ALLO','INTE',KSIOIO,NOCTPA*NOCTPB)
*. Offsets for alpha and beta supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*. sigma needed for MXRESC
      CALL IAIBCM(ISSPC,iWORK(KSIOIO))
      CALL IAIBCM(ICSPC,iWORK(KCIOIO))
*. Arrays giving block type
      CALL GETMEM('CBLTP ','ALLO','INTE',KCBLTP,NSMST)
*. Arrays for additional symmetry operation
c      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
c        CALL MEMMAN(KSVST,NSMST,'ADDL  ',2,'SVST  ')
c        CALL SIGVST(WORK(KSVST),NSMST)
c      ELSE
         KSVST = 1
c      END IF
*
*.scratch space for projected matrices and a CI block
*
*. Scratch space for CJKAIB resolution matrices
*. Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
      IF( ISSPC.GE.ICSPC) THEN
         KSCIOIO = KSIOIO
      ELSE
         KSCIOIO = KCIOIO
      END IF
      CALL MXRESCPH(iWORK(KSCIOIO),IOCTPA, IOCTPB,  NOCTPA,  NOCTPB,
     &                  NSMST,NSTFSMSPGP,MXPNSMST,   NSMOB, MXPNGAS,
     &                   NGAS,   NOBPTS,   IPRCIX,   MAXK, NELFSPGP,
     &                   MXCJ,   MXCIJA,   MXCIJB,MXCIJAB,   MXSXBL,
     &               MXADKBLK,   IPHGAS, NHLFSPGP,   MNHL,  IADVICE,
     &              MXCJ_ALLSYM,MXADKBLK_AS,MX_NSPII)
      IF(IPRCIX.GE.3) THEN
        WRITE(6,*) 'SBLOCK : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXCJ_ALLSYM',
     &                       MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXCJ_ALLSYM
         WRITE(6,*) 'SBLOCK : MXADKBLK ', MXADKBLK
         WRITE(6,*) ' MX_NSPII = ', MX_NSPII
      END IF
*. For hardwired routines MXCIJAB is also used
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB,MX_NSPII)
      IF(IPRCIX.GE.3)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR2
*
      IF(IPRCIX.GE.3)  WRITE(6,*) ' LSCR2 = ', LSCR2
C  I assume memory was allocated for blocks, so
      KC2 = KVEC3
*
      KCJRES = KC2
      KSIRES = KC2 + LSCR2
*
      KSSCR = KSIRES
      KCSCR = KCJRES
*
*.vectors able to hold strings of given sym and type
      MAXIK = MAX(MAXI,MAXK)
*. I1 and Xi1s must also be able to hold largest st block
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB*MXTSOB,MXSTBL0)
      CALL GETMEM('I1    ','ALLO','INTE',KI1  ,LSCR3)
      CALL GETMEM('XI1S  ','ALLO','REAL',KXI1S,LSCR3)
      CALL GETMEM('I2    ','ALLO','INTE',KI2  ,LSCR3)
      CALL GETMEM('XI2S  ','ALLO','REAL',KXI2S,LSCR3)
      CALL GETMEM('I3    ','ALLO','INTE',KI3  ,LSCR3)
      CALL GETMEM('XI3S  ','ALLO','REAL',KXI3S,LSCR3)
      CALL GETMEM('I4    ','ALLO','INTE',KI4  ,LSCR3)
      CALL GETMEM('XI4S  ','ALLO','REAL',KXI4S,LSCR3)
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,iWORK(KCBLTP),iWORK(KSVST))
*.Some TTS arrays
      NTTS = MXNTTS
*
*. for partitioning of vector
      CALL GETMEM('LBTC  ','ALLO','INTE',KLLBT ,NTTS  )
      CALL GETMEM('LECTC ','ALLO','INTE',KLLEBT,NTTS  )
      CALL GETMEM('I1BTC ','ALLO','INTE',KLI1BT,NTTS  )
      CALL GETMEM('IBTC  ','ALLO','INTE',KLIBT ,8*NTTS)
*. For scaling for each TTS block
      CALL GETMEM('SCLFAC','ALLO','REAL',KLSCLFAC ,NTTS)

*. Space for four blocks of string occupations and arrays of
*. reordering arrays
*. Also used to hold an NORB*NORB matrix
      LZSCR = (MAX(NAEL,NBEL)+3)*(NOCOB+1) + 2 * NOCOB + NOCOB*NOCOB
      LZ    = (MAX(NAEL,NBEL)+2) * NOCOB
*. Set up to two blocks for orbital conserving operator
C     DO I1234 = 1, 2
      CALL GETMEM('KLOCS ','ALLO','INTE',KLOCSTR(1),MAX_STR_OC_BLK)
      DO I1234 = 1, 2
        CALL GETMEM('KLREO ','ALLO','INTE',KLREO(I1234),MAX_STR_SPGP)
        CALL GETMEM('KLZ   ','ALLO','INTE',KLZ(I1234),LZ)
      END DO
      CALL GETMEM('KLZSCR','ALLO','INTE',KLZSCR,LZSCR)
* 4 arrays containing all strings of given sym. Dimension can  be
*   reduced to largest number of strings in alpha or beta.
C?    WRITE(6,*) ' SBLOCKS : MAX_STR_SPGP = ', MAX_STR_SPGP
*
      IF(I12.EQ.2) THEN
        IDOH2 = 1
      ELSE
        IDOH2 = 0
      END IF
*. Place perturbation integrals over one body integrals
CINSERT_START
      IF(I12.EQ.2) THEN
        IDOH2 = 1
      ELSE
        IDOH2 = 0
      END IF
*
*. Prepare for perturbation calculation
*
C     IF(IPERTOP.NE.0) THEN
*. Matrix specifying partiotioned spaces
      CALL GETMEM('H0SPC ','ALLO','INTE',KLH0SPC,NOCTPA*NOCTPB)
      CALL H0INTSPC(   IH0SPC,   NPTSPC, IOCPTSPC,   NOCTPA,   NOCTPB,
     &                ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &                   NGAS,MXPNGAS,iWORK(KLH0SPC),NELFGP)
C       IF(IH0SPC.EQ.0) THEN
*. Form of perturbation in subspace has not been defined,
*. Use current IPART
          IH0INSPC(1) = IPART
C       END IF
C     END IF
*
* Jesper: Initializing ksvst
c      KCJPA = 1 ! jwk-cleanup
c      KSIPA = 1 ! jwk-cleanup
      CALL SBLOCKS(NBLOCK,IBLOCK(1,IBOFF),CB,HCB,WORK(KC2),
     &             iWORK(KCIOIO),ISMOST(1,ICSM),iWORK(KCBLTP),
     &             iWORK(KNSTSO(IATP)),iWORK(KNSTSO(IBTP)),
     &             NAEL,IATP,NBEL,IBTP,
     &             IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &             NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,MXPNGAS,
     &             ITSOB,MAXK,MAXI,LSCR1,
     &             WORK(KINSCR),WORK(KCSCR),WORK(KSSCR),
     &             iWORK(KSTSTS),iWORK(KSTSTD),SXDXSX,
     &             ADSXA,NGAS,NELFSPGP,IDC,
     &             iWORK(KI1),WORK(KXI1S),iWORK(KI2),WORK(KXI2S),
     &             IDOH2,MXPOBS,iWORK(KSVST),
     &             PSSIGN,IPRDIA,LUC,ICJKAIB,
     &             WORK(KCJRES),WORK(KSIRES),
     &             iWORK(KI3),WORK(KXI3S),iWORK(KI4),WORK(KXI4S),
     &             MXSXST,MXSXBL,MOCAA,
     &             iWORK(KLLBT),iWORK(KLLEBT),
     &             iWORK(KLI1BT),iWORK(KLIBT),
     &             IRESTRICT,
     &             iWORK(KCONSPA),iWORK(KCONSPB),WORK(KLSCLFAC),
     &             IPERTOP,IH0INSPC,iWORK(KLH0SPC),
     &             ICBAT_RES,ICBAT_INI,ICBAT_END,IUSE_PH,IPHGAS,
     &             I_RES_AB,ISIMSYM,WORK(KINSCR2))
*
* CALL SBLOCKS --> 91
*
*
      IF(IDC.EQ.2) THEN
*. reform
        CALL RFTTS(HCB,CB,IBLOCK(1,IBOFF),NBLOCK,1,NSMST,
     &             iWORK(KNSTSO(IATP)),iWORK(KNSTSO(IBTP)),
     &             IDC,PSSIGN,1,NTEST)
*. scale
        CALL SCDTTS(HCB,IBLOCK(1,IBOFF),NBLOCK, NSMST,
     &              iWORK(KNSTSO(IATP)),iWORK(KNSTSO(IBTP)),
     &              IDC,1,NTEST)
      END IF
*
      IF(LUCBLK.GT.0) THEN
        CALL ITODS([-1],1,-1,LUCBLK)
      END IF
*. Eliminate local memory
*     IDUM = 0
*     CALL MEMMAN(IDUM ,IDUM,'FLUSM ',2,'SBLOCK')
      CALL GETMEM('CONSPA','FREE','INTE',KCONSPA,NOCTPA**2)
      CALL GETMEM('CONSPB','FREE','INTE',KCONSPB,NOCTPB**2)
      CALL GETMEM('KSTSTS','FREE','INTE',KSTSTS,NSMST ** 2)
      CALL GETMEM('KSTSTD','FREE','INTE',KSTSTD,NSMST ** 2)
      CALL GETMEM('INSCR ','FREE','REAL',KINSCR,INTSCR)
      CALL GETMEM('INSCR2','FREE','REAL',KINSCR2,INTSCR)
      CALL GETMEM('CIOIO ','FREE','INTE',KCIOIO,NOCTPA*NOCTPB)
      CALL GETMEM('SIOIO ','FREE','INTE',KSIOIO,NOCTPA*NOCTPB)
      CALL GETMEM('CBLTP ','FREE','INTE',KCBLTP,NSMST)
      CALL GETMEM('I1    ','FREE','INTE',KI1  ,LSCR3)
      CALL GETMEM('XI1S  ','FREE','REAL',KXI1S,LSCR3)
      CALL GETMEM('I2    ','FREE','INTE',KI2  ,LSCR3)
      CALL GETMEM('XI2S  ','FREE','REAL',KXI2S,LSCR3)
      CALL GETMEM('I3    ','FREE','INTE',KI3  ,LSCR3)
      CALL GETMEM('XI3S  ','FREE','REAL',KXI3S,LSCR3)
      CALL GETMEM('I4    ','FREE','INTE',KI4  ,LSCR3)
      CALL GETMEM('XI4S  ','FREE','REAL',KXI4S,LSCR3)
      CALL GETMEM('SCLFAC','FREE','REAL',KLSCLFAC ,NTTS)
      CALL GETMEM('LBTC  ','FREE','INTE',KLLBT ,NTTS  )
      CALL GETMEM('LECTC ','FREE','INTE',KLLEBT,NTTS  )
      CALL GETMEM('I1BTC ','FREE','INTE',KLI1BT,NTTS  )
      CALL GETMEM('IBTC  ','FREE','INTE',KLIBT ,8*NTTS)
      CALL GETMEM('KLOCS ','FREE','INTE',KLOCSTR(1),MAX_STR_OC_BLK)
      DO I1234 = 1, 2
        CALL GETMEM('KLREO ','FREE','INTE',KLREO(I1234),MAX_STR_SPGP)
        CALL GETMEM('KLZ   ','FREE','INTE',KLZ(I1234),LZ)
      END DO
      CALL GETMEM('KLZSCR','FREE','INTE',KLZSCR,LZSCR)
      CALL GETMEM('H0SPC ','FREE','INTE',KLH0SPC,NOCTPA*NOCTPB)
      RETURN
      END
