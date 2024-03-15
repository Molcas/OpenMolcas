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
      use stdalloc, only: mma_allocate, mma_deallocate
      use GLBBAS, only: VEC3
      use hidscr, only: ZSCR, ZOCSTR => OCSTR, REO, Z
      use Local_Arrays, only: CLBT, CLEBT, CI1BT, CIBT, CBLTP,
     &                        Allocate_Local_Arrays,
     &                      Deallocate_Local_Arrays
      use strbas
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
* CONSPA,CONSPB  added October 1996
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
#include "cands.fh"
*. Sigma blocks require
      INTEGER IBLOCK(8,*)
*
*./ORBINP/ : NACOB used
#include "orbinp.fh"
#include "cicisp.fh"
#include "cstate.fh"
#include "strinp.fh"
#include "stinf.fh"
#include "csm.fh"
#include "crun.fh"
#include "gasstr.fh"
#include "cgas.fh"
#include "lucinp.fh"
#include "cprnt.fh"
#include "oper.fh"
#include "io_util.fh"
*
#include "csmprd.fh"
#include "cintfo.fh"
      DIMENSION CB(*),HCB(*)
      Integer, Allocatable:: CONSPA(:), CONSPB(:)
      Real*8, Allocatable:: INSCR(:), INSCR2(:)
      Integer, Allocatable:: STSTS(:), STSTD(:)
      Integer, Allocatable, Target:: CIOIO(:), SIOIO(:)
      Integer, Pointer:: SCIOIO(:)
      Integer, Allocatable:: I1(:), I2(:), I3(:), I4(:)
      Real*8, Allocatable:: XI1S(:), XI2S(:), XI3S(:), XI4S(:)
      Real*8, Allocatable:: LSCLFAC(:)
      Integer, Allocatable:: SVST(:)
      Integer, Allocatable:: H0SPC(:)
*
*     IDUM = 0
*     CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SBLOCK')
*
      NTEST = 00
#ifdef _DEBUGPRINT_
      IF(NTEST.GE.5) WRITE(6,*) ' SBLOCK : ISSPC,ICSPC ', ISSPC,ICSPC
#endif
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
      Call mma_allocate(CONSPA,NOCTPA**2,Label='CONSPA')
      Call mma_allocate(CONSPB,NOCTPB**2,Label='CONSPB')
C     SPGRPCON(IOFSPGRP,NSPGRP,NGAS,MXPNGAS,IELFSPGRP,ISPGRPCON,IPRNT)
      CALL SPGRPCON(   IOCTPA,   NOCTPA,     NGAS,  MXPNGAS, NELFSPGP,
     &              CONSPA,IPRCIX)
      CALL SPGRPCON(   IOCTPB,   NOCTPB,     NGAS,  MXPNGAS, NELFSPGP,
     &              CONSPB,IPRCIX)
*
* string sym, string sym => sx sym
* string sym, string sym => dx sym
      Call mma_allocate(STSTS,NSMST**2,Label='STSTS')
      Call mma_allocate(STSTD,NSMST**2,Label='STSTD')
      CALL STSTSM(STSTS,STSTD,NSMST)
*. Largest block of strings in zero order space
      MXSTBL0 = MXNSTR
*. Largest number of strings of given symmetry and type
      MAXA = 0
      MAXA0 = IMNMX(NSTSO(IATP)%I,NSMST*NOCTYP(IATP),2)
      MAXA = MAX(MAXA,MAXA0)
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(NSTSO(IATPM1)%I,NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(NSTSO(IATPM2)%I,NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
*
      MAXB = 0
      MAXB0 = IMNMX(NSTSO(IBTP)%I,NSMST*NOCTYP(IBTP),2)
      MAXB = MAX(MAXB,MAXB0)
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(NSTSO(IBTPM1)%I,NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(NSTSO(IBTPM2)%I,NSMST*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MXSTBL = MAX(MAXA,MAXB)
#ifdef _DEBUGPRINT_
      IF(IPRCIX.GE.3 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
#endif
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
#ifdef _DEBUGPRINT_
      IF(IPRCIX.GE.3)
     &WRITE(6,*) ' ICISTR,LSCR1 ',ICISTR,LSCR1
#endif
*.SCRATCH space for integrals
* A 4 index integral block with four indices belonging OS class
      INTSCR = MAX(MXTSOB ** 4, NTOOB**2)
#ifdef _DEBUGPRINT_
      IF(IPRCIX.GE.3)
     &WRITE(6,*) ' Integral scratch space ',INTSCR
#endif
      Call mma_allocate(INSCR,INTSCR,Label='INSCR')
      Call mma_allocate(INSCR2,INTSCR,Label='INSCR2')
*. Arrays giving allowed type combinations '
      Call mma_allocate(CIOIO,NOCTPA*NOCTPB,Label='CIOIO')
      Call mma_allocate(SIOIO,NOCTPA*NOCTPB,Label='SIOIO')
*. Offsets for alpha and beta supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*. sigma needed for MXRESC
      CALL IAIBCM(ISSPC,SIOIO)
      CALL IAIBCM(ICSPC,CIOIO)
*. Arrays for additional symmetry operation
c      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
c        Call mma_allocate(SVST,NSMST,Label='SVST')
c        CALL SIGVST(SVST,NSMST)
c      ELSE
         Call mma_allocate(SVST,1,Label='SVST')
c      END IF
*
*.scratch space for projected matrices and a CI block
*
*. Scratch space for CJKAIB resolution matrices
*. Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
      IF( ISSPC.GE.ICSPC) THEN
         SCIOIO => SIOIO
      ELSE
         SCIOIO => CIOIO
      END IF
      CALL MXRESCPH(SCIOIO,IOCTPA, IOCTPB,  NOCTPA,  NOCTPB,
     &                  NSMST,NSTFSMSPGP,MXPNSMST,   NSMOB, MXPNGAS,
     &                   NGAS,   NOBPTS,   IPRCIX,   MAXK, NELFSPGP,
     &                   MXCJ,   MXCIJA,   MXCIJB,MXCIJAB,   MXSXBL,
     &               MXADKBLK,   IPHGAS, NHLFSPGP,   MNHL,  IADVICE,
     &              MXCJ_ALLSYM,MXADKBLK_AS,MX_NSPII)
#ifdef _DEBUGPRINT_
      IF(IPRCIX.GE.3) THEN
        WRITE(6,*) 'SBLOCK : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXCJ_ALLSYM',
     &                       MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXCJ_ALLSYM
         WRITE(6,*) 'SBLOCK : MXADKBLK ', MXADKBLK
         WRITE(6,*) ' MX_NSPII = ', MX_NSPII
      END IF
#endif
*. For hardwired routines MXCIJAB is also used
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB,MX_NSPII)
#ifdef _DEBUGPRINT_
      IF(IPRCIX.GE.3)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR2
*
      IF(IPRCIX.GE.3)  WRITE(6,*) ' LSCR2 = ', LSCR2
#endif
C  I assume memory was allocated for blocks, so
*
*.vectors able to hold strings of given sym and type
      MAXIK = MAX(MAXI,MAXK)
*. I1 and Xi1s must also be able to hold largest st block
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB*MXTSOB,MXSTBL0)
      Call mma_allocate(I1,LSCR3,Label='I1')
      Call mma_allocate(I2,LSCR3,Label='I2')
      Call mma_allocate(I3,LSCR3,Label='I3')
      Call mma_allocate(I4,LSCR3,Label='I4')
      Call mma_allocate(XI1S,LSCR3,Label='XI1S')
      Call mma_allocate(XI2S,LSCR3,Label='XI2S')
      Call mma_allocate(XI3S,LSCR3,Label='XI3S')
      Call mma_allocate(XI4S,LSCR3,Label='XI4S')

*.Some TTS arrays
      NTTS = MXNTTS
*
*. for partitioning of vector
      Call Allocate_Local_Arrays(NTTS,NSMST)
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,CBLTP,SVST)
*. For scaling for each TTS block
      Call mma_allocate(LSCLFAC,8*NTTS,Label='LSCLFAC')

*. Space for four blocks of string occupations and arrays of
*. reordering arrays
*. Also used to hold an NORB*NORB matrix
      LZSCR = (MAX(NAEL,NBEL)+3)*(NOCOB+1) + 2 * NOCOB + NOCOB*NOCOB
      LZ    = (MAX(NAEL,NBEL)+2) * NOCOB
*. Set up to two blocks for orbital conserving operator
      K12=1
      Call mma_allocate(ZOCSTR,MAX_STR_OC_BLK,K12,Label='ZOCSTR')
      I1234=2
      Call mma_allocate(REO,MAX_STR_SPGP,I1234,Label='REO')
      call mma_allocate(Z,LZ,I1234,Label='Z')
      Call mma_allocate(ZSCR,LZSCR,Label='ZSCR')
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
      Call mma_allocate(H0SPC,NOCTPA*NOCTPB,Label='H0SPC')
      CALL H0INTSPC(   IH0SPC,   NPTSPC, IOCPTSPC,   NOCTPA,   NOCTPB,
     &                ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &                   NGAS,MXPNGAS,H0SPC,NELFGP)
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
      CALL SBLOCKS(NBLOCK,IBLOCK(1,IBOFF),CB,HCB,VEC3,
     &             CIOIO,ISMOST(1,ICSM),CBLTP,
     &             NSTSO(IATP)%I,NSTSO(IBTP)%I,
     &             NAEL,IATP,NBEL,IBTP,
     &             IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &             NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,MXPNGAS,
     &             ITSOB,MAXK,MAXI,LSCR1,
     &             INSCR,VEC3,VEC3(1+LSCR2),
     &             STSTS,STSTD,SXDXSX,
     &             ADSXA,NGAS,NELFSPGP,IDC,
     &             I1,XI1S,I2,XI2S,
     &             IDOH2,MXPOBS,SVST,
     &             PSSIGN,IPRDIA,LUC,ICJKAIB,
     &             VEC3,VEC3(1+LSCR2),
     &             I3,XI3S,I4,XI4S,
     &             MXSXST,MXSXBL,MOCAA,
     &             CLBT,CLEBT,
     &             CI1BT,CIBT,
     &             IRESTRICT,
     &             CONSPA,CONSPB,LSCLFAC,
     &             IPERTOP,IH0INSPC,H0SPC,
     &             ICBAT_RES,ICBAT_INI,ICBAT_END,IUSE_PH,IPHGAS,
     &             I_RES_AB,ISIMSYM,INSCR2)
*
* CALL SBLOCKS --> 91
*
*
      IF(IDC.EQ.2) THEN
*. reform
        CALL RFTTS(HCB,CB,IBLOCK(1,IBOFF),NBLOCK,1,NSMST,
     &             NSTSO(IATP)%I,NSTSO(IBTP)%I,
     &             IDC,PSSIGN,1,NTEST)
*. scale
        CALL SCDTTS(HCB,IBLOCK(1,IBOFF),NBLOCK, NSMST,
     &              NSTSO(IATP)%I,NSTSO(IBTP)%I,
     &              IDC,1,NTEST)
      END IF
*
      IF(LUCBLK.GT.0) THEN
        CALL ITODS([-1],1,-1,LUCBLK)
      END IF
*. Eliminate local memory
      call mma_deallocate(CONSPA)
      call mma_deallocate(CONSPB)
      call mma_deallocate(STSTS)
      call mma_deallocate(STSTD)
      call mma_deallocate(INSCR)
      call mma_deallocate(INSCR2)
      call mma_deallocate(CIOIO)
      call mma_deallocate(SIOIO)
      SCIOIO => Null()
      Call Deallocate_Local_Arrays()
      call mma_deallocate(I1)
      call mma_deallocate(I2)
      call mma_deallocate(I3)
      call mma_deallocate(I4)
      call mma_deallocate(XI1S)
      call mma_deallocate(XI2S)
      call mma_deallocate(XI3S)
      call mma_deallocate(XI4S)
      call mma_deallocate(LSCLFAC)
      call mma_deallocate(ZOCSTR)
      call mma_deallocate(REO)
      call mma_deallocate(Z)
      Call mma_deallocate(ZSCR)
      Call mma_deallocate(SVST)
      Call mma_deallocate(H0SPC)
      END
