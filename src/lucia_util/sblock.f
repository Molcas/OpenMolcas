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
      use GLBBAS
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
#include "oper.fh"
#include "io_util.fh"
*
#include "csmprd.fh"
#include "hidscr.fh"
#include "cintfo.fh"
      DIMENSION CB(*),HCB(*)
      Integer, Allocatable:: CONSPA(:), CONSPB(:)
      Real*8, Allocatable:: INSCR(:), INSCR2(:)
      Integer, Allocatable:: STSTS(:), STSTD(:)
      Integer, Allocatable, Target:: CIOIO(:), SIOIO(:)
      Integer, Pointer:: SCIOIO(:)
      Integer, Allocatable, Target:: CBLTP(:)
      Integer, Allocatable:: I1(:), I2(:), I3(:), I4(:)
      Real*8, Allocatable:: XI1S(:), XI2S(:), XI3S(:), XI4S(:)
      Integer, Allocatable:: LLBT(:)
      Integer, Allocatable:: LLEBT(:)
      Integer, Allocatable:: LI1BT(:)
      Integer, Allocatable:: LIBT(:)
      Real*8, Allocatable:: LSCLFAC(:)
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
* A 4 index integral block with four indices belonging OS class
      INTSCR = MAX(MXTSOB ** 4, NTOOB**2)
      IF(IPRCIX.GE.3)
     &WRITE(6,*) ' Integral scratch space ',INTSCR
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
*. Arrays giving block type
      Call mma_allocate(CBLTP,NSMST,Label='CBLTP')
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
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,CBLTP,iWORK(KSVST))
*.Some TTS arrays
      NTTS = MXNTTS
*
*. for partitioning of vector
      Call mma_allocate(LLBT,NTTS,Label='LLBT')
      Call mma_allocate(LLEBT,NTTS,Label='LLEBT')
      Call mma_allocate(LI1BT,NTTS,Label='LI1BT')
      Call mma_allocate(LIBT,8*NTTS,Label='LIBT')
*. For scaling for each TTS block
      Call mma_allocate(LSCLFAC,8*NTTS,Label='LSCLFAC')

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
      CALL SBLOCKS(NBLOCK,IBLOCK(1,IBOFF),CB,HCB,VEC3,
     &             CIOIO,ISMOST(1,ICSM),CBLTP,
     &             iWORK(KNSTSO(IATP)),iWORK(KNSTSO(IBTP)),
     &             NAEL,IATP,NBEL,IBTP,
     &             IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &             NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,MXPNGAS,
     &             ITSOB,MAXK,MAXI,LSCR1,
     &             INSCR,VEC3,VEC3(1+LSCR2),
     &             STSTS,STSTD,SXDXSX,
     &             ADSXA,NGAS,NELFSPGP,IDC,
     &             I1,XI1S,I2,XI2S,
     &             IDOH2,MXPOBS,iWORK(KSVST),
     &             PSSIGN,IPRDIA,LUC,ICJKAIB,
     &             VEC3,VEC3(1+LSCR2),
     &             I3,XI3S,I4,XI4S,
     &             MXSXST,MXSXBL,MOCAA,
     &             LLBT,LLEBT,
     &             LI1BT,LIBT,
     &             IRESTRICT,
     &             CONSPA,CONSPB,LSCLFAC,
     &             IPERTOP,IH0INSPC,iWORK(KLH0SPC),
     &             ICBAT_RES,ICBAT_INI,ICBAT_END,IUSE_PH,IPHGAS,
     &             I_RES_AB,ISIMSYM,INSCR2)
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
      call mma_deallocate(CONSPA)
      call mma_deallocate(CONSPB)
      call mma_deallocate(STSTS)
      call mma_deallocate(STSTD)
      call mma_deallocate(INSCR)
      call mma_deallocate(INSCR2)
      call mma_deallocate(CIOIO)
      call mma_deallocate(SIOIO)
      SCIOIO => Null()
      call mma_deallocate(CBLTP)
      call mma_deallocate(I1)
      call mma_deallocate(I2)
      call mma_deallocate(I3)
      call mma_deallocate(I4)
      call mma_deallocate(XI1S)
      call mma_deallocate(XI2S)
      call mma_deallocate(XI3S)
      call mma_deallocate(XI4S)
      call mma_deallocate(LSCLFAC)
      call mma_deallocate(LLBT)
      call mma_deallocate(LLEBT)
      call mma_deallocate(LI1BT)
      call mma_deallocate(LIBT)
      CALL GETMEM('KLOCS ','FREE','INTE',KLOCSTR(1),MAX_STR_OC_BLK)
      DO I1234 = 1, 2
        CALL GETMEM('KLREO ','FREE','INTE',KLREO(I1234),MAX_STR_SPGP)
        CALL GETMEM('KLZ   ','FREE','INTE',KLZ(I1234),LZ)
      END DO
      CALL GETMEM('KLZSCR','FREE','INTE',KLZSCR,LZSCR)
      CALL GETMEM('H0SPC ','FREE','INTE',KLH0SPC,NOCTPA*NOCTPB)
      RETURN
      END
