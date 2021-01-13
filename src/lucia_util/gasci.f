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
* Copyright (C) 1995, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE GASCI(     ISM,    ISPC,   IPRNT,    EREF,IIUSEH0P,
     &                 MPORENP_E)
*
* CI optimization in GAS space number ISPC for symmetry ISM
*
*
* Jeppe Olsen, Winter of 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL MV7
#include "mxpdim.fh"
#include "cicisp.fh"
#include "WrkSpc.fh"
#include "orbinp.fh"
#include "clunit.fh"
#include "csm.fh"
#include "cstate.fh"
#include "crun.fh"
#include "strinp.fh"
#include "stinf.fh"
#include "strbas.fh"
#include "glbbas.fh"
#include "cprnt.fh"
#include "oper.fh"
#include "gasstr.fh"
#include "cgas.fh"
#include "lucinp.fh"
#include "intform.fh"


#include "cintfo.fh"
#include "spinfo_lucia.fh"
#include "rasscf_lucia.fh"
#include "io_util.fh"
*
*. Common block for communicating with sigma
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*
#include "cecore.fh"
      COMMON/CMXCJ/MXCJ,MAXK1_MX,LSCMAX_MX
*
      COMMON/H_OCC_CONS/IH_OCC_CONS
*
      INTEGER IOCCLS_ARR(1), ZERO_ARR(1)
*
*. Should all parameters be tranfered to Molcas?
c      PARAMETER (IALL = 0)
*
      NTEST = 1
      NTEST = MAX(NTEST,IPRNT)
c      MXACJ = 0
c      MXACIJ = 0
c      MXAADST = 0
*. Normal integrals accessed
      IH1FORM = 1
      I_RES_AB = 0
      IH2FORM = 1
*. CI not CC
*. Not just number conserving part
c      IH_OCC_CONS_TEST = 0
c      IF(IH_OCC_CONS_TEST.EQ.1) THEN
c         WRITE(6,*) ' IH_OCC_CONS set to one in GASCI '
c         WRITE(6,*) ' IH_OCC_CONS set to one in GASCI '
c         WRITE(6,*) ' IH_OCC_CONS set to one in GASCI '
c         WRITE(6,*) ' IH_OCC_CONS set to one in GASCI '
c         WRITE(6,*) ' IH_OCC_CONS set to one in GASCI '
c         WRITE(6,*) ' IH_OCC_CONS set to one in GASCI '
c         WRITE(6,*) ' IH_OCC_CONS set to one in GASCI '
c         WRITE(6,*) ' IH_OCC_CONS set to one in GASCI '
c         WRITE(6,*) ' IH_OCC_CONS set to one in GASCI '
c         IH_OCC_CONS = 1
c      END IF
*
      IF(NTEST .GE. 20) THEN
         WRITE(6,*)
         WRITE(6,*) ' ====================================='
         WRITE(6,*) ' Control has been transferred to GASCI'
         WRITE(6,*) ' ====================================='
         WRITE(6,*)
         WRITE(6,*) ' IIUSEH0P = ', IIUSEH0P
         WRITE(6,*) ' MPORENP_E = ', MPORENP_E
      END IF
      IF(NTEST.GE.5) THEN
         WRITE(6,'(A)') '  A few pertinent data : '
         WRITE(6,*)
         WRITE(6,'(A,I2)') '  CI space         ',ISPC
         WRITE(6,*)
         WRITE(6,*) ' Number of GAS spaces included ',LCMBSPC(ISPC)
         WRITE(6,'(A,10I3)') '  GAS spaces included           ',
     &        (ICMBSPC(II,ISPC),II=1,LCMBSPC(ISPC))
         WRITE(6,*)
         WRITE(6,*) ' Occupation constraints : '
         WRITE(6,*) '========================= '
         WRITE(6,*)
         WRITE(6,*)
         DO JJGASSPC = 1, LCMBSPC(ISPC)
            JGASSPC = ICMBSPC(JJGASSPC,ISPC)
            WRITE(6,*)
     &           ' Gas space  Min acc. occupation Max acc. occupation '
            WRITE(6,*)
     &           ' ================================================== '
            DO IGAS = 1, NGAS
               WRITE(6,'(3X,I2,13X,I3,16X,I3)') IGAS,
     &              IGSOCCX(IGAS,1,JGASSPC),IGSOCCX(IGAS,2,JGASSPC)
            END DO
         END DO
*
      END IF
*
      NDET = INT(XISPSM(ISM,ISPC))
      NEL = NELCI(ISPC)
        IF (NTEST .GE. 20)
     &       WRITE(6,*) ' Number of determinants/combinations  ',NDET
      IF(NDET.EQ.0) THEN
         WRITE(6,*) ' The number of determinants/combinations is zero.'
         WRITE(6,*) ' I am sure that fascinating discussions about '
         WRITE(6,*) ' the energy of such a wave function exists, '
         WRITE(6,*) ' but I am just a dumb program, so I will stop'
         WRITE(6,*)
         WRITE(6,*) ' GASCI : Vanishing number of parameters '
*         STOP       ' GASCI : Vanishing number of parameters '
         CALL SYSABENDMSG('lucia_util/gasci','User error',' ')
      END IF
*.Transfer to CANDS
      ICSM = ISM
      ISSM = ISM
      ICSPC = ISPC
      ISSPC = ISPC
*. Complete operator
      I12 = 2
*. Class info
C     IF(ICLSSEL.EQ.1) THEN
*. Number of occupation classes
      IATP = 1
      IBTP = 2
      NEL = NELFTP(IATP)+NELFTP(IBTP)
      ZERO_ARR(1)=0
      CALL OCCLS(         1,    NOCCLS,IOCCLS_ARR,       NEL,      NGAS,
     &           IGSOCC(1,1),IGSOCC(1,2),       0,ZERO_ARR,   NOBPT)
*. and then the occupation classes
*
      IF(NOCSF.EQ.1) THEN
         NVAR = NDET
      ELSE
C*JESPER : Addition start
C        NVAR = NCSASM(ISM)
         NVAR = NCSF_PER_SYM(ISM)
C*JESPER : Addition end
      END  IF
      IF(IPRNT.GE.5) WRITE(6,*) '  NVAR in GASCI ', NVAR
*. Allocate memory for diagonalization
c      IF(ISIMSYM.EQ.0) THEN
         LBLOCK = MXSOOB
c      ELSE
c         LBLOCK = MXSOOB_AS
c      END IF
      LBLOCK = MAX(LBLOCK,LCSBLK)
* JESPER : Should reduce I/O
         LBLOCK = MAX(INT(XISPSM(IREFSM,1)),MXSOOB)
         IF(PSSIGN.NE.0.0D0) LBLOCK = INT(2.0D0*XISPSM(IREFSM,1))
*
*. Information about block structure- needed by new PICO2 routine.
*. Memory for partitioning of C vector
      IATP = 1
      IBTP = 2
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
      NTTS = MXNTTS
      CALL GETMEM('CLBT  ','ALLO','INTE',KLCLBT ,NTTS  )
      CALL GETMEM('CLEBT ','ALLO','INTE',KLCLEBT ,NTTS  )
      CALL GETMEM('CI1BT ','ALLO','INTE',KLCI1BT,NTTS  )
      CALL GETMEM('CIBT  ','ALLO','INTE',KLCIBT ,8*NTTS)
*. Additional info required to construct partitioning
      CALL GETMEM('CIOIO ','ALLO','INTE',KLCIOIO,NOCTPA*NOCTPB)
      CALL GETMEM('CBLTP ','ALLO','INTE',KLCBLTP,NSMST)
*
      CALL IAIBCM(ISPC,IWORK(KLCIOIO))
*. option KSVST not active so
      KSVST = 1
      CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,IWORK(KLCBLTP),IWORK(KSVST))
*
*. Batches  of C vector
      CALL PART_CIV2(      IDC,
     &               IWORK(KLCBLTP),
     &               IWORK(KNSTSO(IATP)),
     &               IWORK(KNSTSO(IBTP)),NOCTPA,NOCTPB, NSMST,LBLOCK,
     &               IWORK(KLCIOIO),
*
     &               ISMOST(1,ISM),
     &                  NBATCH,
     &               IWORK(KLCLBT),
     &               IWORK(KLCLEBT),IWORK(KLCI1BT),IWORK(KLCIBT),
     &               0,ISIMSYM)
*. Number of BLOCKS
      NBLOCK = IFRMR(IWORK(KLCI1BT),1,NBATCH)
     &       + IFRMR(IWORK(KLCLBT),1,NBATCH) - 1
*.
*. Enabling the calculation of excited states in a new way. Lasse
*. This can be realized for GAS1 to GASN (for the GAS version)
*. Here we find which type should be eliminated in the diagonal and
*. in the sigma vector calculation. This to satisfy RASSI.
*. Insert if statement
      IF(I_ELIMINATE_GAS.GE.1) THEN
        CALL I_AM_SO_EXCITED(NBATCH,IWORK(KLCIBT),IWORK(KLCLBT),
     &                       IWORK(KLCI1BT))
      END IF
*. End of story for Lasse
*.
      CALL GETMEM('CLBT  ','FREE','INTE',KLCLBT ,NTTS  )
      CALL GETMEM('CLEBT ','FREE','INTE',KLCLEBT ,NTTS  )
*. Length of each block
      CALL EXTRROW(IWORK(KLCIBT),8,8,NBLOCK,IWORK(KLCI1BT))
      CALL GETMEM('CI1BT ','FREE','INTE',KLCI1BT,NTTS  )
*.
*. Class divisions of  dets
*. ( Well, in principle I am against class division, but I
*    realize that it is a fact so ...)
*
*. If PICO2/SBLOCK are used, three blocks are used in PICO2, so
*...
*. Largest block of strings in zero order space
      MXSTBL0 = MXNSTR
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
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*. Largest number of strings of given symmetry and type
      MAXA = 0
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(IWORK(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
C?          WRITE(6,*) ' MAXA1 1', MAXA1
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(IWORK(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(IWORK(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(IWORK(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MXSTBL = MAX(MAXA,MAXB,MXSTBL0)
      IF(IPRCIX.GE.2 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
      MAXK = MIN( MXINKA,MXSTBL)
*.scratch space for projected matrices and a CI block
*
*. Scratch space for CJKAIB resolution matrices
*. Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
      CALL MXRESCPH(IWORK(KLCIOIO),IOCTPA,IOCTPB, NOCTPA, NOCTPB,
     &                    NSMST,NSTFSMSPGP,MXPNSMST, NSMOB,MXPNGAS,
     &                     NGAS,  NOBPTS,  IPRCIX,    MAXK,NELFSPGP,
     &                     MXCJ,  MXCIJA,  MXCIJB, MXCIJAB,  MXSXBL,
     &                 MXADKBLK,  IPHGAS,NHLFSPGP,    MNHL, IADVICE,
*
     &                 MXCJ_ALLSYM,MXADKBLK_AS,MX_NSPII)
      IF(IPRCIX.GE.2) THEN
        WRITE(6,*) 'GASCI  : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &           MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
        WRITE(6,*) ' MXADKBLK ,MXADKBLK_AS', MXADKBLK, MXADKBLK_AS
      END IF
c         IF(ISIMSYM.EQ.1) THEN
c            MXCJ = MAX(MXCJ_ALLSYM,MX_NSPII)
c            MXADKBLK_AS = MXADKBLK
c         END IF
*. Using hardwired routines, MXCIJAB also used
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB,MX_NSPII)
      IF(IPRCIX.GE.2)
     &     WRITE(6,*) ' Space for two resolution matrices ',2*LSCR2
      LSCR12 = MAX(LBLOCK,2*LSCR2)
CSVC: is KVEC3 used at all before it is deallocated again?
      CALL GETMEM('KC2   ','ALLO','REAL',KVEC3,LSCR12)
      KVEC1 = KCI_POINTER
c     KVEC2 = KSIGMA_POINTER
      KVEC3_LENGTH = MAX(LSCR12,2*LBLOCK,KVEC3_LENGTH)
*
*. CI diagonal - if required
*
      IF(IDIAG.EQ.2) THEN
         LUDIA = LUSC1
      END IF
      IF(.NOT.(IDIAG.EQ.2.AND.IRESTR.EQ.1)) THEN
         IF(ICISTR.GE.2) IDISK(LUDIA)=0
         I12 = 2
         SHIFT = ECORE_ORIG-ECORE
         CALL GASDIAT(WORK(KVEC1),  LUDIA,  SHIFT, ICISTR,    I12,
     &                IWORK(KLCBLTP),NBLOCK,IWORK(KLCIBT))
*
c         IF(IIUSEH0P.EQ.1) THEN
c*. Diagonal with F
c            CALL SWAPVE(WORK(KFI),WORK(KINT1O),NINT1)
c            CALL GASDIAT(WORK(KVEC1),LUSC52,SHIFT,ICISTR,1,
c     &           WORK(KLCBLTP),NBLOCK,WORK(KLCIBT))
c            CALL SWAPVE(WORK(KFI),WORK(KINT1O),NINT1)
c*. diag of (1-Lambda) F + Lambda H
c            FAC1 = 1.0D0 - XLAMBDA
c            FAC2 = XLAMBDA
c            CALL VECSMD(WORK(KVEC1),WORK(KVEC2),FAC1,FAC2,
c     &           LUSC52,LUDIA,LUSC53,1,LBLK)
c            CALL COPVCD(LUSC53,LUDIA,WORK(KVEC1),1,LBLK)
c         END IF
         IF(NOCSF.EQ.1.AND.ICISTR.EQ.1) THEN
            IDISK(LUDIA)=0
            CALL TODSC(WORK(KVEC1),NVAR,-1,LUDIA)
         END IF
         IF(IPRCIX.GE.2) WRITE(6,*) ' Diagonal constructed  '
      ELSE
         WRITE(6,*) ' Diagonal not calculated '
      END IF

      IDUMMY=1
      CALL GETMEM('CIOIO ','FREE','INTE',KLCIOIO,NOCTPA*NOCTPB)
      CALL GETMEM('CBLTP ','FREE','INTE',KLCBLTP,NSMST)
      CALL GETMEM('CIBT  ','FREE','INTE',KLCIBT ,8*NTTS)
      CALL GETMEM('KC2   ','FREE','REAL',KVEC3,LSCR12)
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_real(EREF)
      END
