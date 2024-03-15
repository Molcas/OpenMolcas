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
      SUBROUTINE DETCTL_GAS()
      use stdalloc, only: mma_allocate, mma_deallocate
      use GLBBAS, only: SDREO_I, CONF_OCC, CONF_REO
      use Local_Arrays, only: CLBT, CLEBT, CI1BT, CIBT, CBLTP,
     &                        Allocate_Local_Arrays,
     &                      Deallocate_Local_Arrays
      use strbas
      use rasscf_lucia, only: kvec3_length, Memory_Needed_Lucia

*
      IMPLICIT REAL*8 (A-H, O-Z)
#include "mxpdim.fh"
#include "gasstr.fh"
#include "cgas.fh"
#include "orbinp.fh"
#include "crun.fh"
#include "cstate.fh"
#include "cands.fh"
#include "cicisp.fh"
#include "cprnt.fh"
#include "stinf.fh"
#include "csm.fh"
#include "spinfo_lucia.fh"
#include "strinp.fh"
#include "lucinp.fh"

      INTEGER IOCCLS(1),IBASSPC(1)
      Integer, Allocatable:: LCIOIO(:)
      Integer, Allocatable:: SVST(:)
      Integer, Allocatable:: BASSPC(:)
      Integer, Allocatable:: KLOCCLS(:)

*. Set variables in cands.fh
      JSYM = IREFSM
      ICSM  = JSYM
      ISSM  = JSYM
      ICSPC = 1
      ISSPC = 1
*. Set NDET
      NDET = INT(XISPSM(JSYM,1))
*. Number of occupation classes
      IATP = 1
      IBTP = 2
      NEL = NELFTP(IATP)+NELFTP(IBTP)
      IBASSPC(1)=0
      CALL OCCLS(1,NOCCLS,IOCCLS,NEL,NGAS,
     &     IGSOCC(1,1),IGSOCC(1,2),0,IBASSPC,NOBPT)
*. and then the occupation classes
      CALL mma_allocate(KLOCCLS,NGAS*NOCCLS,Label='KLOCCLS')
      CALL mma_allocate(BASSPC,NOCCLS,Label='BASSPC')
      CALL OCCLS(2,NOCCLS,KLOCCLS,NEL,NGAS,
     &     IGSOCC(1,1),IGSOCC(1,2),1,BASSPC,NOBPT)
      CALL mma_deallocate(BASSPC)
C     END IF
      IF(NOCSF.EQ.0) THEN
*. Initial information on CSF expansion
         CALL CSFDIM_GAS(KLOCCLS,NOCCLS,JSYM,IPRCIX)
*. Prototype dets and csf's and CSF'SD matrices
         CALL CSDTMT_GAS(IPRCIX)
      END  IF
*. Allocate memory for diagonalization
c      IF(ISIMSYM.EQ.0) THEN
         LBLOCK = MXSOOB
c      ELSE
c         LBLOCK = MXSOOB_AS
c      END IF
      LBLOCK = MAX(LBLOCK,LCSBLK)
* JESPER : Should reduce I/O
*PAM06      LBLOCK = MAX(XISPSM(IREFSM,1),DBLE(MXSOOB))
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
      Call Allocate_Local_Arrays(NTTS,NSMST)
*. Additional info required to construct partitioning
*
      Call mma_allocate(LCIOIO,NOCTPA*NOCTPB,Label='LCIOIO')
*
      CALL IAIBCM(ICSPC,LCIOIO)
      Call mma_allocate(SVST,1,Label='SVST')
      CALL ZBLTP(ISMOST(1,jsym),NSMST,IDC,CBLTP,SVST)
      Call mma_deallocate(SVST)
*
*. Batches  of C vector
      CALL PART_CIV2(IDC,CBLTP,NSTSO(IATP)%I,NSTSO(IBTP)%I,
     &     NOCTPA,NOCTPB,NSMST,LBLOCK,LCIOIO,
     &     ISMOST(1,jsym),
     &     NBATCH,CLBT,CLEBT,
     &     CI1BT,CIBT,0,ISIMSYM)
*. Number of BLOCKS
      NBLOCK = IFRMR(CI1BT,1,NBATCH)
     &     + IFRMR(CLBT,1,NBATCH) - 1
*. Length of each block
      CALL EXTRROW(CIBT,8,8,NBLOCK,CI1BT)
      IF (NEL .GT. 0)
     &   CALL CNFORD_GAS(KLOCCLS, NOCCLS, jsym, PSSIGN, IPRCIX,
     &       CONF_OCC(JSYM)%I, CONF_REO(jsym)%I,
     &       SDREO_I(jsym)%I,
     &       CIBT, NBLOCK)
*
      Call Deallocate_Local_Arrays()
*. If PICO2/SBLOCK are used, three blocks are used in PICO2, so
*...
*. Sblock is used in general nowadays so,
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
         MAXA1 = IMNMX(NSTSO(IATPM1)%I,NSMST*NOCTYP(IATPM1),2)
         MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
         MAXA1 = IMNMX(NSTSO(IATPM2)%I,NSMST*NOCTYP(IATPM2),2)
         MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
         MAXB1 = IMNMX(NSTSO(IBTPM1)%I,NSMST*NOCTYP(IBTPM1),2)
         MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
         MAXB1 = IMNMX(NSTSO(IBTPM2)%I,NSMST*NOCTYP(IBTPM2),2)
         MAXB = MAX(MAXB,MAXB1)
      END IF
*
      MXSTBL = MAX(MAXA,MAXB,MXSTBL0)
*
      IF(IPRCIX.GE.2 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
      MAXI = MIN(MXINKA,MXSTBL)
      MAXK = MIN(MXINKA,MXSTBL)
*.scratch space for projected matrices and a CI block
*
*. Scratch space for CJKAIB resolution matrices
*. Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
      CALL MXRESCPH(LCIOIO,IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &        NSMST,NSTFSMSPGP,MXPNSMST,
     &        NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,
     &        NELFSPGP,
     &        MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,
     &        IPHGAS,NHLFSPGP,MNHL,IADVICE,MXCJ_ALLSYM,
     &        MXADKBLK_AS,MX_NSPII)

      Call mma_deallocate(LCIOIO)

      IF(IPRCIX.GE.2) THEN
         WRITE(6,*) 'DETCTL : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &           MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
         WRITE(6,*) ' MXADKBLK ,MXADKBLK_AS', MXADKBLK, MXADKBLK_AS
      END IF
c      IF(ISIMSYM.EQ.1) THEN
c         MXCJ = MAX(MXCJ_ALLSYM,MX_NSPII)
c         MXADKBLK_AS = MXADKBLK
c      END IF
*. Using hardwired routines, MXCIJAB also used
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB,MX_NSPII)
      IF(IPRCIX.GE.2)
     &     WRITE(6,*) ' Space for two resolution matrices ',2*LSCR2
      LSCR12 = MAX(LBLOCK,2*LSCR2)
*
      KVEC3_LENGTH = MAX(LSCR12,2*LBLOCK)
*
* Calculate how much memory the sigma routines needs.
*
* Memory needed in Sigma_Master
      MEMORY_NEEDED_LUCIA = KVEC3_LENGTH
* Memory needed in MV7
      MEMORY_NEEDED_LUCIA = MEMORY_NEEDED_LUCIA
     &                    + 2*1 + NOCTPA*NOCTPB + NSMST + 11*MXNTTS
c      IF (IDC .EQ. 3 .OR. IDC .EQ. 4) THEN
c         MEMORY_NEEDED_LUCIA = MEMORY_NEEDED_LUCIA + NSMST
c      END IF
* Memory needed in SBLOCK
*Largest active orbital block belonging to given type and symmetry
      MXTSOB = 0
      DO IOBTP = 1, NGAS
      DO IOBSM = 1, NSMOB
       MXTSOB = MAX(MXTSOB,NOBPTS(IOBTP,IOBSM))
      END DO
      END DO
      INTSCR = MAX(MXTSOB**4,NTOOB**2)
*.vectors able to hold strings of given sym and type
      MAXIK = MAX(MAXI,MAXK)
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB*MXTSOB,MXNSTR)
*. Space for four blocks of string occupations and arrays of
*. reordering arrays
*. Also used to hold an NORB*NORB matrix
      LZSCR = (MAX(NAEL,NBEL)+3)*(NOCOB+1) + 2 * NOCOB + NOCOB*NOCOB
      LZ    = (MAX(NAEL,NBEL)+2) * NOCOB
      MEMORY_NEEDED_LUCIA = MEMORY_NEEDED_LUCIA
     &                    + NOCTPA**2 + NOCTPB**2 + 2*NSMST**2
     &                    + 2*INTSCR + 3*NOCTPA*NOCTPB + NSMST
     &                    + 8*LSCR3 + 12*NTTS + MAX_STR_OC_BLK
     &                    + 2*MAX_STR_SPGP + 2*LZ + LZSCR
c      IF (IDC .EQ. 3 .OR. IDC .EQ. 4) THEN
c         MEMORY_NEEDED_LUCIA = MEMORY_NEEDED_LUCIA + NSMST
c      END IF
c      IF(ISIMSYM.EQ.0) THEN
         MEMORY_NEEDED_LUCIA = MEMORY_NEEDED_LUCIA + 4*MAX_STR_SPGP
c      END IF
*
      CALL LUCIA2MOLCAS(
     &     CONF_OCC(jsym)%I,SDREO_I(jsym)%I,
     &     ndet, ncsf_per_sym, nsd_per_sym, nconf_per_sym, mxpcsm,
     &     mxporb, nconf_per_open, npdtcnf, npcscnf, mults,
     &     nCSF_HEXS)

      CALL mma_deallocate(KLOCCLS)

      END SUBROUTINE DETCTL_GAS
*
      SUBROUTINE DETCTL_FREE()
      use strbas
      IMPLICIT REAL*8 (A-H, O-Z)
#include "mxpdim.fh"
#include "gasstr.fh"
#include "cgas.fh"
#include "orbinp.fh"
#include "crun.fh"
#include "cstate.fh"
#include "cands.fh"
#include "cicisp.fh"
#include "cprnt.fh"
#include "stinf.fh"
#include "csm.fh"
#include "spinfo_lucia.fh"
#include "strinp.fh"
#include "lucinp.fh"


      JSYM = IREFSM
      CALL CSFDIM_FREE(JSYM)

      CALL LUCIA2MOLCAS_FREE()

      END SUBROUTINE DETCTL_FREE
