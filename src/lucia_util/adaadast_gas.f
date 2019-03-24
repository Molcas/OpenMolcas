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
* Copyright (C) 1995,1997, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE ADAADAST_GAS(    IOB,  IOBSM,  IOBTP,   NIOB,    IAC,
     &                            JOB,  JOBSM,  JOBTP,   NJOB,    JAC,
     &                          ISPGP,    ISM,    ITP,   KMIN,   KMAX,
     &                             I1,   XI1S,    LI1,     NK,   IEND,
     &                          IFRST,  KFRST,    I12,    K12, SCLFAC)
*
*
* Obtain two-operator mappings
* a+/a IORB a+/a JORB !KSTR> = +/-!ISTR>
*
* Whether creation- or annihilation operators are in use depends
* upon IAC, JAC : 1=> Annihilation,
*                 2=> Creation
*
* In the form
* I1(KSTR) =  ISTR if a+/a IORB a+/a JORB !KSTR> = +/-!ISTR> , ISTR is in
* ISPGP,ISM,IGRP.
* (numbering relative to TS start)
*. Only excitations IOB. GE. JOB are included
* The orbitals are in GROUP-SYM IOBTP,IOBSM, JOBTP,JOBSM respectively,
* and IOB (JOB) is the first orbital to be used, and the number of orbitals
* to be checked is NIOB ( NJOB).
*
* Only orbital pairs IOB .gt. JOB are included (if the types are identical)
*
* The output is given in I1(KSTR,I,J) = I1 ((KSTR,(J-1)*NIOB + I)
*
* Above +/- is stored in XI1S
* Number of K strings checked is returned in NK
* Only Kstrings with relative numbers from KMIN to KMAX are included
*
* If IEND .ne. 0 last string has been checked
*
* Jeppe Olsen , August of 95   ( adadst)
*               November 1997 : annihilation added
*
*
* ======
*. Input
* ======
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.fh"
*./BIGGY
#include "WrkSpc.fh"
!      COMMON/BIGGY/WORK(MXPWRD)
*./ORBINP/
#include "orbinp.fh"
#include "strinp.fh"
#include "strbas.fh"
#include "cgas.fh"
#include "gasstr.fh"
*. Local scratch
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
      COMMON/SSAVE/NELIS(4), NSTRKS(4)
      COMMON/UMMAGUMMA/NSTRIA(4)
C-jwk-cleanup      INTEGER KELFGRP(MXPNGAS)
      INTEGER KGRP(MXPNGAS)
      DIMENSION IDUM(1)
c      COMMON/COMJEP/MXACJ,MXACIJ,MXAADST
*
* =======
*. Output
* =======
*
      INTEGER I1(*)
      DIMENSION XI1S(*)
* Some dummy initializations
      IIGRP = 0 ! jwk-cleanup
      JJGRP = 0 ! jwk-cleanup
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ====================== '
        WRITE(6,*) ' ADAADST_GAS in service '
        WRITE(6,*) ' ====================== '
        WRITE(6,*)
        WRITE(6,*) ' IOB,IOBSM,IOBTP,IAC ', IOB,IOBSM,IOBTP,IAC
        WRITE(6,*) ' JOB,JOBSM,JOBTP,JAC ', JOB,JOBSM,JOBTP,JAC
        WRITE(6,*) ' I12, K12 ', I12, K12
        WRITE(6,*) ' IFRST,KFRST', IFRST,KFRST
      END IF
*
*
*. Internal affairs
*
      IF(I12.LE.4.AND.K12.LE.1) THEN
        KLLOC = KLOCSTR(K12)
        KLLZ = KLZ(I12)
        KLLREO = KLREO(I12)
      ELSE
        WRITE(6,*)
     &  ' ADST_GAS : Illegal value of I12 or K12 ', I12, K12
*        STOP' ADST_GAS : Illegal value of I12 or K12  '
        CALL SYSABENDMSG('lucia_util/adst_gas',
     &                    'Internal error',' ')
        RETURN
      END IF
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' KLLOC KLLREO',KLLOC,KLLREO
      END IF

*
*. Supergroup and symmetry of K strings
*
      CALL SYMCOM(2,0,IOBSM,K1SM,ISM)
      CALL SYMCOM(2,0,JOBSM,KSM,K1SM)
      IF(NTEST.GE.100) WRITE(6,*) ' K1SM,KSM : ',  K1SM,KSM
      ISPGPABS = IBSPGPFTP(ITP)-1+ISPGP
      IACADJ = 2
      IDELTA =-1
      IF(IAC.EQ.2) THEN
        IACADJ = 1
        IDELTA = 1
      END IF
      JACADJ = 2
      JDELTA =-1
      IF(JAC.EQ.2) THEN
        JACADJ = 1
        JDELTA = 1
      END IF
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' IACADJ, JACADJ', IACADJ,JACADJ
       WRITE(6,*) ' IDELTA, JDELTA', IDELTA, JDELTA
      END IF
*. Occupation of K-strings
      IF(IOBTP.EQ.JOBTP) THEN
        IEL = NELFSPGP(IOBTP,ISPGPABS)-IDELTA-JDELTA
        JEL = IEL
      ELSE
        IEL = NELFSPGP(IOBTP,ISPGPABS)-IDELTA
        JEL = NELFSPGP(JOBTP,ISPGPABS)-JDELTA
      END IF
      IF(NTEST.GE.100) WRITE(6,*) ' IEL, JEL', IEL,JEL
*. Trivial zero ? (Nice, then mission is complete )
      ITRIVIAL = 0
      IF(IEL.LT.0.OR.JEL.LT.0.OR.
     &   IEL.GT.NOBPT(IOBTP).OR.JEL.GT.NOBPT(JOBTP)) THEN
*. No strings with this number of elecs - be happy : No work
        NK = 0
        KACT = 0
        KACGRP = 0
        IF(NTEST.GE.100) WRITE(6,*) ' Trivial zero excitations'
        ITRIVIAL = 1
C       GOTO 9999
      ELSE
*. Find group with IEL electrons in IOBTP, JEL in JOBTP
        IIGRP = 0
        DO IGRP = IBGPSTR(IOBTP),IBGPSTR(IOBTP)+NGPSTR(IOBTP)-1
          IF(NELFGP(IGRP).EQ.IEL) IIGRP = IGRP
        END DO
        JJGRP = 0
        DO JGRP = IBGPSTR(JOBTP),IBGPSTR(JOBTP)+NGPSTR(JOBTP)-1
          IF(NELFGP(JGRP).EQ.JEL) JJGRP = JGRP
        END DO
C?      WRITE(6,*) ' ADAADA : IIGRP, JJGRP', IIGRP,JJGRP
*
        IF(IIGRP.EQ.0.OR.JJGRP.EQ.0) THEN
          WRITE(6,*)' ADAADAST : cul de sac, active K groups not found'
          WRITE(6,*)' Active GAS spaces  ' ,IOBTP, JOBTP
          WRITE(6,*)' Number of electrons', IEL, JEL
*          STOP      ' ADAADAST : cul de sac, active K groups not found'
        CALL SYSABENDMSG('lucia_util/adaadast_gas',
     &                    'Internal error',' ')
        END IF
*
      END IF
*. Groups defining Kstrings
      IF(ITRIVIAL.NE.1) THEN
      CALL ICOPVE(ISPGPFTP(1,ISPGPABS),KGRP,NGAS)
      KGRP(IOBTP) = IIGRP
      KGRP(JOBTP) = JJGRP
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Groups in KGRP '
        CALL IWRTMA(KGRP,1,NGAS,1,NGAS)
      END IF
      END IF
*
* In ADADS1_GAS we need : Occupation of KSTRINGS
*                         lexical => Actual order for I strings
* Generate if required
*
      IF(IFRST.NE.0) THEN
*.. Generate information about I strings
*. Arc weights for ISPGP
        NTEST2 = NTEST
        CALL WEIGHT_SPGP(iWORK(KLLZ),
     &                       NGAS,
     &                   NELFSPGP(1,ISPGPABS),
     &                      NOBPT,
     &                   iWORK(KLZSCR),
*
     &                     NTEST2)
        NELI = NELFTP(ITP)
        NELIS(I12) = NELI
*. Reorder array for I strings
        CALL GETSTR_TOTSM_SPGP(    ITP,  ISPGP,    ISM,   NELI,  NSTRI,
     &                         IWORK(KLLOC),
     &                           NOCOB,
     &                               1,
     &                         IWORK(KLLZ),
*
     &                         IWORK(KLLREO))
        IF(NTEST.GE.1000) THEN
         write(6,*) ' Info on I strings generated '
         write(6,*) ' NSTRI = ', NSTRI
         WRITE(6,*) ' REORDER array '
         CALL IWRTMA(iWORK(KLLREO),1,NSTRI,1,NSTRI)
       END IF
       NSTRIA(I12) = NSTRI
*
      END IF
      IF(NTEST.GE.1000) THEN
       WRITE(6,*) ' REORDER array for I STRINGS'
       CALL IWRTMA(iWORK(KLLREO),1,NSTRI,1,NSTRI)
      END IF
*
      IF(ITRIVIAL.EQ.1) GOTO 9999
      NELK = NELIS(I12)
      IF(IAC.EQ.1) THEN
        NELK = NELK + 1
      ELSE
        NELK = NELK - 1
      END IF
      IF(JAC.EQ.1) THEN
        NELK = NELK + 1
      ELSE
        NELK = NELK - 1
      END IF
      IF(NTEST.GE.100) WRITE(6,*) ' NELK = ' , NELK
      IF(KFRST.NE.0) THEN
*. Generate occupation of K STRINGS
       IDUM(1)=0
       CALL GETSTR2_TOTSM_SPGP(   KGRP,   NGAS,    KSM,   NELK,  NSTRK,
     &                         iWORK(KLLOC),NOCOB,    0, IDUM, IDUM)
C     GETSTR2_TOTSM_SPGP(IGRP,NIGRP,ISPGRPSM,NEL,NSTR,ISTR,
C    &                              NORBT,IDOREO,IZ,IREO)
       NSTRKS(K12) = NSTRK
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' K strings generated '
         WRITE(6,*) ' Reorder array after generation of K strings'
         CALL IWRTMA(iWORK(KLLREO),1,NSTRI,1,NSTRI)
       END IF
      END IF
*
      NSTRK = NSTRKS(K12)
*
      IIOB = IOBPTS(IOBTP,IOBSM) + IOB - 1
      JJOB = IOBPTS(JOBTP,JOBSM) + JOB - 1
*
      IZERO = 0
      ZERO = 0.0D0
      CALL ISETVC(I1  ,IZERO,LI1*NIOB*NJOB)
COLD  CALL SETVEC(XI1S,ZERO ,LI1*NIOB*NJOB)
*
      CALL ADAADAS1_GAS(      NK,      I1,    XI1S,     LI1,    IIOB,
     &                      NIOB,     IAC,    JJOB,    NJOB,     JAC,
     &                  iWORK(KLLOC), NELK,NSTRK,
     &                  iWORK(KLLREO),iWORK(KLLZ),
     &                     NOCOB,    KMAX,    KMIN,    IEND,  SCLFAC,
     &                  NSTRIA(I12))
*
*
       IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Reorder array after ADAADAS1'
         CALL IWRTMA(WORK(KLLREO),1,NSTRI,1,NSTRI)
       END IF
 9999 CONTINUE
*
C     WRITE(6,*) ' Memcheck at end of ADAADAS1 '
C     CALL MEMCHK
      RETURN
      END
