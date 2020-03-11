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
* Copyright (C) 1988,1991,1992,1998, Per Ake Malmqvist                 *
************************************************************************
*--------------------------------------------*
* 1998  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE POLY3(IFF)
      IMPLICIT NONE
C  IBM TEST VERSION 0, 1988-06-23.
C  NEW VERSION 1991-02-23, FOR USE WITH RASSCF IN MOLCAS PACKAGE.
C  NEW VERSION 1992-12-05, FOR MOLCAS-3 VERSION.
C  NEW VERSION 1998-10-02
C  AUTHOR: PER-AAKE MALMQUIST

C THIS PROGRAM CALCULATES 1-EL, 2-EL, AND 3-EL
C DENSITY MATRICES FOR A CASSCF WAVE FUNCTION.
C IF THE INTEGER KEY IFF.EQ.1, THEN
C IT ALSO PRODUCES THE CONTRACTIONS OF 1-EL -- 4-EL
C DENSITY MATRICES WITH THE FOCK OPERATOR USED IN
C THE CASSCF-MP2 PROGRAM. THE RESULTS ARE WRITTEN
C TO FILE IN SEVERAL FORMS, TO SUPPORT BOTH KERSTINS
C PRESENT PROGRAM AND ALSO SUCH NEW PROCEDURES WHICH
C MIGHT TAKE ADVANTAGE OF ALL INDEX PERMUTATION SYMMETRIES.
C THE RDSTAT AND THE GUGA ROUTINES USED IN THIS
C PROGRAM ASSUMES THE JOBIPH IS PRODUCED BY THE RASSCF PROGRAM.

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "pt2_guga.fh"
#include "SysDef.fh"

      INTEGER IFF

      INTEGER ILEV
      INTEGER NG3MAX,IPAD,LIDXG3
      INTEGER ILUID

      INTEGER IDCI
      INTEGER J

      INTEGER IPARDIV

      CALL QENTER('POLY3')

      IF (IFF.EQ.1) THEN
C ORBITAL ENERGIES IN CI-COUPLING ORDER:
        DO ILEV=1,NLEV
          ETA(ILEV)=EPSA(L2ACT(ILEV))
        END DO
      END IF

      CALL GETMEM('G1','ALLO','REAL',LG1,NG1)
      CALL GETMEM('G2','ALLO','REAL',LG2,NG2)

C-SVC20100831: recompute approximate max NG3 size needed
      NG3MAX=iPARDIV(NG3TOT,NG2)

C-SVC20100831: allocate local G3 matrices
      CALL GETMEM('G3','ALLO','REAL',LG3,NG3MAX)
      iPad=ItoB-MOD(6*NG3MAX,ItoB)
      CALL GETMEM('idxG3','ALLO','CHAR',LidxG3,6*NG3MAX+iPad)

      WORK(LG1)=0.0D0
      WORK(LG2)=0.0D0
      WORK(LG3)=0.0D0

C ALLOCATE SPACE FOR CORRESPONDING COMBINATIONS WITH H0:
      IF (IFF.EQ.1) THEN
        CALL GETMEM('LF1','ALLO','REAL',LF1,NG1)
        CALL GETMEM('LF2','ALLO','REAL',LF2,NG2)
        CALL GETMEM('LF3','ALLO','REAL',LF3,NG3MAX)
      ELSE
        LF1=LG1
        LF2=LG2
        LF3=LG3
      END IF

* NG3 will change inside subroutine MKFG3 to the actual
* number of nonzero elements, that is why here we allocate
* with NG3MAX, but we only store (PT2_PUT) the first NG3
* elements of the G3 and F3
      NG3=NG3MAX

      CALL GETMEM('LCI','ALLO','REAL',LCI,NCONF)

      IF (.NOT.DoCumulant.AND.ISCF.EQ.0) THEN
        IDCI=IDTCEX
        DO J=1,JSTATE-1
          CALL DDAFILE(LUCIEX,0,WORK(LCI),NCONF,IDCI)
        END DO
        CALL DDAFILE(LUCIEX,2,WORK(LCI),NCONF,IDCI)
        IF (IPRGLB.GE.VERBOSE) THEN
          WRITE(6,*)
          IF (NSTATE.GT.1) THEN
            WRITE(6,'(A,I4)')
     &      ' With new orbitals, the CI array of state ',MSTATE(JSTATE)
          ELSE
            WRITE(6,*)' With new orbitals, the CI array is:'
          END IF
          CALL PRWF_CP2(LSYM,NCONF,WORK(LCI),CITHR)
        END IF
      ELSE
        WORK(LCI)=1.0D0
      END IF

      IF (ISCF.NE.0.AND.NACTEL.NE.0) THEN
        CALL SPECIAL( WORK(LG1),WORK(LG2),WORK(LG3),
     &                WORK(LF1),WORK(LF2),WORK(LF3),
     &                i1WORK(LidxG3))
      ELSE IF (ISCF.EQ.0) THEN
C-SVC20100903: during mkfg3, NG3 is set to the actual value
#if defined _ENABLE_BLOCK_DMRG_ || defined _ENABLE_CHEMPS2_DMRG_
        IF (.NOT.DoCumulant) THEN
#endif
          CALL MKFG3(IFF,WORK(LCI),WORK(LG1),WORK(LF1),WORK(LG2),
     &               WORK(LF2),WORK(LG3),WORK(LF3),i1WORK(LidxG3))
#if defined _ENABLE_BLOCK_DMRG_ || defined _ENABLE_CHEMPS2_DMRG_
        ELSE
          CALL MKFG3DM(IFF,WORK(LG1),WORK(LF1),WORK(LG2),WORK(LF2),
     &                       WORK(LG3),WORK(LF3),i1WORK(LidxG3))
        END IF
#endif
      END IF

      CALL GETMEM('LCI','FREE','REAL',LCI,NCONF)

      IF(NLEV.GT.0) THEN
        CALL PT2_PUT(NG1,' GAMMA1',WORK(LG1))
        CALL PT2_PUT(NG2,' GAMMA2',WORK(LG2))
        CALL PT2_PUT(NG3,' GAMMA3',WORK(LG3))
        iLUID=0
        iPad=ItoB-MOD(6*NG3,ItoB)
        CALL CDAFILE(LUSOLV,1,cWORK(LidxG3),6*NG3+iPad,iLUID)
        IF(IFF.EQ.1) THEN
          CALL PT2_PUT(NG1,' DELTA1',WORK(LF1))
          CALL PT2_PUT(NG2,' DELTA2',WORK(LF2))
          CALL PT2_PUT(NG3,' DELTA3',WORK(LF3))
        END IF
      END IF

      IF(NLEV.GT.0) THEN
        CALL GETMEM('LG1','FREE','REAL',LG1,NG1)
        CALL GETMEM('LG2','FREE','REAL',LG2,NG2)
        CALL GETMEM('LG3','FREE','REAL',LG3,NG3MAX)
        iPad=ItoB-MOD(6*NG3MAX,ItoB)
        CALL GETMEM('idxG3','FREE','CHAR',LidxG3,6*NG3MAX+iPad)
        IF(IFF.EQ.1) THEN
          CALL GETMEM('LF1','FREE','REAL',LF1,NG1)
          CALL GETMEM('LF2','FREE','REAL',LF2,NG2)
          CALL GETMEM('LF3','FREE','REAL',LF3,NG3MAX)
        END IF
      END IF

      CALL QEXIT('POLY3')

      RETURN
      END
