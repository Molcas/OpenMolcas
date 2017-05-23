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
* Copyright (C) 1989, Jeppe Olsen                                      *
*               1989, Markus P. Fuelscher                              *
************************************************************************
      SUBROUTINE PHPCSF(PHP,IPCSF,IPCNF,MXPDIM,
     *                  DTOC,IPRODT,ICONF,
     *                  IREFSM,ONEBOD,ECORE,NACTOB,
     *                  SCR,NCONF,NEL,NAEL,NBEL,NPCSF,NPCNF,
     *                  DIAG,TUVX,NTEST,ExFac,IREOTS)
C
C     Obtain primary subspace and obtain
C     explicit representation of hamilton matrix in subspace
C
C     ARGUMENTS :
C     ===========
C     PHP    : hamilton matrix in subspace ( output)
C     IPCSF  : CSF's defining subspace (output)
C     IPCNF  : Configurations defining subspace (output )
C     MXPDIM : Largest allowed dimension of subspace (Input)
C     DTOC   : Transformation matrix between CSF's and DET's (input)
C     IPRODT : Prototype determinants (input)
C     ICONF  : List of configurations  (input)
C     IREFSM : symmetry of considered CI space (input)
C     Onebod : one body hamilton matrix in rectangular form (input)
C     ECORE  : Core energy (input)
C     NACTOB : Number of active orbitals (input)
C     SCR    : Scratch array of length ????
C     NCONF  : Number of configurations of symmetry IREFSM
C     NPCNF  : Number of primary configurations obtained (output)
C     NPCSF  : Number of primary CSF's obtained  (OUTPUT)
C     TUVX   : Two-electron integrals (MO space)
C     DIAG   : Hamilton diagonal over CSF's
*
*     IREOTS : Type => symmetry reordering array
C
C     Jeppe Olsen , Summer of '89
C     adapted to DETRAS by M.P. Fuelscher, Oktober 1989
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "spinfo.fh"
#include "warnings.fh"
C
      DIMENSION PHP(*),IPCSF(*),IPCNF(*),DIAG(*)
      DIMENSION DTOC(*),IPRODT(*),ICONF(*)
      DIMENSION ONEBOD(*),SCR(*)
      DIMENSION TUVX(*), IREOTS(*)
*
*     Assumed machine accuray (give and take)
*
      Acc=1.0D-13

* construct the diagonal of the Hamilton matrix in CNF basis
      ICSFMN = 0
      IICNF  = 1
      IICSF  = 1
      DO ITYP = 1, NTYP
        NJCNF = NCNFTP(ITYP,IREFSM)
        NIRREP = NCSFTP(ITYP)
        DO ICNF = 1, NJCNF
          SCR(IICNF) = DIAG(IICSF)
          IICNF = IICNF + 1
          IICSF = IICSF + NIRREP
        END DO
      END DO
      If (NTEST.GE.30) Call RecPrt('SCR',' ',SCR,1,IICNF-1)

* find the elements of the subspace
* (MXPDIM elements smallest in energy)

      XMAX = FNDMNX(SCR,NCONF,2)
      NPCSF = 0
      NPCNF = 0
      IFINIT = 0
400   CONTINUE
      XMIN = XMAX + 1.0D0
      IMIN = 0
      IICNF = 1
      IICSF = 1
      DO ITYP = 1, NTYP
        NJCNF = NCNFTP(ITYP,IREFSM)
        NIRREP = NCSFTP(ITYP)
        DO ICNF = 1,NJCNF
          IF( SCR(IICNF)+Acc.LT.XMIN ) THEN
            XMIN = SCR(IICNF)
            IMIN = IICNF
            ICSFMN = IICSF
            NCSFMN = NIRREP
          END IF
          IICNF = IICNF + 1
          IICSF = IICSF + NIRREP
        END DO
      END DO
      IF( (NPCSF+NCSFMN).LE. MXPDIM ) THEN
*       add new configuration
        NPCNF = NPCNF + 1
        IPCNF(NPCNF) = IMIN
        CALL ISTVC2(IPCSF(NPCSF+1),ICSFMN-1,1,NCSFMN)
        NPCSF = NPCSF + NCSFMN
        SCR(IMIN) = XMAX + 1.0D0
      ELSE
* No space for this configuration , remove previous
* configurations with the same diagonal value
* PAM 2011: Do not remove any previous configurations.
        IFINIT = 1
*        IICNF = NPCNF+1
*600     CONTINUE
*        IICNF = IICNF - 1
*        DIAVAL = SCR(IPCNF(IICNF))
*        IF( ABS(DIAVAL-XMIN).LE.1.0D-10) THEN
*          NPCNF = NPCNF -1
*          CALL GETCNF_LUCIA(SCR(NCONF+1),ITYP,IPCNF(IICNF),
*     &                ICONF,IREFSM,NEL)
*          NPCSF = NPCSF - NCSFTP(ITYP)
*          GOTO 600
*        END IF
      END IF
      IF( (IFINIT.EQ.0) .AND. (NPCNF.LT.NCONF) ) GOTO 400

      IF(NPCNF.eq.0) THEN
        Call WarningMessage(2,'Making explicit Hamiltonian failed.')
        WRITE(6,*)' An unforeseen catastrophic failure occurred'
        WRITE(6,*)' in the CI solver. The size of the explicit'
        WRITE(6,*)' part of the CI Hamiltonian matrix was not'
        WRITE(6,*)' sufficient. Suggested fix: Change the size'
        WRITE(6,*)' by adding ''SDAV=XXXXX'' to the rasscf input.'
        WRITE(6,*)'  XXXXX is some integer at least ',MXPDIM+NCSFMN
        WRITE(6,*)' Sorry about this. Consider telling the'
        WRITE(6,*)' Molcas group about this failure.'
        Call Quit(_RC_INTERNAL_ERROR_)
      END IF
*
      IF( NTEST.GE.30 ) THEN
        WRITE(6,*) ' Output from PHPCSF '
        WRITE(6,*) ' ================== '
        WRITE(6,*)
     &  ' Number of Configurations in primary subspace ',NPCNF
        WRITE(6,*)
     &  ' Number of CSFs in primary subspace ',NPCSF
        WRITE(6,*) ' Configurations included : '
        CALL IWRTMA(IPCNF,1,NPCNF,1,NPCNF)
        WRITE(6,*) ' CSFs included : '
        CALL IWRTMA(IPCSF,1,NPCSF,1,NPCSF)
      END IF

* construct Hamiltonian matrix in subspace

      MXCSFC = 0
      DO ITYP = 1, NTYP
        MXCSFC = MAX(MXCSFC,NCSFTP(ITYP))
      END DO
*
      KLFREE = 1
      KLCONF = KLFREE
      KLFREE = KLFREE+NEL
      KRCONF = KLFREE
      KLFREE = KLFREE+NEL
      KLPHPS = KLFREE
      KLFREE = KLFREE+MXCSFC*MXCSFC
*
      IILB = 1
      DO ICNL = 1, NPCNF
        CALL GETCNF_LUCIA(SCR(KLCONF),ILTYP,IPCNF(ICNL),ICONF,IREFSM,
     &                    NEL)
        NCSFL = NCSFTP(ILTYP)
        IIRB = 1
        DO ICNR = 1, ICNL
          CALL GETCNF_LUCIA(SCR(KRCONF),IRTYP,IPCNF(ICNR),ICONF,IREFSM,
     &                      NEL)
          NCSFR = NCSFTP(IRTYP)
          CALL CNHCN(SCR(KLCONF),ILTYP,SCR(KRCONF),IRTYP,SCR(KLPHPS),
     &               SCR(KLFREE),NAEL,NBEL,ECORE,
     &               ONEBOD,IPRODT,DTOC,NACTOB,TUVX,NTEST,ExFac,IREOTS)
          DO IIL = 1, NCSFL
            IF(IILB.EQ.IIRB) THEN
              IIRMAX = IIL
            ELSE
              IIRMAX = NCSFR
            END IF
            DO IIR = 1, IIRMAX
              IIRACT = IIRB-1+IIR
              IILACT = IILB-1+IIL
              ILRI = (IIR-1)*NCSFL+IIL
              ILRO = ((IILACT*IILACT-IILACT)/2)+IIRACT
              PHP(ILRO) = SCR(KLPHPS-1+ILRI)
            END DO
          END DO
          IIRB = IIRB + NCSFR
        END DO
      IILB = IILB + NCSFL
      END DO
*
      RETURN
      END
