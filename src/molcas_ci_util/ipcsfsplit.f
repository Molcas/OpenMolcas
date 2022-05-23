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
*               Giovanni Li Manni                                      *
************************************************************************
      SUBROUTINE ipcsfsplit(PHPCSF,IPCSF,IPCNF,MXPDIM,MXSPLI,
     *                  DTOC,IPRODT,ICONF,
     *                  IREFSM,ONEBOD,ECORE,NACTOB,
     *                  SCR,NCONF,NEL,NAEL,NBEL,
     *                  DIAG,TUVX,NTEST,ExFac,IREOTS)

************** Author : GLMJ ****************
C
C     Obtain primary subspace and obtain
C     explicit representation of hamilton matrix in subspace
C
C     ARGUMENTS :
C     ===========
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
C     DIAG   : Hamilton diagonal over CSF's (INPUT)
*
*     IREOTS : Type => symmetry reordering array
C
C     Jeppe Olsen , Summer of '89
C     adapted to DETRAS by M.P. Fuelscher, Oktober 1989
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "spinfo.fh"
#include "splitcas.fh"
C
      DIMENSION IPCSF(*),IPCNF(*),DIAG(*), PHPCSF(*)
      DIMENSION DTOC(*),IPRODT(*),ICONF(*)
      DIMENSION ONEBOD(*),SCR(*)
      DIMENSION TUVX(*), IREOTS(*)
*
      CALL IPCSFSPLIT_INTERNAL(SCR)
*
*     This is to allow type punning without an explicit interface
      CONTAINS
      SUBROUTINE IPCSFSPLIT_INTERNAL(SCR)
      USE ISO_C_BINDING
      REAL*8, TARGET :: SCR(*)
      INTEGER, POINTER :: iSCR(:)
*
*     Assumed machine accuracy (give and take)
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
*      Call RecPrt('DIAG',' ',DIAG,1,IICSF-1)
*      Call RecPrt('SCR',' ',SCR,1,IICNF-1)


* find the elements of the subspace
* (MXPDIM elements smallest in energy)

      XMAX = FNDMNX(SCR,NCONF,2)
      NPCSF = 0
      NPCNF = 0
C      IFINIT = 0
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
C      IF( (NPCSF+NCSFMN).LE. MXPDIM ) THEN
*       add new configuration
        NPCNF = NPCNF + 1
        IPCNF(NPCNF) = IMIN
        CALL ISTVC2(IPCSF(NPCSF+1),ICSFMN-1,1,NCSFMN)
        if( (NPCSF+NCSFMN).LE.MXSPLI) then
          iDimBlockA = NPCSF+NCSFMN
          iDimBlockACNF = NPCNF
        end if
        NPCSF = NPCSF + NCSFMN
        SCR(IMIN) = XMAX + 1.0D0
C      ELSE
*       No space for this configuration , remove previous
*       configurations with the same diagonal value
C        IFINIT = 1
C        IICNF = NPCNF+1
C  600     CONTINUE
C        IICNF = IICNF - 1
C        DIAVAL = SCR(IPCNF(IICNF))
C        IF( ABS(DIAVAL-XMIN).LE.1.0D-10) THEN
C          NPCNF = NPCNF -1
C          CALL GETCNF_LUCIA(SCR(NCONF+1),ITYP,IPCNF(IICNF),
C     &                ICONF,IREFSM,NEL)
C          NPCSF = NPCSF - NCSFTP(ITYP)
C          GOTO 600
C        END IF
C      END IF
C      IF( (IFINIT.EQ.0) .AND. (NPCNF.LT.NCONF) ) GOTO 400
      IF( NPCNF.LT.NCONF ) GOTO 400
*
      IF( NTEST.GE.30 ) THEN
        WRITE(6,*) ' Output from ipCSFSplit '
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

* construct the diagonal array out of the Hamiltonian matrix

      MXCSFC = 0
      DO ITYP = 1, NTYP
        MXCSFC = MAX(MXCSFC,NCSFTP(ITYP))
      END DO
*
      KLFREE = 1
      KLCONF = KLFREE
      KLFREE = KLFREE+NEL
      KLFREE = KLFREE+NEL
      KLPHPS = KLFREE
      KLFREE = KLFREE+MXCSFC*MXCSFC
C      KLSCRS = KLFREE
C      KLFREE = KLFREE+MXCSFC
*
      IILB = 1
      DO ICNL = 1, NCONF
*       write(*,*) 'IILB', IILB
        CALL C_F_POINTER(C_LOC(SCR(KLCONF)),iSCR,[1])
        CALL GETCNF_LUCIA(iSCR,ILTYP,IPCNF(ICNL),ICONF,IREFSM,
     &                    NEL)
        NULLIFY(iSCR)
        NCSFL = NCSFTP(ILTYP)
*       write(*,*)'NCSFL = ', NCSFL
        CALL C_F_POINTER(C_LOC(SCR(KLCONF)),iSCR,[1])
        CALL CNHCN(iSCR,ILTYP,iSCR,ILTYP,SCR(KLPHPS),
     &             SCR(KLFREE),NAEL,NBEL,ECORE,
     &             ONEBOD,IPRODT,DTOC,NACTOB,TUVX,NTEST,ExFac,IREOTS)
        NULLIFY(iSCR)
        DO IIL = 1, NCSFL
          IILACT = IILB-1+IIL
          ILRI = IIL*IIL
          PHPCSF(IILACT) = SCR(KLPHPS-1+ILRI)
*          write(*,*) 'IILACT, ILRI = ',IILACT, ILRI
*          write(*,*)'PHPCSF(IILACT)', PHPCSF(IILACT)
C          SCR(KLSCRS-1+IIL)= SCR(KLPHPS-1+ILRI)
        END DO
C        XMAX = -FNDMNX(SCR(KLSCRS),NCSFL,2)
C        PHPCNF(ICNL)= XMAX
C        write(*,*)'PHPCNF(IILACT)', PHPCNF(ICNL)
        IILB = IILB + NCSFL
      END DO

      RETURN
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(MXPDIM)
      END SUBROUTINE IPCSFSPLIT_INTERNAL
*
      END
