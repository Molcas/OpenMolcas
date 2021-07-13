!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1989, Jeppe Olsen                                      *
!               1989, Markus P. Fuelscher                              *
!               Giovanni Li Manni                                      *
!***********************************************************************
      SUBROUTINE ipcsfsplit(PHPCSF,IPCSF,IPCNF,MXPDIM,MXSPLI,           &
     &                  DTOC,IPRODT,ICONF,                              &
     &                  IREFSM,ONEBOD,ECORE,NACTOB,                     &
     &                  SCR,NCONF,NEL,NAEL,NBEL,                        &
     &                  DIAG,TUVX,NTEST,ExFac,IREOTS)

!************* Author : GLMJ ****************
!
!     Obtain primary subspace and obtain
!     explicit representation of hamilton matrix in subspace
!
!     ARGUMENTS :
!     ===========
!     IPCSF  : CSF's defining subspace (output)
!     IPCNF  : Configurations defining subspace (output )
!     MXPDIM : Largest allowed dimension of subspace (Input)
!     DTOC   : Transformation matrix between CSF's and DET's (input)
!     IPRODT : Prototype determinants (input)
!     ICONF  : List of configurations  (input)
!     IREFSM : symmetry of considered CI space (input)
!     Onebod : one body hamilton matrix in rectangular form (input)
!     ECORE  : Core energy (input)
!     NACTOB : Number of active orbitals (input)
!     SCR    : Scratch array of length ????
!     NCONF  : Number of configurations of symmetry IREFSM
!     NPCNF  : Number of primary configurations obtained (output)
!     NPCSF  : Number of primary CSF's obtained  (OUTPUT)
!     TUVX   : Two-electron integrals (MO space)
!     DIAG   : Hamilton diagonal over CSF's (INPUT)
!
!     IREOTS : Type => symmetry reordering array
!
!     Jeppe Olsen , Summer of '89
!     adapted to DETRAS by M.P. Fuelscher, Oktober 1989
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
#include "spinfo.fh"
#include "splitcas.fh"
!
      DIMENSION IPCSF(*),IPCNF(*),DIAG(*), PHPCSF(*)
      DIMENSION DTOC(*),IPRODT(*),ICONF(*)
      DIMENSION ONEBOD(*),SCR(*)
      DIMENSION TUVX(*), IREOTS(*)
!
      CALL IPCSFSPLIT_INTERNAL(SCR)
!
!     This is to allow type punning without an explicit interface
      CONTAINS
      SUBROUTINE IPCSFSPLIT_INTERNAL(SCR)
      USE ISO_C_BINDING
      REAL*8, TARGET :: SCR(*)
      INTEGER, POINTER :: iSCR(:)
!
!     Assumed machine accuracy (give and take)
!
      Acc=1.0D-13
! construct the diagonal of the Hamilton matrix in CNF basis
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
!      Call RecPrt('DIAG',' ',DIAG,1,IICSF-1)
!      Call RecPrt('SCR',' ',SCR,1,IICNF-1)


! find the elements of the subspace
! (MXPDIM elements smallest in energy)

      XMAX = FNDMNX(SCR,NCONF,2)
      NPCSF = 0
      NPCNF = 0
!      IFINIT = 0
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
!      IF( (NPCSF+NCSFMN).LE. MXPDIM ) THEN
!       add new configuration
        NPCNF = NPCNF + 1
        IPCNF(NPCNF) = IMIN
        CALL ISTVC2(IPCSF(NPCSF+1),ICSFMN-1,1,NCSFMN)
        if( (NPCSF+NCSFMN).LE.MXSPLI) then
          iDimBlockA = NPCSF+NCSFMN
          iDimBlockACNF = NPCNF
        end if
        NPCSF = NPCSF + NCSFMN
        SCR(IMIN) = XMAX + 1.0D0
!      ELSE
!       No space for this configuration , remove previous
!       configurations with the same diagonal value
!        IFINIT = 1
!        IICNF = NPCNF+1
!  600     CONTINUE
!        IICNF = IICNF - 1
!        DIAVAL = SCR(IPCNF(IICNF))
!        IF( ABS(DIAVAL-XMIN).LE.1.0D-10) THEN
!          NPCNF = NPCNF -1
!          CALL GETCNF_LUCIA(SCR(NCONF+1),ITYP,IPCNF(IICNF),
!     &                ICONF,IREFSM,NEL)
!          NPCSF = NPCSF - NCSFTP(ITYP)
!          GOTO 600
!        END IF
!      END IF
!      IF( (IFINIT.EQ.0) .AND. (NPCNF.LT.NCONF) ) GOTO 400
      IF( NPCNF.LT.NCONF ) GOTO 400
!
      IF( NTEST.GE.30 ) THEN
        WRITE(6,*) ' Output from ipCSFSplit '
        WRITE(6,*) ' ================== '
        WRITE(6,*)                                                      &
     &  ' Number of Configurations in primary subspace ',NPCNF
        WRITE(6,*)                                                      &
     &  ' Number of CSFs in primary subspace ',NPCSF
        WRITE(6,*) ' Configurations included : '
        CALL IWRTMA(IPCNF,1,NPCNF,1,NPCNF)
        WRITE(6,*) ' CSFs included : '
        CALL IWRTMA(IPCSF,1,NPCSF,1,NPCSF)
      END IF

! construct the diagonal array out of the Hamiltonian matrix

      MXCSFC = 0
      DO ITYP = 1, NTYP
        MXCSFC = MAX(MXCSFC,NCSFTP(ITYP))
      END DO
!
      KLFREE = 1
      KLCONF = KLFREE
      KLFREE = KLFREE+NEL
      KLFREE = KLFREE+NEL
      KLPHPS = KLFREE
      KLFREE = KLFREE+MXCSFC*MXCSFC
!      KLSCRS = KLFREE
!      KLFREE = KLFREE+MXCSFC
!
      IILB = 1
      DO ICNL = 1, NCONF
!       write(*,*) 'IILB', IILB
        CALL C_F_POINTER(C_LOC(SCR(KLCONF)),iSCR,[1])
        CALL GETCNF_LUCIA(iSCR,ILTYP,IPCNF(ICNL),ICONF,IREFSM,          &
     &                    NEL)
        NULLIFY(iSCR)
        NCSFL = NCSFTP(ILTYP)
!       write(*,*)'NCSFL = ', NCSFL
        CALL C_F_POINTER(C_LOC(SCR(KLCONF)),iSCR,[1])
        CALL CNHCN(iSCR,ILTYP,iSCR,ILTYP,SCR(KLPHPS),                   &
     &             SCR(KLFREE),NAEL,NBEL,ECORE,                         &
     &             ONEBOD,IPRODT,DTOC,NACTOB,TUVX,NTEST,ExFac,IREOTS)
        NULLIFY(iSCR)
        DO IIL = 1, NCSFL
          IILACT = IILB-1+IIL
          ILRI = IIL*IIL
          PHPCSF(IILACT) = SCR(KLPHPS-1+ILRI)
!          write(*,*) 'IILACT, ILRI = ',IILACT, ILRI
!          write(*,*)'PHPCSF(IILACT)', PHPCSF(IILACT)
!          SCR(KLSCRS-1+IIL)= SCR(KLPHPS-1+ILRI)
        END DO
!        XMAX = -FNDMNX(SCR(KLSCRS),NCSFL,2)
!        PHPCNF(ICNL)= XMAX
!        write(*,*)'PHPCNF(IILACT)', PHPCNF(ICNL)
        IILB = IILB + NCSFL
      END DO

      RETURN
! Avoid unused argument warnings
      If (.False.) Call Unused_integer(MXPDIM)
      END SUBROUTINE IPCSFSPLIT_INTERNAL
!
      END
