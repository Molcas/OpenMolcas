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
      SUBROUTINE DiagOrd(PHPCSF,PHPCNF,IPORDCSF,IPORDCNF,MXPDIM,        &
     &                  condition,iter,                                 &
     &                  DTOC,IPRODT,ICONF,                              &
     &                  IREFSM,ONEBOD,ECORE,NACTOB,                     &
     &                  SCR,NCONF,NEL,NAEL,NBEL,                        &
     &                  TUVX,NTEST,ExFac,IREOTS)

!     Obtain primary subspace and obtain
!     explicit representation of hamilton matrix in subspace
!
!     ARGUMENTS :
!     ===========
!     PHPCSF     : Diagonal hamilton matrix elements in CSF (output)
!     PHPCNF     : Diagonal hamilton matrix elements in CNF (output)
!     IPORDCSF: index array containing energetic order in CSF (output)
!     IPORDCNF: index array containing energetic order in CNF (output)
!     IPCSF   : CSF's defining subspace (input)
!     IPCNF   : Configurations defining subspace (input)
!     MXPDIM  : Largest allowed dimension of subspace (Input)
!     DTOC    : Transformation matrix between CSF's and DET's (input)
!     IPRODT  : Prototype determinants (input)
!     ICONF   : List of configurations  (input)
!     IREFSM  : symmetry of considered CI space (input)
!     Onebod  : one body hamilton matrix in rectangular form (input)
!     ECORE   : Core energy (input)
!     NACTOB  : Number of active orbitals (input)
!     SCR     : Scratch array of length ????
!     NCONF   : Number of configurations of symmetry IREFSM
!     NPCNF   : = NCONF
!     NPCSF   : = MXPDIM
!     TUVX    : Two-electron integrals (MO space)
!
!     IREOTS : Type => symmetry reordering array
!

!*********** Author: GLMJ ****************
!   history:
!     Jeppe Olsen , Summer of '89
!     adapted to DETRAS by M.P. Fuelscher, Oktober 1989
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
#include "spinfo.fh"
#include "splitcas.fh"
!
      DIMENSION PHPCSF(*),PHPCNF(*)
!      DIMENSION IPCSF(*),IPCNF(*)
      DIMENSION IPORDCSF(*),IPORDCNF(*)
      DIMENSION DTOC(*),IPRODT(*),ICONF(*)
      DIMENSION ONEBOD(*),SCR(*)
      DIMENSION TUVX(*), IREOTS(*)

      CALL DIAGORD_INTERNAL(SCR)
!
!     This is to allow type punning without an explicit interface
      CONTAINS
      SUBROUTINE DIAGORD_INTERNAL(SCR)
      USE ISO_C_BINDING
      REAL*8, TARGET :: SCR(*)
      INTEGER, POINTER :: iSCR(:)
      ICSFMN = 0 ! dummy initialize to avoid warning

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
      KLSCRS = KLFREE
      KLFREE = KLFREE+MXCSFC
      KLCSFO = KLFREE
      KLFREE = KLFREE+MXPDIM
      KLCNFO = KLFREE
      KLFREE = KLFREE+NCONF
!
      IILB = 1
      DO ICNL = 1, NCONF
!       write(*,*) 'IILB', IILB
        CALL C_F_POINTER(C_LOC(SCR(KLCONF)),iSCR,[1])
        CALL GETCNF_LUCIA(iSCR,ILTYP,ICNL,ICONF,IREFSM,                 &
     &                    NEL)
        NULLIFY(iSCR)
!        CALL GETCNF_LUCIA(SCR(KLCONF),ILTYP,IPCNF(ICNL),ICONF,IREFSM,
!     &                    NEL)
        NCSFL = NCSFTP(ILTYP)
!       write(*,*)'NCSFL = ', NCSFL
!       call xflush(6)
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
!          call xflush(6)
          SCR(KLSCRS-1+IIL)= SCR(KLPHPS-1+ILRI)
        END DO
        XMAX = -FNDMNX(SCR(KLSCRS),NCSFL,2)
        PHPCNF(ICNL)= XMAX
!        write(*,*)'PHPCNF(IILACT)', PHPCNF(ICNL)
!          call xflush(6)
        IILB = IILB + NCSFL
      END DO
!*********************** ORDERING ***************************
      Acc = 1.0d-13
      XMAX = FNDMNX(PHPCNF,NCONF,2)
      RefSplit = -XMAX
      NPCSF = 0
      NPCNF = 0
! To avoid warnings:
        iTmpDimBlockACNF = 0
        iTmpDimBlockA    = 0

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
          IF( PHPCNF(IICNF)+Acc.LT.XMIN ) THEN
            XMIN = PHPCNF(IICNF)
            IMIN = IICNF
!            write(*,*)'IMIN         :', IMIN
            ICSFMN = IICSF
            NCSFMN = NIRREP
          END IF
          IICNF = IICNF + 1
          IICSF = IICSF + NIRREP
        END DO
      END DO
!      write(*,*)'XMIN         :', XMIN
!      write(*,*)'NPCSF+NCSFMN :', NPCSF+NCSFMN
!      IF( (NPCSF+NCSFMN).LE. MXPDIM ) THEN
        SCR(KLCNFO+NPCNF) = PHPCNF(IMIN)
        NPCNF = NPCNF + 1
        IPORDCNF(NPCNF) = IMIN
!        write(*,*)'NPCNF, IPORDCNF(NPCNF) :', NPCNF,IPORDCNF(NPCNF)
        CALL ISTVC2(IPORDCSF(NPCSF+1),ICSFMN-1,1,NCSFMN)
!        write(*,*) 'IPORDCSF(NPCSF) :',(IPORDCSF(NPCSF+i),i=1,NCSFMN)
!        call IVCPRT('IPORDCSF(NPCSF) :',' ', IPORDCSF(NPCSF+1),NCSFMN)
        NPCSF = NPCSF + NCSFMN

        if (iter.eq.1) then
         if (EnerSplit) then
           if ((PHPCNF(IMIN)-RefSplit).le.condition) then
             iDimBlockACNF = NPCNF
             iDimBlockA = NPCSF
           end if
         else if ( PerSplit ) then
           if ((REAL(NPCSF)/REAL(MXPDIM))*1.0d2.le.condition+Acc) then
             iDimBlockACNF = NPCNF
             iDimBlockA = NPCSF
           end if
         end if
        else
           if (NPCNF.le.iDimBlockACNF) then
             iTmpDimBlockACNF = NPCNF
             iTmpDimBlockA = NPCSF
           end if
        end if

        PHPCNF(IMIN) = XMAX + 1.0D0
        IF( (NPCNF.LT.NCONF) ) GOTO 400

        if (iter.ne.1) then
          iDimBlockACNF = iTmpDimBlockACNF
          iDimBlockA = iTmpDimBlockA
!          write(6,*) 'iDimBlockACNF, iDimBlockACSF',
!     &              iDimBlockACNF,iDimBlockA
        end if
!      ELSE
!       No space for this configuration , remove previous
!       configurations with the same diagonal value
!        IFINIT = 1
!        IICNF = NPCNF+1
!  600     CONTINUE
!        IICNF = IICNF - 1
!        DIAVAL = PHPCNF(IPORDCNF(IICNF))
!        IF( ABS(DIAVAL-XMIN).LE.1.0D-10) THEN
!          NPCNF = NPCNF -1
!          CALL GETCNF_LUCIA(SCR(KLFREE),ITYP,IPCNF(IICNF),
!     &                ICONF,IREFSM,NEL)
!          CALL GETCNF_LUCIA(PHPCNF(NCONF+1),ITYP,IPCNF(IICNF),
!     &                ICONF,IREFSM,NEL)
!          NPCSF = NPCSF - NCSFTP(ITYP)
!          GOTO 600
!        END IF
!      END IF
!      IF( (IFINIT.EQ.0) .AND. (NPCNF.LT.NCONF) ) GOTO 400

      do i=0,MXPDIM-1
        SCR(KLCSFO+i) = PHPCSF(IPORDCSF(i+1))
      end do

      call dcopy_(NCONF, SCR(KLCNFO),1,PHPCNF,1)
      call dcopy_(MXPDIM,SCR(KLCSFO),1,PHPCSF,1)

!        WRITE(6,*) ' OUTPUT from DiagOrd '
!        WRITE(6,*) ' ================== '
!        WRITE(6,*)
!     &  ' Number of Configurations in primary subspace ',NCONF
!        WRITE(6,*)
!     &  ' Number of CSFs in primary subspace ',MXPDIM
!        WRITE(6,*) ' Configurations included : '
!        CALL IWRTMA(IPORDCNF,1,NCONF,1,NCONF)
!        WRITE(6,*) ' CSFs included : '
!        CALL IWRTMA(IPORDCSF,1,MXPDIM,1,MXPDIM)
!        write(6,*) 'Ordered Diagonal array in CSF basis :'
!        do i = 1, MXPDIM
!          write(6,*) PHPCSF(IPORDCSF(i))
!          write(6,*) PHPCSF(i)
!        end do
!        write(6,*) 'Ordered Diagonal array in CNF basis :'
!        do i = 1, NCONF
!          write(6,*) PHPCNF(i)
!        end do

      RETURN
      END SUBROUTINE DIAGORD_INTERNAL
!
      END
