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
* Copyright (C) 2006, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 2006  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE MKRPTORB(FIFA,NFIFA,TORB,NTORB,CMO,NCMO)
      use fciqmc_interface, only: DoFCIQMC, NonDiagonal
      use caspt2_global, only: LUCIEX, IDCIEX, IDTCEX
      use stdalloc, only: mma_allocate, mma_deallocate
#if defined(_DMRG_)
      use qcmaquis_interface, only: c_bool, c_int,
     &                              qcmaquis_interface_rotate_rdms
      use caspt2_module, only: DMRG
#endif
#if defined (_ENABLE_BLOCK_DMRG_) || defined (_DMRG_)
      use caspt2_module, only: jState, nAshT
#endif
      use caspt2_module, only: iSCF, nConf, nOMx,
     &                         nState, nSym, STSym, nIsh, nAsh, nRas1,
     &                         nRas2, nRas3, nSsh, nOrb, nBas, nFro,
     &                         EPS, EPSI, EPSA, nDel, EPSE
#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_)
      use caspt2_module, only: DoCumulant
#endif
      use constants, only: zero
      use definitions, only: iwp, wp
      IMPLICIT NONE
C Transform to orbitals that diagonalize the diagonal blocks of FIFA.
C Affected data sets are CMO, C EPS, EPSI, EPSA, and EPSE. Also, the
C CI arrays are transformed on file LUCIEX.
C Note: FIFA is unchanged  and is not valid for the new orbitals.
C It will be recomputed later.
C The transformation matrices are returned in TORB.
      INTEGER(kind=iwp), INTENT(IN) :: NFIFA,NTORB,NCMO
      REAL(kind=wp), INTENT(IN) :: FIFA(NFIFA)
* -------------------------------------------
      REAL(kind=wp), INTENT(OUT) :: TORB(NTORB)
      REAL(kind=wp), INTENT(INOUT) :: CMO(NCMO)

C     indices
      INTEGER(kind=iwp) I,II,IST,ISYM
      INTEGER(kind=iwp) ITOSTA,ITOEND
      INTEGER(kind=iwp) ICMOSTA,ICMOEND
      INTEGER(kind=iwp) IDR,IDW
      INTEGER(kind=iwp) IEPS,IEPSI,IEPSA,IEPSE
      INTEGER(kind=iwp) IOSTA,IOEND
      INTEGER(kind=iwp) NFOCK,NFES
#if defined(_ENABLE_BLOCK_DMRG_) || defined(_DMRG_)
      INTEGER(kind=iwp) NXMAT
      REAL(kind=wp), ALLOCATABLE:: XMAT(:)
#endif
C     #orbitals per symmetry
      INTEGER(kind=iwp) NI,NA,NR1,NR2,NR3,NS,NO,NB
      INTEGER(kind=iwp) NSCT,NCMOSCT
C     work-arrays
      REAL(kind=wp), ALLOCATABLE:: FOCK(:), CMO2(:), CI(:)


* Allocate space for temporary square Fock matrix in each symmetry:
* NBMX=Max number of basis functions in any symmetry, in common in caspt2_module.F90
      NFOCK=NOMX**2
      CALL mma_allocate(FOCK,NFOCK,LABEL='FOCK')
* Allocate space for new CMO coefficients:
      CALL mma_allocate(CMO2,NCMO,LABEL='CMO2')

* In the loop over symmetries, NFES is the nr of Fock matrix
* elements processed in earlier symmetries.
      NFES=0
      IEPS=0
      IEPSI=0
      IEPSA=0
      IEPSE=0
* The transformation matrices for each symmetry
* will be collected into TORB and returned. This is
* necessary for later backtransformation to original
* MO basis.
* ITOSTA,ITOEND: Section of TORB for each subspace.
      ITOEND=0
* ICMOSTA,ICMOEND: Section of CMO for each subspace.
      ICMOEND=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NA=NASH(ISYM)
        NR1=NRAS1(ISYM)
        NR2=NRAS2(ISYM)
        NR3=NRAS3(ISYM)
        NS=NSSH(ISYM)
        NO=NORB(ISYM)
        NB=NBAS(ISYM)
* Put Fock matrix in square format in FOCK
        IF(NO.GT.0) THEN
        CALL SQUARE(FIFA(NFES+1),FOCK,NO,1,NO)
        END IF
        ! the zero-ing out of the off-diagonal blocks happens here
        ! diafck is not structured like rasscf/fckpt2.f
* Number of orbitals processed so far in this symmetry:
        IOEND=0
* Frozen orbitals: Just copy frozen CMO coefficients.
        NSCT=NFRO(ISYM)
        IF(NSCT.GT.0) THEN
         NCMOSCT=NBAS(ISYM)*NSCT
         ICMOSTA=ICMOEND+1
         ICMOEND=ICMOEND+NCMOSCT
         CALL DCOPY_(NCMOSCT,CMO(ICMOSTA),1,CMO2(ICMOSTA),1)
        END IF
* Inactive block: Section length NSCT=NISH(ISYM)
        NSCT=NISH(ISYM)
        IF(NSCT.GT.0) THEN
         IOSTA=IOEND+1
         IOEND=IOEND+NSCT
         NCMOSCT=NBAS(ISYM)*NSCT
         ICMOSTA=ICMOEND+1
         ICMOEND=ICMOEND+NCMOSCT
         ITOSTA=ITOEND+1
         ITOEND=ITOEND+NSCT**2
         CALL DIAFCK(NO,FOCK,IOSTA,IOEND,TORB(ITOSTA),
     &               NB,CMO(ICMOSTA),CMO2(ICMOSTA))
         DO I=1,NSCT
          II=IOSTA-1+I
          IEPS=IEPS+1
          EPS(IEPS)=FOCK(II+NO*(II-1))
          IEPSI=IEPSI+1
          EPSI(IEPSI)=EPS(IEPS)
         END DO
        END IF
* RAS1 block: Section length NSCT=NRAS1(ISYM)
        NSCT=NRAS1(ISYM)
        IF(NSCT.GT.0) THEN
         IOSTA=IOEND+1
         IOEND=IOEND+NSCT
         NCMOSCT=NBAS(ISYM)*NSCT
         ICMOSTA=ICMOEND+1
         ICMOEND=ICMOEND+NCMOSCT
         ITOSTA=ITOEND+1
         ITOEND=ITOEND+NSCT**2
         CALL DIAFCK(NO,FOCK,IOSTA,IOEND,TORB(ITOSTA),
     &               NB,CMO(ICMOSTA),CMO2(ICMOSTA))
         DO I=1,NSCT
          II=IOSTA-1+I
          IEPS=IEPS+1
          EPS(IEPS)=FOCK(II+NO*(II-1))
          IEPSA=IEPSA+1
          EPSA(IEPSA)=EPS(IEPS)
         END DO
        END IF
* RAS2 block: Section length NSCT=NRAS2(ISYM)
        NSCT=NRAS2(ISYM)
        IF(NSCT.GT.0) THEN
         IOSTA=IOEND+1
         IOEND=IOEND+NSCT
         NCMOSCT=NBAS(ISYM)*NSCT
         ICMOSTA=ICMOEND+1
         ICMOEND=ICMOEND+NCMOSCT
         ITOSTA=ITOEND+1
         ITOEND=ITOEND+NSCT**2
         CALL DIAFCK(NO,FOCK,IOSTA,IOEND,TORB(ITOSTA),
     &               NB,CMO(ICMOSTA),CMO2(ICMOSTA))
         DO I=1,NSCT
          II=IOSTA-1+I
          IEPS=IEPS+1
          EPS(IEPS)=FOCK(II+NO*(II-1))
          IEPSA=IEPSA+1
          EPSA(IEPSA)=EPS(IEPS)
         END DO
        END IF
* RAS3 block: Section length NSCT=NRAS3(ISYM)
        NSCT=NRAS3(ISYM)
        IF(NSCT.GT.0) THEN
         IOSTA=IOEND+1
         IOEND=IOEND+NSCT
         NCMOSCT=NBAS(ISYM)*NSCT
         ICMOSTA=ICMOEND+1
         ICMOEND=ICMOEND+NCMOSCT
         ITOSTA=ITOEND+1
         ITOEND=ITOEND+NSCT**2
         CALL DIAFCK(NO,FOCK,IOSTA,IOEND,TORB(ITOSTA),
     &               NB,CMO(ICMOSTA),CMO2(ICMOSTA))
         DO I=1,NSCT
          II=IOSTA-1+I
          IEPS=IEPS+1
          EPS(IEPS)=FOCK(II+NO*(II-1))
          IEPSA=IEPSA+1
          EPSA(IEPSA)=EPS(IEPS)
         END DO
        END IF
* Secondary (virtual) block: Section length NSCT=NSSH(ISYM)
        NSCT=NSSH(ISYM)
        IF(NSCT.GT.0) THEN
         IOSTA=IOEND+1
         IOEND=IOEND+NSCT
         NCMOSCT=NBAS(ISYM)*NSCT
         ICMOSTA=ICMOEND+1
         ICMOEND=ICMOEND+NCMOSCT
         ITOSTA=ITOEND+1
         ITOEND=ITOEND+NSCT**2
         CALL DIAFCK(NO,FOCK,IOSTA,IOEND,TORB(ITOSTA),
     &               NB,CMO(ICMOSTA),CMO2(ICMOSTA))
         DO I=1,NSCT
          II=IOSTA-1+I
          IEPS=IEPS+1
          EPS(IEPS)=FOCK(II+NO*(II-1))
          IEPSE=IEPSE+1
          EPSE(IEPSE)=EPS(IEPS)
         END DO
        END IF
* Deleted orbitals: Just copy deleted CMO coefficients.
        NSCT=NDEL(ISYM)
        IF(NSCT.GT.0) THEN
         NCMOSCT=NBAS(ISYM)*NSCT
         ICMOSTA=ICMOEND+1
         ICMOEND=ICMOEND+NCMOSCT
         CALL DCOPY_(NCMOSCT,CMO(ICMOSTA),1,CMO2(ICMOSTA),1)
        END IF

        NFES=NFES+(NO*(NO+1))/2

      END DO

* Actually, we will not use the old CMO array any more, so just
* overwrite it with the new ones and get rid of the allocated array.
      CMO(:)=CMO2(:)
      CALL mma_deallocate(CMO2)

* We will not use the Fock matrix either. It was just used temporarily
* for each turn of the symmetry loop. Skip it.
      CALL mma_deallocate(FOCK)

C Finally, loop again over symmetries, transforming the CI:

      IF(ISCF.EQ.0) THEN
#ifdef _DMRG_
        if (DMRG) then
          NXMAT=NASHT**2
          CALL mma_allocate(XMAT,NXMAT,LABEL='XMAT')
          XMAT(:)=Zero
          CALL MKXMAT(TORB,XMAT)

          CALL qcmaquis_interface_rotate_rdms(int(JSTATE-1, c_int),
     &      int(JSTATE-1, c_int), int(0, c_int), XMAT,
     &      logical(.false., c_bool))
          do I=1,NSTATE
            if (JSTATE .ne. I) then
              write(6,*) "QCMaquis> Rotating tRDMs", JSTATE-1, I-1
            CALL qcmaquis_interface_rotate_rdms(int(JSTATE-1, c_int),
     &          int(I-1, c_int), int(0, c_int), XMAT,
     &          logical(.false., c_bool))
            end if
          end do
          CALL mma_deallocate(XMAT)
          end if
#endif

        if (DoFCIQMC) then
          if (NonDiagonal) then
           write(6,*)'Transforming CASPT2 intermediates to '//
     &               'pseudo-canonical orbitals.'
          else
            write(6,*)'FCIQMC-CASPT2 assumes pseudo-canonical orbitals.'
          end if
        else
#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_)
          IF(.NOT.DoCumulant) THEN
#endif
            CALL mma_allocate(CI,NCONF,Label='CI')
            IDR=IDCIEX(1)
            IDW=IDTCEX(1)
            DO IST=1,NSTATE
             CALL DDAFILE(LUCIEX,2,CI,NCONF,IDR)

             Call mkTraCI(nTORB,TORB,STSYM,nConf,CI)

             CALL DDAFILE(LUCIEX,1,CI,NCONF,IDW)
            END DO
            CALL mma_deallocate(CI)
#ifdef _ENABLE_BLOCK_DMRG_
          ELSE
* Transforming 2,3-RDMs from Block DMRG (1-RDM is computed from 2-RDM)
* NN.14 : For the time, Block's dump files of RDMs are directly loaded,
*         but those should be stored in JobIph file eventually.
          NXMAT=NASHT**2
* Workspace for transformation matrix
            CALL mma_allocate(XMAT,NXMAT,LABEL='XMAT')
            XMAT(:)=Zero
            CALL MKXMAT(TORB,XMAT)

            CALL block_tran2pdm(NASHT,XMAT,JSTATE,JSTATE)
            CALL block_tran3pdm(NASHT,XMAT,JSTATE,JSTATE)

            CALL mma_deallocate(XMAT)
          END IF
#elif _ENABLE_CHEMPS2_DMRG_
          ELSE
            write(6,*) 'CHEMPS2> MKRPTORB assumes '//
     &    'PSEUDOCANONICAL orbitals!'
          END IF
#endif
        end if
      END IF

      END SUBROUTINE MKRPTORB
