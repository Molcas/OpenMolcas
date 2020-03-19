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
      SUBROUTINE MKRPTORB(FIFA,TORB,CMO)
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
C Transform to orbitals that diagonalize the diagonal
C blocks of FIFA. Affected data sets are CMO,
C EPS, EPSI, EPSA, and EPSE. Also, the CI arrays are
C transformed on file LUONEM. Note: FIFA is unchanged
C and is not valid for the new orbitals. It will be
C recomputed later.
C The transformation matrices are returned in TORB.
      REAL*8, INTENT(IN) :: FIFA(NFIFA)
* PAM Feb 2015: NTORB is in Include/caspt2.fh
*      INTEGER, INTENT(IN) :: NTORB
* -------------------------------------------
      REAL*8, INTENT(OUT) :: TORB(NTORB)
      REAL*8, INTENT(INOUT) :: CMO(NCMO)

C     indices
      INTEGER I,II,IST,ISYM,ISTART
      INTEGER ITO,ITOSTA,ITOEND
      INTEGER ICMOSTA,ICMOEND
      INTEGER IDR,IDW
      INTEGER IEPS,IEPSI,IEPSA,IEPSE
      INTEGER IOSTA,IOEND
C     work-array pointers
      INTEGER LCI,LCMO2,LFOCK
      INTEGER NFOCK,NFES
#ifdef _ENABLE_BLOCK_DMRG_
      INTEGER LXMAT,NXMAT
#endif
C     #orbitals per symmetry
      INTEGER NF,NI,NA,NR1,NR2,NR3,NS,NO,NB
      INTEGER NSCT,NCMOSCT

      CALL QENTER('MKRPTORB')

* Allocate space for temporary square Fock matrix in each symmetry:
* NBMX=Max number of basis functions in any symmetry, in common in caspt2.fh
      NFOCK=NOMX**2
      CALL GETMEM('FOCK','ALLO','REAL',LFOCK,NFOCK)
* Allocate space for new CMO coefficients:
      CALL GETMEM('CMO2','ALLO','REAL',LCMO2,NCMO)

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
        NF=NFRO(ISYM)
        NI=NISH(ISYM)
        NA=NASH(ISYM)
        NR1=NRAS1(ISYM)
        NR2=NRAS2(ISYM)
        NR3=NRAS3(ISYM)
        NS=NSSH(ISYM)
        NO=NORB(ISYM)
        NB=NBAS(ISYM)
* Put Fock matrix in square format in WORK(LFOCK)
        IF(NO.GT.0) THEN
        CALL SQUARE(FIFA(NFES+1),WORK(LFOCK),NO,1,NO)
        END IF
* Number of orbitals processed so far in this symmetry:
        IOEND=0
* Frozen orbitals: Just copy frozen CMO coefficients.
        NSCT=NFRO(ISYM)
        IF(NSCT.GT.0) THEN
         NCMOSCT=NBAS(ISYM)*NSCT
         ICMOSTA=ICMOEND+1
         ICMOEND=ICMOEND+NCMOSCT
         CALL DCOPY_(NCMOSCT,CMO(ICMOSTA),1,WORK(LCMO2-1+ICMOSTA),1)
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
         CALL DIAFCK(NO,WORK(LFOCK),IOSTA,IOEND,TORB(ITOSTA),
     &               NB,CMO(ICMOSTA),WORK(LCMO2-1+ICMOSTA))
         DO I=1,NSCT
          II=IOSTA-1+I
          IEPS=IEPS+1
          EPS(IEPS)=WORK(LFOCK-1+II+NO*(II-1))
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
         CALL DIAFCK(NO,WORK(LFOCK),IOSTA,IOEND,TORB(ITOSTA),
     &               NB,CMO(ICMOSTA),WORK(LCMO2-1+ICMOSTA))
         DO I=1,NSCT
          II=IOSTA-1+I
          IEPS=IEPS+1
          EPS(IEPS)=WORK(LFOCK-1+II+NO*(II-1))
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
         CALL DIAFCK(NO,WORK(LFOCK),IOSTA,IOEND,TORB(ITOSTA),
     &               NB,CMO(ICMOSTA),WORK(LCMO2-1+ICMOSTA))
         DO I=1,NSCT
          II=IOSTA-1+I
          IEPS=IEPS+1
          EPS(IEPS)=WORK(LFOCK-1+II+NO*(II-1))
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
         CALL DIAFCK(NO,WORK(LFOCK),IOSTA,IOEND,TORB(ITOSTA),
     &               NB,CMO(ICMOSTA),WORK(LCMO2-1+ICMOSTA))
         DO I=1,NSCT
          II=IOSTA-1+I
          IEPS=IEPS+1
          EPS(IEPS)=WORK(LFOCK-1+II+NO*(II-1))
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
         CALL DIAFCK(NO,WORK(LFOCK),IOSTA,IOEND,TORB(ITOSTA),
     &               NB,CMO(ICMOSTA),WORK(LCMO2-1+ICMOSTA))
         DO I=1,NSCT
          II=IOSTA-1+I
          IEPS=IEPS+1
          EPS(IEPS)=WORK(LFOCK-1+II+NO*(II-1))
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
         CALL DCOPY_(NCMOSCT,CMO(ICMOSTA),1,WORK(LCMO2-1+ICMOSTA),1)
        END IF

        NFES=NFES+(NO*(NO+1))/2

      END DO

* Actually, we will not use the old CMO array any more, so just
* overwrite it with the new ones and get rid of the allocated array.
      CALL DCOPY_(NCMO,WORK(LCMO2),1,CMO,1)
      CALL GETMEM('CMO2','FREE','REAL',LCMO2,NCMO)
* We will not use the Fock matrix either. It was just used temporarily
* for each turn of the symmetry loop. Skip it.
      CALL GETMEM('FOCK','FREE','REAL',LFOCK,NFOCK)

C Finally, loop again over symmetries, transforming the CI:
      IF(ISCF.EQ.0) THEN
#if defined _ENABLE_BLOCK_DMRG_ || defined _ENABLE_CHEMPS2_DMRG_
        IF(.NOT.DoCumulant) THEN
#endif
          CALL GETMEM('LCI3','ALLO','REAL',LCI,NCONF)
          IDR=IDCIEX
          IDW=IDTCEX
          DO IST=1,NSTATE
           CALL DDAFILE(LUCIEX,2,WORK(LCI),NCONF,IDR)
           ITOEND=0
           DO ISYM=1,NSYM
            NI=NISH(ISYM)
            NA=NASH(ISYM)
            NR1=NRAS1(ISYM)
            NR2=NRAS2(ISYM)
            NR3=NRAS3(ISYM)
            NS=NSSH(ISYM)
            NO=NORB(ISYM)
            NB=NBAS(ISYM)
            ITOSTA=ITOEND+1
            ITOEND=ITOEND+NI**2+NR1**2+NR2**2+NR3**2+NS**2

            ITO=ITOSTA+NI**2
            IF(NA.GT.0) THEN
              IF(NR1.GT.0) THEN
                ISTART=NAES(ISYM)+1
                CALL TRACI_RPT2(ISTART,NR1,TORB(ITO),LSYM,
     &                                         NCONF,WORK(LCI))
              END IF
              ITO=ITO+NR1**2
              IF(NR2.GT.0) THEN
                ISTART=NAES(ISYM)+NR1+1
                CALL TRACI_RPT2(ISTART,NR2,TORB(ITO),LSYM,
     &                                         NCONF,WORK(LCI))
              END IF
              ITO=ITO+NR2**2
              IF(NR3.GT.0) THEN
                ISTART=NAES(ISYM)+NR1+NR2+1
                CALL TRACI_RPT2(ISTART,NR3,TORB(ITO),LSYM,
     &                                         NCONF,WORK(LCI))
              END IF
            END IF
           END DO
           CALL DDAFILE(LUCIEX,1,WORK(LCI),NCONF,IDW)
          END DO
          CALL GETMEM('LCI3','FREE','REAL',LCI,NCONF)
#ifdef _ENABLE_BLOCK_DMRG_
        ELSE
* Transforming 2,3-RDMs from Block DMRG (1-RDM is computed from 2-RDM)
* NN.14 : For the time, Block's dump files of RDMs are directly loaded,
*         but those should be stored in JobIph file eventually.
          NXMAT=NASHT**2
* Workspace for transformation matrix
          CALL GETMEM('XMAT','ALLO','REAL',LXMAT,NXMAT)
          CALL DCOPY_(NXMAT,[0.0D0],0,WORK(LXMAT),1)
          CALL MKXMAT(TORB,WORK(LXMAT))

          CALL block_tran2pdm(NASHT,WORK(LXMAT),JSTATE,JSTATE)
          CALL block_tran3pdm(NASHT,WORK(LXMAT),JSTATE,JSTATE)

          CALL GETMEM('XMAT','FREE','REAL',LXMAT,NXMAT)
        END IF
#elif _ENABLE_CHEMPS2_DMRG_
        ELSE
          write(6,*) 'CHEMPS2> MKRPTORB assumes '//
     & 'PSEUDOCANONICAL orbitals!'
        END IF
#endif
      END IF


      CALL QEXIT('MKRPTORB')

      RETURN
      END
