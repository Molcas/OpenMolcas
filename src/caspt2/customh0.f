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
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
************************************************************************
* This file contains boiler-plate code to allow developers to experiment
* with modifications to the zero-order hamiltonian. To enable them,
* remove the #if 0/#endif blocks, and use 'HZero = Custom' in the input.
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE NEWB
#if 0
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"

      INTEGER ICASE,ISYM,NAS,NIS,NCOEF
      INTEGER LS,IDS,NS,LB,IDB,NB

C Modify B matrices, if requested.

      CALL QENTER('NEWB')

      DO ICASE=1,11
        DO ISYM=1,NSYM
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          NCOEF=NAS*NIS
          IF(NCOEF.EQ.0) CYCLE
          NS=(NAS*(NAS+1))/2
          CALL GETMEM('SMAT','ALLO','REAL',LS,NS)
          IDS=IDSMAT(ISYM,ICASE)
          CALL DDAFILE(LUSBT,2,WORK(LS),NS,IDS)
          NB=NS
          CALL GETMEM('BMAT','ALLO','REAL',LB,NB)
          IDB=IDBMAT(ISYM,ICASE)
          CALL DDAFILE(LUSBT,2,WORK(LB),NB,IDB)
C Modify B matrix, using S matrix and some other data.
          CALL GETMEM('BMAT','FREE','REAL',LB,NB)
          CALL GETMEM('SMAT','FREE','REAL',LS,NS)
        END DO
      END DO

      CALL QEXIT('NEWB')

      RETURN
#endif
      END

      SUBROUTINE NEWDIA
#if 0
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "SysDef.fh"

      INTEGER ICASE,ISYM,NIN,NAS,NIS,I
      INTEGER LBD,LID,LC1,LC2,ID

C Post-diagonalization modification of diagonal energy
C denominator terms for active and for non-active superindex.

      CALL QENTER('NEWDIA')

      DO ICASE=1,13
        DO ISYM=1,NSYM
          NIN=NINDEP(ISYM,ICASE)
          IF(NIN.EQ.0) CYCLE
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF(NIS.EQ.0) CYCLE
C Remember: NIN values in BDIAG, but must read NAS for correct
C positioning.
          CALL GETMEM('LBD','ALLO','REAL',LBD,NAS)
          CALL GETMEM('LID','ALLO','REAL',LID,NIS)
          CALL GETMEM('LC1','ALLO','REAL',LC1,NAS)
          CALL GETMEM('LC2','ALLO','REAL',LC2,NIS)
          ID=IDBMAT(ISYM,ICASE)
C Active, and non-active, energy denominators:
          CALL DDAFILE(LUSBT,2,WORK(LBD),NAS,ID)
          CALL DDAFILE(LUSBT,2,WORK(LID),NIS,ID)
C Active, and non-active, corrections:
C (Replace this strange example with something sensible)
          CALL DCOPY_(NAS,[0.0d0],0,WORK(LC1),1)
          CALL DCOPY_(NIS,[0.0d0],0,WORK(LC2),1)
C Modifications are added to the usual diagonal energies:
          DO I=1,NAS
            WORK(LBD-1+I)=WORK(LBD-1+I)+WORK(LC1-1+I)
          END DO
          DO I=1,NIS
            WORK(LID-1+I)=WORK(LID-1+I)+WORK(LC2-1+I)
          END DO
          ID=IDBMAT(ISYM,ICASE)
          CALL DDAFILE(LUSBT,1,WORK(LBD),NAS,ID)
          CALL DDAFILE(LUSBT,1,WORK(LID),NIS,ID)
C Added modifications are saved on LUSBT.
          CALL DDAFILE(LUSBT,1,WORK(LC1),NAS,ID)
          CALL DDAFILE(LUSBT,1,WORK(LC2),NIS,ID)

          CALL GETMEM('LBD','FREE','REAL',LBD,NAS)
          CALL GETMEM('LID','FREE','REAL',LID,NIS)
          CALL GETMEM('LC1','FREE','REAL',LC1,NAS)
          CALL GETMEM('LC2','FREE','REAL',LC2,NIS)
          ID=IDBMAT(ISYM,ICASE)
        END DO
      END DO

      CALL QEXIT('NEWDIA')

      RETURN
#endif
      END
