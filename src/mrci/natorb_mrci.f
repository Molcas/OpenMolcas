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
      SUBROUTINE NATORB_MRCI(CMO,DMO,CNO,OCC,SCR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CMO(NCMO),DMO(NBTRI),CNO(NCMO),OCC(NBAST)
      DIMENSION SCR( (NBMAX*(NBMAX+1))/2 )

#include "SysDef.fh"

#include "mrci.fh"
      CALL DCOPY_(NBAST,0.0D00,0,OCC,1)
      CALL DCOPY_(NCMO,CMO,1,CNO,1)
C LOOP OVER SYMMETRY LABELS
C Present index of end of processed CMO block:
      IECMO=0
C Present end of index to MO-s, to be translated by ICH array:
      IEO=0
C Present end of index to basis functions:
      IEB=0
      DO 100 ISYM=1,NSYM
        NO=NORB(ISYM)
        NIAV=NISH(ISYM)+NASH(ISYM)+NVIR(ISYM)
        NB=NBAS(ISYM)
        If (NB.eq.0) Go To 100
        ISB=IEB+1
        IEB=IEB+NB
C ORBITALS PRE-FROZEN IN MOTRA, OR FROZEN IN MRCI:
        NF=NFMO(ISYM)+NFRO(ISYM)
        NBF=NB*NF
        ISCMO=IECMO+1
        IECMO=IECMO+NBF
C (DO NOTHING WITH THE FROZEN ORBITALS)
        IF(NF.GT.0) CALL DCOPY_(NF,2.0D00,0,OCC(ISB),1)
        IEO=IEO+NFRO(ISYM)
C ORBITALS EXPLICITLY USED IN CI:
        NBO=NB*NIAV
        ISCMO=IECMO+1
        IECMO=IECMO+NBO
        ISO=IEO+1
        IEO=IEO+NIAV
C TRANSFER SYMMETRY BLOCK OF DMO TO TRIANGULAR SCRATCH MATRIX:
        I12=0
        DO 10 I=ISO,IEO
          IO1=ICH(I)
          DO 10 J=ISO,I
            IO2=ICH(J)
            IO12=(IO1*(IO1-1))/2+IO2
            IF(IO1.LT.IO2) IO12=(IO2*(IO2-1))/2+IO1
            I12=I12+1
            SCR(I12)=DMO(IO12)
10      CONTINUE
C DIAGONALIZE AND TRANSFORM ORBITALS:
        CALL JACOB(SCR,CNO(ISCMO),NIAV,NB)
C PICK OCCUP NR FROM DIAGONAL:
        II=0
        DO 20 I=1,NIAV
          II=II+I
          OCC(ISB+NF-1+I)=SCR(II)
20      CONTINUE
C ORDER BY DECREASING NATURAL OCCUPANCY:
        NN=NO-NFRO(ISYM)
        DO 40 I=1,NN-1
          OMAX=OCC(ISB+NF-1+I)
          IMAX=I
          DO 30 J=I+1,NN
            OC=OCC(ISB+NF-1+J)
            IF(OMAX.GE.OC) GOTO 30
            IMAX=J
            OMAX=OC
30        CONTINUE
          IF(IMAX.EQ.I) GOTO 40
          OCC(ISB+NF-1+IMAX)=OCC(ISB+NF-1+I)
          OCC(ISB+NF-1+I)=OMAX
          ISTA1=ISCMO+NB*(I-1)
          ISTA2=ISCMO+NB*(IMAX-1)
          CALL DCOPY_(NB,CNO(ISTA1),1,SCR,1)
          CALL DCOPY_(NB,CNO(ISTA2),1,CNO(ISTA1),1)
          CALL DCOPY_(NB,SCR,1,CNO(ISTA2),1)
40      CONTINUE
C ORBITALS PRE-DELETED IN MOTRA OR DELETED IN MRCI:
        ND=NDMO(ISYM)+NDEL(ISYM)
        NBD=NB*ND
        ISCMO=IECMO+1
        IECMO=IECMO+NBD
        IEO=IEO+NDEL(ISYM)
C (DO NOTHING WITH THE DELETED ORBITALS)
100   CONTINUE
      RETURN
      END
