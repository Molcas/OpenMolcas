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
      SUBROUTINE PMATEL(ISTATE,JSTATE,PROP,PINT,SMAT,CNO,OCC,
     *                  SFOLD,AFOLD,TDAO)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PINT(NBTRI),SFOLD(NBTRI),AFOLD(NBTRI),CNO(NCMO)
      DIMENSION TDAO(NBAST,NBAST),SMAT(*),OCC(NBAST)
      DIMENSION IDUM(1)
      INTEGER   ISYMLB
      REAL*8    PROP(NRROOT,NRROOT,NPROP)

#include "SysDef.fh"
#include "mrci.fh"
      SAVE ICALL
      DATA ICALL /0/
      IF(ISTATE.EQ.JSTATE) THEN
C READ OVERLAP INTEGRALS FROM TRAONE.
        CALL RDONE(IRTC,6,'MLTPL  0',1,SMAT,IDUMMY)
C CALCULATE AND WRITE MULLIKEN CHARGES.
        WRITE(6,*)
      CALL XFLUSH(6)
        WRITE(6,'(A,I2)')' MULLIKEN CHARGES FOR STATE NR ',ISTATE
      CALL XFLUSH(6)
        CALL CHARGE(NSYM,NBAS,NAME,CNO,OCC,SMAT,2,.True.,.True.)
        WRITE(6,*)' ',('*',I=1,70)
      CALL XFLUSH(6)
      END IF
C FOLD TDAO SYMMETRICALLY (ANTI-SYMM) INTO SFOLD (AFOLD):
C MOLCAS2 UPDATE: SYMMETRY-BLOCKED STORAGE.
      CALL DCOPY_(NBTRI,[0.0D00],0,SFOLD,1)
      CALL DCOPY_(NBTRI,[0.0D00],0,AFOLD,1)
      IJ=0
      IEND=0
      DO ISY=1,NSYM
        ISTA=IEND+1
        IEND=IEND+NBAS(ISY)
        DO I=ISTA,IEND
          DO J=ISTA,I-1
            IJ=IJ+1
            SFOLD(IJ)=TDAO(I,J)+TDAO(J,I)
            AFOLD(IJ)=TDAO(I,J)-TDAO(J,I)
          END DO
          IJ=IJ+1
          SFOLD(IJ)=TDAO(I,I)
          AFOLD(IJ)=0.0D00
        END DO
      END DO
      NSIZ=0
      DO 100 IPROP=1,NPROP
C PICK UP MATRIX ELEMENTS FROM ONE-ELECTRON FILE:
        CALL iRDONE(IRTC,1,PNAME(IPROP),IPCOMP(IPROP),IDUM,ISYMLB)
        IF(IRTC.EQ.0) NSIZ=IDUM(1)
        CALL RDONE(IRTC,0,PNAME(IPROP),IPCOMP(IPROP),PINT,ISYMLB)
C SEPARATE OUT THE OPERATOR GAUGE ORIGIN, AND NUCLEAR CONTRIBUTION:
        IF(ICALL.EQ.0) THEN
          PORIG(1,IPROP)=PINT(NSIZ+1)
          PORIG(2,IPROP)=PINT(NSIZ+2)
          PORIG(3,IPROP)=PINT(NSIZ+3)
          PNUC(IPROP)=PINT(NSIZ+4)
        END IF
        IF(ISYMLB.NE.1) THEN
C NON-DIAGONAL SYMMETRY BLOCKS MUST BE COMPRESSED AWAY:
          IFROM=1
          ITO=1
          DO 40 ISY1=1,NSYM
            NB1=NBAS(ISY1)
            IF(NB1.EQ.0) GOTO 40
            DO 30 ISY2=1,ISY1
              NB2=NBAS(ISY2)
              IF(NB2.EQ.0) GOTO 30
              ISY12=MUL(ISY1,ISY2)
CPAM96              MASK=2**(ISY12-1)
CPAM96              IF(IAND(ISYMLB,MASK).EQ.0) GOTO 30
              IF(MOD(ISYMLB,2**(ISY12)).EQ.0) GOTO 30
              NB12=NB1*NB2
              IF(ISY12.EQ.1) THEN
                NB12=(NB12+NB1)/2
                IF(IFROM.GT.ITO)
     *            CALL DCOPY_(NB12,PINT(IFROM),1,PINT(ITO),1)
                ITO=ITO+NB12
              END IF
              IFROM=IFROM+NB12
30          CONTINUE
40        CONTINUE
          NSIZ=ITO
        END IF
C PUT DDOT OF TR DENS MATRIX AND INTEGRALS INTO PROPER MATRIX ELEMENT
C FOR MULTIPOLES, USE NEGATIVE SIGN OF ELECTRONIC PART.
        SGN=1.0D00
        IF(PNAME(IPROP)(1:5).EQ.'MLTPL') SGN=-SGN
        IF(PTYPE(IPROP).EQ.'HERM') THEN
          X=SGN*DDOT_(NBTRI,SFOLD,1,PINT,1)
          PROP(ISTATE,JSTATE,IPROP)=X
          PROP(JSTATE,ISTATE,IPROP)=X
        ELSE
          X=SGN*DDOT_(NBTRI,AFOLD,1,PINT,1)
          PROP(ISTATE,JSTATE,IPROP)=X
          PROP(JSTATE,ISTATE,IPROP)=-X
        END IF
100   CONTINUE
      ICALL=1
      RETURN
      END
