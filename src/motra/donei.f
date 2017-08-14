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
      SUBROUTINE DONEI(DLT,DSQ,CMO)
*
      IMPLICIT REAL*8 (A-H,O-Z)
*

#include "motra_global.fh"
*
      Real*8 CMO(*)
      DIMENSION DSQ(*),DLT(*)
*
      Call qEnter('Donei')
*
      ISTSQ=0
      ISTLT=0
      DO 100 ISYM=1,NSYM
        NF=NFRO(ISYM)
        NB=NBAS(ISYM)
        If ( NB*NF.GT.0 )
     &  CALL DGEMM_('N','T',
     &              NB,NB,NF,
     &              1.0d0,CMO(ISTSQ+1),NB,
     &              CMO(ISTSQ+1),NB,
     &              0.0d0,DSQ(ISTSQ+1),NB)
        CALL DSCAL_(NB*NB,2.0D0,DSQ(ISTSQ+1),1)
        IJ=ISTLT
        DO 130 IB=1,NB
          DO 140 JB=1,IB
            IJ=IJ+1
            DLT(IJ)=2.0D0*DSQ(ISTSQ+JB+(IB-1)*NB)
140       CONTINUE
          DLT(IJ)=0.5D0*DLT(IJ)
130     CONTINUE
        ISTSQ=ISTSQ+NB*NB
        ISTLT=ISTLT+NB*(NB+1)/2
100   CONTINUE
*
      IF( IPRINT.GE.5 .OR. DEBUG.NE.0 ) THEN
        WRITE(6,'(6X,A)')'Frozen one-body density matrix in AO basis'
        ISTLT=1
        DO ISYM=1,NSYM
          NB=NBAS(ISYM)
          IF ( NB.GT.0 ) THEN
            WRITE(6,'(6X,A,I2)')'symmetry species:',ISYM
            CALL TRIPRT(' ',' ',DLT(ISTLT),NB)
            ISTLT=ISTLT+NB*(NB+1)/2
          END IF
        END DO
      END IF
*
      Call qExit('Donei')
*
      RETURN
      END
