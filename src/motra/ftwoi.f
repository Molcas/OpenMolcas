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
      SUBROUTINE FTWOI(DLT,DSQ,FLT,nFLT,FSQ,LBUF,X1,X2)
*
      IMPLICIT REAL*8 (A-H,O-Z)
*

#include "motra_global.fh"
#include "files_motra.fh"
*
      DIMENSION FSQ(*),FLT(nFLT),DSQ(*),DLT(*),X1(*),X2(*)
      DIMENSION NBSX(8),KEEP(8)
      Logical FoundTwoEls,ISQUAR
*
      Call qEnter('FTWOI')
*
      Call f_Inquire(FnTwoAo,FoundTwoEls)
      If (.not.FoundTwoEls) Then
        Write (6,*) 'FTwoi: OrdInt not found!'
        Call Abend()
      EndIf
*
      CALL OPNORD(IRC,0,FNTWOAO,LUTWOAO)
      CALL GETORD(IRC,ISQUAR,NSYM2,NBSX,KEEP)
*
*     Compare content of 1el and 2el integral file
*
      IF ( NSYM2.NE.NSYM ) THEN
        Write (6,*) 'FTwoi: NSYM2.NE.NSYM'
        Write (6,*) 'NSYM2=',NSYM2
        Write (6,*) 'NSYM=',NSYM
        Call QTrace()
        Call Abend()
      END IF
      DO ISYM=1,NSYM
        NB1=NBAS(ISYM)
        NB2=NBSX(ISYM)
        IF ( NB1.NE.NB2 ) THEN
           Write (6,*) 'FTwoi: NB1.NE.NB2'
           Write (6,*) 'NB1=',NB1
           Write (6,*) 'NB2=',NB2
           Call QTrace()
           Call Abend()
        END IF
      END DO
*                                                                      *
************************************************************************
*                                                                      *
      ExFac=1.0D0
      Call FockTwo(nSym,nBas,nFro,Keep,
     &             DLT,DSQ,FLT,nFLT,FSQ,LBUF,X1,X2,ExFac)
*                                                                      *
************************************************************************
*                                                                      *
*     Closed electron repulsion integral file
*
      CALL CLSORD(IRC,0)
*                                                                      *
************************************************************************
*                                                                      *
*     Print the Fock-matrix
*
      IF( IPRINT.GE.5 .OR. DEBUG.NE.0 ) THEN
        WRITE(6,'(6X,A)')'Fock matrix in AO basis'
        ISTLTT=1
        DO ISYM=1,NSYM
          NB=NBAS(ISYM)
          IF ( NB.GT.0 ) THEN
            WRITE(6,'(6X,A,I2)')'symmetry species:',ISYM
            CALL TRIPRT(' ',' ',FLT(ISTLTT),NB)
            ISTLTT=ISTLTT+NB*(NB+1)/2
          END IF
        END DO
      END IF
*
      Call qExit('FTWOI')
*
      RETURN
      END
