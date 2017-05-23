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
* Copyright (C) 2001, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE REFORM_CONF_OCC(IOCC_EXP,IOCC_PCK,NEL,NOCOB,IWAY)
*
* Reform between two ways of writing occupations
*
* IOCC_EXP : Occupation in expanded form, i.e. the orbital for each
*            electron is given
*
* IOCC_PCK  : Occupation is given in packed form, i.e. each occupied
*             orbitals is given once, and a negative index indicates
*             a double occupattion
*
* IWAY = 1 Expanded to Packed form
* IWAY = 2 Packed to expanded form
*
* Jeppe Olsen, Nov. 2001
*
#include "implicit.fh"
*. Input/Output
      INTEGER IOCC_EXP(NEL),IOCC_PCK(NOCOB)
*
      IF(IWAY.EQ.1) THEN
*
*. Expanded => Packed form
*
*. Loop over electrons
        IEL = 1
        IOCC = 0
 1000   CONTINUE
C         IEL = IEL + 1
          IF(IEL.LT.NEL) THEN
            IF(IOCC_EXP(IEL).EQ.IOCC_EXP(IEL+1)) THEN
              IOCC = IOCC + 1
              IOCC_PCK(IOCC) = -IOCC_EXP(IEL)
              IEL = IEL + 2
            ELSE
              IOCC = IOCC + 1
              IOCC_PCK(IOCC) =  IOCC_EXP(IEL)
              IEL = IEL + 1
            END IF
          ELSE
*. Last occupation was not identical to previous, so single occupied
            IOCC = IOCC + 1
            IOCC_PCK(IOCC) =  IOCC_EXP(IEL)
            IEL = IEL + 1
          END IF
        IF(IEL.LE.NEL) GOTO 1000
*
      ELSE IF( IWAY.EQ.2) THEN
*
* Packed to expanded form
*
        IEL = 0
        DO IORB = 1, NOCOB
          IF(IOCC_PCK(IORB).LT.0) THEN
            JORB = - IOCC_PCK(IORB)
            IEL = IEL +1
            IOCC_EXP(IEL) = JORB
            IEL = IEL + 1
            IOCC_EXP(IEL) = JORB
          END IF
        END DO
      ELSE
        WRITE(6,*) ' REFORM_CONF... in error, IWAY = ', IWAY
*        STOP       ' REFORM_CONF... in error, IWAY '
         CALL SYSABENDMSG('lucia_util/reform_conv',
     &                    'Internal error',' ')
      END IF
*     ^ End of IWAY switch
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Reforming form of configuration '
        IF(IWAY.EQ.1) THEN
           WRITE(6,*) ' Expanded to packed form '
        ELSE
           WRITE(6,*) ' Packed to expanded form '
        END IF
        WRITE(6,*) ' IOCC_EXP : '
        CALL IWRTMA(IOCC_EXP,1,NEL,1,NEL)
        WRITE(6,*) ' IOCC_PCK : '
        CALL IWRTMA(IOCC_PCK,1,NOCOB,1,NOCOB)
      END IF
*
      RETURN
      END
