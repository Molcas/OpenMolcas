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
      FUNCTION ILEX_FOR_CONF(ICONF,NOCC_ORB,NORB,NEL,IARCW,IDOREO,IREO)
*
* A configuration ICONF of NOCC_ORB orbitals are given
* ICONF(I) = IORB implies  IORB is singly occupied
* ICONF(I) = -IORB  implies that IORB is doubly occupied
*
* Find lexical address
*
* IF IDOREO .ne. 0, IREO is used to reorder lexical number
* Jeppe Olsen, November 2001
*
#include "implicit.fh"
*. Arcweights for single and doubly occupied arcs
      INTEGER IARCW(NORB,NEL,2)
*. Reorder array
      INTEGER IREO(*)
*. Configuration
      INTEGER ICONF(NOCC_ORB)
*
      IEL = 0
      ILEX = 1

      DO IOCC = 1, NOCC_ORB
       IF(ICONF(IOCC).GT.0) THEN
         IEL = IEL + 1
         ILEX = ILEX + IARCW(ICONF(IOCC),IEL,1)
       ELSE IF(ICONF(IOCC).LT.0) THEN
         IEL = IEL + 2
         ILEX = ILEX + IARCW(-ICONF(IOCC),IEL,2)
       END IF
      END DO
*
      IF(IDOREO.NE.0) THEN
       ILEX_FOR_CONF = IREO(ILEX)
      ELSE
       ILEX_FOR_CONF = ILEX
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Configuration '
        CALL IWRTMA(ICONF,1,NOCC_ORB,1,NOCC_ORB)
        WRITE(6,*) ' Lexical number = ', ILEX
        IF(IDOREO.NE.0)
     &  WRITE(6,*) ' Reordered number = ', ILEX_FOR_CONF
      END IF
*
      RETURN
      END
C       CALL GEN_CONF_FOR_OCCLS(IB_OCCLS,
C    &     INITIALIZE_CONF_COUNTERS,
C    &     IOCCLS,NGAS,ISYM,MINOP,MAXOP,NSMST,1,NOCOB,
C    &     NOBPT,NCONF_PER_SYM(ISYM),
C    &     NCONF_PER_OPEN(1,ISYM),IBCONF_PER_OPEN(1,ISYM),
C    &     IDUM,IDUM,IDUM,
C    &     IDOREO,IDUMMY,IDUMMY,NCONF_ALL_SYM)
