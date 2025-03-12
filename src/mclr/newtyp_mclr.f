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
* Copyright (C) 1993, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE NEWTYP_MCLR(INCLS,INTP,IACOP,ITPOP,NOP,OUTCLS,OUTTP)
*
* an input string of given group and occupation type are given
* apply an string of elementary operators to this group and
* obtain group and type of new string
*
* Jeppe Olsen, October 1993
*
* ------
* Input
* ------
*
* INGRP : Group of input string
* INTP  : Type of input string ( occupation type, # e in ras1,ras3)
* IACOP(I) = 1 : operator I is an annihilation operator
*          = 2 : operator I is a  creation   operator
* ITPOP(I) : orbitals space of operator I
* NOP : Number of operators
*
* Output
* ------
* OUTCLS : group of resulting string
* OUTTP  : Type of resulting string
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER ITPOP(*),IACOP(*)
*
      INTEGER OUTCLS,OUTTP
*. Number of electrons in RAS1,RAS3
      CALL  GTSTTP(INCLS,IEL1,IEL3,INTP,2)
*
      IDELTA = 0
      DO 100 IOP = 1, NOP
*. Change in number of orbitals
        IF(IACOP(IOP).EQ.1) THEN
          IDELTA = IDELTA - 1
        ELSE
          IDELTA = IDELTA + 1
        END IF
*. Change in RAS1, RAS3
        IF(ITPOP(IOP).EQ.1) THEN
           IF(IACOP(IOP).EQ.1) THEN
             IEL1 = IEL1 - 1
           ELSE
             IEL1 = IEL1 + 1
           END IF
         ELSE  IF(ITPOP(IOP).EQ.3) THEN
           IF(IACOP(IOP).EQ.1) THEN
             IEL3 = IEL3 - 1
           ELSE
             IEL3 = IEL3 + 1
           END IF
         END IF
  100 CONTINUE
*. Out class
      OUTCLS = INCLS - IDELTA
C?    write(6,*) ' OUTCLS,IEL1,IEL3 ',
C?   &             OUTCLS,IEL1,IEL3
*. out type
      CALL  GTSTTP(OUTCLS,IEL1,IEL3,OUTTP,1)
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' NEWTYP , OUTCLS, OUTTP ', OUTCLS,OUTTP
      END IF
*
      RETURN
      END
