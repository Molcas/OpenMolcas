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
* Copyright (C) 1993,1995, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE NEWTYP(INSPGP,IACOP,ITPOP,OUTSPGP)
*
* an input  supergroup is given.
* apply an string of elementary operators to this supergroup and
* obtain supergroup mumber of new string
*
* Jeppe Olsen, October 1993
* GAS-version : July 95
*
* ------
* Input
* ------
*
* INSPGP  : input super group ( given occupation in each AS )
* IACOP = 1 : operator is an annihilation operator
*       = 2 : operator is a  creation   operator
* ITPOP : orbitals space of operator
*
* Output
* ------
* OUTSPGP  : supergroup of resulting string
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.fh"
#include "cgas.fh"
#include "strbas.fh"
#include "WrkSpc.fh"
*. Input
      INTEGER INSPGP,IACOP,ITPOP
*. output
      INTEGER OUTSPGP
*
      IF(IACOP.EQ.1) THEN
        OUTSPGP = IWORK(KSPGPAN+ITPOP-1+NGAS*(INSPGP-1))
      ELSE
        OUTSPGP = IWORK(KSPGPCR+ITPOP-1+NGAS*(INSPGP-1))
      END IF
*
      RETURN
      END
