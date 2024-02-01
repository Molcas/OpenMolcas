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
      SUBROUTINE GT1DIA(H1DIA)
      use GLBBAS
*
* Obtain diagonal of one electron matrix over active
* orbitals
*
*. Dec 97 : obtained from KINT1O
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.fh"
#include "WrkSpc.fh"
      DIMENSION H1DIA(*)

*.GLobal pointers

#include "lucinp.fh"
#include "orbinp.fh"
*
CINA  CALL GT1DIS(H1DIA,IREOTS(1+NINOB),WORK(KPINT1),WORK(KINT1),
CINA &            ISMFTO,IBSO,NACOB)
      CALL GT1DIS(    H1DIA,IREOTS(1),IWORK(KPINT1),WORK(KINT1O),ISMFTO,
     &                 IBSO,    NACOB)
*
      RETURN
      END
