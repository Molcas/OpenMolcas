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
      use GLBBAS, Only: INT1O, PINT1
*
* Obtain diagonal of one electron matrix over active
* orbitals
*
*. Dec 97 : obtained from INT1O
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.fh"
      DIMENSION H1DIA(*)

*.GLobal pointers

#include "lucinp.fh"
#include "orbinp.fh"
*
CINA  CALL GT1DIS(H1DIA,IREOTS(1+NINOB),PINT1,INT1,ISMFTO,IBSO,NACOB)
      CALL GT1DIS(H1DIA,IREOTS(1),PINT1,INT1O,ISMFTO,IBSO,NACOB)
*
      END
