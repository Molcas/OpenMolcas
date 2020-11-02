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
       FUNCTION GETH1I_MCLR(I,J)
       use Arrays, only: FIMO
*
* Obtain one -electron integral H(IORB,JOB)
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "detdim.fh"
#include "Input.fh"
#include "Pointers.fh"
#include "orbinp_mclr.fh"
#include "lbbas1.fh"
*
      ISM = ISMFTO(I)
      JSM = ISMFTO(J)
      I1 = IREOTS(I)-IBSO(ISM)+nISH(ISM)+1
      J1 = IREOTS(J)-IBSO(JSM)+nISH(JSM)+1

      IJ=ipCM(iSM)-1+(J1-1)*NORB(ISM)+I1
*
      GETH1I_MCLR = FIMO(IJ)
*
*
      RETURN
      END
