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
      SUBROUTINE diag_master
*
*  To do in this subroutine:
*
*     - Make sure all calling parameters are accounted for
*     - Make sure mxntts is set in previous subroutine
*     - Make sure nsmst is set in previous subroutine
*
*  Set up the diagonal for the CI calculation
*
      implicit real*8 (a-h,o-z)
#include "mxpdim.fh"
#include "cstate.fh"
#include "cprnt.fh"
#include "WrkSpc.fh"
#include "orbinp.fh"
#include "glbbas.fh"
#include "cands.fh"
#include "rasscf_lucia.fh"
*
      NDIM = NTOOB**2
      CALL COPVEC(WORK(KINT1O_POINTER),WORK(KINT1_POINTER),NDIM)
      CALL GASCI(ISSM, 1, IPRDIA, EREF, 0, 0)
*
      RETURN
      END
*
