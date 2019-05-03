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
      SUBROUTINE SIGMADET_CVB(C,HC,IREFSM,PERMS2,NCI)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL PERMS2
#include "WrkSpc.fh"
#include "rasscf_lucia.fh"
      DIMENSION C(NCI),HC(NCI)
      DIMENSION DUMMY(1)
C
C Export arguments to be used in sigma_master_cvb
C
      C_POINTER = KCI_POINTER
      CALL DCOPY_(NCI,C,1,WORK(C_POINTER),1)
C Call the sigma routine
      CALL LUCIA_UTIL('SIGMA_CVB',IREFSM,IDUMMY,DUMMY)
      CALL DCOPY_(NCI,WORK(KSIGMA_POINTER),1,HC,1)
C
      RETURN
c Avoid unused argument warnings
      If (.False.) Call Unused_logical(PERMS2)
      END
