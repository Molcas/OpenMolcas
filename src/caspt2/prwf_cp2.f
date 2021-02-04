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
* Copyright (C) 1994, Per Ake Malmqvist                                *
**********************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE PRWF_CP2(ISYCI,NCO,CI,THR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ISYCI, NCO
      REAL*8, INTENT(IN) :: CI(NCO), THR
#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"

      INTEGER I

      WRITE(6,'(20A4)')('----',I=1,20)
      WRITE(6,'(a,d9.2)')' CI COEFFICIENTS LARGER THAN ',THR
      CALL PRWF1_CP2(IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),
     &           IWORK(LIOW),ISYCI,CI,THR)
      RETURN
      END
