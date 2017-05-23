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
* Copyright (C) 2000, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 2000  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE GETSGM2(ILEV,JLEV,ISYCI,CI,SGM)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "pt2_guga.fh"
      DIMENSION  CI(MXCI),SGM(MXCI)

C GIVEN CI COUPLING LEVELS ILEV, JLEV, COMPUTE SGM=E(ILEV,JLEV)*CI
C ILEV,JLEV ARE IN PRINCIPLE ACTIVE ORBITAL NUMBERS, BUT POSSIBLY
C IN ANOTHER ORDER THAN THE USUAL ONE -- HERE WE USE THE ORDER
C FOLLOWED BY THE GUGA COUPLING SCHEME.
C
C THIS ROUTINE REPLACES EARLIER GETSGM, TO GET RID OF THE PACKING AND
C STORING USED EARLIER.

C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C NOTE!! THE EARLIER CALL GETSGM(ILEV,JLEV,IDARR,SGM) IS REPLACED BY
C GETSGM2(ILEV,JLEV,CI,SGM)!!
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IS=ISM(ILEV)
      JS=ISM(JLEV)
      IJS=MUL(IS,JS)
      ISSG=MUL(IJS,ISYCI)
      NSGM=NCSF(ISSG)
      IF(NSGM.EQ.0) RETURN
      CALL DCOPY_(NSGM,0.0D0,0,SGM,1)
      CALL SIGMA1_CP2(ILEV,JLEV,1.0D00,ISYCI,CI,SGM,
     &      IWORK(LNOCSF),IWORK(LIOCSF),IWORK(LNOW),IWORK(LIOW),
     &      IWORK(LNOCP),IWORK(LIOCP),IWORK(LICOUP),
     &      WORK(LVTAB),IWORK(LMVL),IWORK(LMVR))
      RETURN
      END
