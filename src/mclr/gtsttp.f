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
* Copyright (C) Jeppe Olsen                                            *
************************************************************************
      SUBROUTINE GTSTTP(ICLS,IEL1,IEL3,ITYPE,IWAY)
      Use Str_Info
*
*
* Relation between number of electrons in RAS1, RAS3 and string type
*
* IWAY = 1 :
* Get ITYPE : type of strings belonging to class ICLS
*             with IEL1,IEL3 electrons
* IWAY = 2 :
* GET IEL1,IEL3 : Number of electrons of class ICLS of type ITYPE
*
*
* Jeppe Olsen, Another lonely night in Lund
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
#include "detdim.fh"
*
      CALL GTSTTPS(IEL1,IEL3,
     &             Str(MIN(ICLS,ITYP_Dummy))%EL1,
     &             Str(MIN(ICLS,ITYP_Dummy))%EL3,
     &             NOCTYP(ICLS),ITYPE,IWAY)
*
      RETURN
      END
