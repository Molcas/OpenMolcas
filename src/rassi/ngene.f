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
* Copyright (C) 1998, Per Ake Malmqvist                                *
************************************************************************
      INTEGER FUNCTION NGENE(NEL,MLTPL)
C Nr of genealogical spin couplings
      NGENE=0
      IF(MLTPL.LE.0) RETURN
      IS2=MLTPL-1
      IF(NEL.LT.IS2) RETURN
      NU=(NEL+IS2)/2
      ND=(NEL-IS2)/2
      IF(NU+ND.NE.NEL) RETURN
      NGENE=NOVERM(NEL,NU)-NOVERM(NEL,NU+1)
      RETURN
      END
