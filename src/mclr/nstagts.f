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
      FUNCTION NSTAGTS(IGRP,ITP,ISM)
      Use Str_Info
*
* Number of strings of group IGRP
*                      type  ITP
*                      sym   ISM
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "detdim.fh"
#include "csm.fh"
#include "stinf_mclr.fh"
* element (ITP,ISM) corresponds to adress
       IADRESS = (ISM-1)*NOCTYP(IGRP)+ ITP
*
       NSTAGTS = Str(IGRP)%NSTSO(IADRESS)
*
       NTEST = 0
       IF(NTEST.NE.0) THEN
         WRITE(6,*) ' Number of strings of group,type,sym',
     &   IGRP,ITP,ISM,' is ', NSTAGTS
       END IF
*
      RETURN
      END
