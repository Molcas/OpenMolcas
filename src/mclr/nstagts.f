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
*
* Number of strings of group IGRP
*                      type  ITP
*                      sym   ISM
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "detdim.fh"
#include "WrkSpc.fh"
#include "csm.fh"
#include "strbas_mclr.fh"
#include "stinf_mclr.fh"
* element (ITP,ISM) corresponds to adress
       IADRESS = (ISM-1)*NOCTYP(IGRP)+ ITP
*
       NSTAGTS = iWORK(KNSTSO(IGRP)+IADRESS-1)
*
       NTEST = 0
       IF(NTEST.NE.0) THEN
         WRITE(6,*) ' Number of strings of group,type,sym',
     &   IGRP,ITP,ISM,' is ', NSTAGTS
       END IF
*
      RETURN
      END
