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
* Copyright (C) 1992,1994, Per Ake Malmqvist                           *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE PCLOSE
      IMPLICIT REAL*8 (A-H,O-Z)
C
C  PER-AAKE MALMQUIST 92-12-07
C  Deallocates everything concerned with SGUGA, incl CI array.


#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "pt2_guga.fh"

      IF(DoCumulant) RETURN
      IF(NACTEL.EQ.0) RETURN
      IF(ISCF.NE.0) RETURN
      CALL GETMEM('MVL','FREE','INTEG',LMVL,NMVL)
      CALL GETMEM('MVR','FREE','INTEG',LMVR,NMVR)
      CALL GETMEM('NOW','FREE','INTEG',LNOW,NNOW)
      CALL GETMEM('IOW','FREE','INTEG',LIOW,NIOW)
      CALL GETMEM('NOCP','FREE','INTEG',LNOCP,NNOCP)
      CALL GETMEM('IOCP','FREE','INTEG',LIOCP,NIOCP)
      CALL GETMEM('NOCSF','FREE','INTEG',LNOCSF,NNOCSF)
      CALL GETMEM('IOCSF','FREE','INTEG',LIOCSF,NIOCSF)
      CALL GETMEM('ICASE','FREE','INTEG',LICASE,NICASE)
      CALL GETMEM('ICOUP','FREE','INTEG',LICOUP,(3*NICOUP+1)/2)
      CALL GETMEM('VTAB','FREE','REAL',LVTAB,NVTAB)
      CALL GETMEM('SGTMP','FREE','REAL',LSGTMP,NSGTMP)
      RETURN
      END
