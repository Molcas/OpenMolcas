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
      SUBROUTINE ZTRNSF(N,UR,UI,AR,AI)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION UR(N,N),UI(N,N)
      DIMENSION AR(N,N),AI(N,N)
#include "WrkSpc.fh"

      CALL GETMEM('TMPREAL','ALLO','REAL',LCR,N**2)
      CALL GETMEM('TMPIMAG','ALLO','REAL',LCI,N**2)
      CALL DGEMM_('N','N',N,N,N, 1.0D0,AR,N,UR,N,0.0D0,WORK(LCR),N)
      CALL DGEMM_('N','N',N,N,N,-1.0D0,AI,N,UI,N,1.0D0,WORK(LCR),N)
      CALL DGEMM_('N','N',N,N,N, 1.0D0,AR,N,UI,N,0.0D0,WORK(LCI),N)
      CALL DGEMM_('N','N',N,N,N, 1.0D0,AI,N,UR,N,1.0D0,WORK(LCI),N)
      CALL DGEMM_('T','N',N,N,N, 1.0D0,UR,N,WORK(LCR),N,0.0D0,AR,N)
      CALL DGEMM_('T','N',N,N,N, 1.0D0,UI,N,WORK(LCI),N,1.0D0,AR,N)
      CALL DGEMM_('T','N',N,N,N, 1.0D0,UR,N,WORK(LCI),N,0.0D0,AI,N)
      CALL DGEMM_('T','N',N,N,N,-1.0D0,UI,N,WORK(LCR),N,1.0D0,AI,N)
      CALL GETMEM('TMPREAL','FREE','REAL',LCR,N**2)
      CALL GETMEM('TMPIMAG','FREE','REAL',LCI,N**2)

      RETURN
      END
