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
* Copyright (C) 2019, Stefano Battaglia                                *
************************************************************************
      SUBROUTINE TRANSMAT(MAT,EVEC,NGRP)
      IMPLICIT REAL*8 (A-H,O-Z)

* This subroutine is a wrapper for the diagonalization of an arbitrary
* matrix.

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "warnings.fh"

      REAL*8 MAT(NSTATE,NSTATE)

      DIMENSION EVEC(NGRP**2)

      IF (NSTATE.NE.NGRP) THEN
        WRITE(6,*) ' Number of states in the group has to correspond'
        WRITE(6,*) ' to the total number of states! Aborting...'
        CALL ABEND
      END IF

      CALL GETMEM('MTMP1','ALLO','REAL',LMTMP1,NGRP**2)
      CALL GETMEM('MTMP2','ALLO','REAL',LMTMP2,NGRP**2)

      DO I=1,NGRP
        DO J=1,NGRP
          WORK(LMTMP1+I-1+NGRP*(J-1))=MAT(I,J)
        END DO
      END DO

      CALL DGEMM_('T','N',NGRP,NGRP,NGRP,
     &            1.0d0,EVEC,NGRP,WORK(LMTMP1),NGRP,
     &            0.0d0,WORK(LMTMP2),NGRP)
      CALL DGEMM_('N','N',NGRP,NGRP,NGRP,
     &            1.0d0,WORK(LMTMP2),NGRP,EVEC,NGRP,
     &            0.0d0,WORK(LMTMP1),NGRP)

      DO I=1,NGRP
        DO J=1,NGRP
          MAT(I,J)=WORK(LMTMP1+I-1+NGRP*(J-1))
        END DO
      END DO

      CALL GETMEM('HTMP1','FREE','REAL',LMTMP1,NGRP**2)
      CALL GETMEM('HTMP2','FREE','REAL',LMTMP2,NGRP**2)

      RETURN
      END
