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
      SUBROUTINE TRANSMAT(MAT,EVEC,NDIM)
      IMPLICIT REAL*8 (A-H,O-Z)

* This subroutine is a wrapper for the diagonalization of an arbitrary
* matrix.

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"

      REAL*8 MAT(NDIM,NDIM)

      DIMENSION EVEC(NDIM**2)

      CALL GETMEM('MTMP1','ALLO','REAL',LMTMP1,NDIM**2)
      CALL GETMEM('MTMP2','ALLO','REAL',LMTMP2,NDIM**2)

      DO I=1,NDIM
        DO J=1,NDIM
          WORK(LMTMP1+I-1+NDIM*(J-1))=MAT(I,J)
        END DO
      END DO

      IF (IPRGLB.GE.DEBUG) THEN
        WRITE(6,*)' MAT in TRANSMAT:'
        DO I=1,NSTATE
          WRITE(6,'(1x,5f16.8)')(MAT(I,J),J=1,NSTATE)
        END DO
      END IF

      CALL DGEMM_('T','N',NDIM,NDIM,NDIM,
     &            1.0d0,EVEC,NDIM,WORK(LMTMP1),NDIM,
     &            0.0d0,WORK(LMTMP2),NDIM)
      CALL DGEMM_('N','N',NDIM,NDIM,NDIM,
     &            1.0d0,WORK(LMTMP2),NDIM,EVEC,NDIM,
     &            0.0d0,WORK(LMTMP1),NDIM)

      DO I=1,NDIM
        DO J=1,NDIM
          MAT(I,J)=WORK(LMTMP1+I-1+NDIM*(J-1))
        END DO
      END DO

      CALL GETMEM('HTMP1','FREE','REAL',LMTMP1,NDIM**2)
      CALL GETMEM('HTMP2','FREE','REAL',LMTMP2,NDIM**2)

      RETURN
      END
