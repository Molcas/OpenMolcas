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
      SUBROUTINE GETDREF(DREF,NDREF)
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: debug
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT REAL*8 (A-H,O-Z)
#include "caspt2.fh"
#include "pt2_guga.fh"
      INTEGER NDREF
      REAL*8 DREF(NDREF)

      REAL*8, ALLOCATABLE:: G1(:)

* Get active 1-el density matrix GAMMA1 and
* construct DREF in a tringular storage.

* Remember: NDREF=1 if NASHT=0.
      DREF(1)=0.0d0
      IF(NASHT.EQ.0) RETURN
* Active 1-el density matrix:
      CALL mma_allocate(G1,NG1,Label='G1')
      CALL PT2_GET(NG1,'GAMMA1',G1)
      DO I=1,NASHT
        DO J=1,I
          IJ=(I*(I-1))/2+J
          DREF(IJ)=G1(I+NASHT*(J-1))
        END DO
      END DO
      CALL mma_deallocate(G1)

      IF(IPRGLB.GE.DEBUG) THEN
        WRITE(6,*)' GETDREF has constructed DREF.'
        CALL XFLUSH(6)
      END IF

      END SUBROUTINE GETDREF

