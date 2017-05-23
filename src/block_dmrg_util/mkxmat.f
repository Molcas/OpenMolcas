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
* Copyright (C) 2014, Naoki Nakatani                                   *
************************************************************************
      SUBROUTINE MKXMAT(TORB,XMAT)
      IMPLICIT REAL*8 (A-H,O-Z)
* Make full transformation matrix for active space
* from that stored for each symmetry in IAD1M(4)
* Written by N. Nakatani, Oct. 2014
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
      DIMENSION TORB(NTORB), XMAT(NASHT,NASHT)
      DIMENSION UMAT(NASHT,NASHT)

      IF(NASHT.GT.0) THEN
        ITOEND=0
        DO ISYM=1,NSYM
          NI=NISH(ISYM)
          NA=NASH(ISYM)
          NR1=NRAS1(ISYM)
          NR2=NRAS2(ISYM)
          NR3=NRAS3(ISYM)
          NS=NSSH(ISYM)
          ITOSTA=ITOEND+1
          ITOEND=ITOEND+NI**2+NR1**2+NR2**2+NR3**2+NS**2

* Normally, NRAS2 should only be non-zero, but NRAS1 and NRAS3 should be zero
* for DMRG-CASSCF calculation
          ITO=ITOSTA+NI**2
C         ITO=ITOSTA+NI**2-1
* RAS1
          IF(NA.GT.0) THEN
            IF(NR1.GT.0) THEN
              ISTART=NAES(ISYM)
              DO JR1=1,NR1
                J=ISTART+JR1
                DO IR1=1,NR1
                  I=ISTART+IR1
                  XMAT(I,J)=TORB(ITO)
                  ITO=ITO+1
                END DO
              END DO
            END IF
* RAS2
            IF(NR2.GT.0) THEN
              ISTART=NAES(ISYM)+NR1
              DO JR2=1,NR2
                J=ISTART+JR2
                DO IR2=1,NR2
                  I=ISTART+IR2
                  XMAT(I,J)=TORB(ITO)
                  ITO=ITO+1
                END DO
              END DO
            END IF
* RAS3
            IF(NR3.GT.0) THEN
              ISTART=NAES(ISYM)+NR1+NR2
              DO JR3=1,NR3
                J=ISTART+JR3
                DO IR3=1,NR3
                  I=ISTART+IR3
                  XMAT(I,J)=TORB(ITO)
                  ITO=ITO+1
                END DO
              END DO
            END IF
          END IF
        END DO
      END IF

      RETURN
      END
