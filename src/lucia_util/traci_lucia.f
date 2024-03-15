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
* Copyright (C) 1988, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE TRACI_LUCIA(      X,  LUCIN, LUCOUT,  IXSPC,   IXSM,
     &                          VEC1,   VEC2)
      use stdalloc, only: mma_allocate, mma_deallocate
*
* A rotation matrix X is defining expansion from
* old to new orbitals
*        PHI(NEW) = PHI(OLD) * X
*
* change CI coefficients(sym IXSM, space IXSPC )
* so they corresponds to PHI(NEW) basis
*
* The input CI vector is on LUCIN and the transformed CI vector
* will be delivered on LUCOUT.
*
* Transformation as conceived by Per-AAke Malmquist
* (I.J.Q.C. vol XXX, p479 ,1986 (OCTOBER ISSUE ))
*
*  Jeppe Olsen 1988
*
* New LUCIA version of Jan 1998
*
* note The transformation matrix X is supposed to be in complete form
* as a matrix over NTOOB orbitals.
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
#include "mxpdim.fh"
#include "orbinp.fh"
#include "lucinp.fh"
#include "clunit.fh"
*. Common block for communicating with sigma
#include "cands.fh"
*
      Real*8 X(*),VEC1(*),VEC2(*)
      Real*8, Allocatable:: SCR(:), LT(:)
* Some dummy initializations
      IOFF = 0 ! jwk-cleanup
*
      NTEST = 0
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' ================'
        WRITE(6,*) ' Welcome to TRACI '
        WRITE(6,*) ' ================'
        WRITE(6,*)
        WRITE(6,*) ' IXSPC,IXSM = ', IXSPC,IXSM
      END IF
*. Memory allocation
* for a matrix T
      Call mma_allocate(LT,NTOOB**2,Label='LT')
*. Scratch in PAMTMT
      Call mma_allocate(SCR,NTOOB**2 + NTOOB*(NTOOB+1)/2, Label='SCR')
*. Obtain T matrix used for transformation, for each symmetry separately
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IOFF = 1
        ELSE
          IOFF = IOFF + NTOOBS(ISM-1)**2
        END IF
        IF(NTOOBS(ISM).GT.0) THEN
         CALL PAMTMT(X(IOFF),LT(IOFF),SCR,NTOOBS(ISM))
        END IF
      END DO
*. Transform CI-vector
      ICSPC = IXSPC
      ICSM  = ICSM
      ISSPC = IXSPC
      ISSM  = IXSM
*
      CALL TRACID(LT,    LUCIN,   LUCOUT,    LUSC1,    LUSC2,
     &                LUSC3,     VEC1,     VEC2)
*
      Call mma_deallocate(SCR)
      Call mma_deallocate(LT)
*
*
      RETURN
      END
