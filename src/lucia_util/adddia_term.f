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
* Copyright (C) 2011, Jeppe Olsen                                      *
*               2011, Giovanni Li Manni                                *
************************************************************************
      SUBROUTINE ADDDIA_TERM(FACTOR,CVEC,SVEC,IASPGP,IBSPGP,IASM,IBSM)
      use stdalloc, only: mma_allocate, mma_deallocate
      use strbas
*. Update Sigma vector with diagonal terms for a given block
*     SVEC(IASPGP,IBSPGP) = SVEC(IASPGP,IBSPGP)
*                         + (FACTOR+DIAG(IASPGP,IBSPGP))CVEC(IASPGP,IBSPGP)
*
* Jeppe Olsen and Giovanni Li Manni, September 2011
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.fh"
#include "strinp.fh"
#include "orbinp.fh"
#include "gasstr.fh"
#include "csm.fh"
#include "stinf.fh"
#include "cecore.fh"
#include "cprnt.fh"
*
*. Input
      DIMENSION CVEC(*)
*. Output
      DIMENSION SVEC(*)
      Integer, Allocatable:: LASTR(:), LBSTR(:)
      Real*8, Allocatable:: LSCR(:), LSCR2(:)
      Real*8, Allocatable:: LJ(:), LK(:), LXA(:), LXB(:), LRJKA(:),
     &                      LH1D(:)
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRDIA)
*
      IATP   = 1
      IBTP   = 2
      NAEL   = NELEC(IATP)
      NBEL   = NELEC(IBTP)
      NOCTPA = NOCTYP(IATP)
*
#ifdef _DEBUGPRINT_
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ========================='
        WRITE(6,*) '   ADDDIA_TERM for BK is speaking '
        WRITE(6,*) ' ========================='
        WRITE(6,*) ' NAEL NBEL =', NAEL,NBEL
        WRITE(6,*) ' IASPGP, IBSPGP = ', IASPGP, IBSPGP
      END IF
#endif
*. A bit of scracth
      CALL mma_allocate(LH1D ,NTOOB,Label='LH1D')
      CALL mma_allocate(LJ   ,NTOOB**2,Label='LJ')
      CALL mma_allocate(LK   ,NTOOB**2,Label='LK')
      Call mma_allocate(LSCR2,2*NTOOB**2,Label='LSCR2')
      CALL mma_allocate(LXA  ,NACOB,Label='LXA')
      CALL mma_allocate(LXB  ,NACOB,Label='LXB')
      Call mma_allocate(LSCR,2*NACOB,Label='LSCR')
*. Space for blocks of strings
      Call mma_allocate(LASTR,MXNSTR*NAEL,Label='LASTR')
      Call mma_allocate(LBSTR,MXNSTR*NBEL,Label='LBSTR')

      MAXA = IMNMX(NSTSO(IATP)%I,NSMST*NOCTPA,2)
      CALL mma_allocate(LRJKA,MAXA,Label='LRJKA')
*. Diagonal of one-body integrals and coulomb and exchange integrals
*. Integrals assumed in place so :
      CALL GT1DIA(LH1D)
      CALL GTJK(LJ,LK,NTOOB,LSCR2,IREOTS,IREOST)
*. Core energy not included
      ECOREP = 0.0D0
      CALL GTJK(LJ,LK,NTOOB,LSCR2,IREOTS,IREOST)
*
      SHIFT = ECORE_ORIG-ECORE
      FACTORX = FACTOR + SHIFT
*
      CALL ADDDIA_TERMS(NAEL,LASTR,NBEL,LBSTR,
     &             NACOB,CVEC,SVEC,NSMST,LH1D,
     &             LXA,LXB,LSCR,LJ,
     &             LK,NSTSO(IATP)%I,NSTSO(IBTP)%I,
     &             ECOREP,
     &             IPRDIA,NTOOB,
     &             LRJKA,
     &             IASPGP,IASM,IBSPGP,IBSM,FACTORX)
*.Flush local memory
      CALL mma_deallocate(LH1D)
      CALL mma_deallocate(LJ)
      CALL mma_deallocate(LK)
      Call mma_deallocate(LSCR2)
      CALL mma_deallocate(LXA)
      CALL mma_deallocate(LXB)
      Call mma_deallocate(LSCR)
      Call mma_deallocate(LASTR)
      Call mma_deallocate(LBSTR)
      CALL mma_deallocate(LRJKA)

      END
