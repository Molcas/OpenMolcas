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
      SUBROUTINE GASDIAT(    DIAG,   LUDIA,   ECORE,  ICISTR,     I12,
     &                      IBLTP,  NBLOCK,  IBLKFO)
      use stdalloc, only: mma_allocate, mma_deallocate
      use strbas
*
* CI diagonal in SD basis for state with symmetry ISM in internal
* space ISPC
*
* GAS version, Winter of 95
*
* Driven by table of TTS blocks, May97
*
      IMPLICIT REAL*8(A-H,O-Z)
* =====
*.Input
* =====
*
*./ORBINP/ : NACOB used
*
#include "mxpdim.fh"
#include "orbinp.fh"
#include "cicisp.fh"
#include "cstate.fh"
#include "strinp.fh"
#include "stinf.fh"
#include "csm.fh"
#include "cprnt.fh"
#include "cgas.fh"
#include "gasstr.fh"
#include "io_util.fh"
*
      DIMENSION IBLTP(*)
      DIMENSION IBLKFO(8,NBLOCK)
*
* ======
*.Output
* ======
      DIMENSION DIAG(*)
      Integer, Allocatable:: LASTR(:), LBSTR(:)
      Real*8, Allocatable:: LSCR2(:)
      Real*8, Allocatable:: LJ(:), LK(:), LXB(:), LH1D(:), LRJKA(:)
*
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRDIA)
*
** Specifications of internal space
*
      IATP = 1
      IBTP = 2
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NOCTPA = NOCTYP(IATP)
*
#ifdef _DEBUGPRINT_
      NOCTPB = NOCTYP(IBTP)
*. Offsets for alpha and beta supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ================'
        WRITE(6,*) ' GASDIA speaking '
        WRITE(6,*) ' ================'
        WRITE(6,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
        write(6,*) ' NOCTPA NOCTPB  : ', NOCTPA,NOCTPB
        write(6,*) ' IOCTPA IOCTPB  : ', IOCTPA,IOCTPB
      END IF
#endif
*
**. Local memory
*
      CALL mma_allocate(LJ   ,NTOOB**2,Label='LJ')
      CALL mma_allocate(LK   ,NTOOB**2,Label='LK')
      Call mma_allocate(LSCR2,2*NTOOB**2,Label='LSCR2')
      CALL mma_allocate(LXB  ,NACOB,Label='LXB')
      CALL mma_allocate(LH1D ,NACOB,Label='LH1D')
*. Space for blocks of strings
      Call mma_allocate(LASTR,MXNSTR*NAEL,Label='LASTR')
      Call mma_allocate(LBSTR,MXNSTR*NBEL,Label='LBSTR')
      MAXA = IMNMX(NSTSO(IATP)%I,NSMST*NOCTPA,2)
      CALL mma_allocate(LRJKA,MAXA,Label='LRJKA')
*
**. Diagonal of one-body integrals and coulomb and exchange integrals
*
      CALL GT1DIA(LH1D)
      CALL GTJK(LJ,LK,NTOOB,LSCR2,IREOTS,IREOST)
      IF( LUDIA .GT. 0 ) IDISK(LUDIA)=0
      CALL GASDIAS(NAEL,LASTR,NBEL,LBSTR,
     &             NACOB,DIAG,NSMST,
     &             LH1D,LXB,LJ,LK,
     &             NSTSO(IATP)%I,NSTSO(IBTP)%I,
     &             LUDIA,ECORE,PSSIGN,IPRDIA,NTOOB,ICISTR,
     &             LRJKA,I12,IBLTP,NBLOCK,IBLKFO,
     &             I_AM_OUT,N_ELIMINATED_BATCHES)
*.Flush local memory
      Call mma_deallocate(LJ)
      Call mma_deallocate(LK)
      Call mma_deallocate(LSCR2)
      Call mma_deallocate(LXB)
      Call mma_deallocate(LH1D)
      Call mma_deallocate(LASTR)
      Call mma_deallocate(LBSTR)
      Call mma_deallocate(LRJKA)
*
      END
