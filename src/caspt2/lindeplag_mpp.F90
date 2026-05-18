!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

#if defined(_MOLCAS_MPP_) && defined(_GA_)
      Subroutine LinDepLag_MPP(lg_BDER,lg_SDER,nAS,nIN,iSym,iCase)

#ifdef _SCALAPACK_
      use scalapack_mod, only: GA_PDSYEVX_
#endif
      use caspt2_global, only: LUSTD, idBoriMat
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp, u6
      use caspt2_module, only: THRSHS
      use Constants, only: Zero, One, Two

      implicit none

#include "global.fh"
#include "mafdecls.fh"

      integer(kind=iwp), intent(in) :: lg_BDER, lg_SDER, nAS, nIN, iSym,&
     &  iCase

#if ! defined (_SCALAPACK_)
      real(kind=wp) :: WGRONK(2)
      integer(kind=iwp) :: NSCRATCH, info
      real(kind=wp), allocatable :: VEC(:),SCRATCH(:)
#endif
      logical(kind=iwp) :: bStat
      integer(kind=iwp) :: lg_S, myRank, lg_Vec, iLo, iHi,              &
     &  jLo, jHi, mV, LDV, I, lg_Lag, lg_B, IDB, mB, LDB, J
      real(kind=wp) :: EVAL, FACT
      real(kind=wp), allocatable :: EIG(:)
!
!     Parallel LinDepLag
!     We always use the canonical orthonormalization.
!
      !! Obtain the X matrix
      !! First, read S
      CALL PSBMAT_GETMEM ('S',lg_S,NAS)
      CALL PSBMAT_READ ('S',iCase,iSym,lg_S,NAS)

      call mma_allocate(EIG,NAS,Label='EIG')
      EIG(:) = Zero

      myRank = GA_NodeID()
#ifdef _SCALAPACK_
      CALL PSBMAT_GETMEM('VMAT',lg_Vec,NAS)
      CALL GA_PDSYEVX_ (lg_S, lg_Vec, EIG, 0)
      bSTAT = GA_Destroy (lg_S)
#else
! here for the non-ScaLAPACK version: copy matrix to master process,
! diagonalize using the serial DSYEV routine, and copy the resulting
! eigenvectors back to a global array.  Then distribute the eigenvalues.
      IF (myRank == 0) THEN
        call mma_allocate(VEC,NAS*NAS,Label='VEC')
        CALL GA_Get (lg_S, 1, NAS, 1, NAS, VEC, NAS)
      END IF
      bSTAT = GA_Destroy (lg_S)
      IF (myRank == 0) THEN
        CALL DSYEV_('V','L',NAS,VEC,NAS,EIG,WGRONK,-1,INFO)
        NSCRATCH=INT(WGRONK(1))
        call mma_allocate(SCRATCH,NSCRATCH,Label='SCRATCH')
        CALL DSYEV_('V','L',NAS,VEC,NAS,EIG,SCRATCH,NSCRATCH,INFO)
        call mma_deallocate(SCRATCH)
      END IF
      CALL PSBMAT_GETMEM('VMAT',lg_Vec,NAS)
      IF (myRank == 0) THEN
        CALL GA_Put (lg_Vec, 1, NAS, 1, NAS, VEC, NAS)
        call mma_deallocate(VEC)
      END IF
      CALL GADGOP(EIG,NAS,'+')
#endif

      !! Scale only the independent vectors to avoid
      !! any numerically unstable computation
! Form orthonormal transformation vectors by scaling the eigenvectors.
      call GA_Sync()
      call GA_Distribution (lg_Vec, myRank, iLo, iHi, jLo, jHi)
      IF (iLo /= 0) THEN
        call GA_Access (lg_Vec, iLo, iHi, jLo, jHi, mV, LDV)
        IF ((jHi-jLo+1) /= NAS) THEN
          WRITE(u6,*) 'SBDIAG_MPP: error in striping of lg_Vec, ABORT'
          CALL ABEND()
        END IF
        DO I=1,jHi-jLo+1 ! NAS
          EVAL=EIG(I)
          IF(EVAL < THRSHS) CYCLE
          FACT=One/SQRT(EVAL)
          Call DScal_(iHi-iLo+1,FACT,DBL_MB(mV+LDV*(I-1)),1)
        END DO
        call GA_Release_Update (lg_Vec, iLo, iHi, jLo, jHi)
      END IF
      call GA_Sync()

      call mma_deallocate(EIG)
!
!     X matrix has been prpeared
!
      CALL GA_CREATE_STRIPED ('H',NAS,NAS,'Lag',lg_Lag)
      CALL PSBMAT_GETMEM ('B',lg_B,NAS)

      call GA_Distribution (lg_B, myRank, iLo, iHi, jLo, jHi)
      call GA_Access (lg_B, iLo, iHi, jLo, jHi, mB, LDB)
      IDB=IDBoriMat(ISYM,ICASE)
      CALL DDAFILE(LUSTD,2,DBL_MB(mB),(iHi-iLo+1)*(jHi-jLo+1),IDB)
      call GA_Release_Update (lg_B, iLo, iHi, jLo, jHi)

      !! Compute the partial derivative
      !! Work(LF)  : B --> lg_B
      !! BDER      : D --> lg_BDER
      !! Work(LVEC): X^0 and X --> lg_Vec
      Call GA_DGEMM ('N','T',NAS,NAS,NAS,                               &
     &               Two,lg_B,lg_BDER,                                  &
     &               Zero,lg_Lag)
      Call GA_DGEMM ('N','N',NAS,NAS,NAS,                               &
     &               One,lg_Lag,lg_VEC,                                 &
     &               Zero,lg_B)
      Call GADupl(lg_B,lg_Lag)

      Call GA_DGEMM ('T','N',NAS,NAS,NAS,                               &
     &               One,lg_Vec,lg_Lag,                                 &
     &               Zero,lg_B)
      !! At this point,
      !! Work(LF) = 2 \mathcal{X}^0 * B * D * \mathcal{X}

      !! remove dependent part
      !! (linearly indep-indep and dep-dep)
      call GA_Distribution (lg_B, myRank, iLo, iHi, jLo, jHi)
      if (iLo /= 0) then
        call GA_Access (lg_B, iLo, iHi, jLo, jHi, mV, LDV)
!       if (ilo <= nas-nin) then
          Do J = 1, nAS-nIN
            Do I = 1, min(iHi-iLo+1,nAS-nIN-iLo+1)
              DBL_MB(mV+i-1+LDV*(j-1)) = Zero
            End Do
          End Do
!       end if
!       if (ilo >= nas-nin+1) then
          Do J = nAS-nIN+1, nAS
            Do I = max(1,nAS-nIN+1-iLo+1), LDV
              DBL_MB(mV+i-1+LDV*(j-1)) = Zero
            End Do
          End Do
!       end if
        call GA_Release_Update (lg_B, iLo, iHi, jLo, jHi)
      end if

      !! orthogonal -> non-orthogonal
      !! Finalize Eq. (62)
      Call GA_DGEMM ('N','N',NAS,NAS,NAS,                               &
     &               One,lg_Vec,lg_B,                                   &
     &               Zero,lg_Lag)
      Call GA_DGEMM ('N','T',NAS,NAS,NAS,                               &
     &               One,lg_Lag,lg_Vec,                                 &
     &               One,lg_SDER)

      CALL PSBMAT_FREEMEM(lg_Vec)
      bSTAT = GA_Destroy (lg_Lag)
      CALL PSBMAT_FREEMEM (lg_B)

#include "macros.fh"
      unused_var(bStat)

      End Subroutine LinDepLag_MPP

#elif defined (NAGFOR)
      ! Some compilers do not like empty files
      subroutine empty_LinDepLag_MPP()
      end subroutine empty_LinDepLag_MPP
#endif
