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
! Copyright (C) Steven Vancoillie                                      *
!***********************************************************************
!SVC: compute RHS elements "on demand". If we have access to all the
! Cholesky vectors, we can just instruct a process to compute its own
! block of RHS elements, computing the integrals directly. This is much
! more computationally intensive, but should scale much better since we
! go from a badly scaling scatter algorithm to no communication at all.
! This also eliminates the need for the GA library in creating the RHS.
! FIXME: optimizations needed, remove double computation of integrals

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE RHSOD(IVEC)
      use definitions, only: iwp
#ifdef _DEBUGPRINT_
      use definitions, only: wp
      use caspt2_module, only: NSYM, NASUP, NISUP
#endif
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: VERBOSE
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      IMPLICIT None
      integer(kind=iwp), intent(in):: iVec
#ifdef _DEBUGPRINT_
      real(kind=wp) DNRM2
      real(kind=wp), external:: RHS_DDot
      integer(kind=iwp) ICASE, ISYM, NAS, NIS, lg_W
#endif

      IF (IPRGLB>=VERBOSE) THEN
        WRITE(6,'(1X,A)') ' Using RHS on-demand algorithm'
      END IF

#ifdef _MOLCAS_MPP_
      IF (.NOT.Is_Real_Par()) THEN
        WRITE(6,'(1X,A)') 'RHSOD: error: fake parallel not supported'
        CALL AbEnd()
      END IF
#endif

      CALL RHSOD_A(IVEC)
      CALL RHSOD_B(IVEC)
      CALL RHSOD_C(IVEC)
      CALL RHSOD_D(IVEC)
      CALL RHSOD_E(IVEC)
      CALL RHSOD_F(IVEC)
      CALL RHSOD_G(IVEC)
      CALL RHSOD_H(IVEC)

#ifdef _DEBUGPRINT_
! compute and print RHS fingerprints
      WRITE(6,'(1X,A4,1X,A3,1X,A18)') 'Case','Sym','Fingerprint'
      WRITE(6,'(1X,A4,1X,A3,1X,A18)') '====','===','==========='
      DO ICASE=1,13
        DO ISYM=1,NSYM
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF (NAS*NIS/=0) THEN
            CALL RHS_ALLO (NAS,NIS,lg_W)
            CALL RHS_READ (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
            DNRM2 = RHS_DDOT(NAS,NIS,lg_W,lg_W)
            WRITE(6,'(1X,I4,1X,I3,1X,F18.11)') ICASE,ISYM,DNRM2
          END IF
        END DO
      END DO
#endif

      END SUBROUTINE RHSOD
