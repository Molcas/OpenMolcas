!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
MODULE fmm_J_builder

   USE fmm_global_paras
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_get_E_from_Vff,       &
             fmm_get_E_from_pkd_Vff,   &
             fmm_get_J_from_Vff,       &
             fmm_get_J_from_pkd_Vff

CONTAINS

!-------------------------------------------------------------------------------
! Trivial check that the number of moments matches the number of potentials

   SUBROUTINE fmm_verify_Vff_input(scheme,LHS_mms,Vff,J_or_E)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN) :: scheme
      TYPE(raw_mm_data),  INTENT(IN) :: LHS_mms
      REAL(REALK),        INTENT(IN) :: Vff(:,:)
      CHARACTER(LEN=1),   INTENT(IN) :: J_or_E
      LOGICAL :: A,B,C

      A = (SIZE(LHS_mms%paras) /= SIZE(Vff,2))
      IF (A) CALL fmm_quit('incompatible SIZE of Vff and LHS moments!')

      IF (J_or_E == 'J') THEN
         A = (scheme%LHS_mm_range == NUCLEAR_ONLY)
         B = (scheme%RHS_mm_range == NUCLEAR_ONLY)
         C = (scheme%LHS_mm_range == ALL_MOMENTS )
         IF ( (A.AND.B) .OR. C ) CALL fmm_quit('mm_ranges invalid')
      END IF

   END SUBROUTINE fmm_verify_Vff_input

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_get_E_with_text(scheme,energy,text)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)    :: scheme
      REAL(REALK),        INTENT(INOUT) :: energy
      CHARACTER(LEN=*),   INTENT(OUT)   :: text
      LOGICAL :: A,B,C,D

      A = (scheme%LHS_mm_range == ELECTRONIC_ONLY)
      B = (scheme%RHS_mm_range == ELECTRONIC_ONLY)
      C = (scheme%LHS_mm_range == NUCLEAR_ONLY)
      D = (scheme%RHS_mm_range == NUCLEAR_ONLY)

      IF (scheme%LHS_mm_range == scheme%RHS_mm_range) THEN
         energy = half*energy
         text = "total classical Coulomb energy"
         IF (A) text = "classical Coulomb electronic energy"
         IF (C) text = "classical Coulomb nuclear repulsion"
      ELSE IF (A.OR.B) THEN
         IF (C.OR.D) THEN
            text = "classical Coulomb nuclear attraction"
         ELSE
            text = "e-n + 2*(e-e) energy"
         END IF
      ELSE
         ! range is different and neither is ELECTRONIC_ONLY
         text = "e-n + 2*(n-n) energy"
      END IF

   END SUBROUTINE fmm_get_E_with_text

!-------------------------------------------------------------------------------
! Get energy from full contraction of LHS moments and the potentials
!-------------------------------------------------------------------
! We assume Vff(lm,i) is the potential at i due to ALL the RHS moments
! (so that a simple contraction with the LHS moments will double count
!  interactions if the LHS and RHS moment ranges overlap).

   SUBROUTINE fmm_get_E_from_Vff(scheme,LHS_mms,Vff,energy,text)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)  :: scheme
      TYPE(raw_mm_data),  INTENT(IN)  :: LHS_mms
      REAL(REALK),        INTENT(IN)  :: Vff(:,:)
      REAL(REALK),        INTENT(OUT) :: energy
      CHARACTER(LEN=*),   INTENT(OUT) :: text

      REAL(REALK)   :: g
      INTEGER(INTK) :: u,v, lm_max

      CALL fmm_verify_Vff_input(scheme,LHS_mms,Vff,'E')

      ! although Vff should be the same size, we test for generality
      IF (SIZE(LHS_mms%qlm_T,1) /= SIZE(Vff,1)) STOP 'mm_get_E_from_Vff:2'
      lm_max = MIN(SIZE(LHS_mms%qlm_T,1),SIZE(Vff,1))
      DO u = 1, SIZE(LHS_mms%paras)
         v = LHS_mms%paras(u)%id
         g = DOT_PRODUCT(LHS_mms%qlm_T(:lm_max,v),Vff(:lm_max,v))
         energy = energy + g
      END DO

      CALL fmm_get_E_with_text(scheme,energy,text)

   END SUBROUTINE fmm_get_E_from_Vff

!-------------------------------------------------------------------------------
! Get energy from full contraction of LHS moments and the potentials
! on the basis that the LHS parameters and Vff may have been packed
! into batches (and must be "expanded").
!-------------------------------------------------------------------
! We assume Vff(lm,i) is the potential at i due to ALL the RHS moments
! (so that a simple contraction with the LHS moments will double count
!  interactions if the LHS and RHS moment ranges overlap).

   SUBROUTINE fmm_get_E_from_pkd_Vff(scheme,LHS_mms,Vff,energy,text)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)  :: scheme
      TYPE(raw_mm_data),  INTENT(IN)  :: LHS_mms
      REAL(REALK),        INTENT(IN)  :: Vff(:,:)
      REAL(REALK),        INTENT(OUT) :: energy
      CHARACTER(LEN=*),   INTENT(OUT) :: text

      REAL(REALK)   :: g
      INTEGER(INTK) :: u,v,w, lm_max
      TYPE(id_node), POINTER :: batch_map

      CALL fmm_verify_Vff_input(scheme,LHS_mms,Vff,'E')

      lm_max = MIN(SIZE(LHS_mms%qlm_T,1),SIZE(Vff,1))
      packed_loop: DO u = 1, SIZE(LHS_mms%paras)

         v = LHS_mms%paras(u)%id  ! LHS packed (batch) parameters
         batch_map => LHS_mms%batch_map(v)%head

         batch_members2: DO ! over batch list until pointer disassociated
            w = batch_map%id   ! raw LHS moment ID
            g = DOT_PRODUCT(LHS_mms%qlm_T(:lm_max,w),Vff(:lm_max,v))
            energy = energy + g
             ! only do next raw item in batch list if it exists
            IF (.NOT.ASSOCIATED(batch_map%next)) EXIT batch_members2
            batch_map => batch_map%next
         END DO batch_members2

      END DO packed_loop

      CALL fmm_get_E_with_text(scheme,energy,text)

   END SUBROUTINE fmm_get_E_from_pkd_Vff

!-------------------------------------------------------------------------------
! Build J_matrix components from contracion of the LHS moments and potentials
!----------------------------------------------------------------------------
! We assume Vff(lm,i) is the potential at i due to ALL the RHS moments
! This choice of Vff is fine here since the e-e interactions require a
! factor of 2 in the J-matrix, but NOT the e-n interactions.

   SUBROUTINE fmm_get_J_from_Vff(scheme,LHS_mms,Vff,J_matrix)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)  :: scheme
      TYPE(raw_mm_data),  INTENT(IN)  :: LHS_mms
      REAL(REALK),        INTENT(IN)  :: Vff(:,:)
      REAL(REALK),        INTENT(OUT) :: J_matrix(:,:)

      REAL(REALK)   :: g
      INTEGER(INTK) :: u,v, i,j, lm_max

      CALL fmm_verify_Vff_input(scheme,LHS_mms,Vff,'J')

      ! although Vff should be the same size, we test for generality
      IF (SIZE(LHS_mms%qlm_T,1) /= SIZE(Vff,1)) STOP 'mm_get_J_from_Vff:2'
      lm_max = MIN(SIZE(LHS_mms%qlm_T,1),SIZE(Vff,1))
      DO u = 1, SIZE(LHS_mms%paras)
         v = LHS_mms%paras(u)%id
         g = DOT_PRODUCT(LHS_mms%qlm_T(:lm_max,v),Vff(:lm_max,v))
         i = LHS_mms%J_indices(v)%i_indx
         j = LHS_mms%J_indices(v)%j_indx
         J_matrix(i,j) = J_matrix(i,j) + g
         IF (i/=j) J_matrix(j,i) = J_matrix(j,i) + g
      END DO

   END SUBROUTINE fmm_get_J_from_Vff

!-------------------------------------------------------------------------------
! Build J_matrix components from contracion of the LHS moments and potentials
! This routine recognises that the LHS parameters may have been packed
! and thus expands the potential over LHS batches.
!----------------------------------------------------------------------------
! We assume Vff(lm,i) is the potential at i due to ALL the RHS moments

   SUBROUTINE fmm_get_J_from_pkd_Vff(scheme,LHS_mms,Vff,J_matrix)

      IMPLICIT NONE
      TYPE(scheme_paras), INTENT(IN)  :: scheme
      TYPE(raw_mm_data),  INTENT(IN)  :: LHS_mms
      REAL(REALK),        INTENT(IN)  :: Vff(:,:)
      REAL(REALK),        INTENT(OUT) :: J_matrix(:,:)

      REAL(REALK)   :: g
      INTEGER(INTK) :: u,v,w, i,j, lm_max
      TYPE(id_node), POINTER :: batch_map

      CALL fmm_verify_Vff_input(scheme,LHS_mms,Vff,'J')

      ! Vff should now be the same size as the raw LHS parameters;
      ! However, the size of LHS %qlm_T will be larger;
      ! We need to "expand" Vff using the batch mapping.

      lm_max = MIN(SIZE(LHS_mms%qlm_T,1),SIZE(Vff,1))
      packed_loop: DO u = 1, SIZE(LHS_mms%paras)

         v = LHS_mms%paras(u)%id  ! LHS packed (batch) moments
         batch_map => LHS_mms%batch_map(v)%head

         batch_members: DO ! over batch list until pointer disassociated
            w = batch_map%id   ! raw LHS moment ID
            g = DOT_PRODUCT(LHS_mms%qlm_T(:lm_max,w),Vff(:lm_max,v))
            i = LHS_mms%J_indices(w)%i_indx
            j = LHS_mms%J_indices(w)%j_indx
            J_matrix(i,j) = J_matrix(i,j) + g
!            IF (i/=j) J_matrix(j,i) = J_matrix(j,i) + g

             ! only do next raw item in batch list if it exists
            IF (.NOT.ASSOCIATED(batch_map%next)) EXIT batch_members
            batch_map => batch_map%next
         END DO batch_members
      END DO packed_loop

   END SUBROUTINE fmm_get_J_from_pkd_Vff

!-------------------------------------------------------------------------------

END MODULE fmm_J_builder

