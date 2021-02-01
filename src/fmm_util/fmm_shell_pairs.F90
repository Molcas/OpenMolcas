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
MODULE fmm_shell_pairs

   USE fmm_global_paras, ONLY: INTK, REALK, LUPRI, fmm_sh_pairs, fmm_basis, Zero, One, Half
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_get_shell_pairs, fmm_free_shell_pairs

   TYPE(fmm_sh_pairs), ALLOCATABLE, TARGET :: sh_pairs(:)

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_get_shell_pairs(basis,sh_pairs_ptr)

      TYPE(fmm_basis), INTENT(IN) :: basis
      TYPE(fmm_sh_pairs), POINTER :: sh_pairs_ptr(:)
      INTEGER(INTK) :: n_pairs

      IF (ALLOCATED(sh_pairs)) THEN
         sh_pairs_ptr => sh_pairs
      ELSE
         ! make list of non-vanishing shell pairs
         CALL fmm_make_shell_pairs(basis,n_pairs)  ! first get n_pairs
         ALLOCATE(sh_pairs(n_pairs))
         CALL fmm_make_shell_pairs(basis,n_pairs)  ! now store list
         sh_pairs_ptr => sh_pairs
!FIXME add to fmm_stats
write(LUPRI,*) 'Number of shell pairs =', SIZE(sh_pairs)
      END IF

   END SUBROUTINE fmm_get_shell_pairs

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_free_shell_pairs

      IF (ALLOCATED(sh_pairs)) DEALLOCATE(sh_pairs)

   END SUBROUTINE fmm_free_shell_pairs

!-------------------------------------------------------------------------------

   FUNCTION fmm_extent(ExpPI)

!      USE fmm_md4_globals, ONLY: fmm_shell_pair_ThrFac, fmm_X0
!      USE fmm_md4_globals, ONLY: fmm_grain_inv, fmm_extent_min
!      USE fmm_box_utils,   ONLY: fmm_branch

      IMPLICIT NONE

      REAL(REALK), INTENT(IN) :: ExpPI

      REAL(REALK), PARAMETER :: root3 = 1.7320508075688772d0
!      REAL(REALK)   :: tmp1, tmp2
      REAL(REALK)   :: fmm_extent
!      INTEGER(INTK) :: bra_min

!      REAL(REALK), PARAMETER :: erfc_inv = 3.4589  ! 1e-6
!      REAL(REALK), PARAMETER :: erfc_inv = 4.0522  ! 1e-8
!      REAL(REALK), PARAMETER :: erfc_inv = 4.5728  ! 1e-10
!      REAL(REALK), PARAMETER :: erfc_inv = 5.6739  ! 1e-15

!      fmm_extent = fmm_X0*ExpPI

      ! Extent based on Classical overlap
!      tmp1 = SQRT(ExpPI)*erfc_inv

!      ! Artificial extent to keep all near-field exact integrals
!      bra_min = fmm_branch(fmm_extent_min,fmm_grain_inv)
!      tmp2 = half*root3*(bra_min+1)/fmm_grain_inv

!      tmp2 = zero
!      fmm_extent = MAX(tmp1,tmp2)

      fmm_extent = SQRT(ExpPI*LOG(1d10))

   END FUNCTION fmm_extent

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_make_shell_pairs(basis,n_pairs)

!      USE fmm_md4_globals, ONLY: fmm_shell_pair_ThrFac

      IMPLICIT NONE

      TYPE(fmm_basis), INTENT(IN)  :: basis
      INTEGER(INTK),   INTENT(OUT) :: n_pairs

!      REAL(REALK), PARAMETER :: ThrFac = fmm_shell_pair_ThrFac
      REAL(REALK), PARAMETER :: ThrFac = 1d-12

      INTEGER(INTK) :: II,JJ, I,J
      INTEGER(INTK) :: IPrim1, JPrim1
      INTEGER(INTK) :: IPrim2, JPrim2, IJPrim, JPTemp2
      REAL(REALK)   :: Acentr(3), Bcentr(3), ABcentr(3), P(3), PM(3), RAB(3)
      REAL(REALK)   :: ExpA, ExpB, ExpP, ExpPI, ExpAR2
      REAL(REALK)   :: R2AB, ExpKAB, tmp_ext, extent

      n_pairs = 0

      Ishel: DO II = 1, basis%nshells

         Acentr(:) = basis%Centr(:,basis%KAtom(II))
         IPrim1 = basis%KStart(II)
         IPrim2 = IPrim1 + basis%KontG(II) - 1

         Jshel: DO JJ = 1, II

            Bcentr(:) = basis%Centr(:,basis%KAtom(JJ))
            JPrim1 = basis%KStart(JJ)
            JPrim2 = JPrim1 + basis%KontG(JJ) - 1

            ABcentr(:) = Half * (Acentr(:) + Bcentr(:))
            RAB(:) = Acentr(:) - Bcentr(:)
            R2AB = DOT_PRODUCT(RAB,RAB)

            extent = zero
            IJPrim = 0
            DO I = IPrim1, IPrim2
               ExpA = basis%Expnt(I)
               ExpAR2 = ExpA * R2AB
               JPTemp2 = JPrim2
               IF (II == JJ) JPTemp2 = I
               DO J = JPrim1, JPTemp2
                  ExpB = basis%Expnt(J)
                  ExpP = ExpA + ExpB
                  ExpPI = One / ExpP
                  ExpKAB = - ExpAR2 * ExpB * ExpPI
                  IF (ExpKAB >= ThrFac) THEN
                     IJPrim = IJPrim +1
                     P(:) = (ExpA * Acentr(:) + ExpB * Bcentr(:)) * ExpPI
                     PM(:) = P(:) - ABcentr(:)
                     tmp_ext = fmm_extent(ExpPI)
                     tmp_ext = tmp_ext + SQRT(DOT_PRODUCT(PM,PM))
                     extent = MAX( extent, tmp_ext )
                  END IF
               END DO
            END DO

            IF (IJPrim > 0) THEN
               n_pairs = n_pairs + 1
               IF (ALLOCATED(sh_pairs)) THEN
                  IF (n_pairs > SIZE(sh_pairs)) CALL fmm_quit('get_sh_pairs')
                  sh_pairs(n_pairs)%I = II
                  sh_pairs(n_pairs)%J = JJ
                  sh_pairs(n_pairs)%extent = extent
                  sh_pairs(n_pairs)%centre(:) = ABcentr(:)
               END IF
            END IF

         END DO Jshel

      END DO Ishel

   END SUBROUTINE fmm_make_shell_pairs

!-------------------------------------------------------------------------------

END MODULE fmm_shell_pairs
