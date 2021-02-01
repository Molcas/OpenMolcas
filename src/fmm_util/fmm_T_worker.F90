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
MODULE fmm_T_worker

   USE fmm_global_paras, ONLY: INTK, REALK, LUPRI, ZERO_VECT_TOL, Zero, One, Two, Half

   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_get_SPLTSQ_T_matrix,   &
             fmm_get_FLTSQ_T_matrix,    &
             fmm_get_boundary_T_matrix, &
             fmm_contract_Tq,           &
             fmm_postfac_Vff,           &
             fmm_scale_vec

CONTAINS

!------------------------------------------------------------------------------
! Get sparse (l+j <= LMAX), lower-triangular, square-based T-matrix
!FIXME: would be nice to allow possibility of building rectangular-based too

   SUBROUTINE fmm_get_SPLTSQ_T_matrix(LMAX,r_ab,T_matrix)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN)  :: LMAX
      REAL(REALK),   INTENT(IN)  :: r_ab(3)
      REAL(REALK),   INTENT(OUT) :: T_matrix(:,:)

      REAL(REALK) :: I_sh((1+LMAX)**2)

      CALL fmm_generate_I(LMAX,r_ab,I_sh)
      CALL fmm_generate_T(LMAX,fmm_SPLTSQ_JMAX,I_sh,T_matrix)

   END SUBROUTINE fmm_get_SPLTSQ_T_matrix

   FUNCTION fmm_SPLTSQ_JMAX(L,LMAX)
      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN) :: L, LMAX
      INTEGER(INTK) :: fmm_SPLTSQ_JMAX
      fmm_SPLTSQ_JMAX = LMAX-L
   END FUNCTION fmm_SPLTSQ_JMAX

!------------------------------------------------------------------------------
! Get full (all elements built), lower-triangular, square T-matrix

   SUBROUTINE fmm_get_FLTSQ_T_matrix(LMAX,r_ab,T_matrix)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN)  :: LMAX
      REAL(REALK),   INTENT(IN)  :: r_ab(3)
      REAL(REALK),   INTENT(OUT) :: T_matrix(:,:)

      REAL(REALK) :: I_sh((1+2*LMAX)**2)

      CALL fmm_generate_I(2_INTK*LMAX,r_ab,I_sh)
      CALL fmm_generate_T(LMAX,fmm_FLTSQ_JMAX,I_sh,T_matrix)

   END SUBROUTINE fmm_get_FLTSQ_T_matrix

   FUNCTION fmm_FLTSQ_JMAX(l,LMAX)
      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN) :: l, LMAX
      INTEGER(INTK) :: fmm_FLTSQ_JMAX
      fmm_FLTSQ_JMAX = LMAX
! Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(l)
   END FUNCTION fmm_FLTSQ_JMAX

!------------------------------------------------------------------------------
! Get single-column T-matrix for building boundary potential only

   SUBROUTINE fmm_get_boundary_T_matrix(LMAX,r_ab,T_matrix)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN)  :: LMAX
      REAL(REALK),   INTENT(IN)  :: r_ab(3)
      REAL(REALK),   INTENT(OUT) :: T_matrix(:,:)

      REAL(REALK) :: I_sh((1+2*LMAX)**2)

      CALL fmm_generate_I(2_INTK*LMAX,r_ab,I_sh)
      CALL fmm_generate_boundary_T(LMAX,I_sh,T_matrix)

   END SUBROUTINE fmm_get_boundary_T_matrix

!------------------------------------------------------------------------------

   SUBROUTINE fmm_generate_I(LMAX,vector,I_sh)

      ! Subroutine to generate Scaled Irregular solid harmonics
      ! See page 416, "Molecular Electronic Structure Theory".
      !--------------------------------------------------------
      ! I_sh(M,L)  == COS component of I_sh(M,L)
      ! I_sh(-M,L) == SIN component of I_sh(-M,L)
      ! where L and M are positive
      ! Note that I_sh(0,L) == COS component of I_sh(0,L) !
      ! In fact we combine (l,m) into a single index and
      ! we assume the order 0, -1,1,+1, -2,-1,0,+1,+2 ...
      !--------------------------------------------------------

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN)  :: LMAX
      REAL(REALK),   INTENT(IN)  :: vector(3)
      REAL(REALK),   INTENT(OUT) :: I_sh((LMAX+1)**2)

      REAL(REALK)   :: tmp1, tmp2, tmp3
      REAL(REALK)   :: x,y,z, r_2, r_minus2
      INTEGER(INTK) :: L,m, i,j,k,p,q,u, sign_L

      x = vector(1)
      y = vector(2)
      z = vector(3)
      r_2 = x**2 + y**2 + z**2
      r_minus2 = one/(r_2)

!FIXME;
      IF (r_2 < ZERO_VECT_TOL) THEN
         WRITE(LUPRI,'(3E25.15)') vector
         CALL fmm_quit('Why do we try to do a zero T_vector?')
      END IF

      IF (LMAX==0) THEN
         I_sh(1) = SQRT(r_minus2)
         RETURN
      END IF

      ! Start with initialization of lowest order terms
        I_sh(1) =  SQRT(r_minus2)
        I_sh(2) = -(y)*(r_minus2)*I_sh(1)
        I_sh(3) =  (z)*(r_minus2)*I_sh(1)
        I_sh(4) = -(x)*(r_minus2)*I_sh(1)

      ! Now iterate higher order terms
      sign_L = -1
      DO L = 2, LMAX
         sign_L = -sign_L
         i = (L+1)*(L+1)
         j = L*L+1
         p = j-2*L+1
         q = j-1
         tmp1 = (2*L-1)*r_minus2
         tmp2 = tmp1*y*sign_L
         tmp3 = tmp1*x
         !
         I_sh(i) = tmp2*I_sh(p) - tmp3*I_sh(q)
         I_sh(j) = tmp2*I_sh(q) + tmp3*I_sh(p)
         !
         p = p+L-1
         q = L*(L-3)+3
         u = q+L-2
         k = i-L
         tmp2 = tmp1*z
         DO m = 0, (L-2)
            tmp3 = (u-m*m)*(r_minus2)
            !
            I_sh(k+m) = tmp2*I_sh(p+m) - tmp3*I_sh(q+m)
            I_sh(k-m) = tmp2*I_sh(p-m) - tmp3*I_sh(q-m)
            ! I_sh(L,0)=0
            !
         END DO
         ! Now do (L, L-1) terms to avoid calling elements that do not exist
         m = L-1
         I_sh(k+m) = tmp2*I_sh(p+m)
         I_sh(k-m) = tmp2*I_sh(p-m)
      END DO

   END SUBROUTINE fmm_generate_I

!------------------------------------------------------------------------------

   SUBROUTINE fmm_generate_T(LMAX,JMAX,I_sh,T_matrix)

   ! Subroutine to generate Real Interaction Matrix for
   ! given point from Irregular Solid Harmonics (I_sh)
   ! cf Helgaker et. al. pp 415
   !-----------------------------------------------------------------------
   !  T_cos_cos(lm,jk) = 2*{I_cos(l+j,m+k) + (-1)^k*I_cos(l+j,m-k)}  etc.
   !
   ! Note we only write the non-zero elements here, and the precise form
   ! depends on the function JMAX passed in
   !-----------------------------------------------------------------------
   !FIXME: can optimise for j=l by noting that if we build just lower
   ! triangular, we need only half the elements there

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN)  :: LMAX
      REAL(REALK),   INTENT(IN)  :: I_sh(:)
      REAL(REALK),   INTENT(OUT) :: T_matrix(:,:)
      INTEGER(INTK), EXTERNAL    :: JMAX

      INTEGER(INTK) :: sign_k, sign_m
      INTEGER(INTK) :: L,J,k,m, q, pp,qq, u

      IF (LMAX==0) THEN
         T_matrix(1,1) = two*I_sh(1)     ! NB scaling to make T symmetric
         RETURN
      END IF

      L_loop: DO L = 0, LMAX
         qq = L*(L+1) +1
         J_loop: DO J = L, JMAX(L,LMAX)
            pp = J*(J+1) +1
            u = 1+(L+J)*(L+J+1)

            sign_m = -1
            positive_m: DO m = 0, L
               q = qq+m
               sign_m = -sign_m

               ! cos-cos terms
               !---------------

               sign_k = -1
               DO k = 0, MIN(m,J)
                  sign_k = -sign_k
                  ! (m+k)>0 and (m-k)>0
                  T_matrix(pp+k,q) = I_sh(u+(m+k)) + sign_k*I_sh(u+(m-k))
               END DO

               DO k = (m+1), J
                  ! (m+k)>0 and (m-k)<0
                  T_matrix(pp+k,q) = I_sh(u+(m+k)) + sign_m*I_sh(u+(-(m-k)))
               END DO

               ! cos-sin terms
               !---------------

               ! only need to consider (m+k)=0 separately as cannot
               ! also have (m-k)=0 unless k=0

               DO k = (-J), MIN((-m-1),-1)
                  ! (m+k)<0 and (m-k)>0
                  T_matrix(pp+k,q) = I_sh(u+(m+k)) + sign_m*I_sh(u+(-(m-k)))
               END DO

               DO k = MAX(-J,-m), MIN(-m,-1)
                  ! (m+k)=0 and (m-k)>0
                  T_matrix(pp+k,q) = sign_m*I_sh(u+(-(m-k)))
               END DO

               sign_k = 1
               DO k = -1, MAX((-m+1),(-J)), -1
                  sign_k = -sign_k
                  ! (m+k)>0 and (m-k)>0
                  T_matrix(pp+k,q) = -sign_m*sign_k*I_sh(u+(-(m+k)))         &
                                     +sign_m*I_sh(u+(-(m-k)))
               END DO

            END DO  positive_m

            sign_m = 1
            negative_m: DO m = -1, -L, -1
               q = qq+m
               sign_m = -sign_m

               ! sin-cos terms
               !---------------

               k = MIN(J,(-m))
               sign_k = 1
               IF (BTEST(k,0)) sign_k = -1
               ! (m+k)=0 and (m-k)<0
               T_matrix(pp+k,q) = sign_k*I_sh(u+(m-k))

               sign_k = -1
               DO k = 0, MIN((-m-1),J)
                  sign_k = -sign_k
                  ! (m+k)<0 and (m-k)<0
                 T_matrix(pp+k,q) = I_sh(u+(m+k)) + sign_k*I_sh(u+(m-k))
               END DO

               sign_k = -1
               IF (BTEST(J,0)) sign_k = 1
               DO k = J, (-m+1), -1
                  sign_k = -sign_k
                  ! (m+k)>0 and (m-k)<0
                  T_matrix(pp+k,q) = -sign_m*sign_k*I_sh(u+(-(m+k)))         &
                                     +sign_k*I_sh(u+(m-k))
               END DO

               ! sin-sin terms
               !---------------

               sign_k = -1
               IF (BTEST(J,0)) sign_k = 1
               DO k = (-J), m
                  sign_k = -sign_k
                  ! (m+k)<0 and (m-k)>0
                  T_matrix(pp+k,q) = -sign_m*sign_k*I_sh(u+(-(m+k)))         &
                                     +sign_k*I_sh(u+(m-k))
               END DO

               sign_k = 1
               DO k = -1, MAX((-J),(m+1)), -1
                  sign_k = -sign_k
                  ! (m+k)<0 and (m-k)<0
                  T_matrix(pp+k,q) = -sign_m*sign_k*I_sh(u+(-(m+k)))         &
                                     +sign_m*I_sh(u+(-(m-k)))
               END DO

            END DO  negative_m
         END DO     J_loop
      END DO        L_loop

   END SUBROUTINE fmm_generate_T

!------------------------------------------------------------------------------

   SUBROUTINE fmm_generate_boundary_T(LMAX,I_sh,T_matrix)

!FIXME comments!
! This is a simplified version of fmm_generate_T_matrix when building
! only the boundary potential, when only one column of the T_matrix
! is required.

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN)  :: LMAX
      REAL(REALK),   INTENT(IN)  :: I_sh(:)
      REAL(REALK),   INTENT(OUT) :: T_matrix(:,:)

      INTEGER(INTK) :: J,k, pp, u

      IF (LMAX==0) THEN
         T_matrix(1,1) = two*I_sh(1)     ! NB scaling to make T symmetric
         RETURN
      END IF

      DO J = 0, LMAX
         pp = J*(J+1) +1
         u = 1+(J)*(J+1)
         DO k = -J, J
            T_matrix(pp+k,1) = two*I_sh(u+k)
         END DO
      END DO

   END SUBROUTINE fmm_generate_boundary_T

!------------------------------------------------------------------------------

   FUNCTION fmm_contract_Tq(LMAX,vect,T_matrix)

     IMPLICIT NONE
     INTEGER(INTK), INTENT(IN) :: LMAX
     REAL(REALK),   INTENT(IN) :: vect(:)
     REAL(REALK),   INTENT(IN) :: T_matrix(:,:)
     REAL(REALK) :: fmm_contract_Tq(SIZE(vect))

     INTEGER(INTK) :: L, p,q,s, u, qmin,qmax

     ! do L=0 terms first
     p = (LMAX+1)*(LMAX+1)
     fmm_contract_Tq(1) = half*DOT_PRODUCT(vect(1:p),T_matrix(1:p,1))
     DO s = 2, p
        fmm_contract_Tq(s) = vect(1)*T_matrix(s,1)
     END DO

     contract: DO L = 1, LMAX

        u = L*(L+1) +1
        p = (LMAX-L+1)*(LMAX-L+1)

        qmin = u-L
        qmax = MIN( u+L, p )

        qloop: DO q = qmin, qmax
           fmm_contract_Tq(q) = fmm_contract_Tq(q)           &
                               + DOT_PRODUCT(vect(q:p),T_matrix(q:p,q))
           DO s = q+1, p
              fmm_contract_Tq(s) = fmm_contract_Tq(s) + T_matrix(s,q)*vect(q)
           END DO
        END DO qloop

        ! add on extra post-factors "after" contraction
        fmm_contract_Tq(u) = half*fmm_contract_Tq(u)    ! m=0

     END DO contract

   END FUNCTION fmm_contract_Tq

!------------------------------------------------------------------------------

   SUBROUTINE fmm_scale_vec(LMAX,scl_fact,scale_vec,prefactor)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN)  :: LMAX
      REAL(REALK),   INTENT(IN)  :: scl_fact
      REAL(REALK),   INTENT(OUT) :: scale_vec(:), prefactor

      INTEGER(INTK) :: L,L2, lo,hi, u
      REAL(REALK)   :: prefactor1, prevec

      prefactor1 = one/scl_fact

      prevec       = one
      scale_vec(1) = prevec
      DO L = 1, LMAX
         prevec = prefactor1*prevec
         L2 = L*L
         lo = L2+1
         hi = L2+2*L+1
         DO u = lo, hi
            scale_vec(u) = prevec
         END DO
      END DO

      IF (scl_fact < zero) THEN
         prefactor = -prefactor1
      ELSE
         prefactor =  prefactor1
      END IF

   END SUBROUTINE fmm_scale_vec

!------------------------------------------------------------------------------

   SUBROUTINE fmm_postfac_Vff(LMAX,Vff_tmp)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN)    :: LMAX
      REAL(REALK),   INTENT(INOUT) :: Vff_tmp(:)

      INTEGER(INTK) :: l, u
      ! add extra postfactors  (for half**(delta(0,m))
      DO l = 0, LMAX
         u = l*(l+1) +1
         Vff_tmp(u) = half*Vff_tmp(u)
      END DO

   END SUBROUTINE fmm_postfac_Vff

!------------------------------------------------------------------------------

END MODULE fmm_T_worker
