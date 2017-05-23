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
MODULE fmm_multiple_T_worker

   USE fmm_global_paras

   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_get_SPLTSQ_T_matrices,   &
             fmm_get_FLTSQ_T_matrices,    &
             fmm_contract_multi_Tq

CONTAINS

!------------------------------------------------------------------------------
! Get sparse (l+j <= LMAX), lower-triangular, square-based T-matrix
!FIXME: would be nice to allow possibility of building rectangular-based too

   SUBROUTINE fmm_get_SPLTSQ_T_matrices(ndim,LMAX,r_ab,T_matrix)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN)  :: ndim, LMAX
      REAL(REALK),   INTENT(IN)  :: r_ab(:,:)
      REAL(REALK),   INTENT(OUT) :: T_matrix(:,:,:)

      REAL(REALK) :: I_sh(ndim,(1+LMAX)**2)

      CALL fmm_generate_I(ndim,LMAX,r_ab,I_sh)
      CALL fmm_generate_T(ndim,LMAX,.FALSE.,I_sh,T_matrix)

   END SUBROUTINE fmm_get_SPLTSQ_T_matrices

!------------------------------------------------------------------------------
! Get full (all elements built), lower-triangular, square T-matrix

   SUBROUTINE fmm_get_FLTSQ_T_matrices(ndim,LMAX,r_ab,T_matrix)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN)  :: ndim, LMAX
      REAL(REALK),   INTENT(IN)  :: r_ab(:,:)
      REAL(REALK),   INTENT(OUT) :: T_matrix(:,:,:)

      REAL(REALK) :: I_sh(ndim,(1+2*LMAX)**2)

      CALL fmm_generate_I(ndim,2_INTK*LMAX,r_ab,I_sh)
      CALL fmm_generate_T(ndim,LMAX,.TRUE.,I_sh,T_matrix)

   END SUBROUTINE fmm_get_FLTSQ_T_matrices

!------------------------------------------------------------------------------

   SUBROUTINE fmm_generate_I(ndim,LMAX,r_ab,I_sh)

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
      INTEGER(INTK), INTENT(IN)  :: LMAX, NDIM
      REAL(REALK),   INTENT(IN)  :: r_ab(:,:)
      REAL(REALK),   INTENT(OUT) :: I_sh(:,:)

      REAL(REALK)   :: r_minus2(NDIM)
      INTEGER(INTK) :: L,m, i,j,k,p,q,u, sign_L , n, p1,q1, m1
      REAL(REALK)   :: var1, var2, var3, rmin2

      r_minus2(:) = one/(r_ab(:,1)**2 + r_ab(:,2)**2 + r_ab(:,3)**2)

      IF (LMAX==0) THEN
         I_sh(:,1) = SQRT(r_minus2(:))
         RETURN
      END IF

      I_sh(:,1) =  SQRT(r_minus2(:))
      I_sh(:,2) = -r_ab(:,2)*r_minus2(:)*I_sh(:,1)
      I_sh(:,3) =  r_ab(:,3)*r_minus2(:)*I_sh(:,1)
      I_sh(:,4) = -r_ab(:,1)*r_minus2(:)*I_sh(:,1)

      ! Now iterate higher order terms
      sign_L = -1
      DO L = 2, LMAX

         sign_L = -sign_L
         i = (L+1)*(L+1)
         j = L*L+1
         p = j-2*L+1
         q = j-1
         p1 = p
         q1 = q
         p = p+L-1
         q = L*(L-3)+3
         u = q+L-2
         k = i-L
         m1 = L-1

         DO n = 1, ndim

            rmin2 = r_minus2(n)
            var1 = (2*L-1)*rmin2
            var2 = var1*r_ab(n,2)*sign_L
            var3 = var1*r_ab(n,1)
            I_sh(n,i) = var2*I_sh(n,p1) - var3*I_sh(n,q1)
            I_sh(n,j) = var2*I_sh(n,q1) + var3*I_sh(n,p1)
            var2 = var1*r_ab(n,3)

            DO m = 0, L-2
               var3 = (u-m*m)*rmin2
               I_sh(n,k+m) = var2*I_sh(n,p+m) - var3*I_sh(n,q+m)
               I_sh(n,k-m) = var2*I_sh(n,p-m) - var3*I_sh(n,q-m)
            END DO

            I_sh(n,k+m1) = var2*I_sh(n,p+m1)
            I_sh(n,k-m1) = var2*I_sh(n,p-m1)

         end do

      END DO

   END SUBROUTINE fmm_generate_I

!------------------------------------------------------------------------------

   SUBROUTINE fmm_generate_T(ndim,LMAX,TOLMAX,I_sh,T_matrix)

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
      INTEGER(INTK), INTENT(IN)  :: LMAX, NDIM
      LOGICAL,       INTENT(IN)  :: TOLMAX
      REAL(REALK),   INTENT(IN)  :: I_sh(:,:)
      REAL(REALK),   INTENT(OUT) :: T_matrix(:,:,:)

      INTEGER(INTK) :: sign_m
      INTEGER(INTK) :: L,J,k,m, q, pp,qq, u, Jlast

      IF (LMAX==0) THEN
         T_matrix(:,1,1) = two*I_sh(:,1)  ! NB scaling to make T symmetric
         RETURN
      END IF

      L_loop: DO L = 0, LMAX
         qq = L*(L+1) +1

         IF (tolmax) THEN
            Jlast = LMAX
         ELSE
            Jlast = LMAX - L
         END IF

         J_loop: DO J = L, Jlast
            pp = J*(J+1) +1
            u = 1+(L+J)*(L+J+1)

            sign_m = -1
            positive_m: DO m = 0, L
               q = qq+m
               sign_m = -sign_m

               ! cos-cos terms
               !---------------

               DO k = 0, MIN(m,J)
                  IF (BTEST(k,0)) THEN
                     T_matrix(:,pp+k,q) = I_sh(:,u+m+k) - I_sh(:,u+m-k)
                  ELSE
                     T_matrix(:,pp+k,q) = I_sh(:,u+m+k) + I_sh(:,u+m-k)
                  END IF
               END DO
               IF (sign_m > 0) THEN
                  DO k = m+1, J
                     T_matrix(:,pp+k,q) = I_sh(:,u+m+k) + I_sh(:,u-m+k)
                  END DO
               ELSE
                  DO k = m+1, J
                     T_matrix(:,pp+k,q) = I_sh(:,u+m+k) - I_sh(:,u-m+k)
                  END DO
               END IF

               ! cos-sin terms
               !---------------

               IF (sign_m > 0) THEN
                  DO k = -J, MIN(-m-1,-1)
                     T_matrix(:,pp+k,q) = I_sh(:,u+(m+k)) + I_sh(:,u+(-(m-k)))
                  END DO
                  DO k = MAX(-J,-m), MIN(-m,-1)
                     T_matrix(:,pp+k,q) = I_sh(:,u+(-(m-k)))
                  END DO
                  DO k = -1, MAX(-m+1,-J), -1
                     IF (BTEST(k,0)) THEN
                        T_matrix(:,pp+k,q) =   I_sh(:,u+(-(m+k)))    &
                                             + I_sh(:,u+(-(m-k)))
                     ELSE
                        T_matrix(:,pp+k,q) = - I_sh(:,u+(-(m+k)))    &
                                             + I_sh(:,u+(-(m-k)))
                     END IF
                  END DO
               ELSE
                  DO k = -J, MIN(-m-1,-1)
                     T_matrix(:,pp+k,q) = I_sh(:,u+(m+k)) - I_sh(:,u+(-(m-k)))
                  END DO
                  DO k = MAX(-J,-m), MIN(-m,-1)
                     T_matrix(:,pp+k,q) = - I_sh(:,u+(-(m-k)))
                  END DO
                  DO k = -1, MAX(-m+1,-J), -1
                     IF (BTEST(k,0)) THEN
                        T_matrix(:,pp+k,q) = - I_sh(:,u+(-(m+k)))   &
                                             - I_sh(:,u+(-(m-k)))
                     ELSE
                        T_matrix(:,pp+k,q) =   I_sh(:,u+(-(m+k)))   &
                                             - I_sh(:,u+(-(m-k)))
                     END IF
                  END DO
               END IF

            END DO  positive_m

            sign_m = 1
            negative_m: DO m = -1, -L, -1
               q = qq+m
               sign_m = -sign_m

               ! sin-cos terms
               !---------------

               k = MIN(J,-m)
               IF (BTEST(k,0)) THEN
                  T_matrix(:,pp+k,q) = - I_sh(:,u+(m-k))
               ELSE
                  T_matrix(:,pp+k,q) =   I_sh(:,u+(m-k))
               END IF

               DO k = 0, MIN((-m-1),J)
                  IF (BTEST(k,0)) THEN
                     T_matrix(:,pp+k,q) = I_sh(:,u+(m+k)) - I_sh(:,u+(m-k))
                  ELSE
                     T_matrix(:,pp+k,q) = I_sh(:,u+(m+k)) + I_sh(:,u+(m-k))
                  END IF
               END DO

               IF (sign_m > 0) THEN
                  DO k = J, (-m+1), -1
                     IF (BTEST(k,0)) THEN
                        T_matrix(:,pp+k,q) =   I_sh(:,u+(-(m+k)))   &
                                             - I_sh(:,u+(m-k))
                     ELSE
                        T_matrix(:,pp+k,q) = - I_sh(:,u+(-(m+k)))   &
                                             + I_sh(:,u+(m-k))
                     END IF
                  END DO
               ELSE
                  DO k = J, (-m+1), -1
                     IF (BTEST(k,0)) THEN
                        T_matrix(:,pp+k,q) = - I_sh(:,u+(-(m+k)))   &
                                             - I_sh(:,u+(m-k))
                     ELSE
                        T_matrix(:,pp+k,q) =   I_sh(:,u+(-(m+k)))   &
                                             + I_sh(:,u+(m-k))
                     END IF
                  END DO
               END IF

               ! sin-sin terms
               !---------------

               IF (sign_m > 0) THEN
                  DO k = -J, m
                     IF (BTEST(k,0)) THEN
                        T_matrix(:,pp+k,q) =   I_sh(:,u+(-(m+k)))   &
                                             - I_sh(:,u+(m-k))
                     ELSE
                        T_matrix(:,pp+k,q) = - I_sh(:,u+(-(m+k)))   &
                                             + I_sh(:,u+(m-k))
                     END IF
                  END DO
                  DO k = -1, MAX(-J,m+1), -1
                     IF (BTEST(k,0)) THEN
                        T_matrix(:,pp+k,q) =   I_sh(:,u+(-(m+k)))   &
                                             + I_sh(:,u+(-(m-k)))
                     ELSE
                        T_matrix(:,pp+k,q) = - I_sh(:,u+(-(m+k)))   &
                                             + I_sh(:,u+(-(m-k)))
                     END IF
                  END DO
               ELSE
                  DO k = -J, m
                     IF (BTEST(k,0)) THEN
                        T_matrix(:,pp+k,q) = - I_sh(:,u+(-(m+k)))   &
                                             - I_sh(:,u+(m-k))
                     ELSE
                        T_matrix(:,pp+k,q) =   I_sh(:,u+(-(m+k)))   &
                                             + I_sh(:,u+(m-k))
                     END IF
                  END DO
                  DO k = -1, MAX(-J,m+1), -1
                     IF (BTEST(k,0)) THEN
                        T_matrix(:,pp+k,q) = - I_sh(:,u+(-(m+k)))   &
                                             - I_sh(:,u+(-(m-k)))
                     ELSE
                        T_matrix(:,pp+k,q) =   I_sh(:,u+(-(m+k)))   &
                                             - I_sh(:,u+(-(m-k)))
                     END IF
                  END DO
               END IF

            END DO  negative_m
         END DO     J_loop
      END DO        L_loop

! Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(ndim)
   END SUBROUTINE fmm_generate_T

!------------------------------------------------------------------------------

   FUNCTION fmm_contract_multi_Tq(LMAX,vect,T_mats,ndim)

     IMPLICIT NONE
     INTEGER(INTK), INTENT(IN) :: LMAX, ndim
     REAL(REALK),   INTENT(IN) :: vect(:)
     REAL(REALK),   INTENT(IN) :: T_mats(:,:,:)
     REAL(REALK) :: fmm_contract_multi_Tq(ndim,(LMAX+1)**2)

     INTEGER(INTK) :: L, p,q,r,s ,u, qmin,qmax
     REAL(REALK) :: fac

     ! first do L=0 terms
     p = (LMAX+1)*(LMAX+1)

     fac = half*vect(1)
     DO r = 1, ndim
        fmm_contract_multi_Tq(r,1) = fac*T_mats(r,1,1)
     END DO

     DO s = 2, p
     fac = half*vect(s)
        DO r = 1, ndim
           fmm_contract_multi_Tq(r,1) = fmm_contract_multi_Tq(r,1)    &
                                       + fac*T_mats(r,s,1)
        END DO
     END DO

     fac = vect(1)
     DO s = 2, p
        DO r = 1, ndim
           fmm_contract_multi_Tq(r,s) = fac*T_mats(r,s,1)
        END DO
     END DO

     contract: DO L = 1, LMAX

        u = L*(L+1) +1
        p = (LMAX-L+1)*(LMAX-L+1)

        qmin = u-L
        qmax = MIN( u+L, p )

        qloop: DO q = qmin, qmax

           DO s = q, p
           CALL DAXPY_(ndim,vect(s),T_mats(1,s,q),1,fmm_contract_multi_Tq(1,q),1)
           END DO

           fac = vect(q)
           DO s = q+1, p
              DO r = 1, ndim
                 fmm_contract_multi_Tq(r,s) = fmm_contract_multi_Tq(r,s)    &
                                             + fac*T_mats(r,s,q)
              END DO
           END DO

        END DO qloop

        fmm_contract_multi_Tq(:,u) = half*fmm_contract_multi_Tq(:,u)    ! m=0

     END DO contract

   END FUNCTION fmm_contract_multi_Tq

!------------------------------------------------------------------------------

END MODULE fmm_multiple_T_worker
