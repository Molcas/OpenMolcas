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
MODULE fmm_car_to_sph

   USE fmm_global_paras
   IMPLICIT NONE
   PRIVATE
   ! public procedures
   PUBLIC :: fmm_init_car_to_sph,      &
             fmm_free_car_to_sph,      &
             fmm_transform_car_to_sph

   REAL(REALK), ALLOCATABLE, SAVE :: SphCoef(:,:,:)

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_init_car_to_sph(LMAX)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN) :: LMAX

      IF ( ALLOCATED(SphCoef) ) THEN
         CALL fmm_quit('mm_car_to_sph not freed correctly!')
      END IF

      ALLOCATE( SphCoef(0:LMax*2+1, (LMax+1)*(LMax+2)/2, 0:LMax) )
      CALL fmm_build_sphcoef(LMAX)

   END SUBROUTINE fmm_init_car_to_sph

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_free_car_to_sph

      DEALLOCATE( SphCoef )

   END SUBROUTINE fmm_free_car_to_sph

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_build_sphcoef(LMAX)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN) :: LMAX

      INTEGER(INTK) :: L, LL, I, M
      INTEGER(INTK) :: IX, IY, IZ
      INTEGER(INTK) :: IPX1, IPY1, IPZ1, IPY2, IPZ2, IPX2
      REAL(REALK)   :: Denom, PreFac, Arg
!
!     --- Initialization ---
!
      SphCoef(:,:,:) = Zero
!
!     --- s coefficients ---
!
      SphCoef(1,1,0) = One
      IF (LMax == 0) RETURN
!
!     --- p coefficients ---
!
      SphCoef(1,2,1) = One
      SphCoef(2,3,1) = One
      SphCoef(3,1,1) = One
      IF (LMax == 1) GO TO 100
!
!      SphCoef(1,2,2) =  Sqrt3
!      SphCoef(2,5,2) =  Sqrt3
!      SphCoef(3,1,2) = -Half
!      SphCoef(3,4,2) = -Half
!      SphCoef(3,6,2) =  One
!      SphCoef(4,3,2) =  Sqrt3
!      SphCoef(5,1,2) =  Half * Sqrt3
!      SphCoef(5,4,2) = -Half * Sqrt3
!
!     --- Higher than d ---
!
      DO LL = 2, LMax
         L = LL - 1
         Arg = DBLE(L+L+1) / DBLE(L+L+2)
         PreFac = SQRT(Arg)
!
         I = 0
         DO IX = L, 0, -1
            DO IY = L - IX, 0, -1
               IZ = L - IX - IY
               I = I + 1
               IPX1 = I
               IPY1 = I + (L - IX + 1)                                                    ! Increase 1 Y
               IPZ1 = I + (L - IX + 2)                                                    ! Increase 1 Z
!
!     --- Diagonal recurrence ---
!
               SphCoef(LL+LL+1,IPX1,LL) = SphCoef(LL+LL+1,IPX1,LL) +    &
                                          PreFac * SphCoef(L+L+1,I,L)
               SphCoef(LL+LL+1,IPY1,LL) = SphCoef(LL+LL+1,IPY1,LL) -    &
                                          PreFac * SphCoef(1,I,L)
!
               SphCoef(1,IPY1,LL) = SphCoef(1,IPY1,LL) +    &
                                    PreFac * SphCoef(L+L+1,I,L)
               SphCoef(1,IPX1,LL) = SphCoef(1,IPX1,LL) +    &
                                    PreFac * SphCoef(1,I,L)
!
!     --- Vertical recurrence (1) ---
!
               DO M = -LL + 1, LL - 1
                  Denom = SQRT(DBLE((L+M+1)*(L-M+1)))
                  SphCoef(LL+M+1,IPZ1,LL) = SphCoef(LL+M+1,IPZ1,LL) +   &
                          (DBLE(L+L+1)/Denom) * SphCoef(L+M+1,I,L)
               END DO
!
            END DO   ! IY
         END DO   ! IX
!
!     --- Vertical recurrence (2) ---
!
         I = 0
         DO IX = L - 1, 0, -1
            DO IY = L - 1 - IX, 0, -1
               IZ = L - 1 - IX - IY
               I = I + 1
               IPX2 = I
               IPY2 = I + (L - 1 - IX + 1) + (L - 1 - IX + 2)   ! Increase 2 Y
               IPZ2 = I + (L - 1 - IX + 2) + (L - 1 - IX + 3)   ! Increase 2 Z
!
               DO M = -LL + 1, LL - 1
                  Arg = DBLE((L+M)*(L-M)) / DBLE((L+M+1)*(L-M+1))
                  PreFac = SQRT(Arg)
                  SphCoef(LL+M+1,IPX2,LL) = SphCoef(LL+M+1,IPX2,LL) -   &
                                            PreFac * SphCoef(L+M,I,L-1)
                  SphCoef(LL+M+1,IPY2,LL) = SphCoef(LL+M+1,IPY2,LL) -   &
                                            PreFac * SphCoef(L+M,I,L-1)
                  SphCoef(LL+M+1,IPZ2,LL) = SphCoef(LL+M+1,IPZ2,LL) -   &
                                            PreFac * SphCoef(L+M,I,L-1)
               END DO   ! M
!
            END DO   ! IY
         END DO   ! IX
!
      END DO   ! LL
!
!     o p-1, p0, p+1 -> px, py, pz
!
  100 CONTINUE
      SphCoef(:,:,1) = Zero
      SphCoef(1,1,1) = One
      SphCoef(2,2,1) = One
      SphCoef(3,3,1) = One
!
   END SUBROUTINE fmm_build_sphcoef

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_transform_car_to_sph(CarMpole,SphMpole,ndim,LMAX)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN)  :: ndim, LMAX
      REAL(REALK),   INTENT(IN)  :: CarMpole(ndim,(LMAX+1)*(LMAX+2)/2,0:LMAX )
      REAL(REALK),   INTENT(OUT) :: SphMpole(ndim,2*LMAX+1,0:LMAX )

      INTEGER(INTK) :: MType, IS, IC, IJAO
      REAL(REALK)  :: Tmp, TempSphCoef, SphCoef_yzx(3,3)

      SphMpole(:,:,:) = Zero
      SphCoef_yzx(:,:) = Zero
      SphCoef_yzx(1,2) = One
      SphCoef_yzx(2,3) = One
      SphCoef_yzx(3,1) = One

      DO MType = 0, LMAX
         DO IS = 1, (MType + MType + 1)
            DO IC = 1, (MType + 1) * (MType + 2) / 2

               IF (MType /= 1) THEN
                  TempSphCoef = SphCoef(IS,IC,MType)
               ELSE
                  TempSphCoef = SphCoef_yzx(IS,IC)
               END IF

               DO IJAO = 1, ndim
                  Tmp =  TempSphCoef * CarMpole(IJAO,IC,MType)
                  SphMpole(IJAO,IS,MType) = SphMpole(IJAO,IS,MType) +Tmp
               END DO

            END DO
         END DO
      END DO

   END SUBROUTINE fmm_transform_car_to_sph

!-------------------------------------------------------------------------------

END MODULE fmm_car_to_sph

