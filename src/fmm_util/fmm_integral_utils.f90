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
MODULE fmm_integral_utils

   USE fmm_global_paras
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_get_prim_batch,    &
             fmm_build_Ecoef1,      &
             fmm_build_Ecoef2

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_get_prim_batch(basis,Ish,Jsh,batch,NPrim)

      IMPLICIT NONE
      TYPE(fmm_basis),      INTENT(IN)  :: basis
      INTEGER(INTK),        INTENT(IN)  :: Ish, Jsh
      TYPE(fmm_prim_batch), INTENT(OUT) :: batch(:)
      INTEGER(INTK),        INTENT(OUT) :: NPrim

!fixme
      REAL(REALK), PARAMETER :: ThrFac = -2.30258d+00 * 13.0d+00
      INTEGER(INTK) :: I,J
      INTEGER(INTK) :: IPrim1, JPrim1
      INTEGER(INTK) :: IPrim2, JPrim2, JPTemp2
      REAL(REALK)   :: Acentr(3), Bcentr(3), P(3), RAB(3)
      REAL(REALK)   :: ExpA, ExpB, ExpP, ExpPI, ExpAR2, R2AB, ExpKAB


      Acentr(:) = basis%Centr(:,basis%KAtom(Ish))
      IPrim1 = basis%KStart(Ish)
      IPrim2 = IPrim1 + basis%KontG(Ish) - 1

      Bcentr(:) = basis%Centr(:,basis%KAtom(Jsh))
      JPrim1 = basis%KStart(Jsh)
      JPrim2 = JPrim1 + basis%KontG(Jsh) - 1

      RAB(:) = Acentr(:) - Bcentr(:)
      R2AB = DOT_PRODUCT(RAB,RAB)

      NPrim = 0
      DO I = IPrim1, IPrim2
         ExpA = basis%Expnt(I)
         ExpAR2 = ExpA * R2AB
         JPTemp2 = JPrim2
         IF (Ish == Jsh) JPTemp2 = I
         DO J = JPrim1, JPTemp2
            ExpB = basis%Expnt(J)
            ExpP = ExpA + ExpB
            ExpPI = One / ExpP
            ExpKAB = - ExpAR2 * ExpB * ExpPI
            IF (ExpKAB >= ThrFac) THEN
               NPrim = NPrim +1
               batch(NPrim)%ExpntP = ExpP
               batch(NPrim)%ExpPHalf = half * ExpPI
               batch(NPrim)%PreFactAB = EXP(ExpKAB)
               batch(NPrim)%CCoefAB = basis%CCoef(I) * basis%CCoef(J)
               IF ((Ish == Jsh) .AND. (I /= J)) THEN
                  batch(NPrim)%CCoefAB = two * batch(NPrim)%CCoefAB
               END IF
               P(:) = (ExpA*Acentr(:) + ExpB*Bcentr(:)) * ExpPI
               batch(NPrim)%P(:)  = P(:)
               batch(NPrim)%PA(:) = P(:) - Acentr(:)
               batch(NPrim)%PB(:) = P(:) - Bcentr(:)
!print *, 'ExpP:', ExpP
!print *, 'PreFactAB:', batch(NPrim)%PreFactAB
!print *, 'CCoefAB:', batch(NPrim)%CCoefAB
!print *, 'CntrP:', batch(NPrim)%P
            END IF
         END DO
      END DO

   END SUBROUTINE fmm_get_prim_batch

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_build_Ecoef1(batch,NPrim,IAngl,JAngl,ECoefX,ECoefY,ECoefZ)

      IMPLICIT NONE
      TYPE(fmm_prim_batch), INTENT(IN)  :: batch(:)
      INTEGER(INTK),        INTENT(IN)  :: NPrim, IAngl, JAngl
      REAL(REALK),          INTENT(OUT) :: ECoefX(0:,0:,0:,:)
      REAL(REALK),          INTENT(OUT) :: ECoefY(0:,0:,0:,:)
      REAL(REALK),          INTENT(OUT) :: ECoefZ(0:,0:,0:,:)

      INTEGER :: I, J, It, IJ, JTemp
      REAL(REALK) :: PXAX(NPrim), PXBX(NPrim)
      REAL(REALK) :: PYAY(NPrim), PYBY(NPrim)
      REAL(REALK) :: PZAZ(NPrim), PZBZ(NPrim)
      REAL(REALK) :: ExpHalf(NPrim)
      REAL(REALK) :: PreFact(NPrim)

!     o Generate E-coefficients (Hermite expansion coefficients) with the three-term recurence relation

      PXAX(:) = batch(1:Nprim)%PA(1)
      PYAY(:) = batch(1:Nprim)%PA(2)
      PZAZ(:) = batch(1:Nprim)%PA(3)
      PXBX(:) = batch(1:Nprim)%PB(1)
      PYBY(:) = batch(1:Nprim)%PB(2)
      PZBZ(:) = batch(1:Nprim)%PB(3)
      ExpHalf(:) = batch(1:Nprim)%ExpPHalf
      PreFact(:) = batch(1:Nprim)%PreFactAB

!     o E(0,0,0,NPrim), The prefactor is incorporated.

      ECoefX(0,0,0,1:NPrim) = One
      ECoefY(0,0,0,1:NPrim) = One
      ECoefZ(0,0,0,1:NPrim) = PreFact(1:NPrim)

!     o Recurence relation
!       E(I,J,t,NPrim), 0 <= t <= I + J, IAngl >= JAngl

      IF (IAngl == 0) RETURN
      DO I = 1, IAngl
         JTemp = I
         IF (I == IAngl) JTemp = JAngl
         DO J = 0, JTemp

!      * Case of It = 0

            It = 0
            IF (I+J <= 1) THEN
               DO IJ = 1, NPrim
                  ECoefX(I,J,It,IJ) = PXAX(IJ) * ECoefX(I-1,J,It,IJ)
                  ECoefY(I,J,It,IJ) = PYAY(IJ) * ECoefY(I-1,J,It,IJ)
                  ECoefZ(I,J,It,IJ) = PZAZ(IJ) * ECoefZ(I-1,J,It,IJ)
               END DO
            ELSE
               DO IJ = 1, NPrim
                  ECoefX(I,J,It,IJ) = PXAX(IJ) * ECoefX(I-1,J,It,IJ)                      &
     &                              +            ECoefX(I-1,J,It+1,IJ)
                  ECoefY(I,J,It,IJ) = PYAY(IJ) * ECoefY(I-1,J,It,IJ)                      &
     &                              +            ECoefY(I-1,J,It+1,IJ)
                  ECoefZ(I,J,It,IJ) = PZAZ(IJ) * ECoefZ(I-1,J,It,IJ)                      &
     &                              +            ECoefZ(I-1,J,It+1,IJ)
               END DO
            END IF
!
!      * Case of It >= 1
!
            DO It = 1, I+J
               IF (It == I+J  ) THEN
                  DO IJ = 1, NPrim
                     ECoefX(I,J,It,IJ) = ExpHalf(IJ) * ECoefX(I-1,J,It-1,IJ)
                     ECoefY(I,J,It,IJ) = ExpHalf(IJ) * ECoefY(I-1,J,It-1,IJ)
                     ECoefZ(I,J,It,IJ) = ExpHalf(IJ) * ECoefZ(I-1,J,It-1,IJ)
                  END DO
               ELSE IF (It == I+J-1) THEN
                  DO IJ = 1, NPrim
                     ECoefX(I,J,It,IJ) = ExpHalf(IJ) * ECoefX(I-1,J,It-1,IJ)              &
     &                                 + PXAX(IJ)    * ECoefX(I-1,J,It  ,IJ)
                     ECoefY(I,J,It,IJ) = ExpHalf(IJ) * ECoefY(I-1,J,It-1,IJ)              &
     &                                 + PYAY(IJ)    * ECoefY(I-1,J,It  ,IJ)
                     ECoefZ(I,J,It,IJ) = ExpHalf(IJ) * ECoefZ(I-1,J,It-1,IJ)              &
     &                                 + PZAZ(IJ)    * ECoefZ(I-1,J,It  ,IJ)
                  END DO
               ELSE
                  DO IJ = 1, NPrim
                     ECoefX(I,J,It,IJ) = ExpHalf(IJ) * ECoefX(I-1,J,It-1,IJ)              &
     &                                 + PXAX(IJ)    * ECoefX(I-1,J,It  ,IJ)              &
     &                                 + (It+1)      * ECoefX(I-1,J,It+1,IJ)
                     ECoefY(I,J,It,IJ) = ExpHalf(IJ) * ECoefY(I-1,J,It-1,IJ)              &
     &                                 + PYAY(IJ)    * ECoefY(I-1,J,It  ,IJ)              &
     &                                 + (It+1)      * ECoefY(I-1,J,It+1,IJ)
                     ECoefZ(I,J,It,IJ) = ExpHalf(IJ) * ECoefZ(I-1,J,It-1,IJ)              &
     &                                 + PZAZ(IJ)    * ECoefZ(I-1,J,It  ,IJ)              &
     &                                 + (It+1)      * ECoefZ(I-1,J,It+1,IJ)
                  END DO
               END IF
            END DO
!
            IF (I /= j) THEN
!
!      * Case of It = 0
!
               It = 0
               IF (I+J <= 1) THEN
                  DO IJ = 1, NPrim
                     ECoefX(J,I,It,IJ) = PXBX(IJ) * ECoefX(J,I-1,It,IJ)
                     ECoefY(J,I,It,IJ) = PYBY(IJ) * ECoefY(J,I-1,It,IJ)
                     ECoefZ(J,I,It,IJ) = PZBZ(IJ) * ECoefZ(J,I-1,It,IJ)
                  END DO
               ELSE
                  DO IJ = 1, NPrim
                     ECoefX(J,I,It,IJ) = PXBX(IJ) * ECoefX(J,I-1,It,IJ)                   &
     &                                 +            ECoefX(J,I-1,It+1,IJ)
                     ECoefY(J,I,It,IJ) = PYBY(IJ) * ECoefY(J,I-1,It,IJ)                   &
     &                                 +            ECoefY(J,I-1,It+1,IJ)
                     ECoefZ(J,I,It,IJ) = PZBZ(IJ) * ECoefZ(J,I-1,It,IJ)                   &
     &                                 +            ECoefZ(J,I-1,It+1,IJ)
                  END DO
               END IF
!
!      * Case of It >= 1
!
               DO It = 1, I+J
                  IF (It == I+J  ) THEN
                     DO IJ = 1, NPrim
                        ECoefX(J,I,It,IJ) = ExpHalf(IJ) * ECoefX(J,I-1,It-1,IJ)
                        ECoefY(J,I,It,IJ) = ExpHalf(IJ) * ECoefY(J,I-1,It-1,IJ)
                        ECoefZ(J,I,It,IJ) = ExpHalf(IJ) * ECoefZ(J,I-1,It-1,IJ)
                     END DO
                  ELSE IF (It == I+J-1) THEN
                     DO IJ = 1, NPrim
                        ECoefX(J,I,It,IJ) = ExpHalf(IJ) * ECoefX(J,I-1,It-1,IJ)           &
     &                                    + PXBX(IJ)    * ECoefX(J,I-1,It  ,IJ)
                        ECoefY(J,I,It,IJ) = ExpHalf(IJ) * ECoefY(J,I-1,It-1,IJ)           &
     &                                    + PYBY(IJ)    * ECoefY(J,I-1,It  ,IJ)
                        ECoefZ(J,I,It,IJ) = ExpHalf(IJ) * ECoefZ(J,I-1,It-1,IJ)           &
     &                                    + PZBZ(IJ)    * ECoefZ(J,I-1,It  ,IJ)
                     END DO
                  ELSE
                     DO IJ = 1, NPrim
                        ECoefX(J,I,It,IJ) = ExpHalf(IJ) * ECoefX(J,I-1,It-1,IJ)           &
     &                                    + PXBX(IJ)    * ECoefX(J,I-1,It  ,IJ)           &
     &                                    + (It+1)      * ECoefX(J,I-1,It+1,IJ)
                        ECoefY(J,I,It,IJ) = ExpHalf(IJ) * ECoefY(J,I-1,It-1,IJ)           &
     &                                    + PYBY(IJ)    * ECoefY(J,I-1,It  ,IJ)           &
     &                                    + (It+1)      * ECoefY(J,I-1,It+1,IJ)
                        ECoefZ(J,I,It,IJ) = ExpHalf(IJ) * ECoefZ(J,I-1,It-1,IJ)           &
     &                                    + PZBZ(IJ)    * ECoefZ(J,I-1,It  ,IJ)           &
     &                                    + (It+1)      * ECoefZ(J,I-1,It+1,IJ)
                  END DO
                     END IF
               END DO
!
            END IF
!
         END DO
      END DO

   END SUBROUTINE fmm_build_Ecoef1

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_build_Ecoef2(batch,NPrim,IAngl,JAngl,ECoefX,ECoefY,ECoefZ)

      IMPLICIT NONE
      TYPE(fmm_prim_batch), INTENT(IN)  :: batch(:)
      INTEGER(INTK),        INTENT(IN)  :: NPrim, IAngl, JAngl
      REAL(REALK),          INTENT(OUT) :: ECoefX(0:,0:,0:,:)
      REAL(REALK),          INTENT(OUT) :: ECoefY(0:,0:,0:,:)
      REAL(REALK),          INTENT(OUT) :: ECoefZ(0:,0:,0:,:)

      INTEGER :: I, J, It, IJ, JTemp
      REAL(REALK) :: PXAX(NPrim), PXBX(NPrim)
      REAL(REALK) :: PYAY(NPrim), PYBY(NPrim)
      REAL(REALK) :: PZAZ(NPrim), PZBZ(NPrim)
      REAL(REALK) :: ExpHalf(NPrim)
      REAL(REALK) :: PreFact(NPrim)

!     o Generate E-coefficients (Hermite expansion coefficients) with the three-term recurence relation

      PXAX(:) = batch(1:Nprim)%PA(1)
      PYAY(:) = batch(1:Nprim)%PA(2)
      PZAZ(:) = batch(1:Nprim)%PA(3)
      PXBX(:) = batch(1:Nprim)%PB(1)
      PYBY(:) = batch(1:Nprim)%PB(2)
      PZBZ(:) = batch(1:Nprim)%PB(3)
      ExpHalf(:) = batch(1:Nprim)%ExpPHalf
      PreFact(:) = batch(1:Nprim)%PreFactAB


!     o E(0,0,0,NPrim), The prefactor is incorporated.

      DO IJ = 1, NPrim
         ECoefX(0,0,0,IJ) = One
         ECoefY(0,0,0,IJ) = One
         ECoefZ(0,0,0,IJ) = PreFact(IJ)
      END DO

!     o Recurence relation
!       E(I,J,t,NPrim), 0 <= t <= I + J, IAngl >= JAngl

      IF (JAngl == 0) RETURN
      DO J = 1, JAngl
         JTemp = J
         IF (J == JAngl) JTemp = IAngl
         DO I = 0, JTemp

!      * Case of It = 0

            It = 0
            IF (I+J <= 1) THEN
               DO IJ = 1, NPrim
                  ECoefX(I,J,It,IJ) = PXBX(IJ) * ECoefX(I,j-1,It,IJ)
                  ECoefY(I,J,It,IJ) = PYBY(IJ) * ECoefY(I,j-1,It,IJ)
                  ECoefZ(I,J,It,IJ) = PZBZ(IJ) * ECoefZ(I,j-1,It,IJ)
               END DO
            ELSE
               DO IJ = 1, NPrim
                  ECoefX(I,J,It,IJ) = PXBX(IJ) * ECoefX(I,j-1,It,IJ)                      &
     &                              +            ECoefX(I,j-1,It+1,IJ)
                  ECoefY(I,J,It,IJ) = PYBY(IJ) * ECoefY(I,j-1,It,IJ)                      &
     &                              +            ECoefY(I,j-1,It+1,IJ)
                  ECoefZ(I,J,It,IJ) = PZBZ(IJ) * ECoefZ(I,j-1,It,IJ)                      &
     &                              +            ECoefZ(I,j-1,It+1,IJ)
               END DO
            END IF
!
!      * Case of It >= 1
!
            DO It = 1, I+J
               IF (It == I+J  ) THEN
                  DO IJ = 1, NPrim
                     ECoefX(I,J,It,IJ) = ExpHalf(IJ) * ECoefX(I,j-1,It-1,IJ)
                     ECoefY(I,J,It,IJ) = ExpHalf(IJ) * ECoefY(I,j-1,It-1,IJ)
                     ECoefZ(I,J,It,IJ) = ExpHalf(IJ) * ECoefZ(I,j-1,It-1,IJ)
                  END DO
               ELSE IF (It == I+J-1) THEN
                  DO IJ = 1, NPrim
                     ECoefX(I,J,It,IJ) = ExpHalf(IJ) * ECoefX(I,j-1,It-1,IJ)              &
     &                                 + PXBX(IJ)    * ECoefX(I,j-1,It  ,IJ)
                     ECoefY(I,J,It,IJ) = ExpHalf(IJ) * ECoefY(I,j-1,It-1,IJ)              &
     &                                 + PYBY(IJ)    * ECoefY(I,j-1,It  ,IJ)
                     ECoefZ(I,J,It,IJ) = ExpHalf(IJ) * ECoefZ(I,j-1,It-1,IJ)              &
     &                                 + PZBZ(IJ)    * ECoefZ(I,j-1,It  ,IJ)
                  END DO
               ELSE
                  DO IJ = 1, NPrim
                     ECoefX(I,J,It,IJ) = ExpHalf(IJ) * ECoefX(I,j-1,It-1,IJ)              &
     &                                 + PXBX(IJ)    * ECoefX(I,j-1,It  ,IJ)              &
     &                                 + (It+1)      * ECoefX(I,j-1,It+1,IJ)
                     ECoefY(I,J,It,IJ) = ExpHalf(IJ) * ECoefY(I,j-1,It-1,IJ)              &
     &                                 + PYBY(IJ)    * ECoefY(I,j-1,It  ,IJ)              &
     &                                 + (It+1)      * ECoefY(I,j-1,It+1,IJ)
                     ECoefZ(I,J,It,IJ) = ExpHalf(IJ) * ECoefZ(I,j-1,It-1,IJ)              &
     &                                 + PZBZ(IJ)    * ECoefZ(I,j-1,It  ,IJ)              &
     &                                 + (It+1)      * ECoefZ(I,j-1,It+1,IJ)
                  END DO
               END IF
            END DO
!
            IF (I /= j) THEN
!
!      * Case of It = 0
!
               It = 0
               IF (I+J <= 1) THEN
                  DO IJ = 1, NPrim
                     ECoefX(J,I,It,IJ) = PXAX(IJ) * ECoefX(J-1,I,It,IJ)
                     ECoefY(J,I,It,IJ) = PYAY(IJ) * ECoefY(J-1,I,It,IJ)
                     ECoefZ(J,I,It,IJ) = PZAZ(IJ) * ECoefZ(J-1,I,It,IJ)
                  END DO
               ELSE
                  DO IJ = 1, NPrim
                     ECoefX(J,I,It,IJ) = PXAX(IJ) * ECoefX(J-1,I,It,IJ)                   &
     &                                 +            ECoefX(J-1,I,It+1,IJ)
                     ECoefY(J,I,It,IJ) = PYAY(IJ) * ECoefY(J-1,I,It,IJ)                   &
     &                                 +            ECoefY(J-1,I,It+1,IJ)
                     ECoefZ(J,I,It,IJ) = PZAZ(IJ) * ECoefZ(J-1,I,It,IJ)                   &
     &                                 +            ECoefZ(J-1,I,It+1,IJ)
                  END DO
               END IF
!
!      * Case of It >= 1
!
               DO It = 1, I+J
                  IF (It == I+J  ) THEN
                     DO IJ = 1, NPrim
                        ECoefX(J,I,It,IJ) = ExpHalf(IJ) * ECoefX(J-1,I,It-1,IJ)
                        ECoefY(J,I,It,IJ) = ExpHalf(IJ) * ECoefY(J-1,I,It-1,IJ)
                        ECoefZ(J,I,It,IJ) = ExpHalf(IJ) * ECoefZ(J-1,I,It-1,IJ)
                     END DO
                  ELSE IF (It == I+J-1) THEN
                     DO IJ = 1, NPrim
                        ECoefX(J,I,It,IJ) = ExpHalf(IJ) * ECoefX(J-1,I,It-1,IJ)           &
     &                                    + PXAX(IJ)    * ECoefX(J-1,I,It  ,IJ)
                        ECoefY(J,I,It,IJ) = ExpHalf(IJ) * ECoefY(J-1,I,It-1,IJ)           &
     &                                    + PYAY(IJ)    * ECoefY(J-1,I,It  ,IJ)
                        ECoefZ(J,I,It,IJ) = ExpHalf(IJ) * ECoefZ(J-1,I,It-1,IJ)           &
     &                                    + PZAZ(IJ)    * ECoefZ(J-1,I,It  ,IJ)
                     END DO
                  ELSE
                     DO IJ = 1, NPrim
                        ECoefX(J,I,It,IJ) = ExpHalf(IJ) * ECoefX(J-1,I,It-1,IJ)           &
     &                                    + PXAX(IJ)    * ECoefX(J-1,I,It  ,IJ)           &
     &                                    + (It+1)      * ECoefX(J-1,I,It+1,IJ)
                        ECoefY(J,I,It,IJ) = ExpHalf(IJ) * ECoefY(J-1,I,It-1,IJ)           &
     &                                    + PYAY(IJ)    * ECoefY(J-1,I,It  ,IJ)           &
     &                                    + (It+1)      * ECoefY(J-1,I,It+1,IJ)
                        ECoefZ(J,I,It,IJ) = ExpHalf(IJ) * ECoefZ(J-1,I,It-1,IJ)           &
     &                                    + PZAZ(IJ)    * ECoefZ(J-1,I,It  ,IJ)           &
     &                                    + (It+1)      * ECoefZ(J-1,I,It+1,IJ)
                     END DO
                  END IF
               END DO
!
            END IF
!
         END DO
      END DO

   END SUBROUTINE fmm_build_Ecoef2

!-------------------------------------------------------------------------------

END MODULE fmm_integral_utils
