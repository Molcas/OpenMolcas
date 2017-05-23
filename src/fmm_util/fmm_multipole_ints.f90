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
MODULE fmm_multipole_ints

   USE fmm_global_paras
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_init_multipole_ints,   &
             fmm_free_multipole_ints,   &
             fmm_build_multipoles

   REAL(REALK), ALLOCATABLE :: ECoefX(:,:,:,:)
   REAL(REALK), ALLOCATABLE :: ECoefY(:,:,:,:)
   REAL(REALK), ALLOCATABLE :: ECoefZ(:,:,:,:)

   REAL(REALK), ALLOCATABLE :: MIntX(:,:,:)
   REAL(REALK), ALLOCATABLE :: MIntY(:,:,:)
   REAL(REALK), ALLOCATABLE :: MIntZ(:,:,:)

   REAL(REALK), ALLOCATABLE :: MpoleX(:,:,:,:)
   REAL(REALK), ALLOCATABLE :: MpoleY(:,:,:,:)
   REAL(REALK), ALLOCATABLE :: MpoleZ(:,:,:,:)

   REAL(REALK), ALLOCATABLE :: CarMpole(:,:,:)
   REAL(REALK), ALLOCATABLE :: SphMpole(:,:,:)

CONTAINS

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_init_multipole_ints(basis,MaxMul)

      IMPLICIT NONE
      TYPE(fmm_basis), INTENT(IN) :: basis
      INTEGER(INTK),   INTENT(IN) :: MaxMul

      INTEGER(INTK) :: MaxAngl, MaxSgm2
      MaxAngl = basis%maxangl
      MaxSgm2 = basis%maxsgm2

      ALLOCATE(ECoefX(0:MaxAngl, 0:MaxAngl, 0:MaxAngl*2, MaxSgm2))
      ALLOCATE(ECoefY(0:MaxAngl, 0:MaxAngl, 0:MaxAngl*2, MaxSgm2))
      ALLOCATE(ECoefZ(0:MaxAngl, 0:MaxAngl, 0:MaxAngl*2, MaxSgm2))

      ALLOCATE(MIntX(-1:MaxMul+1, -1:MaxMul+1, MaxSgm2))
      ALLOCATE(MIntY(-1:MaxMul+1, -1:MaxMul+1, MaxSgm2))
      ALLOCATE(MIntZ(-1:MaxMul+1, -1:MaxMul+1, MaxSgm2))

      ALLOCATE(MpoleX(0:MaxAngl, 0:MaxAngl, 0:MaxMul, MaxSgm2))
      ALLOCATE(MpoleY(0:MaxAngl, 0:MaxAngl, 0:MaxMul, MaxSgm2))
      ALLOCATE(MpoleZ(0:MaxAngl, 0:MaxAngl, 0:MaxMul, MaxSgm2))


   END SUBROUTINE fmm_init_multipole_ints

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_free_multipole_ints

      IMPLICIT NONE
      DEALLOCATE(ECoefX)
      DEALLOCATE(ECoefY)
      DEALLOCATE(ECoefZ)

      DEALLOCATE(MIntX)
      DEALLOCATE(MIntY)
      DEALLOCATE(MIntZ)

      DEALLOCATE(MpoleX)
      DEALLOCATE(MpoleY)
      DEALLOCATE(MpoleZ)

      IF (ALLOCATED(CarMpole)) DEALLOCATE(CarMpole)
      IF (ALLOCATED(SphMpole)) DEALLOCATE(SphMpole)

   END SUBROUTINE fmm_free_multipole_ints

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_build_multipoles(basis,MaxMul,sh_pairs)

      USE fmm_car_to_sph,     ONLY: fmm_init_car_to_sph,      &
                                    fmm_free_car_to_sph
      USE fmm_integral_utils, ONLY: fmm_get_prim_batch,       &
                                    fmm_build_Ecoef1,         &
                                    fmm_build_Ecoef2


      IMPLICIT NONE
      TYPE(fmm_basis),    INTENT(IN) :: basis
      INTEGER(INTK),      INTENT(IN) :: MaxMul
      TYPE(fmm_sh_pairs), INTENT(IN) :: sh_pairs(:)

!fixme
      CHARACTER(LEN=10), PARAMETER :: Name = 'multipoles'

      TYPE(fmm_prim_batch) :: batch(basis%MaxSgm2)
      CHARACTER(LEN=255) :: FBuf
      INTEGER(INTK) :: Ish, Jsh, NPrim, IAnglA, IAnglB
      INTEGER(INTK) :: i, ij, nmoms

      CALL fmm_init_car_to_sph(MaxMul)

      FBuf = TRIM(Name)//".fmm1"
      OPEN(UNIT=LUINTM, FILE=TRIM(FBuf), STATUS='REPLACE',  &
           ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
      REWIND(LUINTM)

      nmoms = 0
      DO ij = 1, SIZE(sh_pairs)
!          print *, SIZE(sh_pairs), ij

          Ish = sh_pairs(ij)%I
          Jsh = sh_pairs(ij)%J

          ! Prepare batch of primitives
          CALL fmm_get_prim_batch(basis,Ish,Jsh,batch,NPrim)
          ! Build Hermite expansion E-coefficients
          IAnglA = basis%KType(Ish)
          IAnglB = basis%KType(Jsh)
          IF (IAnglA >= IAnglB) THEN
             CALL fmm_build_Ecoef1(batch,NPrim,IAnglA,IAnglB,    &
                                   ECoefX,ECoefY,ECoefZ)

          ELSE
             CALL fmm_build_Ecoef2(batch,NPrim,IAnglA,IAnglB,    &
                                   ECoefX,ECoefY,ECoefZ)
          END IF

          ! Set multipole origin
          DO i = 1, NPrim
             batch(i)%PC(:) = batch(i)%P(:) - sh_pairs(ij)%centre(:)
          END DO
          ! Build Hermite multipole moment integrals
          CALL fmm_build_Mints(batch,NPrim,IAnglA+IAnglB,MaxMul)

          ! Build spherical multipole integrals for Cartesian GTOs
          ! assuming final Fock matrix in spherical GTOs
          ! (only difference is we use spherical AO contraction coefs)
          ! When using these multipoles need Cartesian form of
          ! density and J_matrices
          CALL fmm_build_Mpole_ints(IAnglA,IAnglB,NPrim,MaxMul)
          CALL fmm_build_SphMpole(basis,batch,Ish,Jsh,NPrim,MaxMul)
          CALL fmm_store_SphMpole(basis,batch,Ish,Jsh,NPrim,MaxMul,  &
                                  sh_pairs(ij)%centre,nmoms)
      END DO

      ! Mark end of file with negative L-index
      WRITE(LUINTM) 0,-1,0,0,0, 0d0,0d0,0d0, 0d0
      CLOSE(UNIT=LUINTM, STATUS='KEEP')

      FBuf = TRIM(Name)//".fmm1header"
      OPEN(UNIT=LUINTM, FILE=TRIM(FBuf), STATUS='REPLACE',   &
           ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
      WRITE(LUINTM) MaxMul, basis%nbas, nmoms
      CLOSE(UNIT=LUINTM, STATUS='KEEP')

      CALL fmm_free_car_to_sph

   END SUBROUTINE fmm_build_multipoles

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_build_Mints(batch,NPrim,MaxAng,MaxMul)

!     o Evaluate 1D multipole-moment integrals

      IMPLICIT NONE

      TYPE(fmm_prim_batch), INTENT(IN) :: batch(:)
      INTEGER(INTK),        INTENT(IN) :: NPrim, MaxAng, MaxMul

      INTEGER :: Ie, It, IJ

      REAL(REALK) :: PXCX(NPrim), PYCY(NPrim), PZCZ(NPrim)
      REAL(REALK) :: ExpPHalf(NPrim), ExpntP(NPrim)

      PXCX(:) = batch(1:Nprim)%PC(1)
      PYCY(:) = batch(1:Nprim)%PC(2)
      PZCZ(:) = batch(1:Nprim)%PC(3)
      ExpPHalf(:) = batch(1:Nprim)%ExpPHalf
      ExpntP(:) = batch(1:Nprim)%ExpntP

!     o Recurence relation
!       M(t,e,NPrim)

      MIntX(:,:,:) = Zero
      MIntY(:,:,:) = Zero
      MIntZ(:,:,:) = Zero

      MIntX(0,0,1:NPrim) = SQRT(Pi/ExpntP(:))
      MIntY(0,0,1:NPrim) = SQRT(Pi/ExpntP(:))
      MIntZ(0,0,1:NPrim) = SQRT(Pi/ExpntP(:))

      DO Ie = 1, MaxMul
         DO It = 0, Ie
            DO IJ = 1, NPrim
               MIntX(It,Ie,IJ) = DBLE(It) * MIntX(It-1,Ie-1,IJ)
               MIntY(It,Ie,IJ) = DBLE(It) * MIntY(It-1,Ie-1,IJ)
               MIntZ(It,Ie,IJ) = DBLE(It) * MIntZ(It-1,Ie-1,IJ)
            END DO
            DO IJ = 1, NPrim
               MIntX(It,Ie,IJ) = MIntX(It,Ie,IJ) + PXCX(IJ) * MIntX(It,Ie-1,IJ)
               MIntY(It,Ie,IJ) = MIntY(It,Ie,IJ) + PYCY(IJ) * MIntY(It,Ie-1,IJ)
               MIntZ(It,Ie,IJ) = MIntZ(It,Ie,IJ) + PZCZ(IJ) * MIntZ(It,Ie-1,IJ)
            END DO
            DO IJ = 1, NPrim
               MIntX(It,Ie,IJ) = MIntX(It,Ie,IJ) + ExpPHalf(IJ) * MIntX(It+1,Ie-1,IJ)
               MIntY(It,Ie,IJ) = MIntY(It,Ie,IJ) + ExpPHalf(IJ) * MIntY(It+1,Ie-1,IJ)
               MIntZ(It,Ie,IJ) = MIntZ(It,Ie,IJ) + ExpPHalf(IJ) * MIntZ(It+1,Ie-1,IJ)
            END DO
         END DO
      END DO

! Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(MaxAng)
   END SUBROUTINE fmm_build_Mints

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_build_Mpole_ints(IAnglA,IAnglB,NPrim,MaxMul)

      IMPLICIT NONE

      INTEGER(INTK), INTENT(IN) :: IAnglA, IAnglB, NPrim, MaxMul

      INTEGER(INTK) :: It1, It2, Ie, It, IJ
      INTEGER(INTK) :: ItMax

      MpoleX(:,:,:,:) = zero
      MpoleY(:,:,:,:) = zero
      MpoleZ(:,:,:,:) = zero

      DO It1 = 0, IAnglA
         DO It2 = 0, IAnglB
            DO Ie = 0, MaxMul
               ItMax = MIN((It1+It2),Ie)
               DO It = 0, ItMax
                  DO IJ = 1, NPrim
                     MpoleX(It1,It2,Ie,IJ) = MpoleX(It1,It2,Ie,IJ) +    &
                                             ECoefX(It1,It2,It,IJ) *    &
                                             MIntX(It,Ie,IJ)
                     MpoleY(It1,It2,Ie,IJ) = MpoleY(It1,It2,Ie,IJ) +    &
                                             ECoefY(It1,It2,It,IJ) *    &
                                             MIntY(It,Ie,IJ)
                     MpoleZ(It1,It2,Ie,IJ) = MpoleZ(It1,It2,Ie,IJ) +    &
                                             ECoefZ(It1,It2,It,IJ) *    &
                                             MIntZ(It,Ie,IJ)
                  END DO
               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE fmm_build_Mpole_ints

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_build_SphMpole(basis,batch,Ish,Jsh,NPrim,MaxMul)

      USE fmm_car_to_sph, ONLY: fmm_transform_car_to_sph

      IMPLICIT NONE
      TYPE(fmm_basis),      INTENT(IN) :: basis
      TYPE(fmm_prim_batch), INTENT(IN) :: batch(:)
      INTEGER(INTK),        INTENT(IN) :: Ish, Jsh, NPrim, MaxMul

!fixme
      REAL(REALK), PARAMETER :: ThrInt = 1.0d-12

      LOGICAL :: IEqJ
      INTEGER :: IAnglA, IAnglB
      INTEGER :: IL1, IL2, IL3
      INTEGER :: IeX, IeY, IeZ, MType, MCar
      INTEGER :: It1, Iu1, Iv1
      INTEGER :: It2, Iu2, Iv2
      INTEGER :: IJ, IBatch
      INTEGER :: IIBatch
      INTEGER :: IL2Temp
      INTEGER :: Labelp, Labelq, ndim
      REAL(REALK) :: Temp1, Temp2

      IAnglA = basis%KType(Ish)
      IAnglB = basis%KType(Jsh)

      ndim =  (basis%LtuvMax_Car(IAnglA) - basis%LtuvMin_Car(IAnglA) + 1)   &
            * (basis%LtuvMax_Car(IAnglB) - basis%LtuvMin_Car(IAnglB) + 1)
      ndim = ndim * NPrim

      ALLOCATE(CarMpole(ndim, (MaxMul+1)*(MaxMul+2)/2, 0:MaxMul))
      IF (.NOT.ALLOCATED(SphMpole)) THEN
         ALLOCATE(SphMpole(ndim, (2*MaxMul+1), 0:MaxMul))
      END IF

      CarMpole(:,:,:) = zero

      IEqJ = (ISh == JSh)
      IBatch = 0
      IIBatch = 0
      Labelp = basis%KLoc_Car(ISh)
      DO IL1 = basis%LtuvMin_Car(IAnglA), basis%LtuvMax_Car(IAnglA)
         Labelp = Labelp + 1
         It1 = basis%Lt(IL1)
         Iu1 = basis%Lu(IL1)
         Iv1 = basis%Lv(IL1)

         Labelq = basis%KLoc_Car(JSh)
         IL2Temp = basis%LtuvMax_Car(IAnglB)
         IF (IEqJ) IL2Temp = IL1
         DO IL2 = basis%LtuvMin_Car(IAnglB), IL2Temp
            IBatch = IBatch + 1
            Labelq = Labelq + 1
            It2 = basis%Lt(IL2)
            Iu2 = basis%Lu(IL2)
            Iv2 = basis%Lv(IL2)

!contracted functions
            DO MType = 0, MaxMul
               DO MCar = 1, (MType+1)*(MType+2)/2
                  IL3 = MCar + basis%LtuvMin_Car(MType) - 1
                  IeX = basis%Lt(IL3)
                  IeY = basis%Lu(IL3)
                  IeZ = basis%Lv(IL3)

                  Temp1 = Zero
                  DO IJ = 1, NPrim
                     Temp2 = MpoleX(It1,It2,IeX,IJ) *    &
                             MpoleY(Iu1,Iu2,IeY,IJ) *    &
                             MpoleZ(Iv1,Iv2,IeZ,IJ) *    &
                             batch(IJ)%CCoefAB
                     Temp1 = Temp1 + Temp2
                  END DO

                  IF (ABS(Temp1) >= ThrInt) &
                     CarMpole(IBatch,MCar,MType) = Temp1

               END DO
            END DO

!primitive functions
!            DO IJ = 1, NPrim
!               IIbatch = IIbatch +1
!               DO MType = 0, MaxMul
!                  DO MCar = 1, (MType+1)*(MType+2)/2
!
!                     IL3 = MCar + basis%LtuvMin_Car(MType) - 1
!                     IeX = basis%Lt(IL3)
!                     IeY = basis%Lu(IL3)
!                     IeZ = basis%Lv(IL3)
!                     Temp1 = MpoleX(It1,It2,IeX,IJ) *    &
!                             MpoleY(Iu1,Iu2,IeY,IJ) *    &
!                             MpoleZ(Iv1,Iv2,IeZ,IJ) *    &
!                             batch(IJ)%CCoefAB
!                     IF (ABS(Temp1) >= ThrInt) &
!                        CarMpole(IIBatch,MCar,MType) = Temp1
!
!                  END DO
!               END DO
!            END DO

!------------
         END DO
      END DO

!      print *, 'carmpole:'
!      do iex = 1, size(carmpole,1)
!         DO MType = 0, MaxMul
!            DO MCar = 1, (MType+1)*(MType+2)/2
!               print *, mtype,mcar, carmpole(iex,mcar,mtype)
!            end do
!         end do
!      end do

      ! Transform batch of cartesian integrals to spherical moments
!      ndim = MAX(Ibatch,IIbatch)
      ndim = SIZE(CarMpole,1)
      CALL fmm_transform_car_to_sph(CarMpole,SphMpole,ndim,MaxMul)

      DEALLOCATE(CarMpole)

   END SUBROUTINE fmm_build_SphMpole

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_store_SphMpole(basis,batch,Ish,Jsh,NPrim,MaxMul,mcntr,nmoms)

      IMPLICIT NONE
      TYPE(fmm_basis),      INTENT(IN)    :: basis
      TYPE(fmm_prim_batch), INTENT(IN)    :: batch(:)
      INTEGER(INTK),        INTENT(IN)    :: ISh, JSh, MaxMul, NPrim
      REAL(REALK),          INTENT(IN)    :: mcntr(3)
      INTEGER(INTK),        INTENT(INOUT) :: nmoms

      REAL(REALK), PARAMETER :: MomScrn = 1d-15
      LOGICAL :: IEqJ, flag
      INTEGER :: IAnglA, IAnglB
      INTEGER :: IL1, IL2 !, IJ
      INTEGER :: MType, MSph
      INTEGER :: IBatch
      INTEGER :: IIBatch
      INTEGER :: IL2Temp
      INTEGER :: Labelp, Labelq

      IEqJ = (ISh == JSh)
      IAnglA = basis%KType(Ish)
      IAnglB = basis%KType(Jsh)

      IBatch = 0
      IIBatch = 0
      Labelp = basis%KLoc_Car(ISh)
      DO IL1 = basis%LtuvMin_Car(IAnglA), basis%LtuvMax_Car(IAnglA)
         Labelp = Labelp + 1
         Labelq = basis%KLoc_Car(JSh)
         IL2Temp = basis%LtuvMax_Car(IAnglB)
         IF (IEqJ) IL2Temp = IL1
         DO IL2 = basis%LtuvMin_Car(IAnglB), IL2Temp
            IBatch = IBatch + 1
            Labelq = Labelq + 1

!contracted functions
            flag = .FALSE.
            DO MType = 0, MaxMul
               DO MSph = 1, (MType+MType+1)
                  ! Weak screening
                  IF (ABS(SphMpole(IBatch,MSph,MType)) < MomScrn) CYCLE
                  flag = .TRUE.
                  WRITE(LUINTM) nmoms+1, MType, (MSph-MType-1),      &
                                Labelp, Labelq,                      &
                                mcntr(1:3),                          &
                                SphMpole(IBatch,MSph,MType)
               END DO
            END DO
            IF (flag) nmoms = nmoms + 1

!primitive functions

!            DO IJ = 1, NPrim
!               IIBatch = IIBatch + 1
!               flag = .FALSE.
!               DO MType = 0, MaxMul
!                  DO MSph = 1, (MType+MType+1)
!                     ! Weak screening
!                     IF (ABS(SphMpole(IIBatch,MSph,MType)) < MomScrn) CYCLE
!                     WRITE(LUINTM) nmoms, MType, (MSph-MType-1),        &
!                     flag = .TRUE.
!                                   Labelp, Labelq,                      &
!                                   batch(IJ)%P(1:3),                    &
!                                   SphMpole(IIBatch,MSph,MType)
!                  END DO
!               END DO
!               IF (flag) nmoms = nmoms + 1
!            END DO

!------------
         END DO
      END DO

      IF (ALLOCATED(SphMpole)) DEALLOCATE(SphMpole)

! Avoid unused argument warnings
      IF (.FALSE.) THEN
         CALL Unused_real_array(batch) ! not really real, but well...
         CALL Unused_integer(NPrim)
      END IF
   END SUBROUTINE fmm_store_SphMpole

!-------------------------------------------------------------------------------

END MODULE fmm_multipole_ints

