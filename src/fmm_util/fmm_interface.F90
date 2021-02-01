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
MODULE fmm_interface

   USE fmm_global_paras, ONLY: INTK, REALK, LUPRI, LUINTM, fmm_basis, fmm_sh_pairs, scheme_paras, GFC_FMM, Zero, Two, Half
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_initial,                  &
             fmm_final,                    &
             fmm_get_J_matrix,             &
             fmm_get_boundary_potential
!             fmm_get_boundary_potential,   &
!             fefmm_get_J_matrix

   TYPE(fmm_basis), SAVE :: basis

CONTAINS

!-------------------------------------------------------------------------------

! INTEGER :: NAtom    =  number of atoms
! INTEGER :: NShel    =  number of shells
! INTEGER :: NPrim    =  number of primitive exponents
! INTEGER :: NBF_Car  =  dimension of AO Fock matrix

! INTEGER :: MaxAngl  =  maximum angular momentum primitive
! INTEGER :: MaxSgm2  =  maximum number of primitive pairs in any AO-pair

! INTEGER :: KAtom(NShel)     =  map from shells to atoms
! INTEGER :: KType(NShel)     =  angular momentum of shell
! INTEGER :: KStart(NShel)    =  first primitive in this shell
! INTEGER :: KontG(NShel)     =  number of primitives in this shell
! INTEGER :: KLoc_Car(NShel)  =  first (contracted) AO in this shell -1

! REAL    :: Centr(3,NAtom)   =  array of atomic positions
! REAL    :: Expnt(NPrim)     =  primitive exponents
! REAL    :: CCoef(NPrim)     =  contraction coefficients (normalized)

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_initial(NAtom,NShel,NPrim,NBF_Car,MaxAngl,MaxSgm2,    &
                          KAtom,KType,KStart,KontG,KLoc_Car,            &
                          Centr,Expnt,CCoef,mode)

      USE fmm_scheme_builder, ONLY: fmm_init_scheme, fmm_get_scheme
      USE fmm_shell_pairs,    ONLY: fmm_get_shell_pairs

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN) :: NAtom, NShel, NPrim, NBF_Car
      INTEGER(INTK), INTENT(IN) :: MaxSgm2, MaxAngl, mode
      INTEGER(INTK), INTENT(IN) :: KAtom(NShel)
      INTEGER(INTK), INTENT(IN) :: KType(NShel)
      INTEGER(INTK), INTENT(IN) :: KStart(NShel)
      INTEGER(INTK), INTENT(IN) :: KontG(NShel)
      INTEGER(INTK), INTENT(IN) :: KLoc_Car(NShel)
      REAL(REALK),   INTENT(IN) :: Centr(3,NAtom)
      REAL(REALK),   INTENT(IN) :: Expnt(NPrim)
      REAL(REALK),   INTENT(IN) :: CCoef(NPrim)

      TYPE(fmm_sh_pairs), POINTER :: sh_pairs(:)
      TYPE(scheme_paras), POINTER :: scheme
      INTEGER(INTK) :: lmax

      ! Initialise FMM run-type parameters
      CALL fmm_init_scheme(mode)
      CALL fmm_get_scheme(scheme)
      lmax = scheme%raw_lmax

      ! Initialise basis set information
      CALL fmm_init_basis(NAtom,NShel,NPrim,NBF_Car,MaxAngl,lmax,     &
                          maxsgm2,KAtom,KType,KStart,KontG,KLoc_Car,  &
                          Centr,Expnt,CCoef)

      ! Get list of non-vanishing shell pairs
      CALL fmm_get_shell_pairs(basis,sh_pairs)

      ! Build and store on disk density-independent integral components
!      SELECT CASE (mode)
!         CASE ( FE_FMM )
!            CALL fmm_init_fefmm(scheme,basis,lmax,sh_pairs)
!         CASE DEFAULT
            CALL fmm_init_md4fmm(basis,lmax,sh_pairs)
!       END SELECT

      NULLIFY(sh_pairs)

   END SUBROUTINE fmm_initial

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_init_md4fmm(basis,lmax,sh_pairs)

      USE fmm_multipole_ints, ONLY: fmm_init_multipole_ints,          &
                                    fmm_free_multipole_ints,          &
                                    fmm_build_multipoles

      IMPLICIT NONE
      TYPE(fmm_basis),    INTENT(IN) :: basis
      TYPE(fmm_sh_pairs), INTENT(IN) :: sh_pairs(:)
      INTEGER(INTK),      INTENT(IN) :: lmax

      CALL fmm_init_multipole_ints(basis,lmax)
      CALL fmm_build_multipoles(basis,lmax,sh_pairs)
      CALL fmm_free_multipole_ints

   END SUBROUTINE fmm_init_md4fmm

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_get_J_matrix(nBas,dens,fockAO)

      USE fmm_driver, ONLY: fmm_build_J_matrix

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN)    :: nBas
      REAL(REALK),   INTENT(IN)    :: dens(nBas*(nBas+1)/2)
      REAL(REALK),   INTENT(INOUT) :: fockAO(nBas*(nBas+1)/2)

      CHARACTER(LEN=10), PARAMETER :: Name = 'multipoles'
      CHARACTER(LEN=255) :: FBuf
      REAL(REALK)   :: J_matrix(nBas,nBas)
      REAL(REALK)   :: sq_dens(nBas,nBas)
      INTEGER(INTK) :: i,j,ij
      INTEGER(INTK), EXTERNAL :: IsFreeUnit

      ! Write null header file for nuclear moments (not computed!)
      FBuf = TRIM(Name)//".fmm2header"
      LUINTM = IsFreeUnit(LUINTM)
      OPEN(UNIT=LUINTM, FILE=TRIM(FBuf), STATUS='REPLACE',   &
           ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
      WRITE(LUINTM) 0
      CLOSE(UNIT=LUINTM, STATUS='KEEP')

      sq_dens(:,:) = zero
      ij = 0
      DO j = 1, nBas
         DO i = 1, j
            ij = ij+1
            sq_dens(i,j) = two*dens(ij)
            sq_dens(j,i) = sq_dens(i,j)
         END DO
         sq_dens(j,j) = half*sq_dens(j,j)
      END DO

      CALL fmm_build_J_matrix('TWO_EL',sq_dens,J_matrix)

      ij = 0
      DO j = 1, nBas
      DO i = 1, j
         ij = ij+1
         fockAO(ij) = fockAO(ij) + J_matrix(i,j)
      END DO
      END DO

   END SUBROUTINE fmm_get_J_matrix

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_get_boundary_potential(npoints,nBas,coor,dens,potential)

      USE fmm_stats,                ONLY: stat_points
      USE fmm_scheme_builder,       ONLY: fmm_get_scheme
      USE fmm_driver,               ONLY: fmm_get_multipole_potential
!      USE fmm_boxed_multipole_ints, ONLY: fmm_pack_boxed_mpoles

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN)  :: npoints, nBas
      REAL(REALK),   INTENT(IN)  :: coor(3,npoints)
      REAL(REALK),   INTENT(IN)  :: dens(nBas*(nBas+1)/2)
      REAL(REALK),   INTENT(OUT) :: potential(npoints)

      TYPE(scheme_paras), POINTER :: scheme
      REAL(REALK)   :: sq_dens(nBas,nBas)
      REAL(REALK)   :: Vtmp(1,npoints)
!      REAL(REALK)   :: dummy_dens(1,1) = one
      INTEGER(INTK) :: i,j,ij

      stat_points = npoints

      sq_dens(:,:) = zero
      ij = 0
      DO j = 1, nBas
         DO i = 1, j
            ij = ij+1
            sq_dens(i,j) = two*dens(ij)
            sq_dens(j,i) = sq_dens(i,j)
         END DO
         sq_dens(j,j) = half*sq_dens(j,j)
      END DO

      CALL fmm_get_scheme(scheme)
!      ! Pack moments (with density) to reduce memory in FMM readin
!      CALL fmm_pack_boxed_mpoles(scheme%raw_lmax,nBas,sq_dens)

      ! Rewrite files for LHS potential grid points
      CALL fmm_initialise_gfc_grid(npoints,coor)

!      CALL fmm_get_multipole_potential(GFC_FMM,dummy_dens,Vtmp)
      CALL fmm_get_multipole_potential(GFC_FMM,sq_dens,Vtmp)
      potential(:) = Vtmp(1,:)
      CALL fmm_matrix_norm('fmm_V',potential,npoints)

   END SUBROUTINE fmm_get_boundary_potential

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_initialise_gfc_grid(npoints,coor)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN)  :: npoints
      REAL(REALK),   INTENT(IN)  :: coor(3,npoints)

      CHARACTER(LEN=10), PARAMETER :: Name = 'multipoles'
      CHARACTER(LEN=255) :: FBuf
      INTEGER(INTK) :: i
      INTEGER(INTK), EXTERNAL :: IsFreeUnit

      ! Write grid points to disk
      FBuf = TRIM(Name)//".fmm2"
      LUINTM = IsFreeUnit(LUINTM)
      OPEN(UNIT=LUINTM, FILE=TRIM(FBuf), STATUS='REPLACE',  &
           ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
      REWIND(LUINTM)
      DO i = 1, npoints
         WRITE(LUINTM) 0, 0, 0, 0, 0, coor(1,i), coor(2,i), coor(3,i), 1.0D0
      END DO
      CLOSE(UNIT=LUINTM, STATUS='KEEP')

      ! Write header file
      FBuf = TRIM(Name)//".fmm2header"
      LUINTM = IsFreeUnit(LUINTM)
      OPEN(UNIT=LUINTM, FILE=TRIM(FBuf), STATUS='REPLACE',   &
           ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
      WRITE(LUINTM) npoints
      CLOSE(UNIT=LUINTM, STATUS='KEEP')

   END SUBROUTINE fmm_initialise_gfc_grid

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_init_basis(NAtom,NShel,NPrim,NBF_Car,MaxAngl,lmax,      &
                             maxsgm2,KAtom,KType,KStart,KontG,KLoc_Car,   &
                             Centr,Expnt,CCoef)

      IMPLICIT NONE
      INTEGER(INTK), INTENT(IN) :: NAtom,NShel,NPrim,NBF_Car
      INTEGER(INTK), INTENT(IN) :: MaxAngl,lmax,maxsgm2
      INTEGER(INTK), INTENT(IN) :: KAtom(NShel)
      INTEGER(INTK), INTENT(IN) :: KType(NShel)
      INTEGER(INTK), INTENT(IN) :: KStart(NShel)
      INTEGER(INTK), INTENT(IN) :: KontG(NShel)
      INTEGER(INTK), INTENT(IN) :: KLoc_Car(NShel)
      REAL(REALK),   INTENT(IN) :: Centr(3,NAtom)
      REAL(REALK),   INTENT(IN) :: Expnt(NPrim)
      REAL(REALK),   INTENT(IN) :: CCoef(NPrim)

      INTEGER(INTK), ALLOCATABLE :: tmp(:)
      INTEGER(INTK) :: IL, It, Iu, Iv, Ituv, ii
      INTEGER(INTK) :: Maxtuv, MaxLtuv

      Maxtuv  = MAX(MaxAngl,lmax)
      MaxLtuv = (Maxtuv+1) * (Maxtuv+2) * (Maxtuv+3) /6

      ALLOCATE(basis%Centr(3,NAtom))
      ALLOCATE(basis%KAtom(NShel))
      ALLOCATE(basis%KType(NShel))
      ALLOCATE(basis%KStart(NShel))
      ALLOCATE(basis%KontG(NShel))
      ALLOCATE(basis%KLoc_Car(NShel))
      ALLOCATE(basis%Expnt(NPrim))
      ALLOCATE(basis%CCoef(NPrim))
      ALLOCATE(basis%LtuvMin_Car(0:Maxtuv))
      ALLOCATE(basis%LtuvMax_Car(0:Maxtuv))
      ALLOCATE(basis%Lt(MaxLtuv))
      ALLOCATE(basis%Lu(MaxLtuv))
      ALLOCATE(basis%Lv(MaxLtuv))

      basis%nshells     = NShel
      basis%nbas        = NBF_Car
      basis%maxangl     = MaxAngl
      basis%maxsgm2     = MaxSgm2
      basis%Centr(:,:)  = Centr(:,:)
      basis%KAtom(:)    = KAtom(:)
      basis%KType(:)    = KType(:)
      basis%KStart(:)   = KStart(:)
      basis%KontG(:)    = KontG(:)
      basis%KLoc_Car(:) = KLoc_Car(:)
      basis%Expnt(:)    = Expnt(:)
      basis%CCoef(:)    = CCoef(:)

      ! Now initialise LtuvMin_Car, LtuvMax_Car, Lt, Lu, Lv

      ALLOCATE(tmp(0:Maxtuv))
      tmp(0) = 1
      basis%LtuvMin_Car(0) = 1
      basis%LtuvMax_Car(0) = 1
      DO IL = 1, Maxtuv
         tmp(IL) = (IL+1) * (IL+2) / 2
         basis%LtuvMin_Car(IL) = basis%LtuvMin_Car(IL-1) + tmp(IL-1)
         basis%LtuvMax_Car(IL) = basis%LtuvMax_Car(IL-1) + tmp(IL)
      END DO

      Ituv = 0
      DO IL = 0, Maxtuv
         DO It = Maxtuv, 0, -1
            DO Iu = Maxtuv, 0, -1
               DO Iv = Maxtuv, 0, -1
                  IF (It+Iu+Iv == IL) THEN
                     Ituv = Ituv + 1
                     basis%Lt(Ituv) = It
                     basis%Lu(Ituv) = Iu
                     basis%Lv(Ituv) = Iv
                  END IF
               END DO
            END DO
         END DO
      END DO
      DEALLOCATE(tmp)

      RETURN

      ! Print section
      write(LUPRI,*) '------------'
      write(LUPRI,*) 'NAtom',  NAtom
      write(LUPRI,*) 'NShel',  NShel
      write(LUPRI,*) 'NPrim',  NPrim
      write(LUPRI,*) 'NBF_Car',NBF_Car
      write(LUPRI,*) '------------'
      write(LUPRI,*) 'Maxangl', MaxAngl
      write(LUPRI,*) 'Maxsgm2', MaxSgm2
      write(LUPRI,*) '------------'
      write(LUPRI,*) 'katom:'
      do ii = 1, nshel
         write(LUPRI,*) ii, KAtom(ii)
      end do
      write(LUPRI,*) '------------'
      write(LUPRI,*) 'ktype:'
      do ii = 1, nshel
         write(LUPRI,*) ii, KType(ii)
      end do
      write(LUPRI,*) '------------'
      write(LUPRI,*) 'kstart:'
      do ii = 1, nshel
         write(LUPRI,*) ii, KStart(ii)
      end do
      write(LUPRI,*) '------------'
      write(LUPRI,*) 'kontg:'
      do ii = 1, nshel
         write(LUPRI,*) ii, KontG(ii)
      end do
      write(LUPRI,*) '------------'
      write(LUPRI,*) 'kloc_car:'
      do ii = 1, nshel
         write(LUPRI,*) ii, KLoc_Car(ii)
      end do
      write(LUPRI,*) '------------'
      write(LUPRI,*) 'expnt:'
      do ii = 1, nprim
       write(LUPRI,*) ii,Expnt(ii)
      end do
      write(LUPRI,*) '------------'
      write(LUPRI,*) 'ccoef:'
      do ii = 1, nprim
       write(LUPRI,*) ii, CCoef(ii)
      end do
      write(LUPRI,*) '------------'

   END SUBROUTINE fmm_init_basis

!-------------------------------------------------------------------------------

   SUBROUTINE fmm_final

      USE fmm_shell_pairs, ONLY: fmm_free_shell_pairs

      DEALLOCATE(basis%Centr)
      DEALLOCATE(basis%KAtom)
      DEALLOCATE(basis%KType)
      DEALLOCATE(basis%KStart)
      DEALLOCATE(basis%KontG)
      DEALLOCATE(basis%KLoc_Car)
      DEALLOCATE(basis%Expnt)
      DEALLOCATE(basis%CCoef)
      DEALLOCATE(basis%LtuvMin_Car)
      DEALLOCATE(basis%LtuvMax_Car)
      DEALLOCATE(basis%Lt)
      DEALLOCATE(basis%Lu)
      DEALLOCATE(basis%Lv)
      CALL fmm_free_shell_pairs

   END SUBROUTINE fmm_final

!-------------------------------------------------------------------------------

END MODULE fmm_interface
