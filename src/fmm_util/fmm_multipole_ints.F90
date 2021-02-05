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

module fmm_multipole_ints

use fmm_global_paras, only: INTK, REALK, LUINTM, fmm_basis, fmm_sh_pairs, fmm_prim_batch, Zero, Pi
implicit none
private
! Public procedures
public :: fmm_init_multipole_ints, &
          fmm_free_multipole_ints, &
          fmm_build_multipoles

real(REALK), allocatable :: ECoefX(:,:,:,:)
real(REALK), allocatable :: ECoefY(:,:,:,:)
real(REALK), allocatable :: ECoefZ(:,:,:,:)

real(REALK), allocatable :: MIntX(:,:,:)
real(REALK), allocatable :: MIntY(:,:,:)
real(REALK), allocatable :: MIntZ(:,:,:)

real(REALK), allocatable :: MpoleX(:,:,:,:)
real(REALK), allocatable :: MpoleY(:,:,:,:)
real(REALK), allocatable :: MpoleZ(:,:,:,:)

real(REALK), allocatable :: CarMpole(:,:,:)
real(REALK), allocatable :: SphMpole(:,:,:)

contains

!-------------------------------------------------------------------------------

subroutine fmm_init_multipole_ints(basis,MaxMul)

  implicit none
  type(fmm_basis), intent(in) :: basis
  integer(INTK), intent(in)   :: MaxMul

  integer(INTK) :: MaxAngl, MaxSgm2
  MaxAngl = basis%maxangl
  MaxSgm2 = basis%maxsgm2

  allocate(ECoefX(0:MaxAngl,0:MaxAngl,0:MaxAngl*2,MaxSgm2))
  allocate(ECoefY(0:MaxAngl,0:MaxAngl,0:MaxAngl*2,MaxSgm2))
  allocate(ECoefZ(0:MaxAngl,0:MaxAngl,0:MaxAngl*2,MaxSgm2))

  allocate(MIntX(-1:MaxMul+1,-1:MaxMul+1,MaxSgm2))
  allocate(MIntY(-1:MaxMul+1,-1:MaxMul+1,MaxSgm2))
  allocate(MIntZ(-1:MaxMul+1,-1:MaxMul+1,MaxSgm2))

  allocate(MpoleX(0:MaxAngl,0:MaxAngl,0:MaxMul,MaxSgm2))
  allocate(MpoleY(0:MaxAngl,0:MaxAngl,0:MaxMul,MaxSgm2))
  allocate(MpoleZ(0:MaxAngl,0:MaxAngl,0:MaxMul,MaxSgm2))

end subroutine fmm_init_multipole_ints

!-------------------------------------------------------------------------------

subroutine fmm_free_multipole_ints()

  implicit none
  deallocate(ECoefX)
  deallocate(ECoefY)
  deallocate(ECoefZ)

  deallocate(MIntX)
  deallocate(MIntY)
  deallocate(MIntZ)

  deallocate(MpoleX)
  deallocate(MpoleY)
  deallocate(MpoleZ)

  if (allocated(CarMpole)) deallocate(CarMpole)
  if (allocated(SphMpole)) deallocate(SphMpole)

end subroutine fmm_free_multipole_ints

!-------------------------------------------------------------------------------

subroutine fmm_build_multipoles(basis,MaxMul,sh_pairs)

  use fmm_car_to_sph, only: fmm_init_car_to_sph, &
                            fmm_free_car_to_sph
  use fmm_integral_utils, only: fmm_get_prim_batch, &
                                fmm_build_Ecoef1, &
                                fmm_build_Ecoef2

  implicit none
  type(fmm_basis), intent(in)    :: basis
  integer(INTK), intent(in)      :: MaxMul
  type(fmm_sh_pairs), intent(in) :: sh_pairs(:)

  !fixme
  character(len=10), parameter :: FName = 'multipoles'

  type(fmm_prim_batch) :: batch(basis%MaxSgm2)
  character(len=255) :: FBuf
  integer(INTK) :: Ish, Jsh, NPrim, IAnglA, IAnglB
  integer(INTK) :: i, ij, nmoms
  integer(INTK), external :: IsFreeUnit

  call fmm_init_car_to_sph(MaxMul)

  FBuf = trim(FName)//'.fmm1'
  LUINTM = IsFreeUnit(LUINTM)
  open(unit=LUINTM,file=trim(FBuf),status='REPLACE',access='SEQUENTIAL',form='UNFORMATTED')
  rewind(LUINTM)

  nmoms = 0
  do ij=1,size(sh_pairs)
    !print *, SIZE(sh_pairs), ij

    Ish = sh_pairs(ij)%I
    Jsh = sh_pairs(ij)%J

    ! Prepare batch of primitives
    call fmm_get_prim_batch(basis,Ish,Jsh,batch,NPrim)
    ! Build Hermite expansion E-coefficients
    IAnglA = basis%KType(Ish)
    IAnglB = basis%KType(Jsh)
    if (IAnglA >= IAnglB) then
      call fmm_build_Ecoef1(batch,NPrim,IAnglA,IAnglB,ECoefX,ECoefY,ECoefZ)
    else
      call fmm_build_Ecoef2(batch,NPrim,IAnglA,IAnglB,ECoefX,ECoefY,ECoefZ)
    end if

    ! Set multipole origin
    do i=1,NPrim
      batch(i)%PC(:) = batch(i)%P(:)-sh_pairs(ij)%centre(:)
    end do
    ! Build Hermite multipole moment integrals
    call fmm_build_Mints(batch,NPrim,MaxMul)

    ! Build spherical multipole integrals for Cartesian GTOs
    ! assuming final Fock matrix in spherical GTOs
    ! (only difference is we use spherical AO contraction coefs)
    ! When using these multipoles need Cartesian form of
    ! density and J_matrices
    call fmm_build_Mpole_ints(IAnglA,IAnglB,NPrim,MaxMul)
    call fmm_build_SphMpole(basis,batch,Ish,Jsh,NPrim,MaxMul)
    call fmm_store_SphMpole(basis,batch,Ish,Jsh,NPrim,MaxMul,sh_pairs(ij)%centre,nmoms)
  end do

  ! Mark end of file with negative L-index
  write(LUINTM) 0,-1,0,0,0,zero,zero,zero,zero
  close(unit=LUINTM,status='KEEP')

  FBuf = trim(FName)//'.fmm1header'
  LUINTM = IsFreeUnit(LUINTM)
  open(unit=LUINTM,file=trim(FBuf),status='REPLACE',access='SEQUENTIAL',form='UNFORMATTED')
  write(LUINTM) MaxMul,basis%nbas,nmoms
  close(unit=LUINTM,status='KEEP')

  call fmm_free_car_to_sph()

end subroutine fmm_build_multipoles

!-------------------------------------------------------------------------------

subroutine fmm_build_Mints(batch,NPrim,MaxMul)

  ! o Evaluate 1D multipole-moment integrals

  implicit none

  type(fmm_prim_batch), intent(in) :: batch(:)
  integer(INTK), intent(in)        :: NPrim, MaxMul

  integer(INTK) :: Ie, It, IJ

  real(REALK) :: PXCX(NPrim), PYCY(NPrim), PZCZ(NPrim)
  real(REALK) :: ExpPHalf(NPrim), ExpntP(NPrim)

  PXCX(:) = batch(1:Nprim)%PC(1)
  PYCY(:) = batch(1:Nprim)%PC(2)
  PZCZ(:) = batch(1:Nprim)%PC(3)
  ExpPHalf(:) = batch(1:Nprim)%ExpPHalf
  ExpntP(:) = batch(1:Nprim)%ExpntP

  ! o Recurrence relation
  !   M(t,e,NPrim)

  MIntX(:,:,:) = Zero
  MIntY(:,:,:) = Zero
  MIntZ(:,:,:) = Zero

  MIntX(0,0,1:NPrim) = sqrt(Pi/ExpntP(:))
  MIntY(0,0,1:NPrim) = sqrt(Pi/ExpntP(:))
  MIntZ(0,0,1:NPrim) = sqrt(Pi/ExpntP(:))

  do Ie=1,MaxMul
    do It=0,Ie
      do IJ=1,NPrim
        MIntX(It,Ie,IJ) = real(It,kind=REALK)*MIntX(It-1,Ie-1,IJ)
        MIntY(It,Ie,IJ) = real(It,kind=REALK)*MIntY(It-1,Ie-1,IJ)
        MIntZ(It,Ie,IJ) = real(It,kind=REALK)*MIntZ(It-1,Ie-1,IJ)
      end do
      do IJ=1,NPrim
        MIntX(It,Ie,IJ) = MIntX(It,Ie,IJ)+PXCX(IJ)*MIntX(It,Ie-1,IJ)
        MIntY(It,Ie,IJ) = MIntY(It,Ie,IJ)+PYCY(IJ)*MIntY(It,Ie-1,IJ)
        MIntZ(It,Ie,IJ) = MIntZ(It,Ie,IJ)+PZCZ(IJ)*MIntZ(It,Ie-1,IJ)
      end do
      do IJ=1,NPrim
        MIntX(It,Ie,IJ) = MIntX(It,Ie,IJ)+ExpPHalf(IJ)*MIntX(It+1,Ie-1,IJ)
        MIntY(It,Ie,IJ) = MIntY(It,Ie,IJ)+ExpPHalf(IJ)*MIntY(It+1,Ie-1,IJ)
        MIntZ(It,Ie,IJ) = MIntZ(It,Ie,IJ)+ExpPHalf(IJ)*MIntZ(It+1,Ie-1,IJ)
      end do
    end do
  end do

end subroutine fmm_build_Mints

!-------------------------------------------------------------------------------

subroutine fmm_build_Mpole_ints(IAnglA,IAnglB,NPrim,MaxMul)

  implicit none

  integer(INTK), intent(in) :: IAnglA, IAnglB, NPrim, MaxMul

  integer(INTK) :: It1, It2, Ie, It, IJ
  integer(INTK) :: ItMax

  MpoleX(:,:,:,:) = zero
  MpoleY(:,:,:,:) = zero
  MpoleZ(:,:,:,:) = zero

  do It1=0,IAnglA
    do It2=0,IAnglB
      do Ie=0,MaxMul
        ItMax = min((It1+It2),Ie)
        do It=0,ItMax
          do IJ=1,NPrim
            MpoleX(It1,It2,Ie,IJ) = MpoleX(It1,It2,Ie,IJ)+ECoefX(It1,It2,It,IJ)*MIntX(It,Ie,IJ)
            MpoleY(It1,It2,Ie,IJ) = MpoleY(It1,It2,Ie,IJ)+ECoefY(It1,It2,It,IJ)*MIntY(It,Ie,IJ)
            MpoleZ(It1,It2,Ie,IJ) = MpoleZ(It1,It2,Ie,IJ)+ECoefZ(It1,It2,It,IJ)*MIntZ(It,Ie,IJ)
          end do
        end do
      end do
    end do
  end do

end subroutine fmm_build_Mpole_ints

!-------------------------------------------------------------------------------

subroutine fmm_build_SphMpole(basis,batch,Ish,Jsh,NPrim,MaxMul)

  use fmm_car_to_sph, only: fmm_transform_car_to_sph

  implicit none
  type(fmm_basis), intent(in) :: basis
  type(fmm_prim_batch), intent(in) :: batch(:)
  integer(INTK), intent(in) :: Ish, Jsh, NPrim, MaxMul

  !fixme
  real(REALK), parameter :: ThrInt = 1.0e-12_REALK

  logical :: IEqJ
  integer(INTK) :: IAnglA, IAnglB
  integer(INTK) :: IL1, IL2, IL3
  integer(INTK) :: IeX, IeY, IeZ, MType, MCar
  integer(INTK) :: It1, Iu1, Iv1
  integer(INTK) :: It2, Iu2, Iv2
  integer(INTK) :: IJ, IBatch
  !integer(INTK) :: IIBatch
  integer(INTK) :: IL2Temp
  integer(INTK) :: Labelp, Labelq, ndim
  real(REALK) :: Temp1, Temp2

  IAnglA = basis%KType(Ish)
  IAnglB = basis%KType(Jsh)

  ndim = (basis%LtuvMax_Car(IAnglA)-basis%LtuvMin_Car(IAnglA)+1)*(basis%LtuvMax_Car(IAnglB)-basis%LtuvMin_Car(IAnglB)+1)
  ndim = ndim*NPrim

  allocate(CarMpole(ndim,(MaxMul+1)*(MaxMul+2)/2,0:MaxMul))
  if (.not. allocated(SphMpole)) then
    allocate(SphMpole(ndim,(2*MaxMul+1),0:MaxMul))
  end if

  CarMpole(:,:,:) = zero

  IEqJ = (ISh == JSh)
  IBatch = 0
  !IIBatch = 0
  Labelp = basis%KLoc_Car(ISh)
  do IL1=basis%LtuvMin_Car(IAnglA),basis%LtuvMax_Car(IAnglA)
    Labelp = Labelp+1
    It1 = basis%Lt(IL1)
    Iu1 = basis%Lu(IL1)
    Iv1 = basis%Lv(IL1)

    Labelq = basis%KLoc_Car(JSh)
    IL2Temp = basis%LtuvMax_Car(IAnglB)
    if (IEqJ) IL2Temp = IL1
    do IL2=basis%LtuvMin_Car(IAnglB),IL2Temp
      IBatch = IBatch+1
      Labelq = Labelq+1
      It2 = basis%Lt(IL2)
      Iu2 = basis%Lu(IL2)
      Iv2 = basis%Lv(IL2)

      !contracted functions
      do MType=0,MaxMul
        do MCar=1,(MType+1)*(MType+2)/2
          IL3 = MCar+basis%LtuvMin_Car(MType)-1
          IeX = basis%Lt(IL3)
          IeY = basis%Lu(IL3)
          IeZ = basis%Lv(IL3)

          Temp1 = Zero
          do IJ=1,NPrim
            Temp2 = MpoleX(It1,It2,IeX,IJ)*MpoleY(Iu1,Iu2,IeY,IJ)*MpoleZ(Iv1,Iv2,IeZ,IJ)*batch(IJ)%CCoefAB
            Temp1 = Temp1+Temp2
          end do

          if (abs(Temp1) >= ThrInt) CarMpole(IBatch,MCar,MType) = Temp1

        end do
      end do

      !primitive functions
      !do IJ=1, NPrim
      !  IIbatch = IIbatch +1
      !  do MType=0,MaxMul
      !    do MCar=1,(MType+1)*(MType+2)/2
      !
      !      IL3 = MCar+basis%LtuvMin_Car(MType)-1
      !      IeX = basis%Lt(IL3)
      !      IeY = basis%Lu(IL3)
      !      IeZ = basis%Lv(IL3)
      !      Temp1 = MpoleX(It1,It2,IeX,IJ)*MpoleY(Iu1,Iu2,IeY,IJ)*MpoleZ(Iv1,Iv2,IeZ,IJ)*batch(IJ)%CCoefAB
      !      if (abs(Temp1) >= ThrInt) CarMpole(IIBatch,MCar,MType) = Temp1
      !
      !    end do
      !  end do
      !end do
      !------------

    end do
  end do

  !print *, 'carmpole:'
  !do iex=1,size(carmpole,1)
  !  do MType=0,MaxMul
  !    do MCar=1,(MType+1)*(MType+2)/2
  !      print *, mtype,mcar,carmpole(iex,mcar,mtype)
  !    end do
  !  end do
  !end do

  ! Transform batch of cartesian integrals to spherical moments
  !ndim = MAX(Ibatch,IIbatch)
  ndim = size(CarMpole,1)
  call fmm_transform_car_to_sph(CarMpole,SphMpole,ndim,MaxMul)

  deallocate(CarMpole)

end subroutine fmm_build_SphMpole

!-------------------------------------------------------------------------------

subroutine fmm_store_SphMpole(basis,batch,Ish,Jsh,NPrim,MaxMul,mcntr,nmoms)

#include "macros.fh"

  implicit none
  type(fmm_basis), intent(in)      :: basis
  type(fmm_prim_batch), intent(in) :: batch(:)
  integer(INTK), intent(in)        :: ISh, JSh, MaxMul, NPrim
  real(REALK), intent(in)          :: mcntr(3)
  integer(INTK), intent(inout)     :: nmoms

  real(REALK), parameter :: MomScrn = 1.0e-15_REALK
  logical :: IEqJ, flag
  integer(INTK) :: IAnglA, IAnglB
  integer(INTK) :: IL1, IL2 !, IJ
  integer(INTK) :: MType, MSph
  integer(INTK) :: IBatch
  !integer(INTK) :: IIBatch
  integer(INTK) :: IL2Temp
  integer(INTK) :: Labelp, Labelq

  IEqJ = (ISh == JSh)
  IAnglA = basis%KType(Ish)
  IAnglB = basis%KType(Jsh)

  IBatch = 0
  !IIBatch = 0
  Labelp = basis%KLoc_Car(ISh)
  do IL1=basis%LtuvMin_Car(IAnglA),basis%LtuvMax_Car(IAnglA)
    Labelp = Labelp+1
    Labelq = basis%KLoc_Car(JSh)
    IL2Temp = basis%LtuvMax_Car(IAnglB)
    if (IEqJ) IL2Temp = IL1
    do IL2=basis%LtuvMin_Car(IAnglB),IL2Temp
      IBatch = IBatch+1
      Labelq = Labelq+1

      !contracted functions
      flag = .false.
      do MType=0,MaxMul
        do MSph=1,(MType+MType+1)
          ! Weak screening
          if (abs(SphMpole(IBatch,MSph,MType)) < MomScrn) cycle
          flag = .true.
          write(LUINTM) nmoms+1,MType,(MSph-MType-1),Labelp,Labelq,mcntr(1:3),SphMpole(IBatch,MSph,MType)
        end do
      end do
      if (flag) nmoms = nmoms+1

      unused_var(batch)
      unused_var(NPrim)

      !primitive functions
      !do IJ=1,NPrim
      !  IIBatch = IIBatch + 1
      !  flag = .false.
      !  do MType=0,MaxMul
      !    do MSph=1,(MType+MType+1)
      !      ! Weak screening
      !      if (ABS(SphMpole(IIBatch,MSph,MType)) < MomScrn) cycle
      !      flag = .true.
      !      write(LUINTM) nmoms,MType,(MSph-MType-1),Labelp,Labelq,batch(IJ)%P(1:3),SphMpole(IIBatch,MSph,MType)
      !    end do
      !  end do
      !  if (flag) nmoms = nmoms + 1
      !end do
      !------------

    end do
  end do

  if (allocated(SphMpole)) deallocate(SphMpole)

end subroutine fmm_store_SphMpole

!-------------------------------------------------------------------------------

end module fmm_multipole_ints
