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

module fmm_interface

use fmm_global_paras, only: INTK, REALK, LUPRI, LUINTM, fmm_basis, fmm_sh_pairs, scheme_paras, GFC_FMM, Zero, One, Two, Half
implicit none
private
! Public procedures
public :: fmm_initial, &
          fmm_final, &
          fmm_get_J_matrix, &
          fmm_get_boundary_potential
          !fmm_get_boundary_potential,   &
          !fefmm_get_J_matrix

type(fmm_basis), save :: basis

contains

!-------------------------------------------------------------------------------

! integer(INTK) :: NAtom    =  number of atoms
! integer(INTK) :: NShel    =  number of shells
! integer(INTK) :: NPrim    =  number of primitive exponents
! integer(INTK) :: NBF_Car  =  dimension of AO Fock matrix

! integer(INTK) :: MaxAngl  =  maximum angular momentum primitive
! integer(INTK) :: MaxSgm2  =  maximum number of primitive pairs in any AO-pair

! integer(INTK) :: KAtom(NShel)     =  map from shells to atoms
! integer(INTK) :: KType(NShel)     =  angular momentum of shell
! integer(INTK) :: KStart(NShel)    =  first primitive in this shell
! integer(INTK) :: KontG(NShel)     =  number of primitives in this shell
! integer(INTK) :: KLoc_Car(NShel)  =  first (contracted) AO in this shell -1

! real(REALK)   :: Centr(3,NAtom)   =  array of atomic positions
! real(REALK)   :: Expnt(NPrim)     =  primitive exponents
! real(REALK)   :: CCoef(NPrim)     =  contraction coefficients (normalized)

!-------------------------------------------------------------------------------

subroutine fmm_initial(NAtom,NShel,NPrim,NBF_Car,MaxAngl,MaxSgm2,KAtom,KType,KStart,KontG,KLoc_Car,Centr,Expnt,CCoef,mode)

  use fmm_scheme_builder, only: fmm_init_scheme, fmm_get_scheme
  use fmm_shell_pairs, only: fmm_get_shell_pairs

  implicit none
  integer(INTK), intent(in) :: NAtom, NShel, NPrim, NBF_Car
  integer(INTK), intent(in) :: MaxSgm2, MaxAngl, mode
  integer(INTK), intent(in) :: KAtom(NShel)
  integer(INTK), intent(in) :: KType(NShel)
  integer(INTK), intent(in) :: KStart(NShel)
  integer(INTK), intent(in) :: KontG(NShel)
  integer(INTK), intent(in) :: KLoc_Car(NShel)
  real(REALK), intent(in)   :: Centr(3, NAtom)
  real(REALK), intent(in)   :: Expnt(NPrim)
  real(REALK), intent(in)   :: CCoef(NPrim)

  type(fmm_sh_pairs), pointer :: sh_pairs(:)
  type(scheme_paras), pointer :: scheme
  integer(INTK) :: lmax

  ! Initialise FMM run-type parameters
  call fmm_init_scheme(mode)
  call fmm_get_scheme(scheme)
  lmax = scheme%raw_lmax

  ! Initialise basis set information
  call fmm_init_basis(NAtom,NShel,NPrim,NBF_Car,MaxAngl,lmax,maxsgm2,KAtom,KType,KStart,KontG,KLoc_Car,Centr,Expnt,CCoef)

  ! Get list of non-vanishing shell pairs
  call fmm_get_shell_pairs(basis,sh_pairs)

  ! Build and store on disk density-independent integral components
  select case (mode)
    !case (FE_FMM)
    !  call fmm_init_fefmm(scheme,basis,lmax,sh_pairs)
    case default
      call fmm_init_md4fmm(basis,lmax,sh_pairs)
  end select

  nullify(sh_pairs)

end subroutine fmm_initial

!-------------------------------------------------------------------------------

subroutine fmm_init_md4fmm(basis,lmax,sh_pairs)

  use fmm_multipole_ints, only: fmm_init_multipole_ints, &
                                fmm_free_multipole_ints, &
                                fmm_build_multipoles

  implicit none
  type(fmm_basis), intent(in)    :: basis
  type(fmm_sh_pairs), intent(in) :: sh_pairs(:)
  integer(INTK), intent(in)      :: lmax

  call fmm_init_multipole_ints(basis,lmax)
  call fmm_build_multipoles(basis,lmax,sh_pairs)
  call fmm_free_multipole_ints()

end subroutine fmm_init_md4fmm

!-------------------------------------------------------------------------------

subroutine fmm_get_J_matrix(nBas,dens,fockAO)

  use fmm_driver, only: fmm_build_J_matrix

  implicit none
  integer(INTK), intent(in)  :: nBas
  real(REALK), intent(in)    :: dens(nBas*(nBas+1)/2)
  real(REALK), intent(inoUt) :: fockAO(nBas*(nBas+1)/2)

  character(len=10), parameter :: FName='multipoles'
  character(len=255) :: FBuf
  real(REALK) :: J_matrix(nBas,nBas)
  real(REALK) :: sq_dens(nBas,nBas)
  integer(INTK) :: i, j, ij
  integer(INTK), external :: IsFreeUnit

  ! Write null header file for nuclear moments (not computed!)
  FBuf = trim(FName)//'.fmm2header'
  LUINTM = IsFreeUnit(LUINTM)
  open(unit=LUINTM,file=trim(FBuf),status='REPLACE',access='SEQUENTIAL',form='UNFORMATTED')
  write(LUINTM) 0
  close(unit=LUINTM,status='KEEP')

  sq_dens(:,:) = zero
  ij = 0
  do j=1,nBas
    do i=1,j
      ij = ij+1
       sq_dens(i,j) = two*dens(ij)
     sq_dens(j,i) = sq_dens(i,j)
    end do
    sq_dens(j,j) = half*sq_dens(j,j)
  end do

  call fmm_build_J_matrix('TWO_EL',sq_dens,J_matrix)

  ij=0
  do j=1,nBas
    do i=1,j
      ij = ij+1
      fockAO(ij) = fockAO(ij)+J_matrix(i,j)
    end do
  end do

end subroutine fmm_get_J_matrix

!-------------------------------------------------------------------------------

subroutine fmm_get_boundary_potential(npoints,nBas,coor,dens,potential)

  use fmm_stats, only: stat_points
  use fmm_scheme_builder, only: fmm_get_scheme
  use fmm_driver, only: fmm_get_multipole_potential
  !use fmm_boxed_multipole_ints, only: fmm_pack_boxed_mpoles
  use fmm_utils, only: fmm_matrix_norm

  implicit none
  integer(INTK), intent(in) :: npoints, nBas
  real(REALK), intent(in)   :: coor(3,npoints)
  real(REALK), intent(in)   :: dens(nBas*(nBas+1)/2)
  real(REALK), intent(out)  :: potential(npoints)

  type(scheme_paras), pointer :: scheme
  real(REALK) :: sq_dens(nBas,nBas)
  real(REALK) :: Vtmp(1,npoints)
  !real(REALK) :: dummy_dens(1,1) = one
  integer(INTK) :: i, j, ij

  stat_points = npoints

  sq_dens(:,:) = zero
  ij = 0
  do j=1,nBas
    do i=1,j
      ij = ij+1
      sq_dens(i,j) = two*dens(ij)
      sq_dens(j,i) = sq_dens(i,j)
    end do
    sq_dens(j,j) = half*sq_dens(j,j)
  end do

  call fmm_get_scheme(scheme)
  ! Pack moments (with density) to reduce memory in FMM readin
  !call fmm_pack_boxed_mpoles(scheme%raw_lmax,nBas,sq_dens)

  ! Rewrite files for LHS potential grid points
  call fmm_initialise_gfc_grid(npoints,coor)

  !call fmm_get_multipole_potential(GFC_FMM,dummy_dens,Vtmp)
  call fmm_get_multipole_potential(GFC_FMM,sq_dens,Vtmp)
  potential(:) = Vtmp(1,:)
  call fmm_matrix_norm('fmm_V',potential,npoints)

end subroutine fmm_get_boundary_potential

!-------------------------------------------------------------------------------

subroutine fmm_initialise_gfc_grid(npoints,coor)

  implicit none
  integer(INTK), intent(in) :: npoints
  real(REALK), intent(in)   :: coor(3,npoints)

  character(len=10), parameter :: FName = 'multipoles'
  character(len=255) :: FBuf
  integer(INTK) :: i
  integer(INTK), external :: IsFreeUnit

  ! Write grid points to disk
  FBuf = trim(FName)//'.fmm2'
  LUINTM  =IsFreeUnit(LUINTM)
  open(unit=LUINTM,file=trim(FBuf),status='REPLACE',access='SEQUENTIAL',FORM='UNFORMATTED')
  rewind(LUINTM)
  do i=1,npoints
    write(LUINTM) 0,0,0,0,0,coor(1,i),coor(2,i),coor(3,i),one
  enddo
  close(unit=LUINTM,status='KEEP')

  ! Write header file
  FBuf = trim(FName)//'.fmm2header'
  LUINTM = IsFreeUnit(LUINTM)
  open(unit=LUINTM,file=trim(FBuf),status='REPLACE',access='SEQUENTIAL',form='UNFORMATTED')
  write(LUINTM) npoints
  close(unit=LUINTM,status='KEEP')

end subroutine fmm_initialise_gfc_grid

!-------------------------------------------------------------------------------

subroutine fmm_init_basis(NAtom,NShel,NPrim,NBF_Car,MaxAngl,lmax,maxsgm2,KAtom,KType,KStart,KontG,KLoc_Car,Centr,Expnt,CCoef)

  implicit none
  integer(INTK), intent(IN) :: NAtom, NShel, NPrim, NBF_Car
  integer(INTK), intent(IN) :: MaxAngl, lmax, maxsgm2
  integer(INTK), intent(IN) :: KAtom(NShel)
  integer(INTK), intent(IN) :: KType(NShel)
  integer(INTK), intent(IN) :: KStart(NShel)
  integer(INTK), intent(IN) :: KontG(NShel)
  integer(INTK), intent(IN) :: KLoc_Car(NShel)
  real(REALK), intent(IN)   :: Centr(3,NAtom)
  real(REALK), intent(IN)   :: Expnt(NPrim)
  real(REALK), intent(IN)   :: CCoef(NPrim)

  integer(INTK), allocatable :: tmp(:)
  integer(INTK) :: IL, It, Iu, Iv, Ituv, ii
  integer(INTK) :: Maxtuv, MaxLtuv

  Maxtuv = max(MaxAngl,lmax)
  MaxLtuv = (Maxtuv+1)*(Maxtuv+2)*(Maxtuv+3)/6

  allocate(basis%Centr(3,NAtom))
  allocate(basis%KAtom(NShel))
  allocate(basis%KType(NShel))
  allocate(basis%KStart(NShel))
  allocate(basis%KontG(NShel))
  allocate(basis%KLoc_Car(NShel))
  allocate(basis%Expnt(NPrim))
  allocate(basis%CCoef(NPrim))
  allocate(basis%LtuvMin_Car(0:Maxtuv))
  allocate(basis%LtuvMax_Car(0:Maxtuv))
  allocate(basis%Lt(MaxLtuv))
  allocate(basis%Lu(MaxLtuv))
  allocate(basis%Lv(MaxLtuv))

  basis%nshells = NShel
  basis%nbas = NBF_Car
  basis%maxangl = MaxAngl
  basis%maxsgm2 = MaxSgm2
  basis%Centr(:,:) = Centr(:,:)
  basis%KAtom(:) = KAtom(:)
  basis%KType(:) = KType(:)
  basis%KStart(:) = KStart(:)
  basis%KontG(:) = KontG(:)
  basis%KLoc_Car(:) = KLoc_Car(:)
  basis%Expnt(:) = Expnt(:)
  basis%CCoef(:) = CCoef(:)

  ! Now initialise LtuvMin_Car, LtuvMax_Car, Lt, Lu, Lv

  allocate(tmp(0:Maxtuv))
  tmp(0) = 1
  basis%LtuvMin_Car(0) = 1
  basis%LtuvMax_Car(0) = 1
  do IL=1,Maxtuv
    tmp(IL) = (IL+1)*(IL+2)/2
    basis%LtuvMin_Car(IL) = basis%LtuvMin_Car(IL-1)+tmp(IL-1)
    basis%LtuvMax_Car(IL) = basis%LtuvMax_Car(IL-1)+tmp(IL)
  end do

  Ituv = 0
  do IL=0,Maxtuv
    do It=Maxtuv,0,-1
      do Iu=Maxtuv,0,-1
        do Iv=Maxtuv,0,-1
          if (It+Iu+Iv == IL) then
            Ituv = Ituv+1
            basis%Lt(Ituv) = It
            basis%Lu(Ituv) = Iu
            basis%Lv(Ituv) = Iv
          end if
        end do
      end do
    end do
  end do
  deallocate(tmp)

  return

  ! Print section
  write(LUPRI,*) '------------'
  write(LUPRI,*) 'NAtom',NAtom
  write(LUPRI,*) 'NShel',NShel
  write(LUPRI,*) 'NPrim',NPrim
  write(LUPRI,*) 'NBF_Car',NBF_Car
  write(LUPRI,*) '------------'
  write(LUPRI,*) 'Maxangl',MaxAngl
  write(LUPRI,*) 'Maxsgm2',MaxSgm2
  write(LUPRI,*) '------------'
  write(LUPRI,*) 'katom:'
  do ii=1,nshel
    write(LUPRI,*) ii,KAtom(ii)
  end do
  write(LUPRI,*) '------------'
  write(LUPRI,*) 'ktype:'
  do ii=1,nshel
    write(LUPRI,*) ii,KType(ii)
  end do
  write(LUPRI,*) '------------'
  write(LUPRI,*) 'kstart:'
  do ii=1,nshel
    write(LUPRI,*) ii,KStart(ii)
  end do
  write(LUPRI,*) '------------'
  write(LUPRI,*) 'kontg:'
  do ii=1,nshel
    write(LUPRI,*) ii,KontG(ii)
  end do
  write(LUPRI,*) '------------'
  write(LUPRI,*) 'kloc_car:'
  do ii=1,nshel
    write(LUPRI,*) ii,KLoc_Car(ii)
  end do
  write(LUPRI,*) '------------'
  write(LUPRI,*) 'expnt:'
  do ii=1,nprim
    write(LUPRI,*) ii,Expnt(ii)
  end do
  write(LUPRI,*) '------------'
  write(LUPRI,*) 'ccoef:'
  do ii=1,nprim
    write(LUPRI,*) ii,CCoef(ii)
  end do
  write(LUPRI,*) '------------'

end subroutine fmm_init_basis

!-------------------------------------------------------------------------------

subroutine fmm_final()

  use fmm_shell_pairs, only: fmm_free_shell_pairs
  implicit none

  deallocate(basis%Centr)
  deallocate(basis%KAtom)
  deallocate(basis%KType)
  deallocate(basis%KStart)
  deallocate(basis%KontG)
  deallocate(basis%KLoc_Car)
  deallocate(basis%Expnt)
  deallocate(basis%CCoef)
  deallocate(basis%LtuvMin_Car)
  deallocate(basis%LtuvMax_Car)
  deallocate(basis%Lt)
  deallocate(basis%Lu)
  deallocate(basis%Lv)
  call fmm_free_shell_pairs()

end subroutine fmm_final

!-------------------------------------------------------------------------------

end module fmm_interface
