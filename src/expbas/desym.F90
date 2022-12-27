!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2009, Giovanni Li Manni                                *
!               2009, Francesco Aquilante                              *
!               2020, Oskar Weser                                      *
!***********************************************************************

module desymmetrize_mod

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
private
public :: desym

integer(kind=iwp), parameter :: nNumber = 61
character, parameter :: number(nNumber) = ['1','2','3','4','5','6','7','8','9','0', &
                                           'a','b','c','d','e','f','g','h','i','j', &
                                           'k','l','m','n','o','p','q','r','s','t', &
                                           'u','v','w','x','y','z','A','B','C','D', &
                                           'E','F','G','H','I','J','K','L','M','N', &
                                           'O','P','Q','R','S','T','V','W','X','Y', &
                                           'Z']

! NOTE: These global variables are ugly as hell, but we need
!  it to support the shitty SUN compiler.
#ifndef INTERNAL_PROC_ARG
!> This is the orbital kind for each orbital.
integer(kind=iwp), allocatable :: kind_per_orb(:)
real(kind=wp), allocatable :: energy(:), occ(:)
real(kind=wp), allocatable :: CMO(:,:)
#endif

contains

! symmetry-----> C1 INPORB
subroutine desym(ireturn)
!***********************************************************************
!                                                                      *
! Purpose: Convert INPORB file with symmetry to DeSymOrb               *
!          file with C1 symmetry.                                      *
!                                                                      *
!          G. Li Manni, F. Aquilante, University of Geneva, June 2009. *
!          Oskar Weser, MPI FKF Stuttgart,  November 2020.             *
!                                                                      *
! Example of input:                                                    *
!                                                                      *
!   &SEWARD                                                            *
!    coord                                                             *
!     file.xyz                                                         *
!    basis                                                             *
!     ano-rcc-mb                                                       *
!    oneonly                                                           *
!    nodk                                                              *
!                                                                      *
!   >>COPY $CurrDir/scf_sym.ScfOrb INPORB                              *
!                                                                      *
!   &EXPBAS                                                            *
!    NoEx                                                              *
!    Desy                                                              *
!   END                                                                *
!                                                                      *
!***********************************************************************

  use Basis_Info, only: nBas, nCnttp, dbsc, Shells, MolWgh
  use Center_Info, only: dc
  use linalg_mod, only: verify_
  use Symmetry_Info, only: nIrrep
  use info_expbas_mod, only: DoExpbas, EB_FileOrb, LenIn, MxAtom, mxsym, n_orb_kinds
  use Constants, only: Zero

  integer(kind=iwp), intent(out) :: ireturn

  integer(kind=iwp) :: iErr, nOrb(mxsym)
  real(kind=wp), parameter :: EorbThr = 50.0_wp
  character(len=512) :: FilesOrb
# ifdef INTERNAL_PROC_ARG
  !> This is the orbital kind for each orbital.
  integer(kind=iwp), allocatable :: kind_per_orb(:)
  real(kind=wp), allocatable :: energy(:), occ(:), CMO(:,:)
# endif
  !> This is the number of orbitals for every kind.
  integer(kind=iwp) :: n_kinds(n_orb_kinds)

  character(len=50) :: VTitle
  character(len=128) :: SymOrbName
  logical(kind=iwp) :: exists, found
  logical(kind=iwp), parameter :: y_cart = .false.

  integer(kind=iwp) :: nAtom, nData, nDeg, nTot, nTot2, nB, iCnttp, file_id, iWfType, iatom, iDeg, ishell, iIrrep, mdc, kk, i, j, &
                       ik, k, l, kk_Max, ii, iB, ipp, ic, iv, ipc, icontr, nBasisi, icntr
  real(kind=wp) :: dummy(1)
  integer(kind=iwp), allocatable :: Cent(:), Center(:), iBas_Lab(:), nCent(:), Phase(:)
  real(kind=wp), allocatable :: AdCMO(:), AdEor(:), AdOcc(:), CMO2(:), Coor(:,:), Vector(:), Znuc(:)
  character(len=LenIn), allocatable :: AtomLabel(:)
  character(len=LenIn+8), allocatable :: Label(:)
  character(len=LenIn+9), allocatable :: gtolabel(:)
  integer(kind=iwp), parameter :: notSymm = 1, arbitrary_number = 42, noUHF = 0, iWarn = 1
  character(len=*), parameter :: baslab_0(1) = ['01s     '], &
                                 baslab_1(3) = ['02px    ','02py    ','02pz    '], &
                                 baslab_2(5) = ['03d02-  ','03d01-  ','03d00   ','03d01+  ','03d02+  '], &
                                 baslab_3(7) = ['04f03-  ','04f02-  ','04f01-  ','04f00   ','04f01+  ','04f02+  ','04f03+  '], &
                                 baslab_4(9) = ['05g04-  ','05g03-  ','05g02-  ','05g01-  ','05g00   ', &
                                                '05g01+  ','05g02+  ','05g03+  ','05g04+  '], &
                                 baslab_5(11) = ['06h05-  ','06004-  ','06h03-  ','06h02-  ','06h01-  ','06h00   ', &
                                                 '06h01+  ','06h02+  ','06h03+  ','06h04+  ','06h05+  ']

  ireturn = 0

  file_id = arbitrary_number

  call f_Inquire('RUNFILE',exists)
  call verify_(exists,'Error finding RUNFILE')

  ! Read the characteristics of all different basis sets,
  ! provide each atom with a nuclear charge and establish
  ! a link between an atom and its basis set ---
  !
  ! NOTICE!!!
  ! This call will also fill Basis_Info and Center_Info

  call mma_allocate(AtomLabel,MxAtom,label='AtomLabel')
  call mma_allocate(iBas_Lab,MxAtom,label='iBas_Lab')
  call mma_allocate(Coor,3,MxAtom,label='Coor')
  call mma_allocate(Znuc,MxAtom,label='Znuc')
  call Inter1(AtomLabel,iBas_Lab,Coor,Znuc,nAtom)
  call mma_deallocate(iBas_Lab)
  call mma_deallocate(Coor)
  call mma_deallocate(Znuc)
  call Qpg_iArray('nOrb',found,nData)
  if (found) then
    call Get_iArray('nOrb',nOrb,nData)
  else
    nOrb(:nIrrep) = nBas(:nIrrep)
  end if

  ! Compute memory requirements and allocate memory

  nB = sum(nBas(0:nIrrep-1))
  call mma_allocate(Cent,8*nB,label='Cent')
  call mma_allocate(Phase,8*nB,label='Phase')
  call mma_allocate(nCent,nB,label='nCent')
  call mma_allocate(Center,nB,label='Center')
  call mma_allocate(CMO2,nB**2,label='CMO2')
  call mma_allocate(Vector,nB**2,label='Vector')
  call mma_allocate(gtolabel,nB,label='gtolabel')
  CMO2(:) = Zero
  Vector(:) = Zero

  ! Read exponents and contraction coefficients of each unique basis.
  ! Write the present basis set (iCnttp) to the molden.input file for
  ! the appropriate atoms.
  ! Moreover, a list is constructed which contains a label for each
  ! GTO (gtolabel). This list follows the MOLDEN order of GTOs.
  ! Later this list will be used in the transformation of sabf (the
  ! symmetry adapted basis functions).

  iatom = 0
  mdc = 0
  kk = 0
  do iCnttp=1,nCnttp                ! loop over unique basis sets
    if (.not. (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag)) then
      do iCntr=1,dbsc(iCnttp)%nCntr ! loop over symmetry unique centers
        mdc = mdc+1
        nDeg = nIrrep/dc(mdc)%nStab
        do iDeg=1,nDeg              ! loop over centers
          iAtom = iAtom+1

          if (dbsc(iCnttp)%nVal > 6) then
            write(u6,*) 'Desym: too high angular momentum!'
            write(u6,*) 'iCnttp and dbsc(iCnttp)%nVal= ',iCnttp,dbsc(iCnttp)%nVal
            call Abend()
          end if

          do l=0,dbsc(iCnttp)%nVal-1
            ishell = dbsc(iCnttp)%iVal+l
            nBasisi = Shells(iShell)%nBasis
            !write(u6,*) 'nBasisi', Shells(iShell)%nBasis
            if (nBasisi > nNumber) then
              write(u6,*) 'Desym: too many contracted functions!'
              write(u6,*) 'nBasisi=',Shells(iShell)%nBasis
              call Abend()
            end if

            ! Iterate over each contracted GTO

            select case (l)
              case (0)
                do i=1,size(baslab_0)
                  do icontr=1,nBasisi
                    kk = kk+1
                    gtolabel(kk) = AtomLabel(iAtom)//baslab_0(i)//number(icontr)
                    Center(kk) = iAtom
                  end do
                end do
              case (1)
                do i=1,size(baslab_1)
                  do icontr=1,nBasisi
                    kk = kk+1
                    gtolabel(kk) = AtomLabel(iAtom)//baslab_1(i)//number(icontr)
                    Center(kk) = iAtom
                  end do
                end do
              case (2)
                if (.not. y_cart) then
                  do i=1,size(baslab_2)
                    do icontr=1,nBasisi
                      kk = kk+1
                      gtolabel(kk) = AtomLabel(iAtom)//baslab_2(i)//number(icontr)
                      Center(kk) = iAtom
                    end do
                  end do
                end if
              case (3)
                if (.not. y_cart) then
                  do i=1,size(baslab_3)
                    do icontr=1,nBasisi
                      kk = kk+1
                      gtolabel(kk) = AtomLabel(iAtom)//baslab_3(i)//number(icontr)
                      Center(kk) = iAtom
                    end do
                  end do
                end if
              case (4)
                if (.not. y_cart) then
                  do i=1,size(baslab_4)
                    do icontr=1,nBasisi
                      kk = kk+1
                      gtolabel(kk) = AtomLabel(iAtom)//baslab_4(i)//number(icontr)
                      Center(kk) = iAtom
                    end do
                  end do
                end if
              case (5)
                if (.not. y_cart) then
                  do i=1,size(baslab_5)
                    do icontr=1,nBasisi
                      kk = kk+1
                      gtolabel(kk) = AtomLabel(iAtom)//baslab_5(i)//number(icontr)
                      Center(kk) = iAtom
                    end do
                  end do
                end if
            end select
          end do
        end do
      end do
    end if
  end do
  kk_Max = kk
  if (nB > kk_max) then
    write(u6,*) 'nB > kk_max'
    write(u6,*) 'nB,kk_Max=',nB,kk_Max
    call ClsSew()
    return
  end if

  call mma_deallocate(AtomLabel)

  nTot = sum(nBas(0:nIrrep-1))
  nTot2 = sum(nBas(0:nIrrep-1)**2)

  call mma_allocate(kind_per_orb,nTot)
  call mma_allocate(AdOcc,nTot,label='AdOcc')
  call mma_allocate(AdEor,nTot,label='AdEor')
  call mma_allocate(AdCMO,nTot2,label='AdCMO')
  AdOcc(:) = Zero
  AdEor(:) = Zero
  AdCMO(:) = Zero

  ! Read HF CMOs from file
  FilesOrb = EB_FileOrb
  if (len_trim(FilesOrb) == 0) FilesOrb = 'INPORB'
  if (DoExpbas) FilesOrb = 'EXPORB'
  call RdVec_(trim(FilesOrb),file_id,'COEI',noUHF,nIrrep,nBas,nBas,AdCMO,dummy,AdOcc,dummy,AdEor,dummy,kind_per_orb,VTitle,iWarn, &
              iErr,iWfType)

  if (iErr /= 0) then
    ireturn = 1
    return
  end if

  ! Get the coeff. of sym adapted basis functions (CMO2)

  call Dens_IF_SCF(CMO2,AdCMO,'F')
  call mma_deallocate(AdCMO)

  ! Back 'transformation' of the symmetry adapted basis functions.
  ! Probably somewhat clumsy, but it seems to work. If someone
  ! knows a more elegant way to do it, please improve this part!
  !
  ! PART 1: Obtain symmetry information (soout), construct a label
  !         for each sabf, which will be used in part 2 to find the
  !         corresponding GTO in the MOLDEN list by comparing with
  !         gtolabel
  !
  ! nB     --- Total number of contracted basis functions
  ! nCent  --- degeneracy of a basis function
  ! Cent   --- centres over which the basis function is delocalized
  ! Phase  --- phase of the AO in the linear combination

  call mma_allocate(Label,nB,label='Label')
  Phase(:) = 0
  Cent(:) = 0
  call SOout(label,Cent,Phase)
  ipc = 1
  do iContr=1,nB
    nCent(iContr) = 0
    do k=1,8
      if (Cent(ipc) /= 0) nCent(iContr) = nCent(iContr)+1
      ipc = ipc+1
    end do
  end do
  ! Part 2: -Take a MOLCAS symmetry functions (loop i)
  !         -Find the corresponding label in the MOLDEN list (loop j)
  !         -Copy the coeff of the sabf in the MOLDEN MO (vector), multiply
  !          by the appropriate factor (Phase),and divide by the number of
  !          centres over which the sabf is delocalized (Center).
  !         -The vectors are copied by rows!
  !
  ! loop over MOLCAS symmetry functions
  i = 0
  ik = 0
  do iIrrep=0,nIrrep-1

    do iB=1,nBas(iIrrep)
      i = i+1

      if (iB == 1) then
        ik = 1
      else
        if (label(i-1) == label(i)) then
          ik = ik+1
        else
          ik = 1
        end if
      end if

      do j=1,nB

        if (gtolabel(j) == label(i)//number(ik)) then
          do k=1,8
            ipc = (i-1)*8+k
            ipp = ipc
            if (Cent(ipc) == Center(j)) then
              do ii=1,nB
                ic = (ii-1)*nB+i
                iv = (ii-1)*nB+j
                if (MolWgh == 0) then
                  Vector(iv) = Vector(iv)+CMO2(ic)*real(Phase(ipp),kind=wp)/real(nCent(i),kind=wp)
                else
                  Vector(iv) = Vector(iv)+CMO2(ic)*real(Phase(ipp),kind=wp)/sqrt(real(nCent(i),kind=wp))
                end if
              end do
            end if
          end do
        end if
      end do
    end do
  end do
  call mma_deallocate(Label)
  call mma_deallocate(gtolabel)

  !**************************** START SORTING ****************************

  call mma_allocate(CMO,nTot,nTot)
  call mma_allocate(occ,nTot)
  call mma_allocate(energy,nTot)

  energy(:) = AdEor(:)
  occ(:) = AdOcc(:)

  do i=1,nTot
    CMO(:,i) = Vector(nTot*(i-1)+1:nTot*i)
  end do

# ifdef INTERNAL_PROC_ARG
  call reorder_orbitals(nTot,kind_per_orb,CMO,occ,energy)
# else
  call reorder_orbitals()
# endif

  n_kinds(:) = 0
  do i=lbound(kind_per_orb,1),ubound(kind_per_orb,1)
    n_kinds(kind_per_orb(i)) = n_kinds(kind_per_orb(i))+1
  end do

  SymOrbName = 'DESORB'
  VTitle = 'Basis set desymmetrized orbital file DESORB'
  call WrVec_(SymOrbName,file_id,'COEI',noUHF,notSymm,[nTot],[nTot],CMO,dummy,occ,dummy,energy,dummy,n_kinds,VTitle,iWFtype)
  call Add_Info('desym CMO',CMO,999,8)

  call mma_deallocate(occ)
  call mma_deallocate(CMO)
  call mma_deallocate(energy)
  call mma_deallocate(kind_per_orb)
  call mma_deallocate(AdOcc)
  call mma_deallocate(AdEor)
  call mma_deallocate(Cent)
  call mma_deallocate(Phase)
  call mma_deallocate(nCent)
  call mma_deallocate(Center)
  call mma_deallocate(CMO2)
  call mma_deallocate(Vector)

  call ClsSew()

end subroutine desym

#ifdef INTERNAL_PROC_ARG

subroutine reorder_orbitals(nTot,kind_per_orb,CMO,occ,energy)

  use sorting, only: sort

  integer(kind=iwp), intent(in) :: nTot
  integer(kind=iwp), intent(inout) :: kind_per_orb(nTot)
  real(kind=wp), intent(inout) :: CMO(nTot,nTot), occ(nTot), energy(nTot)

  integer(kind=iwp) :: i
  integer(kind=iwp), allocatable :: idx(:)

  call mma_allocate(idx,nTot,label='idx')

  idx(:) = [(i,i=1,nTot)]

  call sort(idx,closure_compare)

  kind_per_orb(:) = kind_per_orb(idx)
  energy(:) = energy(idx)
  CMO(:,:) = CMO(:,idx)
  occ(:) = occ(idx)

  call mma_deallocate(idx)

contains

  !> @brief
  !>  Comparison function for generic sort.
  !>
  !> @details
  !>  Sort non-strict i.e. (compare(i, j) .and. compare(j, i)) can be true)
  !>  Sort first by orbital kind ascendingly (frozen, inactive, RAS1, ...),
  !>      second by occupation number descendingly (2.0, 2.0, 1.x, 0., ...),
  !>      third by energy ascendingly (-3., -2., -2., 0., 1., ...).
  !>  Note, that `sort` uses a stable sorting algorithm and since the
  !>  input orbitals are automatically by irrep the output
  !>  will be sorted last ascendingly by irrep.
  pure function closure_compare(i,j) result(res)
    logical(kind=iwp) :: res
    integer(kind=iwp), intent(in) :: i, j

    if (kind_per_orb(i) /= kind_per_orb(j)) then
      res = kind_per_orb(i) < kind_per_orb(j)
    else if (occ(i) /= occ(j)) then
      res = occ(i) > occ(j)
    else if (energy(i) /= energy(j)) then
      res = energy(i) < energy(j)
    else
      ! All relevant values are equal and our comparison has
      ! to be non-strict.
      res = .true.
    end if
  end function closure_compare

end subroutine reorder_orbitals

#else

subroutine reorder_orbitals()

  use sorting, only: sort

  integer(kind=iwp) :: nTot, i
  integer(kind=iwp), allocatable :: idx(:)

  nTot = size(kind_per_orb)
  call mma_allocate(idx,nTot,label='idx')
  idx(:) = [(i,i=1,nTot)]

  call sort(idx,compare)

  kind_per_orb(:) = kind_per_orb(idx)
  energy(:) = energy(idx)
  CMO(:,:) = CMO(:,idx)
  occ(:) = occ(idx)

  call mma_deallocate(idx)

# ifdef _WARNING_WORKAROUND_
  return
  ! There was a wrong and stupid warning from GFortran 4.8
  ! If we stop supporting this compiler remove this ugly workaround.
  if (.false.) then
    if (compare(1,2)) return
  end if
# endif

end subroutine reorder_orbitals

!> @brief
!>  Comparison function for generic sort.
!>
!> @details
!>  Sort non-strict i.e. (compare(i, j) .and. compare(j, i)) can be true)
!>  Sort first by orbital kind ascendingly (frozen, inactive, RAS1, ...),
!>      second by occupation number descendingly (2.0, 2.0, 1.x, 0., ...),
!>      third by energy ascendingly (-3., -2., -2., 0., 1., ...).
!>  Note, that `sort` uses a stable sorting algorithm and since the
!>  input orbitals are automatically by irrep the output
!>  will be sorted last ascendingly by irrep.
pure function compare(i,j) result(res)

  logical(kind=iwp) :: res
  integer(kind=iwp), intent(in) :: i, j

  if (kind_per_orb(i) /= kind_per_orb(j)) then
    res = kind_per_orb(i) < kind_per_orb(j)
  else if (occ(i) /= occ(j)) then
    res = occ(i) > occ(j)
  else if (energy(i) /= energy(j)) then
    res = energy(i) < energy(j)
  else
    ! All relevant values are equal and our comparison has
    ! to be non-strict.
    res = .true.
  end if

end function compare

#endif

end module desymmetrize_mod
