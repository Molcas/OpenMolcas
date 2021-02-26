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

use Basis_Info, only: nBas, nCnttp, dbsc, Shells, MolWgh
use Center_Info, only: dc
use linalg_mod, only: verify_
use Symmetry_Info, only: nIrrep, lIrrep
use sorting, only: sort
use info_expbas_mod, only: DoExpbas, EB_FileOrb
use stdalloc, only: mma_allocate, mma_deallocate
use definitions, only: wp

implicit none
private
public :: desym

integer, parameter :: nNumber = 61
character(len=1), parameter :: number(nNumber) = ['1','2','3','4','5','6','7','8','9','0', &
                                                  'a','b','c','d','e','f','g','h','i','j', &
                                                  'k','l','m','n','o','p','q','r','s','t', &
                                                  'u','v','w','x','y','z','A','B','C','D', &
                                                  'E','F','G','H','I','J','K','L','M','N', &
                                                  'O','P','Q','R','S','T','V','W','X','Y', &
                                                  'Z']
!> Different kinds of orbitals
!> f, i, 1, 2, 3, s, d
integer, parameter :: n_orb_kinds = 7

! NOTE: These global variables are ugly as hell, but we need
!  it to support the shitty SUN compiler.
#ifndef INTERNAL_PROC_ARG
!> This is the orbital kind for each orbital.
integer, allocatable :: kind_per_orb(:)
real(wp), allocatable :: energy(:), occ(:)
real(wp), allocatable :: CMO(:,:)
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
  integer, intent(out) :: ireturn

  integer :: ierr
# include "Molcas.fh"
# include "WrkSpc.fh"
  real(wp), parameter :: EorbThr = 50._wp
  real(wp) :: Coor(3,MxAtom), Znuc(MxAtom)
  character(len=LENIN) :: AtomLabel(MxAtom)
  character(len=512) :: FilesOrb
  character(len=LENIN8), allocatable :: label(:)
  character(len=8) :: MO_Label(maxbfn)
  integer :: ibas_lab(MxAtom), nOrb(8)
# ifdef INTERNAL_PROC_ARG
  !> This is the orbital kind for each orbital.
  integer, allocatable :: kind_per_orb(:)
  real(wp), allocatable :: energy(:), occ(:)
  real(wp), allocatable :: CMO(:,:)
# endif
  !> This is the number of orbitals for every kind.
  integer :: n_kinds(n_orb_kinds)

  character(len=LENIN8+1) :: gtolabel(maxbfn)
  character(len=50) :: VTitle
  character(len=128) :: SymOrbName
  logical :: exists, found
  logical, parameter :: y_cart = .false.

  integer :: nAtom, nData, nDeg, nTot, nTot2, nB
  integer :: iCnttp
  integer :: ipCent, ipCent2, ipCent3
  integer :: ipPhase, ipC2, ipV
  integer :: mAdOcc, mAdEor, mAdCMO
  integer :: file_id, iWfType
  integer, parameter :: notSymm = 1, arbitrary_number = 42, noUHF = 0, iWarn = 1
  integer :: iatom, iDeg, ishell, iIrrep

  integer :: mdc, kk, i, j, ik, k, l, kk_Max, ii, iB, ipp, ic, iv
  integer :: ipc, icontr, nBasisi, icntr

  ireturn = 0

  file_id = arbitrary_number

  call f_Inquire('RUNFILE',exists)
  call verify_(exists,'Error finding RUNFILE')

  !-----Read the characteristics of all different basis sets,
  !     provide each atom with a nuclear charge and establish
  !     a link between an atom and its basis set ---
  !
  !     NOTICE!!!
  !     This call will also fill info.fh and the dynamic storage in
  !     Work(ipInf)

  call Inter1(AtomLabel,iBas_Lab,Coor,Znuc,nAtom)
  call Qpg_iArray('nOrb',found,nData)
  if (found) then
    call Get_iArray('nOrb',nOrb,nData)
  else
    nOrb(:nIrrep) = nBas(:nIrrep)
  end if

  ! Compute memory requirements and allocate memory

  nB = sum(nBas(0:nIrrep-1))
  call GetMem('ICENT','ALLO','INTE',ipCent,8*nB)
  call GetMem('IPHASE','ALLO','INTE',ipPhase,8*nB)
  call GetMem('nCENT','ALLO','INTE',ipCent2,nB)
  call GetMem('ICENTER','ALLO','INTE',ipCent3,nB)
  call GetMem('CMO2','ALLO','REAL',ipC2,nB**2)
  call GetMem('VECTOR','ALLO','REAL',ipV,nB**2)
  call dcopy_(nB**2,[0._wp],0,Work(ipV),1)
  call FZero(Work(ipC2),nB**2)

  ! Read exponents and contraction coefficients of each unique basis.
  ! Write the present basis set (iCnttp) to the molden.input file for
  ! the appropriate atoms.
  ! Moreover, a list is constructed which contains a label for each
  ! GTO (gtolabel). This list follows the MOLDEN order of GTOs.
  ! Later this list will be used in the transformation of sabf (the
  ! symmetry adapted basis functions).

  iatom = 0; mdc = 0; kk = 0
  do iCnttp=1,nCnttp             ! loop over unique basis sets
    if (.not.(dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag)) then
      do iCntr=1,dbsc(iCnttp)%nCntr! loop over symmetry unique centers
        mdc = mdc+1
        nDeg = nIrrep/dc(mdc)%nStab
        do iDeg=1,nDeg             ! loop over centers
          iAtom = iAtom+1

          if (dbsc(iCnttp)%nVal > 6) then
            write(6,*) 'Desym: too high angular momentum!'
            write(6,*) 'iCnttp and dbsc(iCnttp)%nVal= ',iCnttp,dbsc(iCnttp)%nVal
            call Abend()
          end if

          do l=0,dbsc(iCnttp)%nVal-1
            ishell = dbsc(iCnttp)%iVal+l
            nBasisi = Shells(iShell)%nBasis
            !write(6,*) 'nBasisi', Shells(iShell)%nBasis
            if (nBasisi > nNumber) then
              write(6,*) 'Desym: too many contracted functions!'
              write(6,*) 'nBasisi=',Shells(iShell)%nBasis
              call Abend()
            end if

            ! Iterate over each contracted GTO

            if (l == 0) then
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'01s     '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
            end if
            if (l == 1) then
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'02px    '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'02py    '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'02pz    '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
            end if
            if ((l == 2) .and. (.not. y_cart)) then
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'03d02-  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'03d01-  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'03d00   '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'03d01+  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'03d02+  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
            end if
            if ((l == 3) .and. (.not. y_cart)) then
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'04f03-  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'04f02-  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'04f01-  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'04f00   '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'04f01+  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'04f02+  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'04f03+  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
            end if
            if ((l == 4) .and. (.not. y_cart)) then
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'05g04-  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'05g03-  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'05g02-  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'05g01-  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'05g00   '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'05g01+  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'05g02+  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'05g03+  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'05g04+  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
            end if
            if ((l == 5) .and. (.not. y_cart)) then
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'06h05+  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'06h04-  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'06h03-  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'06h02-  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'06h01-  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'06h00   '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'06h01+  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'06h02+  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'06h03+  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'06h04+  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
              do icontr=1,nBasisi
                kk = kk+1
                gtolabel(kk) = AtomLabel(iAtom)//'06h05+  '//number(icontr)
                iWork(ipCent3+kk-1) = iAtom
              end do
            end if
          end do
        end do
      end do
    end if
  end do
  kk_Max = kk
  if (nB > kk_max) then
    write(6,*) 'nB > kk_max'
    write(6,*) 'nB,kk_Max=',nB,kk_Max
    call ClsSew()
    return
  end if

  nTot = sum(nBas(0:nIrrep-1))
  nTot2 = sum(nBas(0:nIrrep-1)**2)

  call mma_allocate(kind_per_orb,nTot)
  call GetMem('Occ','Allo','Real',mAdOcc,nTot)
  call GetMem('Eor','Allo','Real',mAdEor,nTot)
  call GetMem('CMO','Allo','Real',mAdCMO,nTot2)
  call FZero(Work(mAdOcc),nTot)
  call FZero(Work(mAdEor),nTot)
  call FZero(Work(mAdCMO),nTot2)

  !---- Read HF CMOs from file
  FilesOrb = EB_FileOrb
  if (len_trim(FilesOrb) == 0) FilesOrb = 'INPORB'
  if (DoExpbas) FilesOrb = 'EXPORB'
  call RdVec_(trim(FilesOrb),file_id,'COEI',noUHF,nIrrep,nBas,nBas,Work(mAdCMO),Work(ip_Dummy),Work(mAdOcc),Work(ip_Dummy), &
              Work(mAdEor),Work(ip_Dummy),kind_per_orb,VTitle,iWarn,iErr,iWfType)

  if (iErr /= 0) then
    ireturn = 1
    return
  end if

  ! Get the coeff. of sym adapted basis functions (ipC2)

  call Dens_IF_SCF(Work(ipC2),Work(mAdCMO),'F')
  call GetMem('CMO','Free','Real',mAdCMO,nTot2)

  !  Back 'transformation' of the symmetry adapted basis functions.
  !  Probably somewhat clumsy, but it seems to work.If someone
  !  knows a more elegant way to do it, please improve this part!
  !
  !  PART 1: Obtain symmetry information (soout), construct a label
  !          for each sabf, which will be used in part 2 to find the
  !          corresponding GTO in the MOLDEN list by comparing with
  !          gtolabel
  !
  !  nB       --- Total number of contracted basis functions
  !  ipcent2  --- degeneracy of a basis function
  !  ipCent   --- centres over which the basis function is
  !               delocalized
  !  ipPhase  --- phase of the AO in the linear combination

  call mma_allocate(Label,MaxBfn+MaxBfn_Aux,label='Label')
  call icopy(8*nB,[0],0,iWork(ipPhase),1)
  call icopy(8*nB,[0],0,iWork(ipCent),1)
  call SOout(label,iWork(ipCent),iWork(ipPhase))
  ipc = 0
  do iContr=1,nB
    iWork(ipCent2+iContr-1) = 0
    do k=1,8
      if (iWork(ipCent+ipc) /= 0) then
        iWork(ipcent2+iContr-1) = iWork(ipCent2+iContr-1)+1
      end if
      ipc = ipc+1
    end do
  end do
  ! Part 2: -Take a MOLCAS symmetry functions (loop i)
  !         -Find the corresponding label in the MOLDEN list (loop j)
  !         -Copy the coeff of the sabf in the MOLDEN MO (vector), multiply
  !          by the appropriate factor (ipPhase),and divide by the number of
  !          centres over which the sabf is delocalized (ipCent3).
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
      write(MO_Label(i),'(I5,A3)') iB,lirrep(iIrrep)

      do j=1,nB

        if (gtolabel(j) == label(i)//number(ik)) then
          do k=1,8
            ipc = (i-1)*8+k-1
            ipp = ipc
            if (iWork(ipCent+ipc) == iWork(ipcent3+j-1)) then
              do ii=1,nB
                ic = (ii-1)*nB+(i-1)
                iv = (ii-1)*nB+(j-1)
                if (MolWgh == 0) then
                  Work(ipV+iv) = Work(ipV+iv)+Work(ipC2+ic)*dble(iWork(ipPhase+ipp))/dble(iWork(ipcent2+i-1))
                else
                  Work(ipV+iv) = Work(ipV+iv)+Work(ipC2+ic)*dble(iWork(ipPhase+ipp))/sqrt(dble(iWork(ipcent2+i-1)))
                end if
              end do
            end if
          end do
        end if
      end do
    end do
  end do
  call mma_deallocate(Label)

  !**************************** START SORTING ****************************

  call mma_allocate(CMO,nTot,nTot)
  call mma_allocate(occ,nTot)
  call mma_allocate(energy,nTot)

  energy(:) = Work(mAdEor:mAdEor+nTot-1)
  occ(:) = Work(mAdOcc:mAdocc+nTot-1)

  do i=0,nTot-1
    l = ipV+nTot*i
    CMO(:,i+1) = work(l:l+nTot-1)
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
  call WrVec_(SymOrbName,file_id,'COEI',noUHF,notSymm,[nTot],[nTot],CMO,Work(ip_Dummy),occ,Work(ip_Dummy),energy,Work(ip_Dummy), &
              n_kinds,VTitle,iWFtype)
  call Add_Info('desym CMO',CMO,999,8)

  call mma_deallocate(occ)
  call mma_deallocate(CMO)
  call mma_deallocate(energy)
  call mma_deallocate(kind_per_orb)
  call GetMem('Eor','Free','Real',mAdEor,nTot)
  call GetMem('Occ','Free','Real',mAdOcc,nTot)
  call GetMem('ICENT','FREE','INTE',ipCent,8*nB)
  call GetMem('IPHASE','FREE','INTE',ipPhase,8*nB)
  call GetMem('nCENT','FREE','INTE',ipCent2,nB)
  call GetMem('ICENTER','FREE','INTE',ipCent3,nB)
  call GetMem('CMO2','FREE','REAL',ipC2,nB**2)
  call GetMem('VECTOR','FREE','REAL',ipV,nB**2)

  call ClsSew()

end subroutine desym

#ifdef INTERNAL_PROC_ARG

subroutine reorder_orbitals(nTot,kind_per_orb,CMO,occ,energy)
  integer, intent(in) :: nTot
  integer, intent(inout) :: kind_per_orb(nTot)
  real(wp), intent(inout) :: CMO(nTot,nTot), occ(nTot), energy(nTot)

  integer :: i
  integer, allocatable :: idx(:)

  allocate(idx(nTot))

  idx(:) = [(i,i=1,nTot)]

  call sort(idx,closure_compare)

  kind_per_orb(:) = kind_per_orb(idx)
  energy(:) = energy(idx)
  CMO(:,:) = CMO(:,idx)
  occ(:) = occ(idx)

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
    integer, intent(in) :: i, j
    logical :: res

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
  end function
end subroutine reorder_orbitals

#else

subroutine reorder_orbitals()
  integer :: nTot, i
  integer, allocatable :: idx(:)

  nTot = size(kind_per_orb)
  allocate(idx(nTot))
  idx(:) = [(i,i=1,nTot)]

  call sort(idx,compare)

  kind_per_orb(:) = kind_per_orb(idx)
  energy(:) = energy(idx)
  CMO(:,:) = CMO(:,idx)
  occ(:) = occ(idx)

# ifdef _WARNING_WORKAROUND_
  ! There was a wrong and stupid warning from GFortran 4.8
  ! If we stop supporting this compiler remove this ugly workaround.
  if (.false.) then
    if (compare(1,2)) continue
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
  integer, intent(in) :: i, j
  logical :: res

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
