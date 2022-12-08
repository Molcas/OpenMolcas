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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************

#include "compiler_features.h"

module Basis_Info

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
private

public :: Basis_Info_Dmp, Basis_Info_Free, Basis_Info_Get, Basis_Info_Init, dbsc, Distinct_Basis_set_Centers, Gaussian_Type, &
          iCnttp_Dummy, Max_Shells, mGaussian_Type, MolWgh, nBas, nBas_Aux, nBas_Frag, nCnttp, nFrag_LineWords, Nuclear_Model, &
          PAMExp, Point_Charge, Shells

#include "Molcas.fh"
#include "itmax.fh"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! D E C L A R E   D E R I V E D   T Y P E S
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Work in progress
!
! nCntr     : number of centers associated with a dbsc
! Coor      : the coordinates of a dbsc
!
! nM1       : number of ECP M1 type terms on the ith unique center
! M1xp      : ECP M1-type exponents for ith unique center
! M1cf      : ECP M1 type coefficients for ith unique cntr
! nM2       : number of ECP M2 type terms on the ith unique center
! M2xp      : ECP M2-type exponents for ith unique center
! M2cf      : ECP M2 type coefficients for ith unique cntr
!
! nFragType : number of unique centers in a fragment (0=not a frag)
! nFragCoor : total number of centers in a fragment
! nFragEner : number of orbital energies/occupied orbitals in a given fragment
! nFragDens : size of the fragments density matrix
! FragType  : the data of the fragment's unique centers (associated basis set, size nFragType)
! FragCoor  : the data of all fragment's centers (atom type / relative coordinates, size nFragCoor)
! FragEner  : the fragment's orbital's energies (size nFragEner)
! FragCoef  : the fragment's MO coefficients (size nFragDens*nFragEner)
! ECP       : Flag if dbsc is a ECP basis set
! Frag      : Flag if dbsc is a Fragment basis set
! Aux       : Flag if dbsc is an auxiliary basis set
! FOp       : Flag if dbsc has a Fock Operator
! IsMM      : integer flag to indicate that associated centers are treated as MM centers in QM/MM calculations

type Distinct_Basis_set_centers
  real(kind=wp), pointer :: Coor(:,:) => null()
  real(kind=wp), allocatable :: Coor_Hidden(:,:)
  integer(kind=iwp) :: nCntr = 0
  integer(kind=iwp) :: nM1 = 0
  real(kind=wp), allocatable :: M1xp(:), M1cf(:)
  integer(kind=iwp) :: nM2 = 0
  real(kind=wp), allocatable :: M2xp(:), M2cf(:)
  integer(kind=iwp) :: nFragType = 0, nFragCoor = 0, nFragEner = 0, nFragDens = 0
  real(kind=wp), allocatable :: FragType(:,:), FragCoor(:,:), FragEner(:), FragCoef(:,:)
  logical(kind=iwp) :: lPAM2 = .false.
  integer(kind=iwp) :: nPAM2 = -1
  real(kind=wp), allocatable :: PAM2(:)
  logical(kind=iwp) :: ECP = .false.
  logical(kind=iwp) :: Aux = .false.
  logical(kind=iwp) :: Frag = .false.
  logical(kind=iwp) :: FOp = .false.
  integer(kind=iwp) :: IsMM = 0
  integer(kind=iwp) :: Parent_iCnttp = 0
  integer(kind=iwp) :: lOffAO = 0
  integer(kind=iwp) :: nOpt = 0
  integer(kind=iwp) :: mdci = 0
  integer(kind=iwp) :: iVal = 0, nVal = 0
  integer(kind=iwp) :: iPrj = 0, nPrj = 0
  integer(kind=iwp) :: iSRO = 0, nSRO = 0
  integer(kind=iwp) :: iSOC = 0, nSOC = 0
  integer(kind=iwp) :: kDel(0:iTabMx)
  integer(kind=iwp) :: iPP = 0, nPP = 0
  integer(kind=iwp) :: nShells = 0
  integer(kind=iwp) :: AtmNr = 0
  real(kind=wp) :: Charge = Zero
  logical(kind=iwp) :: NoPair = .false.
  logical(kind=iwp) :: SODK = .false.
  logical(kind=iwp) :: pChrg = .false.
  logical(kind=iwp) :: Fixed = .false.
  real(kind=wp) :: CrRep = Zero
  real(kind=wp) :: FragCharge = Zero
  real(kind=wp) :: aCD_Thr = One
  real(kind=wp) :: fmass = One
  real(kind=wp) :: CntMass = Zero
  real(kind=wp) :: ExpNuc = -One
  real(kind=wp) :: w_mGauss = One
# ifdef CHAR_MEMBER_INIT
  ! Some GCC versions give spurious warning with initialization of character members
  ! if there is more than one allocatable member...
  character(len=80) :: Bsl = '', Bsl_Old = ''
# else
  character(len=80) :: Bsl, Bsl_Old
# endif
end type Distinct_Basis_set_centers

! nExp   : number of exponents of the ith shell
! Exp    : the exponents of the ith shell
! nBasis : number of contracted radial functions of the ith shell
! Cff_c  : Contraction coefficients in processed and raw input form
! Cff_p  : Contraction coefficient in the case of no contraction, processed and raw
! Cff    : copy of Cff_c or Cff_p
! Transf : Cartesian transformed to real sphericals.
! Projct : real sphericals without contaminations (3s, 4d, etc.)
! Bk     : ECP proj shift parameters for ith shell, the number of parameters is given by nBasis
! Occ    : Occupation numbers for core ECP orbitals
! FockOp : the Fock operator matrix
! Aux    : Logical flag for auxiliary basis set shells
! Frag   : Logical flag for fragment shells

type Shell_Info
  integer(kind=iwp) :: nExp = 0
  real(kind=wp), allocatable :: Exp(:)
  integer(kind=iwp) :: nBasis = 0
  integer(kind=iwp) :: nBasis_c = 0
  real(kind=wp), allocatable :: pCff(:,:)
  real(kind=wp), allocatable :: Cff_c(:,:,:), Cff_p(:,:,:)
  logical(kind=iwp) :: Transf = .true.
  logical(kind=iwp) :: Prjct = .true.
  integer(kind=iwp) :: nBk = 0
  real(kind=wp), allocatable :: Bk(:)
  real(kind=wp), allocatable :: Occ(:)
  integer(kind=iwp) :: nAkl = 0
  real(kind=wp), allocatable :: Akl(:,:,:)
  integer(kind=iwp) :: nFockOp = 0
  real(kind=wp), allocatable :: FockOp(:,:)
  logical(kind=iwp) :: Aux = .false.
  logical(kind=iwp) :: Frag = .false.
  integer(kind=iwp) :: kOffAO = 0
end type Shell_Info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! E N D   D E C L A R E   D E R I V E D   T Y P E S
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Actual content of Basis_Info

! MolWgh: integer flag to indicate the normalization of the symmetry transformation
!         0: double coset representative normalization
!         1: as in MOLECULE
!         2: as in MOLPRO

integer(kind=iwp), parameter :: Point_Charge = 0, Gaussian_Type = 1, mGaussian_Type = 2

real(kind=wp), allocatable :: PAMexp(:,:)
integer(kind=iwp) :: iCnttp_Dummy = 0, Max_Shells = 0, mFields = 11, MolWgh = 2, nBas(0:7) = 0, nBas_Aux(0:7) = 0, &
                     nBas_Frag(0:7) = 0, nCnttp = 0, nFields = 33+(1+iTabMx), nFrag_LineWords = 0, Nuclear_Model = Point_Charge
logical(kind=iwp) :: Initiated = .false.

type(Distinct_Basis_set_centers), allocatable, target :: dbsc(:)
type(Shell_Info), allocatable :: Shells(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Private extensions to mma interfaces

interface cptr2loff
  module procedure dbsc_cptr2loff
  module procedure shell_cptr2loff
end interface
interface mma_Allocate
  module procedure dbsc_mma_allo_1D, dbsc_mma_allo_1D_lim
  module procedure shell_mma_allo_1D, shell_mma_allo_1D_lim
end interface
interface mma_Deallocate
  module procedure dbsc_mma_free_1D
  module procedure shell_mma_free_1D
end interface

contains

!***********************************************************************
!***********************************************************************
!
! This to make either the initial allocation of dbsc and Shells according to the default sizes
! as defined by the parameters in Molcas.fh or according to the actual sizes as recorded on the
! run file.

subroutine Basis_Info_Init()

  use Definitions, only: u6

# include "macros.fh"
  unused_proc(mma_allocate(dbsc,[0,0]))
  unused_proc(mma_allocate(Shells,[0,0]))

# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'Enter Basis_Info_Init'
  write(u6,*)
# endif
  if (Initiated) then
    write(u6,*) ' Basis_Info already initiated!'
    write(u6,*) ' Maybe there is missing a Basis_Info_Free call.'
    call Abend()
  end if
  if (nCnttp == 0) then
    call mma_allocate(dbsc,Mxdbsc,label='dbsc')
  else
    call mma_allocate(dbsc,nCnttp,label='dbsc')
  end if
# ifndef CHAR_MEMBER_INIT
  ! See above
  dbsc(:)%Bsl = ''
  dbsc(:)%Bsl_Old = ''
# endif
  if (Max_Shells == 0) then
    call mma_allocate(Shells,MxShll,label='Shells')
  else
    call mma_allocate(Shells,Max_Shells,label='Shells')
  end if
  Initiated = .true.
# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'Exit Basis_Info_Init'
  write(u6,*)
# endif

  return

end subroutine Basis_Info_Init

!***********************************************************************
!***********************************************************************

subroutine Basis_Info_Dmp()

  use stdalloc, only: mma_allocate, mma_deallocate
  use Definitions, only: u6

  integer(kind=iwp) :: i, j, lcDmp, nAkl, nAtoms, nAux, nAux2, nBasis, nBk, nExp, nFockOp, nFragCoor, nM1, nM2
  integer(kind=iwp), allocatable :: iDmp(:,:)
  real(kind=wp), allocatable, target :: rDmp(:,:)
  real(kind=wp), pointer :: qDmp(:,:)
  character(len=2*80), allocatable :: cDmp(:)

# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'Enter Basis_Info_Dmp'
  write(u6,*)
  write(u6,*) 'Coordinates:'
  do i=1,nCnttp
    do j=1,dbsc(i)%nCntr
      write(u6,*) dbsc(i)%Coor(:,j)
    end do
  end do
# endif

  ! Integer dbsc stuff

  call mma_Allocate(iDmp,nFields,nCnttp+1,Label='iDmp')
  nAtoms = 0
  nAux = 0
  do i=1,nCnttp
    iDmp(1,i) = dbsc(i)%nCntr
    iDmp(2,i) = dbsc(i)%nM1
    iDmp(3,i) = dbsc(i)%nM2
    iDmp(4,i) = dbsc(i)%nFragType
    iDmp(5,i) = dbsc(i)%nFragCoor
    iDmp(6,i) = dbsc(i)%nFragEner
    iDmp(7,i) = dbsc(i)%nFragDens
    iDmp(8,i) = 0
    if (dbsc(i)%ECP) iDmp(8,i) = 1
    iDmp(9,i) = 0
    if (dbsc(i)%Frag) iDmp(9,i) = 1
    iDmp(10,i) = 0
    if (dbsc(i)%Aux) iDmp(10,i) = 1
    iDmp(11,i) = 0
    if (dbsc(i)%FOp) iDmp(11,i) = 1
    iDmp(12,i) = dbsc(i)%IsMM
    iDmp(13,i) = dbsc(i)%Parent_iCnttp
    iDmp(14,i) = dbsc(i)%lOffAO
    iDmp(15,i) = dbsc(i)%nOpt
    iDmp(16,i) = dbsc(i)%mdci
    iDmp(17,i) = dbsc(i)%iVal
    iDmp(18,i) = dbsc(i)%nVal
    iDmp(19,i) = dbsc(i)%iPrj
    iDmp(20,i) = dbsc(i)%nPrj
    iDmp(21,i) = dbsc(i)%iSRO
    iDmp(22,i) = dbsc(i)%nSRO
    iDmp(23,i) = dbsc(i)%iSOC
    iDmp(24,i) = dbsc(i)%nSOC
    iDmp(25,i) = dbsc(i)%iPP
    iDmp(26,i) = dbsc(i)%nPP
    iDmp(27,i) = dbsc(i)%nShells
    iDmp(28,i) = dbsc(i)%AtmNr
    iDmp(29,i) = 0
    if (dbsc(i)%NoPair) iDmp(29,i) = 1
    iDmp(30,i) = 0
    if (dbsc(i)%SODk) iDmp(30,i) = 1
    iDmp(31,i) = 0
    if (dbsc(i)%pChrg) iDmp(31,i) = 1
    iDmp(32,i) = 0
    if (dbsc(i)%Fixed) iDmp(32,i) = 1
    iDmp(33,i) = 0
    if (dbsc(i)%lPAM2) iDmp(33,i) = 1
    do j=0,iTabMx
      iDmp(34+j,i) = dbsc(i)%kDel(j)
    end do
    if ((.not. dbsc(i)%Aux) .or. (i == iCnttp_Dummy)) then
      nAtoms = nAtoms+dbsc(i)%nCntr
    end if
    nFragCoor = max(0,dbsc(i)%nFragCoor)  ! Fix the misuse in FragExpand
    nAux = nAux+2*dbsc(i)%nM1+2*dbsc(i)%nM2+nFrag_LineWords*dbsc(i)%nFragType+5*nFragCoor+dbsc(i)%nFragEner+ &
           dbsc(i)%nFragDens*dbsc(i)%nFragEner
#   ifdef _DEBUGPRINT_
    write(u6,'(A,8I4)') 'iCnttp=',i,nFrag_LineWords,dbsc(i)%nFragType,dbsc(i)%nFragCoor,dbsc(i)%nFragEner,dbsc(i)%nFragDens,nAux, &
                        dbsc(i)%iVal
#   endif

    if (dbsc(i)%nPAM2 /= -1) then
      write(u6,*) 'Not yet implemented for PAM2 integrals.'
      call Abend()
    end if
  end do
  iDmp(1,nCnttp+1) = nFrag_LineWords
  iDmp(2,nCnttp+1) = nCnttp
  iDmp(3,nCnttp+1) = iCnttp_Dummy
  iDmp(4,nCnttp+1) = Max_Shells
  iDmp(5,nCnttp+1) = Nuclear_Model
  iDmp(6:13,nCnttp+1) = nBas(0:7)
  iDmp(14:21,nCnttp+1) = nBas_Aux(0:7)
  iDmp(22:29,nCnttp+1) = nBas_Frag(0:7)
  iDmp(30,nCnttp+1) = MolWgh
  call Put_iArray('iDmp',iDmp,nFields*(nCnttp+1))
  call mma_deallocate(iDmp)

  ! Integer shells stuffs

  call mma_Allocate(iDmp,mFields,Max_Shells-1,Label='iDmp')
  nAux2 = 0
  do i=1,Max_Shells-1
    iDmp(1,i) = Shells(i)%nBK
    iDmp(2,i) = Shells(i)%nAkl
    iDmp(3,i) = Shells(i)%nFockOp
    iDmp(4,i) = Shells(i)%nExp
    iDmp(5,i) = Shells(i)%nBasis
    iDmp(6,i) = Shells(i)%nBasis_c
    iDmp(7,i) = 0
    if (Shells(i)%Transf) iDmp(7,i) = 1
    iDmp(8,i) = 0
    if (Shells(i)%Prjct) iDmp(8,i) = 1
    iDmp(9,i) = 0
    if (Shells(i)%Frag) iDmp(9,i) = 1
    iDmp(10,i) = 0
    if (Shells(i)%Aux) iDmp(10,i) = 1
    iDmp(11,i) = Shells(i)%kOffAO
    nAux2 = nAux2+2*Shells(i)%nBK+2*Shells(i)%nAkl**2+Shells(i)%nFockOp**2+Shells(i)%nExp+2*Shells(i)%nExp*Shells(i)%nBasis+ &
            2*Shells(i)%nExp**2
#   ifdef _DEBUGPRINT_
    write(u6,'(A,7I4)') 'iShll=',i,Shells(i)%nBK,Shells(i)%nAkl,Shells(i)%nFockOp,Shells(i)%nExp,Shells(i)%nBasis
#   endif

  end do
  call Put_iArray('iDmp:S',iDmp,mFields*(Max_Shells-1))
  call mma_deallocate(iDmp)

  ! Real Stuff

  call mma_allocate(rDmp,3,nAtoms+3*nCnttp,Label='rDmp')
  nAtoms = 0
  do i=1,nCnttp
    !call RecPrt('dbsc(i)%Coor',' ',dbsc(i)%Coor(1,1),3,dbsc(i)%nCntr)
    if ((.not. dbsc(i)%Aux) .or. (i == iCnttp_Dummy)) then
      do j=1,dbsc(i)%nCntr
        nAtoms = nAtoms+1
        rDmp(1:3,nAtoms) = dbsc(i)%Coor(1:3,j)
      end do
    end if
    nAtoms = nAtoms+1
    rDmp(1,nAtoms) = dbsc(i)%Charge
    rDmp(2,nAtoms) = dbsc(i)%CrRep
    rDmp(3,nAtoms) = dbsc(i)%FragCharge
    nAtoms = nAtoms+1
    rDmp(1,nAtoms) = dbsc(i)%aCD_Thr
    rDmp(2,nAtoms) = dbsc(i)%fMass
    rDmp(3,nAtoms) = dbsc(i)%CntMass
    nAtoms = nAtoms+1
    rDmp(1,nAtoms) = dbsc(i)%ExpNuc
    rDmp(2,nAtoms) = dbsc(i)%w_mGauss
    rDmp(3,nAtoms) = Zero
  end do
  call Put_dArray('rDmp',rDmp,3*nAtoms)
  call mma_deallocate(rDmp)

  if (nAux > 0) then
    !write(u6,*) 'nAux=',nAux
    call mma_allocate(rDmp,nAux,1,Label='rDmp')
    nAux = 0
    do i=1,nCnttp
      nM1 = dbsc(i)%nM1
      if (nM1 > 0) then
        !call RecPrt('M1xp',' ',dbsc(i)%M1xp,1,nM1)
        !call RecPrt('M1cf',' ',dbsc(i)%M1cf,1,nM1)
        rDmp(nAux+1:nAux+nM1,1) = dbsc(i)%M1xp(:)
        nAux = nAux+nM1
        rDmp(nAux+1:nAux+nM1,1) = dbsc(i)%M1cf(:)
        nAux = nAux+nM1
      end if
      nM2 = dbsc(i)%nM2
      if (nM2 > 0) then
        !call RecPrt('M2xp',' ',dbsc(i)%M2xp,1,nM2)
        !call RecPrt('M2cf',' ',dbsc(i)%M2cf,1,nM2)
        rDmp(nAux+1:nAux+nM2,1) = dbsc(i)%M2xp(:)
        nAux = nAux+nM2
        rDmp(nAux+1:nAux+nM2,1) = dbsc(i)%M2cf(:)
        nAux = nAux+nM2
      end if

      !write(u6,*) 'iAux=',nAux
      !write(u6,*) nFrag_LineWords, dbsc(i)%nFragType
      if (dbsc(i)%nFragType > 0) then
        qDmp(1:nFrag_LineWords,1:dbsc(i)%nFragType) => rDmp(nAux+1:nAux+nFrag_LineWords*dbsc(i)%nFragType,1)
        qDmp(:,:) = dbsc(i)%FragType(:,:)
        nAux = nAux+nFrag_LineWords*dbsc(i)%nFragType
        nullify(qDmp)
      end if
      !write(u6,*) dbsc(i)%nFragCoor
      nFragCoor = max(0,dbsc(i)%nFragCoor)
      if (nFragCoor > 0) then
        qDmp(1:5,1:nFragCoor) => rDmp(nAux+1:nAux+5*nFragCoor,1)
        qDmp(:,:) = dbsc(i)%FragCoor(:,:)
        nAux = nAux+5*nFragCoor
        nullify(qDmp)
      end if
      !write(u6,*) dbsc(i)%nFragEner
      if (dbsc(i)%nFragEner > 0) then
        rDmp(nAux+1:nAux+dbsc(i)%nFragEner,1) = dbsc(i)%FragEner(:)
        nAux = nAux+dbsc(i)%nFragEner
      end if
      !write(u6,*) dbsc(i)%nFragDens
      if (dbsc(i)%nFragDens*dbsc(i)%nFragEner > 0) then
        qDmp(1:dbsc(i)%nFragDens,1:dbsc(i)%nFragEner) => rDmp(nAux+1:nAux+dbsc(i)%nFragDens*dbsc(i)%nFragEner,1)
        qDmp(:,:) = dbsc(i)%FragCoef(:,:)
        nAux = nAux+dbsc(i)%nFragDens*dbsc(i)%nFragEner
        nullify(qDmp)
      end if
    end do
    !call RecPrt('rDmp:A',' ',rDmp,1,nAux)
    call Put_dArray('rDmp:A',rDmp,nAux)
    call mma_deallocate(rDmp)
  end if

  if (nAux2 > 0) then
    call mma_allocate(rDmp,nAux2,1,Label='rDmp')
    nAux2 = 0
    do i=1,Max_Shells-1
      nBk = Shells(i)%nBk
      if (nBk > 0) then
        rDmp(nAux2+1:nAux2+nBK,1) = Shells(i)%Bk(:)
        nAux2 = nAux2+nBK
        rDmp(nAux2+1:nAux2+nBK,1) = Shells(i)%Occ(:)
        nAux2 = nAux2+nBK
      end if
      nAkl = Shells(i)%nAkl
      if (nAkl > 0) then
        call DCopy_(2*nAkl**2,Shells(i)%Akl,1,rDmp(nAux2+1,1),1)
        nAux2 = nAux2+2*nAkl**2
      end if
      nFockOp = Shells(i)%nFockOp
      if (nFockOp > 0) then
        call DCopy_(nFockOp**2,Shells(i)%FockOp,1,rDmp(nAux2+1,1),1)
        nAux2 = nAux2+nFockOp**2
      end if
      nExp = Shells(i)%nExp
      if (nExp > 0) then
        call DCopy_(nExp,Shells(i)%Exp,1,rDmp(nAux2+1,1),1)
        nAux2 = nAux2+nExp
      end if
      nBasis = Shells(i)%nBasis
      ! Note: the contraction coefficients are not always there.
      if (nExp*nBasis > 0) then
        call DCopy_(2*nExp**2,Shells(i)%Cff_p,1,rDmp(nAux2+1,1),1)
        nAux2 = nAux2+2*nExp**2
      end if
      if (nExp*nBasis > 0) then
        call DCopy_(2*nExp*nBasis,Shells(i)%Cff_c,1,rDmp(nAux2+1,1),1)
        nAux2 = nAux2+2*nExp*nBasis
      end if
    end do
    call Put_dArray('rDmp:S',rDmp,nAux2)
    call mma_deallocate(rDmp)
  end if

  ! Character Stuff

  call mma_allocate(cDmp,nCnttp,Label='cDmp')
  do i=1,nCnttp
    cDmp(i)(1:80) = dbsc(i)%Bsl
    cDmp(i)(81:160) = dbsc(i)%Bsl_Old
  end do
  lcDmp = 2*80*nCnttp
  call Put_cArray('cDmp',cDmp(1),lcDmp)
  call mma_deallocate(cDmp)

# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'Exit Basis_Info_Dmp'
  write(u6,*)
# endif

  return

end subroutine Basis_Info_Dmp

!***********************************************************************
!***********************************************************************

subroutine Basis_Info_Get()

  use stdalloc, only: mma_allocate, mma_deallocate
  use Definitions, only: u6

  integer(kind=iwp) :: i, j, lcDmp, Len1, Len2, nAkl, nAtoms, nAux, nAux2, nBasis, nBK, nExp, nFockOp, nFragCoor, nFragDens, &
                       nFragEner, nFragType, nM1, nM2
  logical(kind=iwp) :: Found
  integer(kind=iwp), allocatable :: iDmp(:,:)
  real(kind=wp), allocatable, target :: rDmp(:,:)
  real(kind=wp), pointer :: pDmp(:), qDmp(:,:)
  character(len=2*80), allocatable :: cDmp(:)

# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'Enter Basis_Info_Get'
  write(u6,*)
# endif

  call qpg_iArray('iDmp',Found,Len1)
  Len2 = Len1/nFields
  call mma_Allocate(iDmp,nFields,Len2,Label='iDmp')
  if (Found) then
    call Get_iArray('iDmp',iDmp,Len1)
  else
    write(u6,*) 'Basis_Info_Get: iDmp not found!'
    call Abend()
  end if
  nFrag_LineWords = iDmp(1,Len2)
  nCnttp = iDmp(2,Len2)
  iCnttp_Dummy = iDmp(3,Len2)
  Max_Shells = iDmp(4,Len2)
  Nuclear_Model = iDmp(5,Len2)
  nBas(0:7) = iDmp(6:13,Len2)
  nBas_Aux(0:7) = iDmp(14:21,Len2)
  nBas_Frag(0:7) = iDmp(22:29,Len2)
  MolWgh = iDmp(30,Len2)
  nAux = 0

  ! Initiate the memory allocation of dsbc and Shells

  if (.not. Initiated) call Basis_Info_Init()

  do i=1,nCnttp
    dbsc(i)%nCntr = iDmp(1,i)
    dbsc(i)%nM1 = iDmp(2,i)
    dbsc(i)%nM2 = iDmp(3,i)
    dbsc(i)%nFragType = iDmp(4,i)
    dbsc(i)%nFragCoor = iDmp(5,i)
    dbsc(i)%nFragEner = iDmp(6,i)
    dbsc(i)%nFragDens = iDmp(7,i)
    dbsc(i)%ECP = iDmp(8,i) == 1
    dbsc(i)%Frag = iDmp(9,i) == 1
    dbsc(i)%Aux = iDmp(10,i) == 1
    dbsc(i)%FOp = iDmp(11,i) == 1
    dbsc(i)%IsMM = iDmp(12,i)
    dbsc(i)%Parent_iCnttp = iDmp(13,i)
    dbsc(i)%lOffAO = iDmp(14,i)
    dbsc(i)%nOpt = iDmp(15,i)
    dbsc(i)%mdci = iDmp(16,i)
    dbsc(i)%iVal = iDmp(17,i)
    dbsc(i)%nVal = iDmp(18,i)
    dbsc(i)%iPrj = iDmp(19,i)
    dbsc(i)%nPrj = iDmp(20,i)
    dbsc(i)%iSRO = iDmp(21,i)
    dbsc(i)%nSRO = iDmp(22,i)
    dbsc(i)%iSOC = iDmp(23,i)
    dbsc(i)%nSOC = iDmp(24,i)
    dbsc(i)%iPP = iDmp(25,i)
    dbsc(i)%nPP = iDmp(26,i)
    dbsc(i)%nShells = iDmp(27,i)
    dbsc(i)%AtmNr = iDmp(28,i)
    dbsc(i)%NoPair = iDmp(29,i) == 1
    dbsc(i)%SODK = iDmp(30,i) == 1
    dbsc(i)%pChrg = iDmp(31,i) == 1
    dbsc(i)%Fixed = iDmp(32,i) == 1
    dbsc(i)%lPAM2 = iDmp(33,i) == 1
    do j=0,iTabMx
      dbsc(i)%kDel(j) = iDmp(34+j,i)
    end do
    nFragCoor = max(0,dbsc(i)%nFragCoor)
    nAux = nAux+2*dbsc(i)%nM1+2*dbsc(i)%nM2+nFrag_LineWords*dbsc(i)%nFragType+5*nFragCoor+dbsc(i)%nFragEner+ &
           dbsc(i)%nFragDens*dbsc(i)%nFragEner
#   ifdef _DEBUGPRINT_
    write(u6,'(A,8I4)') 'iCnttp=',i,nFrag_LineWords,dbsc(i)%nFragType,dbsc(i)%nFragCoor,dbsc(i)%nFragEner,dbsc(i)%nFragDens,nAux, &
                        dbsc(i)%iVal
#   endif
  end do
  call mma_deallocate(iDmp)

  call qpg_iArray('iDmp:S',Found,Len1)
  call mma_Allocate(iDmp,mFields,Len1/mFields,Label='iDmp')
  if (Found) then
    call get_iArray('iDmp:S',iDmp,Len1)
  else
    write(u6,*) 'Basis_Info_Get: iDmp:S not found!'
    call Abend()
  end if
  nAux2 = 0
  do i=1,Max_Shells-1
    Shells(i)%nBK = iDmp(1,i)
    Shells(i)%nAkl = iDmp(2,i)
    Shells(i)%nFockOp = iDmp(3,i)
    Shells(i)%nExp = iDmp(4,i)
    Shells(i)%nBasis = iDmp(5,i)
    Shells(i)%nBasis_c = iDmp(6,i)
    Shells(i)%Transf = iDmp(7,i) == 1
    Shells(i)%Prjct = iDmp(8,i) == 1
    Shells(i)%Frag = iDmp(9,i) == 1
    Shells(i)%Aux = iDmp(10,i) == 1
    Shells(i)%kOffAO = iDmp(11,i)
    nAux2 = nAux2+2*Shells(i)%nBK+2*Shells(i)%nAkl**2+Shells(i)%nFockOp**2+Shells(i)%nExp
    ! Coefficients only there is nBasis =/=0
    if (Shells(i)%nBasis > 0) then
      nAux2 = nAux2+2*Shells(i)%nExp*Shells(i)%nBasis+2*Shells(i)%nExp**2
    end if
  end do
  call mma_deallocate(iDmp)

  call qpg_dArray('rDmp',Found,Len1)
  if (.not. Found) then
    write(u6,*) 'rDMP not found on the run file.'
    call Abend()
  end if
  nAtoms = Len1/3
  call mma_allocate(rDmp,3,nAtoms,Label='rDmp')
  call Get_dArray('rDmp',rDmp,3*nAtoms)
  nAtoms = 0
  do i=1,nCnttp
    if ((.not. dbsc(i)%Aux) .or. (i == iCnttp_Dummy)) then
      if (dbsc(i)%nCntr > 0) then
        call mma_Allocate(dbsc(i)%Coor_Hidden,3,dbsc(i)%nCntr,Label='dbsc:C')
        dbsc(i)%Coor => dbsc(i)%Coor_Hidden(:,:)
      end if
      do j=1,dbsc(i)%nCntr
        nAtoms = nAtoms+1
        dbsc(i)%Coor(1:3,j) = rDmp(1:3,nAtoms)
      end do
    else
      j = dbsc(i)%Parent_iCnttp
      dbsc(i)%Coor => dbsc(j)%Coor(1:3,1:dbsc(j)%nCntr)
    end if
    nAtoms = nAtoms+1
    dbsc(i)%Charge = rDmp(1,nAtoms)
    dbsc(i)%CrRep = rDmp(2,nAtoms)
    dbsc(i)%FragCharge = rDmp(3,nAtoms)
    nAtoms = nAtoms+1
    dbsc(i)%aCD_Thr = rDmp(1,nAtoms)
    dbsc(i)%fMass = rDmp(2,nAtoms)
    dbsc(i)%CntMass = rDmp(3,nAtoms)
    nAtoms = nAtoms+1
    dbsc(i)%ExpNuc = rDmp(1,nAtoms)
    dbsc(i)%w_mGauss = rDmp(2,nAtoms)
  end do
  call mma_deallocate(rDmp)

  if (nAux > 0) then
    call qpg_dArray('rDmp:A',Found,Len1)
    !write(u6,*) 'nAux=',nAux
    call mma_allocate(rDmp,nAux,1,Label='rDmp')
    call Get_dArray('rDmp:A',rDmp,Len1)
    nAux = 0
    do i=1,nCnttp

      ! ECP stuff

      nM1 = dbsc(i)%nM1
      if (nM1 > 0) then
        if (.not. allocated(dbsc(i)%M1xp)) call mma_allocate(dbsc(i)%M1xp,nM1,Label='dbsc:M1xp')
        dbsc(i)%M1xp(:) = rDmp(nAux+1:nAux+nM1,1)
        nAux = nAux+nM1
        if (.not. allocated(dbsc(i)%M1cf)) call mma_allocate(dbsc(i)%M1cf,nM1,Label='dbsc:M1cf')
        dbsc(i)%M1cf(:) = rDmp(nAux+1:nAux+nM1,1)
        nAux = nAux+nM1
        !call RecPrt('M1xp',' ',dbsc(i)%M1xp,1,nM1)
        !call RecPrt('M1cf',' ',dbsc(i)%M1cf,1,nM1)
      end if
      nM2 = dbsc(i)%nM2
      if (nM2 > 0) then
        if (.not. allocated(dbsc(i)%M2xp)) call mma_allocate(dbsc(i)%M2xp,nM2,Label='dbsc:M2xp')
        dbsc(i)%M2xp(:) = rDmp(nAux+1:nAux+nM2,1)
        nAux = nAux+nM2
        if (.not. allocated(dbsc(i)%M2cf)) call mma_allocate(dbsc(i)%M2cf,nM2,Label='dbsc:M2cf')
        dbsc(i)%M2cf(:) = rDmp(nAux+1:nAux+nM2,1)
        nAux = nAux+nM2
        !call RecPrt('M2xp',' ',dbsc(i)%M2xp,1,nM2)
        !call RecPrt('M2cf',' ',dbsc(i)%M2cf,1,nM2)
      end if

      ! Fragment stuff

      nFragType = dbsc(i)%nFragType
      if (nFragType > 0) then
        if (.not. allocated(dbsc(i)%FragType)) call mma_allocate(dbsc(i)%FragType,nFrag_LineWords,nFragType,Label='FragType')
        qDmp(1:nFrag_LineWords,1:nFragType) => rDmp(nAux+1:nAux+nFrag_LineWords*nFragType,1)
        dbsc(i)%FragType(:,:) = qDmp(:,:)
        nAux = nAux+nFrag_LineWords*nFragType
        nullify(qDmp)
      end if
      nFragCoor = max(0,dbsc(i)%nFragCoor)
      if (nFragCoor > 0) then
        if (.not. allocated(dbsc(i)%FragCoor)) call mma_allocate(dbsc(i)%FragCoor,5,nFragCoor,Label='FragCoor')
        qDmp(1:5,1:nFragCoor) => rDmp(nAux+1:nAux+5*nFragCoor,1)
        dbsc(i)%FragCoor(:,:) = qDmp(:,:)
        nAux = nAux+5*nFragCoor
        nullify(qDmp)
      end if
      nFragEner = dbsc(i)%nFragEner
      if (nFragEner > 0) then
        if (.not. allocated(dbsc(i)%FragEner)) call mma_allocate(dbsc(i)%FragEner,nFragEner,Label='FragEner')
        pDmp(1:nFragEner) => rDmp(nAux+1:nAux+nFragEner,1)
        dbsc(i)%FragEner(:) = pDmp(:)
        nAux = nAux+nFragEner
        nullify(pDmp)
      end if
      nFragDens = dbsc(i)%nFragDens
      if (nFragDens*nFragEner > 0) then
        if (.not. allocated(dbsc(i)%FragCoef)) call mma_allocate(dbsc(i)%FragCoef,nFragDens,nFragEner,Label='FragCoef')
        qDmp(1:nFragDens,1:nFragEner) => rDmp(nAux+1:nAux+nFragDens*nFragEner,1)
        dbsc(i)%FragCoef(:,:) = qDmp(:,:)
        nAux = nAux+nFragDens*nFragEner
        nullify(qDmp)
      end if
    end do
    call mma_deallocate(rDmp)
  end if

  if (nAux2 > 0) then
    call qpg_dArray('rDmp:S',Found,Len1)
    call mma_allocate(rDmp,nAux2,1,Label='rDmp')
    call Get_dArray('rDmp:S',rDmp,Len1)
    nAux2 = 0
    do i=1,Max_Shells-1

      nBk = Shells(i)%nBK
      if (nBk > 0) then
        if (.not. allocated(Shells(i)%Bk)) call mma_allocate(Shells(i)%Bk,nBk,Label='Bk')
        Shells(i)%Bk(:) = rDmp(nAux2+1:nAux2+nBk,1)
        nAux2 = nAux2+nBk
        if (.not. allocated(Shells(i)%Occ)) call mma_allocate(Shells(i)%Occ,nBk,Label='Occ')
        Shells(i)%Occ(:) = rDmp(nAux2+1:nAux2+nBk,1)
        nAux2 = nAux2+nBk
      end if

      nAkl = Shells(i)%nAkl
      if (nAkl > 0) then
        if (.not. allocated(Shells(i)%Akl)) call mma_allocate(Shells(i)%Akl,nAkl,nAkl,2,Label='Akl')
        call DCopy_(2*nAkl**2,rDmp(nAux2+1,1),1,Shells(i)%Akl,1)
        nAux2 = nAux2+2*nAkl**2
      end if

      nFockOp = Shells(i)%nFockOp
      if (nFockOp > 0) then
        if (.not. allocated(Shells(i)%FockOp)) call mma_allocate(Shells(i)%FockOp,nFockOp,nFockOp,Label='FockOp')
        call DCopy_(nFockOp**2,rDmp(nAux2+1,1),1,Shells(i)%FockOp,1)
        nAux2 = nAux2+nFockOp**2
      end if

      nExp = Shells(i)%nExp
      if (nExp > 0) then
        if (.not. allocated(Shells(i)%Exp)) call mma_allocate(Shells(i)%Exp,nExp,Label='Exp')
        call DCopy_(nExp,rDmp(nAux2+1,1),1,Shells(i)%Exp,1)
        nAux2 = nAux2+nExp
      end if

      nBasis = Shells(i)%nBasis
      if (nExp*nBasis > 0) then
        if (.not. allocated(Shells(i)%Cff_p)) call mma_allocate(Shells(i)%Cff_p,nExp,nExp,2,Label='Cff_p')
        call DCopy_(2*nExp**2,rDmp(nAux2+1,1),1,Shells(i)%Cff_p,1)
        nAux2 = nAux2+2*nExp**2
      end if

      if (nExp*nBasis > 0) then
        if (.not. allocated(Shells(i)%Cff_c)) call mma_allocate(Shells(i)%Cff_c,nExp,nBasis,2,Label='Cff_c')
        call DCopy_(2*nExp*nBasis,rDmp(nAux2+1,1),1,Shells(i)%Cff_c,1)
        nAux2 = nAux2+2*nExp*nBasis

        if (.not. allocated(Shells(i)%pCff)) call mma_allocate(Shells(i)%pCff,nExp,nBasis,Label='Cff')
        Shells(i)%pCff(:,:) = Shells(i)%Cff_c(:,:,1)
      end if
    end do
    call mma_deallocate(rDmp)
  end if

  ! Character stuff

  call mma_allocate(cDmp,nCnttp,Label='cDmp')
  lcDmp = 2*80*nCnttp
  call Get_cArray('cDmp',cDmp,lcDmp)
  do i=1,nCnttp
    dbsc(i)%Bsl = cDmp(i)(1:80)
    dbsc(i)%Bsl_Old = cDmp(i)(81:160)
  end do
  call mma_deallocate(cDmp)
# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'Coordinates:'
  do i=1,nCnttp
    do j=1,dbsc(i)%nCntr
      write(u6,*) dbsc(i)%Coor(:,j)
    end do
  end do
  write(u6,*)
  write(u6,*) 'Exit Basis_Info_Get'
  write(u6,*)
# endif

  return

end subroutine Basis_Info_Get

!***********************************************************************
!***********************************************************************

subroutine Basis_Info_Free()

  use stdalloc, only: mma_deallocate
# ifdef _DEBUGPRINT_
  use Definitions, only: u6
# endif

  integer(kind=iwp) :: i

# ifdef _DEBUGPRINT_
  write(u6,*) 'Basis_Info_Free()'
# endif

  ! Deallocate all allocatable parts of dbsc.

  do i=1,nCnttp

    ! Molecular Coordinates

    if (dbsc(i)%nCntr > 0) then
      if ((.not. dbsc(i)%Aux) .or. (i == iCnttp_Dummy)) then
        call mma_deallocate(dbsc(i)%Coor_Hidden)
      end if
      nullify(dbsc(i)%Coor)
      dbsc(i)%nCntr = 0
    end if

    ! ECP stuff

    if (allocated(dbsc(i)%M1xp)) call mma_deallocate(dbsc(i)%M1xp)
    if (allocated(dbsc(i)%M1cf)) call mma_deallocate(dbsc(i)%M1cf)
    dbsc(i)%nM1 = 0
    if (allocated(dbsc(i)%M2xp)) call mma_deallocate(dbsc(i)%M2xp)
    if (allocated(dbsc(i)%M2cf)) call mma_deallocate(dbsc(i)%M2cf)
    dbsc(i)%nM2 = 0

    ! Fragment stuff

    if (allocated(dbsc(i)%FragType)) call mma_deallocate(dbsc(i)%FragType)
    dbsc(i)%nFragType = 0
    if (allocated(dbsc(i)%FragCoor)) call mma_deallocate(dbsc(i)%FragCoor)
    dbsc(i)%nFragCoor = 0
    if (allocated(dbsc(i)%FragEner)) call mma_deallocate(dbsc(i)%FragEner)
    dbsc(i)%nFragEner = 0
    if (allocated(dbsc(i)%FragCoef)) call mma_deallocate(dbsc(i)%FragCoef)
    dbsc(i)%nFragDens = 0

    ! PAM2 stuff

    if (allocated(dbsc(i)%PAM2)) call mma_deallocate(dbsc(i)%PAM2)
    dbsc(i)%nPAM2 = -1
  end do
  nCnttp = 0
  iCnttp_Dummy = 0

  ! Stuff on unique basis set shells

  do i=1,Max_Shells-1
    if (allocated(Shells(i)%Bk)) call mma_deallocate(Shells(i)%Bk)
    if (allocated(Shells(i)%Occ)) call mma_deallocate(Shells(i)%Occ)
    Shells(i)%nBk = 0
    if (allocated(Shells(i)%Akl)) call mma_deallocate(Shells(i)%Akl)
    Shells(i)%nAkl = 0
    if (allocated(Shells(i)%FockOp)) call mma_deallocate(Shells(i)%FockOp)
    Shells(i)%nFockOp = 0
    if (allocated(Shells(i)%Exp)) call mma_deallocate(Shells(i)%Exp)
    Shells(i)%nExp = 0
    if (allocated(Shells(i)%pCff)) call mma_deallocate(Shells(i)%pCff)
    if (allocated(Shells(i)%Cff_c)) call mma_deallocate(Shells(i)%Cff_c)
    if (allocated(Shells(i)%Cff_p)) call mma_deallocate(Shells(i)%Cff_p)
    Shells(i)%nBasis = 0
    Shells(i)%Transf = .true.
  end do
  Max_Shells = 0

  if (allocated(dbsc)) call mma_deallocate(dbsc)
  if (allocated(Shells)) call mma_deallocate(Shells)
  Initiated = .false.

  return

end subroutine Basis_Info_Free

!***********************************************************************
!***********************************************************************

! Private extensions to mma_interfaces, using preprocessor templates
! (see src/mma_util/stdalloc.f)

! Define dbsc_cptr2loff, dbsc_mma_allo_1D, dbsc_mma_allo_1D_lim, dbsc_mma_free_1D
! (using _NO_GARBLE_ because all members are initialized)
#define _TYPE_ type(Distinct_Basis_set_centers)
#  define _NO_GARBLE_
#  define _FUNC_NAME_ dbsc_cptr2loff
#  include "cptr2loff_template.fh"
#  undef _FUNC_NAME_
#  define _SUBR_NAME_ dbsc_mma
#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'dbsc_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#  undef _NO_GARBLE_
#undef _TYPE_

! Define shell_cptr2loff, shell_mma_allo_1D, shell_mma_allo_1D_lim, shell_mma_free_1D
! (using _NO_GARBLE_ because all members are initialized)
#define _TYPE_ type(Shell_Info)
#  define _NO_GARBLE_
#  define _FUNC_NAME_ shell_cptr2loff
#  include "cptr2loff_template.fh"
#  undef _FUNC_NAME_
#  define _SUBR_NAME_ shell_mma
#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'shell_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#  undef _NO_GARBLE_
#undef _TYPE_

end module Basis_Info
