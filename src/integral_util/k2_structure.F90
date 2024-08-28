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
! Copyright (C) 1999,2023, Roland Lindh                                *
!               2012, Andy May                                         *
!***********************************************************************
!***********************************************************************
!                                                                      *
!      These functions should eventually be used globally for any      *
!      seward utility, be it integrals, gradients, direct Fock,        *
!      direct integrals, and second order derivatives.                 *
!                                                                      *
!                          W A R N I N G !                             *
!                                                                      *
!      Observe that the index array (IndZ) should always be placed     *
!      after all arrays with real!                                     *
!                                                                      *
!      Real Arrays:                                                    *
!      alpha  : alpha Gaussian exponent                                *
!      beta   : beta  Gaussian exponent                                *
!      Z      : alpha+beta Gaussian exponent                           *
!      Kappa  : Common factor from Gaussian product theorem            *
!      Pcoor  : Coordinates of Gaussian product                        *
!      ab     : Max(|(ab|ab)|^{1/2})    (primitive basis)              *
!               over all angular components                            *
!      abCon  : Max(|(ab|ab)|^{1/2})* C ("contracted basis")           *
!               over all angular components                            *
!                                                                      *
!      Integer Arrays:                                                 *
!      IndZ   : Pair index (rectangular indexation), the last          *
!               entry in the array contains the effective number       *
!               of elements                                            *
!                                                                      *
!      Scalars:                                                        *
!      EstI   : Estimated largest contracted integral |(ab|ab)|^{1/2}  *
!      ZtMax  : Z of the largest abCon                                 *
!      abMax  : largest abCon                                          *
!      ZetaM  : largest Zeta value                                     *
!      ZtMaxD : Z of the largest ab * D                                *
!      abMaxD : largest ab * D                                         *
!                                                                      *
!      Auxiliary arrays:                                               *
!      HrrMtrx: matrices to use for transformation from intermediate   *
!               integrals to real spherical harmonics                  *
!                                                                      *
!      Author: R. Lindh                                                *
!              Dept. of Chemical Physics                               *
!              University of Lund, Sweden                              *
!              April 1999                                              *
!                                                                      *
!      Converted from statement functions by A. May June 2012          *
!      Converted to a user defined type bt R.L. October 2023           *
!***********************************************************************

module k2_structure

use Constants, only: Zero

implicit none
private

type k2_type
  integer :: nZeta = 0
  integer :: ijCmp = 0
  integer :: nHm = 0
  real*8, pointer :: Zeta(:)
  real*8, pointer :: Kappa(:)
  real*8, pointer :: PCoor(:,:)
  real*8, pointer :: ZInv(:)
  real*8, pointer :: ab(:)
  real*8, pointer :: abG(:,:)
  real*8, pointer :: abCon(:)
  real*8, pointer :: Alpha(:)
  real*8, pointer :: Beta(:)
  real*8, pointer :: HrrMtrx(:,:)
  integer, pointer :: IndZ(:)
  real*8 :: EstI = Zero
  real*8 :: ZtMax = Zero
  real*8 :: ZtMaxD = Zero
  real*8 :: ZetaM = Zero
  real*8 :: abMax = Zero
  real*8 :: abMaxD = Zero
end type k2_type

type(k2_type), allocatable, target :: k2data(:,:)

integer, parameter :: nDArray = 11, nDScalar = 9

real*8, allocatable, target :: ZZZ_r(:)
integer, allocatable, target :: ZZZ_i(:)
integer :: iZZZ_r = 0, iZZZ_i = 0

logical :: k2_processed = .false.
integer, allocatable :: IndK2(:,:)
integer nIndk2

public :: k2_type, k2data, Allocate_k2data, Free_k2data
public :: nDArray, nDScalar, k2_processed, IndK2, nIndk2
public :: ZZZ_i, ZZZ_r

contains

subroutine Allocate_k2data(k2data,nZeta,ijCmp,nHm)

  use Symmetry_Info, only: nIrrep
  use Definitions, only: u6

  integer nZeta, ijCmp, nHm
  type(k2_type), target :: k2data
  integer iS, iE

  k2Data%nZeta = nZeta
  k2Data%nHm = nHm
  k2Data%ijCmp = ijCmp

  iE = iZZZ_r

  iS = iE+1
  iE = iE+nZeta
  k2Data%Zeta(1:nZeta) => ZZZ_r(iS:iE)
  iS = 1+iE
  iE = iE+nZeta
  k2Data%Kappa(1:nZeta) => ZZZ_r(iS:iE)
  iS = 1+iE
  iE = iE+nZeta*3
  k2Data%PCoor(1:nZeta,1:3) => ZZZ_r(iS:iE)
  iS = 1+iE
  iE = iE+nZeta
  k2Data%ZInv(1:nZeta) => ZZZ_r(iS:iE)
  iS = 1+iE
  iE = iE+nZeta
  k2Data%ab(1:nZeta) => ZZZ_r(iS:iE)
  iS = 1+iE
  iE = iE+nZeta
  k2Data%abCon(1:nZeta) => ZZZ_r(iS:iE)
  iS = 1+iE
  iE = iE+nZeta
  k2Data%Alpha(1:nZeta) => ZZZ_r(iS:iE)
  iS = 1+iE
  iE = iE+nZeta
  k2Data%Beta(1:nZeta) => ZZZ_r(iS:iE)
  if (nHm /= 0) then
    iS = 1+iE
    iE = iE+nHm*nIrrep
    k2Data%HrrMtrx(1:nHm,1:nIrrep) => ZZZ_r(iS:iE)
  end if
  if (ijCmp /= 0) then
    iS = 1+iE
    iE = iE+nZeta*ijCmp*2
    k2Data%abG(1:nZeta*ijCmp,1:2) => ZZZ_r(iS:iE)
  end if
  iZZZ_r = iE
  if (iZZZ_r > size(ZZZ_r)) then
    write(u6,*) 'iZZZ_r out for range'
    call Abend()
  end if

  iE = iZZZ_i

  iS = iE+1
  iE = iE+nZeta+1
  k2Data%IndZ(1:nZeta+1) => ZZZ_i(iS:iE)
  iZZZ_i = iE
  if (iZZZ_i > size(ZZZ_i)) then
    write(u6,*) 'iZZZ_i out for range'
    call Abend()
  end if

end subroutine Allocate_k2data

subroutine Free_k2data()

  use stdalloc, only: mma_deallocate

  integer nIrrep, ik2, iIrrep, i

  nIrrep = size(k2data,1)
  ik2 = size(k2data,2)

  do i=1,ik2
    do iIrrep=1,nIrrep
      call Free_k2data_Internal(k2data(iIrrep,i))
    end do
  end do

  call mma_deallocate(ZZZ_r)
  iZZZ_r = 0
  call mma_deallocate(ZZZ_i)
  iZZZ_i = 0

  deallocate(k2data)

end subroutine Free_k2data

subroutine Free_k2data_Internal(k2data_1D)

  type(k2_type) :: k2data_1D

  k2Data%nZeta = 0
  k2Data%nHm = 0
  k2Data%ijCmp = 0

  nullify(k2Data_1D%Zeta,k2Data_1D%Kappa,k2Data_1D%Pcoor,k2Data_1D%ZInv,k2Data_1D%ab,k2Data_1D%abCon,k2Data_1D%Alpha, &
          k2Data_1D%Beta,k2Data_1D%HRRMtrx,k2Data_1D%abG,k2Data_1D%IndZ)

end subroutine Free_k2data_Internal

end module k2_structure
