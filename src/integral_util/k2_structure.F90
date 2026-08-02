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
!      abMax  : largest abCon                                          *
!      ZetaM  : largest Zeta value                                     *
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
!      Converted to a user defined type by R.L. October 2023           *
!***********************************************************************

module k2_structure

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

type k2_type
  integer(kind=iwp) :: nZeta = 0, ijCmp = 0, nHm = 0
  integer(kind=iwp), pointer :: IndZ(:) => null()
  real(kind=wp) :: EstI = Zero, ZetaM = Zero, abMax = Zero, abMaxD = Zero, abConMax = Zero
  real(kind=wp), pointer :: Zeta(:) => null(), Kappa(:) => null(), PCoor(:,:) => null(), ZInv(:) => null(), ab(:) => null(), &
                            abG(:,:) => null(), abCon(:) => null(), Alpha(:) => null(), Beta(:) => null(), HrrMtrx(:,:) => null()
end type k2_type

integer(kind=iwp) :: iZZZ_i = 0, iZZZ_r = 0
logical(kind=iwp) :: k2_processed = .false.
integer(kind=iwp), allocatable :: IndK2(:,:)
integer(kind=iwp), allocatable, target :: ZZZ_i(:)
real(kind=wp), allocatable, target :: ZZZ_r(:)
type(k2_type), allocatable, target :: k2data(:,:)

public :: Allocate_k2data, Allocate_k2data_in, Free_k2data, IndK2, k2_processed, k2_type, k2data, ZZZ_i, ZZZ_r

! Private extensions to mma interfaces

interface mma_allocate
  module procedure :: k2d_mma_allo_2D, k2d_mma_allo_2D_lim
end interface
interface mma_deallocate
  module procedure :: k2d_mma_free_2D
end interface

contains

subroutine Allocate_k2data(nIrrep,nk2)

  integer(kind=iwp), intent(in) :: nIrrep, nk2

  call mma_allocate(k2Data,nIrrep,nk2,label='k2Data')

# include "macros.fh"
  unused_proc(mma_allocate(k2Data,[0,0],[0,0]))

end subroutine Allocate_k2data

subroutine Allocate_k2data_in(k2Data_1D,nZeta,ijCmp,nHm)

  use Symmetry_Info, only: nIrrep
  use Definitions, only: u6

  type(k2_type), target, intent(inout) :: k2Data_1D
  integer(kind=iwp), intent(in) :: nZeta, ijCmp, nHm
  integer(kind=iwp) :: iE, iS

  k2Data_1D%nZeta = nZeta
  k2Data_1D%nHm = nHm
  k2Data_1D%ijCmp = ijCmp

  iE = iZZZ_r

  iS = iE+1
  iE = iE+nZeta
  k2Data_1D%Zeta(1:nZeta) => ZZZ_r(iS:iE)
  iS = 1+iE
  iE = iE+nZeta
  k2Data_1D%Kappa(1:nZeta) => ZZZ_r(iS:iE)
  iS = 1+iE
  iE = iE+nZeta*3
  k2Data_1D%PCoor(1:nZeta,1:3) => ZZZ_r(iS:iE)
  iS = 1+iE
  iE = iE+nZeta
  k2Data_1D%ZInv(1:nZeta) => ZZZ_r(iS:iE)
  iS = 1+iE
  iE = iE+nZeta
  k2Data_1D%ab(1:nZeta) => ZZZ_r(iS:iE)
  iS = 1+iE
  iE = iE+nZeta
  k2Data_1D%abCon(1:nZeta) => ZZZ_r(iS:iE)
  iS = 1+iE
  iE = iE+nZeta
  k2Data_1D%Alpha(1:nZeta) => ZZZ_r(iS:iE)
  iS = 1+iE
  iE = iE+nZeta
  k2Data_1D%Beta(1:nZeta) => ZZZ_r(iS:iE)
  if (nHm /= 0) then
    iS = 1+iE
    iE = iE+nHm*nIrrep
    k2Data_1D%HrrMtrx(1:nHm,1:nIrrep) => ZZZ_r(iS:iE)
  end if
  if (ijCmp /= 0) then
    iS = 1+iE
    iE = iE+nZeta*ijCmp*2
    k2Data_1D%abG(1:nZeta*ijCmp,1:2) => ZZZ_r(iS:iE)
  end if
  iZZZ_r = iE
  if (iZZZ_r > size(ZZZ_r)) then
    write(u6,*) 'iZZZ_r out for range'
    call Abend()
  end if

  iE = iZZZ_i

  iS = iE+1
  iE = iE+nZeta+1
  k2Data_1D%IndZ(1:nZeta+1) => ZZZ_i(iS:iE)
  iZZZ_i = iE
  if (iZZZ_i > size(ZZZ_i)) then
    write(u6,*) 'iZZZ_i out for range'
    call Abend()
  end if

end subroutine Allocate_k2data_in

subroutine Free_k2data()

  integer(kind=iwp) :: i, iIrrep

  do i=1,size(k2data,2)
    do iIrrep=1,size(k2data,1)
      call Free_k2data_Internal(k2data(iIrrep,i))
    end do
  end do

  call mma_deallocate(ZZZ_r)
  iZZZ_r = 0
  call mma_deallocate(ZZZ_i)
  iZZZ_i = 0

  call mma_deallocate(k2data)

end subroutine Free_k2data

subroutine Free_k2data_Internal(k2Data_1D)

  type(k2_type), intent(inout) :: k2Data_1D

  nullify(k2Data_1D%Zeta,k2Data_1D%Kappa,k2Data_1D%Pcoor,k2Data_1D%ZInv,k2Data_1D%ab,k2Data_1D%abG,k2Data_1D%abCon, &
          k2Data_1D%Alpha,k2Data_1D%Beta,k2Data_1D%HRRMtrx,k2Data_1D%IndZ)

end subroutine Free_k2data_Internal

! Private extensions to mma_interfaces, using preprocessor templates
! (see src/mma_util/stdalloc.f)

! Define k2d_mma_allo_2D, k2d_mma_allo_2D_lim, k2d_mma_free_2D
! (using _NO_GARBLE_ because all members are initialized)
#define _TYPE_ type(k2_type)
#  define _NO_GARBLE_
#  define _SUBR_NAME_ k2d_mma
#  define _DIMENSIONS_ 2
#  define _DEF_LABEL_ 'k2_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

end module k2_structure
