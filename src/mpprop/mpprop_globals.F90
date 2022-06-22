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

module MPProp_globals

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
private

#include "LenIn.fh"

! A type for multipole array data, where each element
! is an array of different size (number of components)
! (Different from blockdiagonal_matrices, because here each "block" is not square)
type MltPlArr
  real(kind=wp), allocatable :: M(:,:)
end type MltPlArr

real(kind=wp) :: EneV
character(len=180) :: Title
character(len=8) :: Method
integer(kind=iwp), allocatable :: iAtomType(:), iAtomPar(:), nAtomPBas(:), iAtPrTab(:,:)
real(kind=wp), allocatable :: Cor(:,:,:), CordMltPl(:,:), Frac(:,:), AtPol(:,:), AtBoPol(:,:), Qnuc(:)
character(len=LenIn), allocatable :: Labe(:)
character(len=LenIn*2+2), allocatable :: Cen_Lab(:)
logical(kind=iwp), allocatable :: BondMat(:,:)
type(MltPlArr), allocatable :: AtBoMltPl(:), AtBoMltPlCopy(:), AtBoMltPlTot(:), AtMltPl(:), AtMltPlTot(:), MltPl(:)

public :: AtBoMltPl, AtBoMltPlCopy, AtBoMltPlTot, AtBoPol, AtMltPl, AtMltPlTot, AtPol, BondMat, Cen_Lab, Cor, CordMltPl, EneV, &
          Frac, iAtomPar, iAtomType, iAtPrTab, Labe, Method, MltPl, nAtomPBas, Qnuc, Title

! Private extensions to mma interfaces
interface cptr2loff
  module procedure :: mltpl_cptr2loff
end interface
interface mma_allocate
  module procedure :: mltpl_mma_allo_1D, mltpl_mma_allo_1D_lim
end interface
interface mma_deallocate
  module procedure :: mltpl_mma_free_1D
end interface

public :: Alloc_MltPlArr, Free_MltPlArr

contains

! Private extensions to mma_interfaces, using preprocessor templates
! (see src/mma_util/stdalloc.f)

! Define mltpl_cptr2loff, mltpl_mma_allo_1D, mltpl_mma_allo_1D_lim, mltpl_mma_free_1D
#define _TYPE_ type(MltPlArr)
#  define _FUNC_NAME_ mltpl_cptr2loff
#  include "cptr2loff_template.fh"
#  undef _FUNC_NAME_
#  define _SUBR_NAME_ mltpl_mma
#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'mltpl_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

subroutine Alloc_MltPlArr(Array,N,Label)
  type(MltPlArr), allocatable, intent(inout) :: Array(:)
  integer(kind=iwp), intent(in) :: N(2)
  character(len=*), intent(in) :: Label
# ifdef _GARBLE_
  interface
    subroutine c_null_alloc(A)
      import wp
      real(kind=wp), allocatable :: A(:,:)
    end subroutine c_null_alloc
  end interface
  integer(kind=iwp) :: i
# endif
  call mma_allocate(Array,N,label=Label)
# ifdef _GARBLE_
  ! Garbling corrupts the allocation status of allocatable components, use a hack to reset it
  do i=N(1),N(2)
    call c_null_alloc(Array(i)%M)
  end do
# endif
end subroutine

subroutine Free_MltPlArr(Array)
  type(MltPlArr), allocatable, intent(inout) :: Array(:)
  integer(kind=iwp) :: i
  do i=lbound(Array,1),ubound(Array,1)
    if (allocated(Array(i)%M)) call mma_deallocate(Array(i)%M)
  end do
  call mma_deallocate(Array)
# ifdef _WARNING_WORKAROUND_
  ! unused subroutine
  if (.false.) then
    call mma_allocate(Array,0)
  end if
# endif
end subroutine Free_MltPlArr

end module MPProp_globals
