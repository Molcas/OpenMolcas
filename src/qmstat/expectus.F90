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

subroutine Expectus(QMMethod,HmatOld,Vmat,VpolMat,Smat,VEC,nDim,lEig,iEig,ExpVal)

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
character(len=5), intent(in) :: QMMethod
integer(kind=iwp), intent(in) :: nDIM, iEig
real(kind=wp), intent(in) :: HmatOld(nTri_Elem(nDim)), Vmat(nTri_Elem(nDim)), VpolMat(nTri_Elem(nDim)), Smat(nTri_Elem(nDim)), &
                             VEC(nDim,nDim)
logical(kind=iwp), intent(in) :: lEig
real(kind=wp), intent(_OUT_) :: ExpVal(4,*)
integer(kind=iwp) :: iRoot, nDTri, nRoots
real(kind=wp), allocatable :: DTmp(:)
real(kind=wp), external :: Ddot_
#include "warnings.h"

! Take different path for different QM-method.

if (QMMethod(1:5) == 'RASSI') then

  ! For how many roots are the eigenvalues to be computed.

  if (lEig) then
    nRoots = iEig
  else
    nRoots = nDim
  end if

  ! Loop over roots and compute expectation values according to
  ! well-known formulas.

  nDTri = nTri_Elem(nDim)
  call mma_allocate(DTmp,nDTri,label='DenTemp')
  do iRoot=1,nRoots

    ! Generate density matrix for relevant root.

    call DensiSt(DTmp,VEC,iRoot,nDim,nDim)

    ! Expectation values.

    ExpVal(1,iRoot) = Ddot_(nDTri,DTmp,1,HmatOld,1)
    ExpVal(2,iRoot) = Ddot_(nDTri,DTmp,1,Vmat,1)
    ExpVal(3,iRoot) = Ddot_(nDTri,DTmp,1,Vpolmat,1)
    ExpVal(4,iRoot) = Ddot_(nDTri,DTmp,1,Smat,1)
  end do
  call mma_deallocate(DTmp)

else if (QMMethod(1:5) == 'SCF  ') then
  ! If it's SCF we are running.

  nDTri = nTri_Elem(nDim)
  call mma_allocate(DTmp,nDTri,label='DenTemp')
  call Densi_MO(DTmp,VEC,1,iEig,nDim,nDim)

  ! Expectation values.

  ExpVal(1,1) = Ddot_(nDTri,DTmp,1,HmatOld,1)
  ExpVal(2,1) = Ddot_(nDTri,DTmp,1,Vmat,1)
  ExpVal(3,1) = Ddot_(nDTri,DTmp,1,Vpolmat,1)
  ExpVal(4,1) = Ddot_(nDTri,DTmp,1,Smat,1)
  call mma_deallocate(DTmp)

else
  ! Shit happens.

  write(u6,*)
  write(u6,*) ' Now how did this happen, says Expectus!'
  call Quit(_RC_INTERNAL_ERROR_)
end if

! What's you major malfunction, numb nuts!

return

end subroutine Expectus
