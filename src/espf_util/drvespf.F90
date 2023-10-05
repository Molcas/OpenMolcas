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

subroutine Drvespf(Grad,Temp,nGrad,CCoor)
! Driver to compute the ESPF B*dV contribution to the gradient
! This is a hack of the alaska/drvh1 subroutine with a little
! piece of (extinct) integral_util/drvprop subroutine

use Index_Functions, only: nTri_Elem, nTri_Elem1
use Basis_Info, only: nBas
use Grd_interface, only: grd_kernel, grd_mem
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nGrad
real(kind=wp), intent(inout) :: Grad(nGrad)
real(kind=wp), intent(out) :: Temp(nGrad)
real(kind=wp), intent(in) :: CCoor(*)
#include "print.fh"
integer(kind=iwp) :: ii, iIrrep, iPL, iPrint, nComp, nDens, nOrdOp
logical(kind=iwp) :: DiffOp
character(len=80) :: Label
integer(kind=iwp), allocatable :: lOper(:)
real(kind=wp), allocatable :: D_Var(:)
integer(kind=iwp), external :: iPL_espf
procedure(grd_kernel) :: BdVGrd
procedure(grd_mem) :: NAMmG

! Prologue
iPrint = 1

! Allocate memory for density matrix

nDens = 0
do iIrrep=0,nIrrep-1
  nDens = nDens+nTri_Elem(nBas(iIrrep))
end do

! Read the variational 1st order density matrix
! density matrix in AO/SO basis

call mma_allocate(D_Var,nDens,Label='D_Var')
call Get_D1ao_Var(D_Var,nDens)
if (iPrint >= 99) then
  write(u6,*) 'variational 1st order density matrix'
  ii = 1
  do iIrrep=0,nIrrep-1
    write(u6,*) 'symmetry block',iIrrep
    call TriPrt(' ',' ',D_Var(ii),nBas(iIrrep))
    ii = ii+nTri_Elem(nBas(iIrrep))
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! nOrdOp: order/rank of the operator
! lOper: lOper of each component of the operator

nPrint(112) = 5
iPL = iPL_espf()
if (iPL >= 3) nPrint(112) = 15
nOrdOp = 0
nComp = nTri_Elem1(nOrdOp)
call mma_allocate(lOper,nComp,Label='lOper')
lOper(:) = 1
!
!***********************************************************************
!                                                                      *
!     Trace the "variational" first order density matrix with the      *
!     gradient of the external field integrals.                        *
!                                                                      *
!***********************************************************************
DiffOp = .true.
Label = ' The ESPF BdV contribution'
call OneEl_g(BdVGrd,NAMmG,Temp,nGrad,DiffOp,CCoor,D_Var,nDens,lOper,nComp,nOrdOp,Label)
Grad(:) = Grad+Temp
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(lOper)
call mma_deallocate(D_Var)

return

end subroutine Drvespf
