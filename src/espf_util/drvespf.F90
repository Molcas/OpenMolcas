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

use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep

implicit real*8(A-H,O-Z)
#include "espf.fh"
external BdVGrd
external NAMmG
#include "print.fh"
#include "disp.fh"
#include "nsd.fh"
#include "setup.fh"
#include "wldata.fh"
#include "stdalloc.fh"
character Label*80
real*8 Grad(nGrad), Temp(nGrad)
real*8 Ccoor(*)
!logical DiffOp, DoRys
logical DiffOp
real*8, allocatable :: D_Var(:)
integer, allocatable :: lOper(:)
! Statement function
nElem(i) = (i+1)*(i+2)/2

! Prologue
iPrint = 1

! Allocate memory for density matrix

nDens = 0
do iIrrep=0,nIrrep-1
  nDens = nDens+nBas(iIrrep)*(nBas(iIrrep)+1)/2
end do

! Read the variational 1st order density matrix
! density matrix in AO/SO basis

call mma_allocate(D_Var,nDens,Label='D_Var')
call Get_D1ao_Var(D_Var,nDens)
if (iPrint >= 99) then
  write(6,*) 'variational 1st order density matrix'
  ii = 1
  do iIrrep=0,nIrrep-1
    write(6,*) 'symmetry block',iIrrep
    call TriPrt(' ',' ',D_Var(ii),nBas(iIrrep))
    ii = ii+nBas(iIrrep)*(nBas(iIrrep)+1)/2
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
nComp = nElem(nOrdOp)
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
call DaXpY_(nGrad,One,Temp,1,Grad,1)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(lOper)
call mma_deallocate(D_Var)

return

end subroutine Drvespf
