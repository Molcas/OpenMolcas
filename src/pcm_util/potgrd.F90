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

subroutine PotGrd(Grad,nGrad)

use Basis_Info, only: nBas
use Symmetry_Info, only: nIrrep
use Index_Functions, only: nTri_Elem1
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nGrad
real(kind=wp), intent(out) :: Grad(nGrad)
integer(kind=iwp) :: ii, iIrrep, iPrint, iRout, lOper(1), nComp, nDens, nOrdOp
real(kind=wp) :: C(3), TCpu1, TCpu2, TWall1, TWall2
logical(kind=iwp) :: DiffOp
character(len=80) :: Label
character(len=8) :: Method
real(kind=wp), allocatable :: D_Var(:)
external :: PCMGrd1, PCMMmg
#include "Molcas.fh"
#include "disp.fh"
#include "print.fh"

! Prologue
iRout = 131
iPrint = nPrint(iRout)
call CWTime(TCpu1,TWall1)

! Allocate memory for density

nDens = 0
do iIrrep=0,nIrrep-1
  nDens = nDens+nBas(iIrrep)*(nBas(iIrrep)+1)/2
end do

! Get the method label
!write(u6,*) ' Read Method label'
call Get_cArray('Relax Method',Method,8)

! Read the variational 1st order density matrix
! density matrix in AO/SO basis
!write(u6,*) ' Read density matrix'
call mma_allocate(D_Var,nDens,Label='D_Var')
call Get_D1ao_Var(D_var,nDens)

if (iPrint >= 99) then
  write(u6,*) 'variational 1st order density matrix'
  ii = 1
  do iIrrep=0,nIrrep-1
    write(u6,*) 'symmetry block',iIrrep
    call TriPrt(' ',' ',D_Var(ii),nBas(iIrrep))
    ii = ii+nBas(iIrrep)*(nBas(iIrrep)+1)/2
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! nOrdOp: order/rank of the operator
! lOper: lOper of each component of the operator

nOrdOp = 0
nComp = nTri_Elem1(nOrdOp)
C(:) = Zero
lOper(1) = 1
DiffOp = .true.
call OneEl_g_pcm(PCMGrd1,PCMMmG,Grad,nGrad,DiffOp,C,D_Var,nDens,lOper,nComp,nOrdOp,Label)
call PrGrad_pcm(' TEST (PCM) contribution',Grad,nGrad,ChDisp,5)
!                                                                      *
!***********************************************************************
!                                                                      *
! Epilogue, end

call mma_deallocate(D_Var)

call CWTime(TCpu2,TWall2)

return

end subroutine PotGrd
