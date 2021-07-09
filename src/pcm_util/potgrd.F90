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
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nGrad
real(kind=wp), intent(out) :: Grad(nGrad)
integer(kind=iwp) :: ii, iIrrep, ip1, ipC, iPrint, iRout, nComp, nDens, nOrdOp
real(kind=wp) :: TCpu1, TCpu2, TWall1, TWall2
logical(kind=iwp) :: DiffOp
character(len=80) :: Label
character(len=8) :: Method
real(kind=wp), allocatable :: D_Var(:)
external :: PCMGrd1, PCMMmg
#include "Molcas.fh"
#include "disp.fh"
#include "print.fh"
#include "WrkSpc.fh"

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
! Work(ip1): lOper of each component of the operator

nOrdOp = 0
nComp = (nOrdOp+1)*(nOrdOp+2)/2
call GetMem('Coor','Allo','Real',ipC,3)
call GetMem('lOper','Allo','Inte',ip1,nComp)
call dcopy_(3,[Zero],0,Work(ipC),1)
iWork(ip1) = 1
DiffOp = .true.
call OneEl_g_mck(PCMGrd1,PCMMmG,Grad,nGrad,DiffOp,Work(ipC),D_Var,nDens,iWork(ip1),nComp,nOrdOp,Label)
call PrGrad_mck(' TEST (PCM) contribution',Grad,nGrad,ChDisp,5)
!                                                                      *
!***********************************************************************
!                                                                      *
call GetMem('lOper','Free','Inte',ip1,nComp)
call GetMem('Coor','Free','Real',ipC,3*nComp)
!                                                                      *
!***********************************************************************
!                                                                      *
! Epilogue, end

call mma_deallocate(D_Var)

call CWTime(TCpu2,TWall2)
call SavTim(3,TCpu2-TCpu1,TWall2-TWall1)

return

end subroutine PotGrd
