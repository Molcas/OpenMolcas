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

subroutine Drvh1_EMB(Grad,Temp,nGrad)

use Basis_Info, only: dbsc, nCnttp, nBas
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nGrad
real(kind=wp), intent(inout) :: Grad(nGrad)
real(kind=wp), intent(out) :: Temp(nGrad)
#include "print.fh"
integer(kind=iwp) :: i, ii, iIrrep, iPrint, iRout, nComp, nDens, nOrdOp
real(kind=wp) :: TCpu1, TCpu2, TWall1, TWall2
logical(kind=iwp) :: DiffOp, lECP, lPP, lFAIEMP
character(len=80) :: Label
integer(kind=iwp), allocatable :: lOper(:)
real(kind=wp), allocatable :: Coor(:,:), D_Var(:)
external :: FragPGrd, FragPMmG, M1Grd, M1MmG, M2Grd, M2MmG, NAGrd, NAMmG, PPGrd, PPMmG, PrjGrd, PrjMmG, SROGrd, SROMmG

!...  Prologue
iRout = 131
iPrint = nPrint(iRout)
call CWTime(TCpu1,TWall1)
call StatusLine(' Alaska:',' Computing 1-el OFE gradients')

call Set_Basis_Mode('Valence')
call Setup_iSD()

lECP = .false.
lPP = .false.
lFAIEMP = .false.
do i=1,nCnttp
  lECP = lECP .or. dbsc(i)%ECP
  lPP = lPP .or. (dbsc(i)%nPP /= 0)
  lFAIEMP = LFAIEMP .or. dbsc(i)%Frag
end do

! Allocate memory for density matrices

nDens = 0
do iIrrep=0,nIrrep-1
  nDens = nDens+nBas(iIrrep)*(nBas(iIrrep)+1)/2
end do

! Read the variational 1st order density matrix
! density matrix in AO/SO basis

call NameRun('AUXRFIL') ! switch RUNFILE name

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

! Annihilate all the components of rho_B in the bsfs of the A subsystem

!SVC: fixed according to instructions from Francesco,
!      as embedding should not deal with symmetry
call Annihil_rho(D_var,nBas(0))

call NameRun('#Pop')    ! switch RUNFILE name
!                                                                      *
!***********************************************************************
!                                                                      *
! nOrdOp: order/rank of the operator
! lOper(:): lOper of each component of the operator

nOrdOp = 0
nComp = (nOrdOp+1)*(nOrdOp+2)/2
call mma_allocate(Coor,3,nComp)
call mma_allocate(lOper,nComp,Label='lOper')
Coor(:,:) = Zero
lOper(:) = 1

!***********************************************************************
!3)                                                                    *
!     Trace the "variational" first order density matrix with the      *
!     gradient of the nuclear attraction integrals.                    *
!                                                                      *
!***********************************************************************

DiffOp = .true.
Label = ' The Nuclear Attraction Contribution'
call OneEl_g(NAGrd,NAMmG,Temp,nGrad,DiffOp,Coor,D_Var,nDens,lOper,nComp,nOrdOp,Label)

call DaXpY_(nGrad,One,Temp,1,Grad,1)

!***********************************************************************
!4)                                                                    *
!     Trace the "variational" first order density matrix with the      *
!     gradient of the ECP integrals.                                   *
!                                                                      *
!***********************************************************************

if (lECP) then
  DiffOp = .true.
  Label = ' The Projection Operator contribution'
  call OneEl_g(PrjGrd,PrjMmG,Temp,nGrad,DiffOp,Coor,D_Var,nDens,lOper,nComp,nOrdOp,Label)
  call DaXpY_(nGrad,One,Temp,1,Grad,1)

  Label = ' The M1 Operator contribution'
  call OneEl_g(M1Grd,M1MmG,Temp,nGrad,DiffOp,Coor,D_Var,nDens,lOper,nComp,nOrdOp,Label)
  call DaXpY_(nGrad,One,Temp,1,Grad,1)

  Label = ' The M2 Operator contribution'
  call OneEl_g(M2Grd,M2MmG,Temp,nGrad,DiffOp,Coor,D_Var,nDens,lOper,nComp,nOrdOp,Label)
  call DaXpY_(nGrad,One,Temp,1,Grad,1)

  Label = ' The SR Operator contribution'
  call OneEl_g(SROGrd,SROMmG,Temp,nGrad,DiffOp,Coor,D_Var,nDens,lOper,nComp,nOrdOp,Label)
  call DaXpY_(nGrad,One,Temp,1,Grad,1)
end if
if (lPP) then
  Label = ' The Pseudo Potential contribution'
  call OneEl_g(PPGrd,PPMmG,Temp,nGrad,DiffOp,Coor,D_Var,nDens,lOper,nComp,nOrdOp,Label)
  call DaXpY_(nGrad,One,Temp,1,Grad,1)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (lFAIEMP) then
  DiffOp = .true.
  Label = ' The FAIEMP Projection Operator Contribution'
  call OneEl_g(FragPGrd,FragPMmG,Temp,nGrad,DiffOp,Coor,D_Var,nDens,lOper,nComp,nOrdOp,Label)
  call DaXpY_(nGrad,One,Temp,1,Grad,1)
  call DrvG_FAIEMP(Grad,Temp,nGrad)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(lOper)
call mma_deallocate(Coor)
!                                                                      *
!***********************************************************************
!                                                                      *
! Epilogue, end

call mma_deallocate(D_Var)

call Free_iSD()
call CWTime(TCpu2,TWall2)

return

end subroutine Drvh1_EMB
