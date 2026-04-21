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
! Copyright (C) 2026, Lila Zapp                                        *
!                                                                      *
!***********************************************************************

subroutine GetNumGrad_PM(CMO,nOrb2Loc,nBasis,fsdim,NumGrad)

! computes the numerical Gradient of the Pipek-Mezey functional

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp,iwp,u6
use Constants, only: Zero
use Localisation_globals, only: nAtoms

implicit none
real(kind=wp), intent(inout) :: NumGrad(fsdim)
integer(kind=iwp), intent(in) :: fsdim,nBasis,nOrb2Loc
real(kind=wp), intent(in) :: CMO(nBasis,nOrb2Loc)

real(kind=wp),allocatable :: infDisp(:), diff(:),DispMat(:,:),rotated_CMO(:,:),gref(:),infUmat(:,:)
real(kind=wp) :: dx,fref,fpdx,fmdx,fp2dx,fm2dx,GradNorm,PA(nOrb2Loc,nOrb2Loc,nAtoms)
integer(kind=iwp) :: i,NumGradMeth
logical:: debug=.false., debug2=.true.
integer(kind=iwp), parameter ::  fourpoint=1,symm=2,asymm=3

call mma_allocate(infDisp,fsdim,label="infDisp")
call mma_allocate(diff,fsdim,label="diff")
call mma_allocate(DispMat, fsdim,fsdim, Label ="DispMat")
call mma_allocate(infUmat, nOrb2Loc,nOrb2Loc, Label ="infUmat")
call mma_allocate(rotated_CMO, nBasis,nOrb2Loc, Label ="rotated_CMO")
call mma_allocate(gref,fsdim,label="gref")


!choose method based on cost and accuracy + choose adequate dx
NumGradMeth = fourpoint
dx = 1.0e-4_wp ! 1e-4 is good for fourpoint; decrease dx for the other methods


! get Func and analytical grad at x=0
call GenerateP(CMO,nBasis,nOrb2Loc,nAtoms,PA)
call ComputeFunc(nAtoms,nOrb2Loc,PA,fref,.false.)
call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,gref)

NumGrad(:) = Zero

do i = 1,fsdim

    select case(NumGradMeth)

    case (asymm)
        ! get Func at x + dx
        infDisp(:) = Zero
        infDisp(i) = dx
        call vec2upper_triag(DispMat(:,:),nOrb2Loc,infDisp(:),fsdim,.true.)
        call RotateNxN(CMO,DispMat,nOrb2Loc,nBasis,infUmat,rotated_CMO)
        call GenerateP(rotated_CMO,nBasis,nOrb2Loc,nAtoms,PA)
        call ComputeFunc(nAtoms,nOrb2Loc,PA,fpdx,.false.)

        ! compute numerical Gradient
        NumGrad(i)=(fpdx-fref)/dx

   case (symm)
        ! get Func at x + dx
        infDisp(:) = Zero
        infDisp(i) = dx
        call vec2upper_triag(DispMat(:,:),nOrb2Loc,infDisp(:),fsdim,.true.)
        call RotateNxN(CMO,DispMat,nOrb2Loc,nBasis,infUmat,rotated_CMO)
        call GenerateP(rotated_CMO,nBasis,nOrb2Loc,nAtoms,PA)
        call ComputeFunc(nAtoms,nOrb2Loc,PA,fpdx,.false.)

        ! get Func at x - dx
        infDisp(:) = Zero
        infDisp(i) = -dx
        call vec2upper_triag(DispMat(:,:),nOrb2Loc,infDisp(:),fsdim,.true.)
        call RotateNxN(CMO,DispMat,nOrb2Loc,nBasis,infUmat,rotated_CMO)
        call GenerateP(rotated_CMO,nBasis,nOrb2Loc,nAtoms,PA)
        call ComputeFunc(nAtoms,nOrb2Loc,PA,fmdx,.false.)

        ! compute numerical Gradient
        NumGrad(i)=(fpdx-fmdx)/(2*dx)

   case (fourpoint)
        ! get Func at x + dx
        infDisp(:) = Zero
        infDisp(i) = dx
        call vec2upper_triag(DispMat(:,:),nOrb2Loc,infDisp(:),fsdim,.true.)
        call RotateNxN(CMO,DispMat,nOrb2Loc,nBasis,infUmat,rotated_CMO)
        call GenerateP(rotated_CMO,nBasis,nOrb2Loc,nAtoms,PA)
        call ComputeFunc(nAtoms,nOrb2Loc,PA,fpdx,.false.)

        ! get Func at x - dx
        infDisp(:) = Zero
        infDisp(i) = -dx
        call vec2upper_triag(DispMat(:,:),nOrb2Loc,infDisp(:),fsdim,.true.)
        call RotateNxN(CMO,DispMat,nOrb2Loc,nBasis,infUmat,rotated_CMO)
        call GenerateP(rotated_CMO,nBasis,nOrb2Loc,nAtoms,PA)
        call ComputeFunc(nAtoms,nOrb2Loc,PA,fmdx,.false.)

        ! get Func at x + 2dx
        infDisp(:) = Zero
        infDisp(i) = 2*dx
        call vec2upper_triag(DispMat(:,:),nOrb2Loc,infDisp(:),fsdim,.true.)
        call RotateNxN(CMO,DispMat,nOrb2Loc,nBasis,infUmat,rotated_CMO)
        call GenerateP(rotated_CMO,nBasis,nOrb2Loc,nAtoms,PA)
        call ComputeFunc(nAtoms,nOrb2Loc,PA,fp2dx,.false.)

        ! get Func at x - 2dx
        infDisp(:) = Zero
        infDisp(i) = -2*dx
        call vec2upper_triag(DispMat(:,:),nOrb2Loc,infDisp(:),fsdim,.true.)
        call RotateNxN(CMO,DispMat,nOrb2Loc,nBasis,infUmat,rotated_CMO)
        call GenerateP(rotated_CMO,nBasis,nOrb2Loc,nAtoms,PA)
        call ComputeFunc(nAtoms,nOrb2Loc,PA,fm2dx,.false.)

        ! compute numerical Gradient
        NumGrad(i)=(8*(fpdx-fmdx)-fp2dx+fm2dx)/(12*dx)

    end select

    if (debug) then
        call RecPrt('infDisp',' ',DispMat,nOrb2Loc,nOrb2Loc)
        call RecPrt('infUmat',' ',infUmat,nOrb2Loc,nOrb2Loc)
        call RecPrt('rotated_CMO-CMO',' ',rotated_CMO-CMO,nBasis,nOrb2Loc)
        write(u6,*) "i, NumGrad(i) =",i, NumGrad(i)
    end if

end do

if (debug2) then
    ! print results
    write(u6,*) "Grad dx =", dx
    call RecPrt('Analytical Gradient',' ',gref(:),fsdim,1)
    call RecPrt('Numerical Gradient',' ',NumGrad(:),fsdim,1)
    diff(:) = Zero
    diff(:) = NumGrad(:) - gref(:)
    !call RecPrt('Difference',' ',NumGrad(:)-gref(:),fsdim,1)
    write(u6,*) "Grad diff norm", sqrt(dot_product(diff,diff))
end if

call mma_Deallocate(infDisp)
call mma_Deallocate(diff)
call mma_Deallocate(DispMat)
call mma_Deallocate(infUmat)
call mma_Deallocate(rotated_CMO)
call mma_Deallocate(gref)



end subroutine GetNumGrad_PM

