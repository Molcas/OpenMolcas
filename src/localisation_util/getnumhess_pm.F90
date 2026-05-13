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

subroutine GetNumHess_PM(CMO,nOrb2Loc,nBasis,fsdim,NumHessSymm,debug2)

! computes the numerical Hessian of the Pipek-Mezey functional

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp,iwp, u6
use Constants, only: Zero,Half
use Localisation_globals, only: nAtoms

implicit none

integer(kind=iwp), intent(in) :: nBasis, nOrb2Loc,fsdim
real(kind=wp), intent(in) :: CMO(nBasis,nOrb2Loc)
real(kind=wp),intent(inout) :: NumHessSymm(fsdim,fsdim)


real(kind=wp),allocatable :: infDisp(:), NumHess(:,:),diff(:),gref(:),gpdx(:),rotated_CMO(:,:),oldPA(:,:,:),&
                             gmdx(:),gp2dx(:),gm2dx(:), NumHdiag(:),href(:),DispMat(:,:), Umat(:,:)
real(kind=wp) :: dx,GradNorm,PA(nOrb2Loc,nOrb2Loc,nAtoms)
integer(kind=iwp) :: i,k,l,NumHessMeth
logical:: debug=.false.
logical, intent(in) :: debug2
integer(kind=iwp), parameter ::  fourpoint=1,symm=2,asymm=3

call mma_allocate(NumHess, fsdim, fsdim,Label = "NumHess")
call mma_allocate(infDisp, fsdim,Label ="infDisp")
call mma_allocate(diff, fsdim,Label ="diff")
call mma_allocate(gref, fsdim,Label ="gref")
call mma_allocate(gpdx, fsdim,Label ="gpdx")
call mma_allocate(gmdx, fsdim,Label ="gmdx")
call mma_allocate(gp2dx, fsdim,Label ="gp2dx")
call mma_allocate(gm2dx, fsdim,Label ="gm2dx")
call mma_allocate(NumHdiag, fsdim,Label ="NumHdiag")
call mma_allocate(DispMat, fsdim,fsdim, Label ="DispMat")
call mma_allocate(href, fsdim, Label ="href")
call mma_allocate(Umat, nOrb2Loc,nOrb2Loc, Label ="Umat")
call mma_allocate(oldPA, nOrb2Loc,nOrb2Loc,nAtoms, Label ="oldPA")
call mma_allocate(rotated_CMO, nBasis,nOrb2Loc, Label ="rotated_CMO")


!choose method based on cost and accuracy + choose adequate dx
NumHessMeth = asymm
dx = 1.0e-8_wp ! 1e-4 is good for fourpoint; decrease dx for the other methods


! get grad and analytical Hdiag at x=0
call generateP(CMO,nBasis,nOrb2Loc,nAtoms,PA)
oldPA(:,:,:) = PA(:,:,:)
call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,gref)
call GetHdiag_PM(nAtoms,nOrb2Loc,PA, href(:),.false.) !dont modify it

NumHess(:,:) = Zero
NumHdiag(:) = Zero
NumHessSymm(:,:) = Zero
diff(:) = Zero

do i = 1,fsdim

    select case(NumHessMeth)

    case (asymm)
        ! get Grad at x + dx
        infDisp(:) = Zero
        infDisp(i) = dx
        call takestep()
        call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,gpdx(:))

        ! compute numerical hessian columnwise
        NumHess(:,i)=(gpdx(:)-gref(:))/dx

   case (symm)
        ! get Grad at x + dx
        infDisp(:) = Zero
        infDisp(i) = dx
        call takestep()
        call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,gpdx(:))

        ! get Func at x - dx
        infDisp(:) = Zero
        infDisp(i) = -dx
        call takestep()
        call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,gmdx(:))

        ! compute numerical hessian columnwise
        NumHess(:,i)=(gpdx(:)-gmdx(:))/(2*dx)

    case (fourpoint)
        ! get Grad at x + dx
        infDisp(:) = Zero
        infDisp(i) = dx
        call takestep()
        call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,gpdx(:))

        ! get Func at x - dx
        infDisp(:) = Zero
        infDisp(i) = -dx
        call takestep()
        call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,gmdx(:))

        ! get Grad at x + 2dx
        infDisp(:) = Zero
        infDisp(i) = 2*dx
        call takestep()
        call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,gp2dx(:))

        ! get Func at x - 2dx
        infDisp(:) = Zero
        infDisp(i) = -2*dx
        call takestep()
        call GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,gm2dx(:))

        ! compute numerical hessian columnwise
        NumHess(:,i)=(8*(gpdx-gmdx)-gp2dx+gm2dx)/(12*dx)

    end select

    if (debug) then
        call RecPrt('infDisp',' ',DispMat,nOrb2Loc,nOrb2Loc)
        call RecPrt('Umat',' ',Umat,nOrb2Loc,nOrb2Loc)
        call RecPrt('rotated_CMO-CMO',' ',rotated_CMO-CMO,nBasis,nOrb2Loc)
        write(u6,*) "i, NumHess(:,i) =",i, NumHess(:,i)
    end if

end do

! retrieve the symmetry of the hessian
do k=1,fsdim
    do l=1,fsdim
        NumHessSymm(k,l) = NumHess(k,l)+NumHess(l,k)
    end do
end do
NumHessSymm(:,:) = Half * NumHessSymm(:,:)

do k=1,fsdim
    NumHdiag(k) = NumHessSymm(k,k)
end do

diff(:) = NumHdiag(:) - href(:)

if (debug2) then
    ! print results
    write(u6,*) "Hess dx =", dx
    call RecPrt('Analytical Hdiag',' ',href(:),fsdim,1)
    call RecPrt('Numerical Hdiag',' ',NumHdiag(:),fsdim,1)
    !call RecPrt('Numerical Hessian',' ',NumHess(:,:),fsdim,fsdim)
    !call RecPrt('Numerical Hessian symm.',' ',NumHessSymm(:,:),fsdim,fsdim)
    call RecPrt('Difference',' ',diff,fsdim,1)
    write(u6,*) "Hess diff norm", sqrt(dot_product(diff,diff))
end if

href(:) = -abs(NumHdiag(:))

call mma_Deallocate(NumHess)
call mma_Deallocate(infDisp)
call mma_Deallocate(diff)
call mma_Deallocate(gref)
call mma_Deallocate(gpdx)
call mma_Deallocate(gmdx)
call mma_Deallocate(gp2dx)
call mma_Deallocate(gm2dx)
call mma_Deallocate(NumHdiag)
call mma_Deallocate(DispMat)
call mma_Deallocate(href)
call mma_Deallocate(Umat)
call mma_Deallocate(rotated_CMO)
call mma_Deallocate(oldPA)

contains

subroutine takestep()
    PA = oldPA
    call vec2upper_triag(DispMat,nOrb2Loc,infDisp(:),fsdim,.true.)
    call expkap_localisation(DispMat,nOrb2Loc,Umat)
    call RotateNxN(CMO,nOrb2Loc,nBasis,Umat,rotated_CMO)
    call transformPA(PA,nOrb2Loc,Umat)
end subroutine takestep

end subroutine GetNumHess_PM


