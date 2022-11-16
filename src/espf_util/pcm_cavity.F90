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

subroutine PCM_Cavity(iPrint,ICharg,NAtm,AtmC,IAtm,IsAtMM,LcAtmC,LcIAtm,JSurf)

use PCM_arrays

implicit real*8(a-h,o-z)
#include "espf.fh"
#include "rctfld.fh"
#include "status.fh"
#include "stdalloc.fh"
real*8 AtmC(3,NAtm), LcAtmC(3,NAtm)
integer IAtm(NAtm), IsAtMM(NAtm), LcIAtm(NAtm)
save Rad_Cor, Surf_Inc
data Rad_Cor/2.0d0/,Surf_Inc/0.5d0/
real*8, allocatable :: Xs(:), Ys(:), Zs(:), Rs(:)
integer, allocatable :: pNs(:)

! Build the cavity.

call PCMDef(ISlPar,RSlPar,iPrint)
RSlPar(3) = 5.0d-1
RSlPar(7) = 2.0d0
RSlPar(9) = Rad_Cor+dble(JSurf)*Surf_Inc

! Possibly print parameter values

if (iPrint >= 99) then
  write(6,'(''PCM parameters'')')
  do I=1,100
    write(6,'(''ISlpar('',i3,'') ='',i6)') I,ISlPar(I)
  end do
  do I=1,100
    write(6,'(''RSlpar('',i3,'') ='',F8.3)') I,RSlPar(I)
  end do
end if

! Recover solvent data

call DataSol(ISlPar(15))

! It is necessary to avoid spurious "atoms" which can cause errors
! in UATM: then let's copy only the real atoms on local arrays
! of coordinates and atomic numbers.

LcI = 0
do I=1,NAtm
  if ((IAtm(I) > 0) .and. (IsAtMM(I) == 0)) then
    LcI = LcI+1
    LcAtmC(1,LcI) = AtmC(1,I)
    LcAtmC(2,LcI) = AtmC(2,I)
    LcAtmC(3,LcI) = AtmC(3,I)
    LcIAtm(LcI) = IAtm(I)
  end if
end do
LcNAtm = LcI
ISlPar(42) = LcNAtm

! Define atomic/group spheres
! Allocate space for X, Y, Z, Radius and NOrd for MxSph spheres

call mma_allocate(Xs,MxSph,Label='Xs')
call mma_allocate(Ys,MxSph,Label='Ys')
call mma_allocate(Zs,MxSph,Label='Zs')
call mma_allocate(Rs,MxSph,Label='Rs')
call mma_allocate(pNs,MxSph,Label='pNs')

NSinit = 0
call FndSph(LcNAtm,ICharg,LcAtmC,LcIAtm,ISlPar(9),ISlPar(14),RSlPar(9),Xs,Ys,Zs,Rs,pNs,MxSph,iPrint)

! Define surface tesserae

call FndTess(iPrint,Xs,Ys,Zs,Rs,pNs,MxSph)

call mma_deallocate(pNs)
call mma_deallocate(Rs)
call mma_deallocate(Zs)
call mma_deallocate(Ys)
call mma_deallocate(Xs)

! If needed compute the geometrical derivatives

if (DoDeriv) then
  RSolv = RSlPar(19)
  LcNAtm = ISlPar(42)
  call mma_allocate(dTes,nTs,LcNAtm,3,Label='dTes')
  call mma_allocate(dPnt,nTs,LcNAtm,3,3,Label='dPnt')
  call mma_allocate(dRad,nS,LcNAtm,3,Label='dRad')
  call mma_allocate(dCntr,nS,LcNAtm,3,3,Label='dCntr')
  call mma_allocate(PCM_SQ,2,nTs,Label='PCM_SQ')
  call Deriva(0,LcNAtm,nTs,nS,nSInit,RSolv,PCMTess,Vert,Centr,PCMSph,PCMiSph,IntSph,PCM_N,NVert,NewSph,dTes,dPnt,dRad,dCntr)
  if (nPCM_info == 0) then
    write(6,'(A)') ' GEPOL failed to compute the grid deriv.'
    write(6,'(A)') ' Reduce the number of surfaces.'
    call Quit_OnUserError()
  end if
end if

return

end subroutine PCM_Cavity
