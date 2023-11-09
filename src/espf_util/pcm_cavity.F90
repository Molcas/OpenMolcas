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

use PCM_arrays, only: Centr, dCntr, dPnt, dRad, dTes, IntSph, MxSph, NewSph, NVert, PCM_N, PCM_SQ, PCMiSph, PCMSph, PCMTess, Vert
use rctfld_module, only: DoDeriv, ISlPar, nPCM_Info, nS, nSInit, nTS, RSlPar, RSolv
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iPrint, ICharg, NAtm, IAtm(NAtm), IsAtMM(NAtm), JSurf
real(kind=wp), intent(in) :: AtmC(3,NAtm)
real(kind=wp), intent(out) :: LcAtmC(3,NAtm)
integer(kind=iwp), intent(out) :: LcIAtm(NAtm)
integer(kind=iwp) :: I, LcI, LcNAtm
integer(kind=iwp), allocatable :: pNs(:)
real(kind=wp), allocatable :: Rs(:), Xs(:), Ys(:), Zs(:)
real(kind=wp), parameter :: Rad_Cor = Two, Surf_Inc = Half

! Build the cavity.

call PCMDef(ISlPar,RSlPar,iPrint)
RSlPar(3) = Half
RSlPar(7) = Two
RSlPar(9) = Rad_Cor+real(JSurf,kind=wp)*Surf_Inc

! Possibly print parameter values

if (iPrint >= 99) then
  write(u6,'(''PCM parameters'')')
  do I=1,100
    write(u6,'(''ISlpar('',i3,'') ='',i6)') I,ISlPar(I)
  end do
  do I=1,100
    write(u6,'(''RSlpar('',i3,'') ='',F8.3)') I,RSlPar(I)
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
    LcAtmC(:,LcI) = AtmC(:,I)
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
    write(u6,'(A)') ' GEPOL failed to compute the grid deriv.'
    write(u6,'(A)') ' Reduce the number of surfaces.'
    call Quit_OnUserError()
  end if
end if

return

end subroutine PCM_Cavity
