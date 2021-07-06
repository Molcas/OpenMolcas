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

subroutine PCM_Init(iPrint,ICharg,NAtm,ToAng,AtmC,IAtm,LcAtmC,LcIAtm,nIrrep,NonEq)

use PCM_arrays, only: Centr, dCntr, dPnt, dRad, dTes, IntSph, NewSph, nVert, PCM_N, PCMDM, PCMiSph, PCMSph, PCMTess, Vert
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iPrint, ICharg, NAtm, IAtm(NAtm), nIrrep
real(kind=wp), intent(in) :: ToAng, AtmC(3,NAtm)
real(kind=wp), intent(out) :: LcAtmC(3,NAtm)
integer(kind=iwp), intent(out) :: LcIAtm(NAtm)
logical(kind=iwp), intent(in) :: NonEq
integer(kind=iwp) :: I, ip_RM, ip_SDM, ip_SM, ip_TM, LcI, LcNAtm, nTs2
real(kind=wp) :: Eps_, RJunk(1), TAbs
integer(kind=iwp), allocatable :: pNs(:), VTS(:)
real(kind=wp), allocatable :: Xs(:), Ys(:), Zs(:), Rs(:)
#include "WrkSpc.fh"
#include "rctfld.fh"

! Build the cavity.
! Write the input file for GeomView.
! Form the PCM matrix.

! Possibly print parameter values
if (iPrint >= 99) then
  write(u6,'(a)') 'PCM parameters'
  do I=1,100
    write(u6,'("ISlpar(",i3,") =",i6)') I,ISlPar(I)
  end do
  do I=1,100
    write(u6,'("RSlpar(",i3,") =",F8.3)') I,RSlPar(I)
  end do
end if

! Recover solvent data

call DataSol(ISlPar(15))

! It is necessary to avoid spurious "atoms" which can cause errors
! in UATM: then let's copy only the real atoms on local arrays
! of coordinates and atomic numbers.

LcI = 0
do I=1,NAtm
  if (IAtm(I) > 0) then
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
pNs(:) = 0

NSinit = 0
call FndSph(LcNAtm,ICharg,ToAng,LcAtmC,LcIAtm,ISlPar(9),ISlPar(14),RSlPar(9),Xs,Ys,Zs,Rs,pNs,iPrint)

! Define surface tesserae

call FndTess(iPrint,ToAng,LcNAtm,Xs,Ys,Zs,Rs,pNs,MxSph)

call mma_deallocate(pNs)
call mma_deallocate(Rs)
call mma_deallocate(Zs)
call mma_deallocate(Ys)
call mma_deallocate(Xs)

! Prepare an input file for GeomView visualization tool

call mma_allocate(VTS,MxVert*nTs,Label='VTS')
call GVWrite(1,nTs,NSinit,LcNAtm,LcAtmC,LcIAtm,PCMSph,PCMTess,NVert,Vert,PCMiSph,RJunk,VTS,MxVert)
call mma_deallocate(VTS)

! If needed compute the geometrical derivatives

if (DoDeriv) then
  RSolv = RSlPar(19)
  call Deriva(0,ToAng,LcNAtm,nTs,nS,nSInit,RSolv,PCMTess,Vert,Centr,PCMSph,PCMiSph,IntSph,PCM_N,NVert,NewSph,DTes,dPnt,dRad,dCntr)
end if

! Compute cavitation energy

TAbs = RSlPar(16)
call Cavitation(DoDeriv,ToAng,LcNAtm,NS,nTs,RSlPar(46),VMol,TAbs,RSolv,PCMSph,PCMTess,PCMiSph)

! Define PCM matrix: the inverse is stored in PCMDM

nTs2 = nTs*nTs
call GetMem('SMat','Allo','Real',ip_SM,nTs2)
call GetMem('SDMat','Allo','Real',ip_SDM,nTs2)
call GetMem('TMat','Allo','Real',ip_TM,nTs2)
call GetMem('RMat','Allo','Real',ip_RM,nTs2)
if (NonEq) then
  Eps_ = EpsInf
else
  Eps_ = Eps
end if
call MatPCM(nTs,Eps_,Conductor,PCMiSph,PCMSph,PCMTess,PCMDM,Work(ip_SM),Work(ip_SDM),Work(ip_TM),Work(ip_RM))
call GetMem('RMat','Free','Real',ip_RM,nTs2)
call GetMem('TMat','Free','Real',ip_TM,nTs2)
call GetMem('SDMat','Free','Real',ip_SDM,nTs2)
call GetMem('SMat','Free','Real',ip_SM,nTs2)

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(nIrrep)
end if

end subroutine PCM_Init
