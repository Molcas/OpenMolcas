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

subroutine GF_on_the_Fly(iDo_dDipM)

use Symmetry_Info, only: nIrrep
use Slapaf_Info, only: Coor, nDimBC, mTROld
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iDo_dDipM
integer(kind=iwp) :: iCtl, iEl, iNeg, iOff, jSym, lModes, Lu_10, mDisp(8), mSym, mTR, nAtom, nDoF, nInter, nModes, nX
real(kind=wp) :: DipM(3)
real(kind=wp), allocatable :: dDipM(:), EVal(:), EVec(:), IRInt(:), RedMas(:), Temp(:)
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *

nX = 3*size(Coor,2)
nAtom = size(Coor,2)
nInter = nDimBC-mTROld
mTR = mTROld
nDoF = nDimBC

call mma_allocate(EVec,nX**2,Label='EVec')
call mma_allocate(EVal,nX,Label='EVal')
call mma_allocate(RedMas,nX,Label='RedMas')
call mma_allocate(dDipM,3*nDoF,Label='dDipM')
dDipM(:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
DipM(:) = Zero
call GF(nX,nDoF,nInter,EVec,EVal,RedMas,iNeg,dDipM,mTR,nAtom,DipM)

#ifdef _DEBUGPRINT_
call RecPrt('EVec',' ',EVec,nX,nX)
call RecPrt('EVal',' ',EVal,1,nX)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Write out frequencies and IR intensities

write(u6,10)
write(u6,10) ' Observe that the harmonic oscillator analysis is only valid at stationary points!'
write(u6,10)
write(u6,10) ' Note that rotational and translational degrees have been automatically removed,'
write(u6,10) ' if the energy is invariant to these degrees of freedom.'
write(u6,10)
write(u6,10)
write(u6,10) ' Harmonic frequencies in cm-1'
write(u6,10)
if (iDo_dDipM == 1) then
  write(u6,10) ' IR Intensities in km/mol'
  write(u6,10)
end if
write(u6,10) ' Normal modes in gf_on_the_fly.f '

iOff = 0
iCtl = iDo_dDipM
iEl = 3
jSym = 1

call mma_allocate(Temp,3*nDoF,Label='Temp')

call DGeTMO(dDipM,3,3,nInter,Temp,nInter)

call mma_deallocate(dDipM)

call mma_allocate(IRInt,nDoF,Label='IRInt')

Lu_10 = 10
Lu_10 = IsFreeUnit(Lu_10)
call Molcas_Open(lu_10,'UNSYM')

write(Lu_10,'(A,I1)') '*NORMAL MODES SYMMETRY: ',jsym
call GF_Print(EVal,EVec,Temp,iEl,nDoF,nInter,iCtl,IRInt,RedMas,Lu_10,iOff)

close(Lu_10)
call mma_deallocate(Temp)

call Add_Info('Approx. Freq.',EVal,nInter,1)
!                                                                      *
!***********************************************************************
!                                                                      *
! Save normal modes for later generation of Molden input.

lModes = 0
nModes = 0
lModes = lModes+nInter*nDoF
nModes = nModes+nInter
#ifdef _DEBUGPRINT_
call RecPrt('NModes',' ',EVec,nDoF,nInter)
#endif

if (nIrrep == 1) call Print_Mode_Components(EVec,EVal,nModes,lModes,mDisp)
!                                                                      *
!***********************************************************************
!                                                                      *
! Write stuff on Molden input file

mSym = 1
mDisp(:) = 0
mDisp(1) = nInter
call Freq_Molden(EVal,nModes,EVec,lModes,mSym,IRInt,mDisp,RedMas)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(IRInt)
call mma_deallocate(RedMas)
call mma_deallocate(EVal)
call mma_deallocate(EVec)
!                                                                      *
!***********************************************************************
!                                                                      *
return

10 format(a)

end subroutine GF_on_the_Fly
