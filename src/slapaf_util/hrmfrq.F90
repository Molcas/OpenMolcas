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

subroutine HrmFrq(nAtom,nInter,iNeg,dDipM,mTR,DipM,IRInt)

use Slapaf_Info, only: lDoubleIso, lTherm, nsRot, nUserPT, UserP, UserT
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtom, nInter, mTR
integer(kind=iwp), intent(out) :: iNeg
real(kind=wp), intent(inout) :: dDipM(3,nInter+mTR)
real(kind=wp), intent(in) :: DipM(3)
real(kind=wp), intent(out) :: IRInt(nInter+mTR)
integer(kind=iwp) :: iCtl, iEl, iOff, iSym, lModes, Lu_10, mDisp(8), nDoF, nIrrep, nModes, nSym, nX
real(kind=wp), allocatable :: EVal(:), EVec(:), RedMas(:), Temp(:)
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
nDoF = nInter+mTR
nX = 3*nAtom

call mma_allocate(EVec,nDoF**2,Label='EVec')
call mma_allocate(EVal,nDoF,Label='EVal')
call mma_allocate(RedMas,nDoF,Label='RedMas')
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute harmonic frequencies and dipole derivatives

call GF(nX,nDoF,nInter,EVec,EVal,RedMas,iNeg,dDipM,mTR,nAtom,DipM)
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
write(u6,10) ' IR Intensities in km/mol'
write(u6,10)

iOff = 0
iCtl = 1
iEl = 3
iSym = 1

call mma_allocate(Temp,3*nDoF,Label='Temp')
call DGeTMO(dDipM,3,3,nInter,Temp,nInter)

Lu_10 = 10
Lu_10 = IsFreeUnit(Lu_10)
call Molcas_Open(lu_10,'UNSYM')

write(Lu_10,'(A,I1)') '*NORMAL MODES SYMMETRY: ',isym
call GF_Print(EVal,EVec,Temp,iEl,nDoF,nInter,iCtl,IRInt,RedMas,Lu_10,iOff)

close(Lu_10)
call mma_deallocate(Temp)

if (lTherm) call Thermo_Driver(UserT,UserP,nUserPT,nsRot,EVal,nInter,lTherm)
!                                                                      *
!***********************************************************************
!                                                                      *
! Do single and double isotope substitutions

call Get_iScalar('NSYM',nIrrep)
if (nIrrep == 1) call IsoLoop(lDoubleIso)
!                                                                      *
!***********************************************************************
!                                                                      *
! Save normal modes for later generation of Molden input.

lModes = 0
nModes = 0
lModes = lModes+nInter*nDoF
nModes = nModes+nInter
!                                                                      *
!***********************************************************************
!                                                                      *
! Write stuff on Molden input file

nSym = 1
mDisp(:) = 0
mDisp(1) = nInter
call Print_Mode_Components(EVec,EVal,nModes,lModes,mDisp)
call Freq_Molden(EVal,nModes,EVec,lModes,nSym,IRInt,mDisp,RedMas)
!                                                                      *
!***********************************************************************
!                                                                      *
! Deallocate memory

call mma_deallocate(RedMas)
call mma_deallocate(EVal)
call mma_deallocate(EVec)
!                                                                      *
!***********************************************************************
!                                                                      *
return

10 format(a)

end subroutine HrmFrq
