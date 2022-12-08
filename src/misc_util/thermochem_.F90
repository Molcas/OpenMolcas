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

subroutine ThermoChem_(UserT,UserP,TotalM,TRotA,TRotB,TRotC,nUserPT,nsRot,iMult,EVal,in_nFreq,lSlapaf)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half, auTocm, auTokcalmol, auTokJ, kBoltzmann
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: UserT(64), UserP, TRotA, TRotB, TRotC
real(kind=wp), intent(in) :: TotalM, EVal(*)
integer(kind=iwp), intent(inout) :: nUserPT
integer(kind=iwp), intent(in) :: nsRot, iMult, in_nFreq
logical(kind=iwp), intent(in) :: lSlapaf
#include "Molcas.fh"
integer(kind=iwp) :: i, nFreq, nTR, nTr2
real(kind=wp) :: dFreqI, Energy, ZPVE
real(kind=wp), allocatable :: Freq(:), VibT(:)
real(kind=wp), parameter :: rk = kBoltzmann*1.0e-3_wp/auTokJ ! Boltzmann constant in a.u./ K

! If nUserPT == 0 (no User-defined Pressure and Temperatures)
! Compute thermodynamic data at 1.00 atm & 298.15 Kelvin.

if (nUserPT == 0) then
  nUserPT = 1
  UserP = One
  UserT(1) = 298.15_wp
end if
call Get_dScalar('Last energy',Energy)

! Is the system linear?

nTR = 6
if (TRotA > 1.0e99_wp) then
  TRotA = Zero
  nTR = nTR-1
end if
if (TRotB > 1.0e99_wp) then
  TRotB = Zero
  nTR = nTR-1
end if
if (TRotC > 1.0e99_wp) then
  TRotC = Zero
  nTR = nTR-1
end if

! Remove translational and rotational frequencies

call mma_allocate(Freq,in_nFreq,label='Freq')
nFreq = 0
nTr2 = nTr
if (lSlapaf) nTr2 = 0
do i=1,in_nFreq
  dFreqI = EVal(i)
  if (dFreqI > 20.0_wp) then
    nFreq = nFreq+1
    Freq(nFreq) = dFreqI
  end if
end do
Freq(nFreq+1:) = Zero
if ((in_nFreq-nFreq-nTR2) > 0) write(u6,*) ' *** Warning: ',(in_nFreq-nFreq-nTR2),' vibrational contributions removed.'

! Convert frequencies from cm-1 to hartree

call mma_allocate(VibT,nFreq,label='VibT')
ZPVE = Zero
do i=1,nFreq
  Freq(i) = Freq(i)/auTocm ! = 219474.625
  VibT(i) = Freq(i)/rk
  ZPVE = ZPVE+Freq(i)*Half
end do
write(u6,'(A)') ' Vibrational temperature (K): '
do i=1,nFreq,5
  if (i <= (nFreq-4)) write(u6,'(1X,5F9.2)') VibT(i:i+4)
  if (i == (nFreq-3)) write(u6,'(1X,4F9.2)') VibT(i:i+3)
  if (i == (nFreq-2)) write(u6,'(1X,3F9.2)') VibT(i:i+2)
  if (i == (nFreq-1)) write(u6,'(1X,2F9.2)') VibT(i:i+1)
  if (i == (nFreq)) write(u6,'(1X, F9.2)') VibT(i)
end do
write(u6,'(A,I2)') ' Number of trans. and rot. degrees of freedom: ',nTR
write(u6,'(A,F9.3,A,F9.6,A)') ' ZPVE             ',ZPVE*auTokcalmol,' kcal/mol     ',ZPVE,' au.'
write(u6,'(A,13X,F15.6,A)') ' ZPVE corrected energy',ZPVE+Energy,' au.'

do i=1,nUserPT
  call Thermo_VibG(nFreq,Freq,UserT(i),UserP,TotalM,nTR,nsRot,TRotA,TRotB,TRotC,iMult,Energy)
end do
call mma_deallocate(Freq)
call mma_deallocate(VibT)

return

end subroutine ThermoChem_
