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

subroutine ThermoChem_(UserT,UserP,TotalM,TRotA,TRotB,TRotC,nUserPT,nsRot,iMult,nAtom,EVal,in_nFreq,lSlapaf)

implicit real*8(a-h,o-z)
#include "Molcas.fh"
#include "constants.fh"
#include "WrkSpc.fh"
#include "real.fh"
integer in_nFreq, nFreq, iMult, nUserPT, nsRot, nAtom
real*8 UserT(64), UserP, EVal(*)
real*8 TotalM, TRotA, TRotB, TRotC
real*8 Freq(MxAtom*3-6), VibT(MxAtom*3-6), dFreqI
real*8 Energy
logical lSlapaf

! If nUserPT == 0 (no User-defined Pressure and Temperatures)
! Compute thermodynamic data at  1.00 atm &  298.15 Kelvin.

if (nUserPT == 0) then
  nUserPT = 1
  UserP = 1.0d0
  UserT(1) = 298.15d0
end if
call Get_dScalar('Last energy',Energy)

! Is the system linear?

nTR = 6
if (TRotA > 1.0d99) then
  TRotA = 0.0d0
  nTR = nTR-1
end if
if (TRotB > 1.0d99) then
  TRotB = 0.0d0
  nTR = nTR-1
end if
if (TRotC > 1.0d99) then
  TRotC = 0.0d0
  nTR = nTR-1
end if

! Remove translational and rotational frequencies

nFreq = 0
nTr2 = nTr
if (lSlapaf) nTr2 = 0
do i=1,in_nFreq
  dFreqI = EVal(i)
  if (dFreqI > 20.0d0) then
    nFreq = nFreq+1
    Freq(nFreq) = dFreqI
  end if
end do
if ((in_nFreq-nFreq-nTR2) > 0) write(6,*) ' *** Warning: ',(in_nFreq-nFreq-nTR2),' vibrational contributions removed.'

! Convert frequencies from cm-1 to hartree

r_k = CONST_BOLTZMANN_
r_J2au = 1.0D-3/CONV_AU_TO_KJ_
rk = r_k*r_J2au ! Bolzmann constant in a.u./ K
ZPVE = 0.0d0
do i=1,nFreq
  Freq(i) = Freq(i)/CONV_AU_TO_CM1_ ! = 219474.625
  VibT(i) = Freq(i)/rk
  ZPVE = ZPVE+Freq(i)/2.0d0
end do
write(6,'(A)') ' Vibrational temperature (K): '
do i=1,nFreq,5
  if (i <= (nFreq-4)) write(6,'(1X,5F9.2)') (VibT(j),j=i,i+4)
  if (i == (nFreq-3)) write(6,'(1X,4F9.2)') (VibT(j),j=i,i+3)
  if (i == (nFreq-2)) write(6,'(1X,3F9.2)') (VibT(j),j=i,i+2)
  if (i == (nFreq-1)) write(6,'(1X,2F9.2)') (VibT(j),j=i,i+1)
  if (i == (nFreq)) write(6,'(1X, F9.2)') VibT(i)
end do
write(6,'(A,I2)') ' Number of trans. and rot. degrees of freedom: ',nTR
write(6,'(A,F9.3,A,F9.6,A)') ' ZPVE             ',ZPVE*6.27509541D+2,' kcal/mol     ',ZPVE,' au.'
write(6,'(A,13X,F15.6,A)') ' ZPVE corrected energy',ZPVE+Energy,' au.'

do i=1,nUserPT
  call Thermo_VibG(nFreq,Freq,UserT(i),UserP,TotalM,nTR,nsRot,TRotA,TRotB,TRotC,iMult,Energy)
end do

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nAtom)

end subroutine ThermoChem_
