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

subroutine Thermo_VibG(nFreq,Freq,T,P,TotalM,nTR,nsRot,TRotA,TRotB,TRotC,iMult,Energy)

use Constants, only: Zero, One, Two, Half, OneHalf, Pi, auTokcalmol, auTokJ, atmToPa, cal_to_J, kBoltzmann, Rgas, rNAVO, rPlanck
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nFreq, nTR, nsRot, iMult
real(kind=wp), intent(in) :: Freq(nFreq), T, P, TotalM, TRotA, TRotB, TRotC, Energy
integer(kind=iwp) :: i
real(kind=wp) :: beta, Const, CQ_e, dFact, dG_TOT, dH_TOT, dMT, dS_e, dS_rot, dS_TOT, dS_tr, dS_vib, dS_vib_Tot, dU_e, dU_rot, &
                 dU_TOT, dU_tr, dU_vib, dU_vib_Tot, dVM, eta, q_e, q_rot, q_TOT, q_tr, q_vib, q_vib_Tot, ZPVE
real(kind=wp), parameter :: R_gas_kcal = 1.0e-3_wp*Rgas/cal_to_J, rk = kBoltzmann*1.0e-3_wp/auTokJ ! Boltzmann constant in a.u./ K

q_e = One
q_tr = One
q_rot = One
q_vib_Tot = One
dS_e = Zero
dS_tr = Zero
dS_rot = Zero
dS_vib_Tot = Zero
dU_e = Zero
dU_tr = Zero
dU_rot = Zero
dU_vib_Tot = Zero
ZPVE = Zero
q_TOT = One
dS_TOT = Zero
dU_TOT = Zero
dH_TOT = Zero
dG_TOT = Zero

! Electronic Contributions

! Molecular Partition Function
q_e = iMult
! Canonical Partition Function
CQ_e = q_e
! Entropy
dS_e = R_gas_kcal*log(CQ_e)
! Thermal
dU_e = Zero

! Translational Contributions

if (T > Zero) then
  ! Molecular Partition Function (q/V in 1/m^3)
  dFact = (Two*Pi*kBoltzmann/(rPlanck**2*rNAVO*1.0e3_wp))**OneHalf ! ((2*PI*k_B)/(h^2*N_A*1000))^(3/2)
  dMT = (TotalM*1.0e-3_wp)*T  ! m = PM/1000
  q_tr = dFact*dMT
  q_tr = q_tr*sqrt(dMT)
  ! Entropy
  dVM = Rgas/atmToPa*T/P ! Gas Molar Volume / m**3
  Const = log(dFact/rNAVO*1.0e3_wp**OneHalf)+2.5_wp ! ~ 18.605
  dS_tr = 1.0e3_wp*R_gas_kcal*(log(dVM)+OneHalf*(log(T)+log(TotalM*1.0e-3_wp))+Const)
  ! Thermal
  dU_tr = OneHalf*R_gas_kcal*T
end if

! Rotational Contributions

if (T > Zero) then
  ! Molecular Partition Function
  if (nTR == 5) then
    q_rot = T/(nsRot*TRotC)
  else if (nTR == 3) then
    q_rot = One
  else
    q_rot = T*T*T
    q_rot = q_rot/TRotA
    q_rot = q_rot/TRotB
    q_rot = q_rot/TRotC
    q_rot = PI*q_rot
    q_rot = sqrt(q_rot)
    q_rot = q_rot/nsRot
  end if
end if
! Entropy
if (nTR == 5) then
  dS_rot = 1.0e3_wp*R_gas_kcal*(log(q_rot)+One)
else if (nTR == 3) then
  dS_rot = Zero
else
  dS_rot = 1.0e3_wp*R_gas_kcal*(log(q_rot)+OneHalf)
end if
! Thermal
if (nTR == 5) then
  dU_rot = R_gas_kcal*T
else if (nTr == 3) then
  dU_rot = Zero
else
  dU_rot = OneHalf*R_gas_kcal*T
end if

! Vibrational Contributions

if (T == Zero) then
  beta = 1.0e99_wp
else
  beta = One/(rk*T)
end if
do i=1,nFreq
  q_vib = One
  dU_vib = Zero
  dS_vib = Zero
  eta = Freq(i)
  if (eta > Zero) then
    ZPVE = ZPVE+eta*Half
    if (T == Zero) then
      q_vib = One
      dU_vib = eta*Half
      dS_vib = Zero
    else
      ! Eq. (6-20)
      q_vib = exp(-eta*beta*Half)/(One-exp(-eta*beta))
      dU_vib = (eta*Half)+(eta/(exp(eta*beta)-One))
      dS_vib = (eta*beta)/(exp(eta*beta)-One)
      dS_vib = dS_vib-log(One-exp(-eta*beta))
    end if
    q_vib_Tot = q_vib_Tot*q_vib
    dU_vib_Tot = dU_vib_Tot+dU_vib
    dS_vib_Tot = dS_vib_Tot+dS_vib
  end if
end do
dU_vib_Tot = dU_vib_Tot*auTokcalmol
dS_vib_Tot = dS_vib_Tot*R_gas_kcal*1.0e3_wp

q_TOT = q_e*q_tr*q_rot*q_vib_Tot
dS_TOT = dS_e+dS_tr+dS_rot+dS_vib_Tot
dU_TOT = dU_e+dU_tr+dU_rot+dU_vib_Tot
dH_TOT = dU_TOT+R_gas_kcal*T
dG_TOT = dH_TOT-dS_TOT*T*1.0e-3_wp

! Print Results

write(u6,*)
write(u6,'(A)') ' *****************************************************'
write(u6,'(A,F8.2,A,F7.2,A)') ' Temperature = ',T,' Kelvin, Pressure =',P,' atm'
write(u6,'(A)') ' -----------------------------------------------------'
write(u6,'(A)') ' Molecular Partition Function and Molar Entropy:'
write(u6,'(A)') '                        q/V (M**-3)    S(kcal/mol*K)'
write(u6,'(A,D17.6,F13.3)') ' Electronic       ',q_e,dS_e
write(u6,'(A,D17.6,F13.3)') ' Translational    ',q_tr,dS_tr
write(u6,'(A,D17.6,F13.3)') ' Rotational       ',q_rot,dS_rot
write(u6,'(A,D17.6,F13.3)') ' Vibrational      ',q_vib_Tot,dS_vib_Tot
write(u6,'(A,D17.6,F13.3)') ' TOTAL            ',q_TOT,dS_TOT

write(u6,*)
write(u6,'(A)') ' Thermal contributions to INTERNAL ENERGY:'
write(u6,'(A,F9.3,A,F9.6,A)') ' Electronic       ',dU_e,' kcal/mol     ',dU_e/auTokcalmol,' au.'
write(u6,'(A,F9.3,A,F9.6,A)') ' Translational    ',dU_tr,' kcal/mol     ',dU_tr/auTokcalmol,' au.'
write(u6,'(A,F9.3,A,F9.6,A)') ' Rotational       ',dU_rot,' kcal/mol     ',dU_rot/auTokcalmol,' au.'
write(u6,'(A,F9.3,A,F9.6,A)') ' Vibrational      ',dU_vib_Tot,' kcal/mol     ',dU_vib_Tot/auTokcalmol,' au.'
write(u6,'(A,F9.3,A,F9.6,A)') ' TOTAL            ',dU_TOT,' kcal/mol     ',dU_TOT/auTokcalmol,' au.'
write(u6,*)
write(u6,'(A)') ' Thermal contributions to'
write(u6,'(A,F9.3,A,F9.6,A)') ' ENTHALPY         ',dH_TOT,' kcal/mol     ',dH_TOT/auTokcalmol,' au.'
write(u6,'(A,F9.3,A,F9.6,A)') ' GIBBS FREE ENERGY',dG_TOT,' kcal/mol     ',dG_TOT/auTokcalmol,' au.'
write(u6,*)
write(u6,'(A)') ' Sum of energy and thermal contributions'
write(u6,'(A,17X,F15.6,A)') ' INTERNAL ENERGY  ',Energy+dU_TOT/auTokcalmol,' au.'
write(u6,'(A,17X,F15.6,A)') ' ENTHALPY         ',Energy+dH_TOT/auTokcalmol,' au.'
write(u6,'(A,17X,F15.6,A)') ' GIBBS FREE ENERGY',Energy+dG_TOT/auTokcalmol,' au.'
write(u6,'(A)') ' -----------------------------------------------------'

return

end subroutine Thermo_VibG
