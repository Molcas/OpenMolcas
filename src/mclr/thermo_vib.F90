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

subroutine Thermo_Vib(nFreq,Freq,T,nTR,iter)

use Constants, only: Zero, One, Half, auTokcalmol, auTokJ, cal_to_j, kBoltzmann, rNAVO
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nFreq, nTR, iter
real(kind=wp), intent(in) :: Freq(nFreq), T
integer(kind=iwp) :: i
real(kind=wp) :: beta, DeltaG, DeltaS, DeltaU, eta, q_vib, q_vib_Tot, U_TR, U_vib, U_vib_Tot, ZPVE
integer(kind=iwp), parameter :: First = 1
real(kind=wp), parameter :: rk = kBoltzmann*1.0e-3_wp/auTokJ, rk_kcalmol = kBoltzmann*rNAVO/(1.0e3_wp*cal_to_J)

!r_k = 1.38065800e-23_wp
!r_J2au = 2.29371049e17_wp ! Convert joules to atomic units
! Bolzmann's constant in a.u./ K
!rk = kBoltzmann*r_J2au
!write(u6,*) "Bolzmann's constant=",rK
if (T == Zero) then
  beta = 1.0e99_wp
else
  beta = One/(rk*T)
end if
!***********************************************************************
!                                                                      *
!      The Canonical partition function Q for indistinguishable        *
!      molecules                                                       *
!                                                                      *
!      is Q = q^N /N!, where N is the number of indistinguishable      *
!                                                                      *
!      molecules and q is the function with                            *
!                                                                      *
!      q = q_e * q_v * q_t * q_r                                       *
!                                                                      *
!      N = n * N_A (N_A is Avogrado's number) and the molecular        *
!                                                                      *
!      partition function is defined as q_m=q/n.                       *
!                                                                      *
!      We approximate q as follows                                     *
!                                                                      *
!      electronic:                                                     *
!      q_e=g_e (g_e is the degeneracy factor, for a singlet            *
!               ground state g_e=1, for a triplet ground state         *
!               g_e to a good approximation is 3, etc).                *
!                                                                      *
!      vibrational:                                                    *
!      q_v=q_v(1)*q_v(2)*...                                           *
!                                                                      *
!            q_v(i)= exp(-e(i)*beta/2) / (1-exp(-e(i)*beta))           *
!                                                                      *
!      translational:                                                  *
!      q_t=V/L^3, where L=(h^2 beta/2 pi m)^(1/2)                      *
!                                                                      *
!      d lnq_t/d beta = - 3kT/2                                        *
!                                                                      *
!      rotational: (high temperature approximations!)                  *
!      for linear systems                                              *
!      q_r=1/(beta sigma hc B) sigma=symmetry number                   *
!                                                                      *
!      d lnq_r/d beta = - 2kT/2                                        *
!                                                                      *
!      for nonlinear systems                                           *
!      q_r=(1/sigma)(1/(beta hc))^(3/2) (pi/ABC)^(1/2)                 *
!                                                                      *
!      d lnq_r/d beta = - 3kT/2                                        *
!                                                                      *
!      Some useful relations!                                          *
!                                                                      *
!      The internal energy:                                            *
!      U-U(0)=-(d lnQ/d beta)_V=-N(d lnq/d beta)_V                     *
!                                                                      *
!      Gibbs free energy:                                              *
!      G-G(0)=-nRT ln(q_m/N_A)                                         *
!                                                                      *
!      The enthalpy:                                                   *
!      H-H(0)=G-G(0)+T(S-S(0))                                         *
!                                                                      *
!            =-(d lnQ/d beta)_V + kTV(d lnQ/d V)_T                     *
!                                                                      *
!            =-N(d lnq/d beta)_V + NkTV(d lnQ/d V)_T                   *
!                                                                      *
!            =U-U(0)+pV                                                *
!                                                                      *
!      p=(NkT/q)(dq/dV)_T                                              *
!                                                                      *
!      since the translational part is the only which depends on the   *
!      volume we have that                                             *
!                                                                      *
!      H-H(0)=U-U(0)+NkT                                               *
!                                                                      *
!***********************************************************************

write(u6,*)
write(u6,*)
write(u6,'(A,F6.2,A)') ' Temperature = ',T,' kelvin'
write(u6,'(A)') ' ---------------------------'
write(u6,*)

! Iterate over vibrations, and compute the molecular
! partition functions of the vibrations

q_vib_Tot = One
ZPVE = Zero
U_vib_Tot = Zero
do i=1,nFreq
  eta = Freq(i)
  if (iter == first) write(u6,*) ' Vibrational temperature =',eta/rk
  if (eta > Zero) then
    ZPVE = ZPVE+eta*Half
    if (T == Zero) then
      q_vib = Zero
      U_vib = eta*Half
    else
      ! Eq. (6-20)
      q_vib = exp(-eta*beta*Half)/(One-exp(-eta*beta))
      U_vib = (eta*Half)+(eta/(exp(eta*beta)-One))
    end if
    q_vib_Tot = q_vib_Tot*q_vib
    U_vib_Tot = U_vib_Tot+U_vib
  end if
end do

! Add translational and rotational energy, kT*nTR/2,
! the classiscal limit is used.
!
! U-U(0)=-(d lnQ/d beta)_V, U(0)=0
!
! lnQ = N lnq - lnN!
!
! U-U(0)=-N(d lnq/d beta)_V

!r_J2kcalmol = 1.43932522e20_wp ! conversion J to kcal/mol
U_TR = real(nTR,kind=wp)*(T*rk_kcalmol*Half)

! G-G(0)=-kT lnQ + kTV(d lnQ /dV)_T, G(0)=0

if (T == Zero) then
  DeltaG = Zero
else
  DeltaG = -log(q_vib_Tot)*rk*T
end if

!auTokcalmol = 6.27509541e2_wp ! Conversion au to kcal/mol
DeltaG = DeltaG*auTokcalmol

DeltaU = U_vib_Tot*auTokcalmol
ZPVE = ZPVE*auTokcalmol
write(u6,'(A,F6.2,A)') '         DeltaG =',DeltaG,' kcal/mol'
write(u6,'(A,F6.2,A)') '           ZPVE =',ZPVE,' kcal/mol'
write(u6,'(A,F6.2,A)') '      TotDeltaU =',DeltaU,' kcal/mol'
write(u6,'(A,F6.2,A)') ' TotDeltaU-ZPVE =',DeltaU-ZPVE,' kcal/mol'

if (T > Zero) then
  DeltaS = 1.0e3_wp*(DeltaU-DeltaG)/T
else
  DeltaS = Zero
end if
write(u6,'(A,F8.4,A)') '      Entropy S =',DeltaS,' cal/(mol*K)'
write(u6,'(A,F8.4,A)') '         U(T&R) =',U_TR,' kcal/mol'
write(u6,'(A,F8.4,A)') '       Tot(temp)=',U_TR+DeltaU-ZPVE,' kcal/mol'

return

end subroutine Thermo_Vib
