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

! Numerical solution of the vibrational Schroedinger equation
! for a diatomic molecule.
!
! PotR is the electronic energy given as input
! Redm is the reduced mass
! Jrot is the rotational quantum number
! G is the radial wave function in a logaritmic scale U=ln(R)
! Umin to Umax is the integration range
! del is the step length
! R=exp(U) is the radial coordinates
! E is the energy
! ndim is the number of integration steps
! E0 is the starting value for an eigenvalue search. E0 is simply
! put equal to the minimum value of the potential.
! nvib is the number of eigenvalues wanted.
! dE0 the original step length in the energy interpolation procedure
!
! ********** MOLCAS-2 Release 90 05 01 **********
subroutine Vibrot(ngrid,nvib,Umin,Umax,R,PotR,E0,dE0,Redm,Req,sc,Temp)

use Vibrot_globals, only: Atom1, Atom2, EoutO, iad12, iadrsp, iadvib, iobs, ispc, J1A, J2A, lambda, npoint, Titobs, Vibwvs
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Four, Five, Six, Twelve
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ngrid, nvib
real(kind=wp), intent(in) :: Umin, Umax, R(ngrid+2), PotR(ngrid+2), dE0, Redm, Req, sc, Temp
real(kind=wp), intent(inout) :: E0
integer(kind=iwp) :: i, ii, ist, iter, ivib, J, ngrid1
real(kind=wp) :: Adif, AGi, Akk, AVi, const, D, D1, D1t, D2, dE, del, del12, del2, del56, delRN, diffm, DMem, E, E0x, E1, Emem, &
                 Est, Evib, fac, fact, GAM0, GAMn, Gmax, H11, H11a, H11c, Hnna, Hnnc, PotR0, PotRn1, R0, Redm2, Rn1, S11, S11a, &
                 S11c, Snna, Snnc, Svec0, SvecN1, term2, term4, V0, V1, Vmax, Vn, VN1, Wii, X2, Xrot
real(kind=wp), allocatable :: V(:), G(:), Svec(:), alpha(:), beta(:), B(:), W(:), Vec(:), Ener(:), S(:,:), H(:,:)
real(kind=wp), parameter :: Thre=1.0e-8_wp, THR=1.0e-3_wp
#include "warnings.fh"

D = Zero    ! dummy initialize
DMem = Zero ! dummy initialize

if (ngrid > npoint) then
  write(u6,999) ngrid,npoint
  call Abend()
end if

call mma_allocate(V,ngrid,label='V')
call mma_allocate(G,ngrid,label='G')
call mma_allocate(Svec,ngrid,label='Svec')
call mma_allocate(alpha,ngrid,label='alpha')
call mma_allocate(beta,ngrid,label='beta')
call mma_allocate(B,ngrid,label='B')
call mma_allocate(W,ngrid-1,label='W')
call mma_allocate(Vec,ngrid+1,label='Vec')
call mma_allocate(Ener,nvib,label='Ener')
call mma_allocate(S,3,ngrid,label='S')
call mma_allocate(H,3,ngrid,label='H')

del = (Umax-Umin)/(ngrid-1)
del2 = del**2
del12 = del2/Twelve
del56 = del2*Five/Six
PotR0 = PotR(ngrid+1)
PotRn1 = PotR(ngrid+2)
R0 = R(ngrid+1)
Rn1 = R(ngrid+2)

write(u6,2700)

! Start loop over rotational quantum numbers
Est = E0
do J=J1A,J2A
  iad12(J-J1A+1) = iadvib
  Xrot = -One/Two+sqrt(One/Four+J*(J+1)-lambda**2)

  ! Construct effective potential V and overlap vector Svec
  const = One/Four+Xrot*(Xrot+1)
  Redm2 = Two*Redm
  Svec0 = Redm2*R0**2
  V0 = const+Svec0*PotR0
  SvecN1 = Redm2*Rn1**2
  VN1 = const+SvecN1*PotRn1
  do i=1,ngrid
    fact = Redm2*R(i)**2
    Svec(i) = fact
    V(i) = const+fact*PotR(i)
  end do

  ! Construct Hamilton and Overlap matrices
  ngrid1 = ngrid-1
  do i=2,ngrid1
    H(1,i) = One-del12*V(i-1)
    S(1,i) = -del12*Svec(i-1)
    H(2,i) = -Two-del56*V(i)
    S(2,i) = -del56*Svec(i)
    H(3,i) = One-del12*V(i+1)
    S(3,i) = -del12*Svec(i+1)
  end do
  H(3,1) = One-del12*V(2)
  S(3,1) = -del12*Svec(2)
  H(1,ngrid) = One-del12*V(ngrid1)
  S(1,ngrid) = -del12*Svec(ngrid1)

  ! H(1,1) and S(1,1) depends explicitly on energy and are
  ! therefore constructed during the iteration process
  !
  ! set up and solve secular problem
  ! search for eigenvalues starting from E0.
  fac = exp(-del*(Xrot+One/Two))
  X2 = Redm/(Two*Xrot+Three)
  term2 = fac*X2*R0**2
  term4 = X2*R(1)**2
  H11a = -Two-del56*V(1)
  Hnna = -Two-del56*V(ngrid)
  H11c = One-del12*V0
  Hnnc = One-del12*VN1
  S11a = -del56*Svec(1)
  Snna = -del56*Svec(ngrid)
  S11c = -del12*Svec0
  Snnc = -del12*SvecN1
  delRN = Rn1-R(ngrid)

  do ivib=1,nvib
   iter = 0
   E = E0-dE0
   ist = 0
   Emem = Zero
   dE = dE0
   do
     ! IFG Stop if above the dissociation limit
     !     (assuming energy at infinity is set to zero)
     if (E > zero) then
       call SysQuitMsg(_RC_NOT_CONVERGED_,'VibRot','Failed to find a state.','Try decreasing STEP.')
     end if
     E = E+dE
     ! Compute h(1,1) and s(1,1)
     GAM0 = (fac+(PotR0-E)*term2)/(One+(PotR0-E)*term4)
     H11 = H11a+GAM0*H11c
     S11 = S11a+GAM0*S11c
     H(2,1) = H11
     S(2,1) = S11
     ! Compute H(N,N) and S(N,N)
     GAMn = exp(-sqrt(-Redm2*E)*delRN-del/Two)
     H(2,ngrid) = HnnA+GAMn*HnnC
     S(2,ngrid) = SnnA+GAMn*SnnC
     ! Calculate secular determinant
     D1t = H11-E*S11
     D1 = D1t/abs(D1t)
     D2 = One
     do i=2,ngrid
       Akk = H(2,i)-E*S(2,i)
       D = D1*Akk/abs(Akk)-D2*(H(3,i-1)-E*S(3,i-1))*(H(1,i)-E*S(1,i))/abs((H(2,i-1)-E*S(2,i-1))*Akk)
       D2 = D1
       D1 = D
     end do

     if (ist /= 0) then
       if (D/Dmem < Zero) then
         ! an eigenvalue has been bracketed
         E1 = E
         E = E1-DE*D/(D-Dmem)
         if (abs(E-Emem) < Thre) exit
         Emem = E
         E = E-DE
         DE = DE/Three
         E = E-DE
         ist = 0
       else
         Dmem = D
       end if
     else
       ist = 1
       Dmem = D
     end if
   end do
   ! One eigenvalue found. Remove scaling
   Ener(ivib) = (E-Est)/sc
   E0 = E+dE0
   alpha(ngrid) = H(2,ngrid)-E*S(2,ngrid)
   beta(ngrid) = H(1,ngrid)-E*S(1,ngrid)
   do i=1,ngrid1
     ii = ngrid-i
     if (ii /= 1) beta(ii) = H(1,ii)-E*S(1,ii)
     Wii = (H(3,ii)-E*S(3,ii))/alpha(ii+1)
     alpha(ii) = H(2,ii)-E*S(2,ii)-Wii*beta(ii+1)
     W(ii) = Wii
   end do
   ! step 2: unit right hand side
   if (abs(alpha(1)) == Zero) alpha(1) = 1.0e-10_wp
   G(1) = One/alpha(1)
   do i=2,ngrid
     G(i) = (One-beta(i)*G(i-1))/alpha(i)
   end do
   ! step 3: g on the right hand side
   do
     B(ngrid) = G(ngrid)
     do i=1,ngrid1
       ii = ngrid-i
       B(ii) = G(ii)-W(ii)*B(ii+1)
     end do
     Vec(1) = B(1)/alpha(1)
     do i=2,ngrid
       Vec(i) = (B(i)-beta(i)*Vec(i-1))/alpha(i)
     end do
     ! normalize g and vec
     Gmax = Zero
     Vmax = Zero
     do i=1,ngrid
       AGi = abs(G(i))
       if (AGi > Gmax) Gmax = AGi
       AVi = abs(Vec(i))
       if (AVi > Vmax) Vmax = AVi
     end do
     diffm = Zero
     do i=1,ngrid
       G(i) = G(i)/Gmax
       Vec(i) = Vec(i)/Vmax
       Adif = abs(Vec(i))-abs(G(i))
       if (Adif > diffm) diffm = Adif
     end do
     if ((diffm < 1.E-06).and.(iter /= 0)) exit
     iter = iter+1
     do i=1,ngrid
       G(i) = Vec(i)
     end do
   end do
   Evib = (E-Est)/sc
   !GG: Remove scaling
   Vec(ngrid+1) = E/sc
   call DDafile(Vibwvs,1,Vec(1),ngrid+1,iadvib)

   ! Print energies and notify on end point values
   V1 = abs(Vec(1))
   Vn = abs(Vec(ngrid))
   if ((V1 < THR).and.(Vn < THR)) write(u6,3000) ivib-1,J,Evib
   if ((V1 < THR).and.(Vn >= THR)) write(u6,3001) ivib-1,J,Evib,Vn
   if ((V1 >= THR).and.(Vn < THR)) write(u6,3002) ivib-1,J,Evib,V1
   if ((V1 >= THR).and.(Vn >= THR)) write(u6,3003) ivib-1,J,Evib,V1,Vn

  end do

  iadrsp(J-J1A+1) = iadvib
  call DDafile(Vibwvs,1,Ener(1),nvib,iadvib)
  E0 = Est
end do

! write address record to Vibwvs
iadvib = 0
call iDafile(Vibwvs,1,iad12,100,iadvib)

! Compute spectroscopic parameters
E0x = Est/sc
if (ispc /= 0) call Spectc(Req,E0x,Atom1,Atom2,nvib)

! Compute matrix elements of observables
if (iobs /= 0) then
  do i=1,iobs
    call Vibmat(nvib,ngrid,Umin,Umax,Titobs(i),R,EoutO(1,i),Temp)
  end do
end if

call mma_deallocate(V)
call mma_deallocate(G)
call mma_deallocate(Svec)
call mma_deallocate(alpha)
call mma_deallocate(beta)
call mma_deallocate(B)
call mma_deallocate(W)
call mma_deallocate(Vec)
call mma_deallocate(Ener)
call mma_deallocate(S)
call mma_deallocate(H)

return

999 format(/1x,'Number of gridpoints',i4,' is larger than' &
           /1x,' dimension npoint=',i3,' in VIBROT.' &
           /1x,' Program cannot continue.')
2700 format(//1x,'Eigenstates'//1x,' Vib. q.n.   Rot. q.n.    Energy')
3000 format(4x,i3,9x,i3,2x,f12.6)
3001 format(4x,i3,9x,i3,2x,f12.6,5x,'Warning:function value at Rmax is',f12.6)
3002 format(4x,i3,9x,i3,2x,f12.6,5x,'Warning:function value at Rmin is',f12.6)
3003 format(4x,i3,9x,i3,2x,f12.6,5x,'Warning:function value at Rmin is',f12.6,' and at Rmax',f12.6)

end subroutine Vibrot
