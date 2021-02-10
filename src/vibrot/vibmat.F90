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

! Calculation of matrix elements over vibrational wave functions
! of a given series of observables by a simpson quadrature.
! The vibrational wave functions have been computed in vibrot
! and stored on unit Vibwvs=12 for each rotational quantum number.
! The observables are read as input for a sequence of R-values
! and fitted to an analytical form in function pot.
! This routine is called when the parameter iobs is not zero.
! integration is performed on a logaritmic scale (u=ln(r))
! between Umin and Umax with a grid size del corresponding to
! ndim integration steps. ndim has to be odd.
!
! ********** MOLCAS-2 Release 91 05 01 **********

subroutine Vibmat(ne,ndim,Umin,Umax,Title,R,PotR,Temp)

use Vibrot_globals, only: J1A, J2A, Vibwvs, iad12
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, kBoltzmann, auTokJ
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ne, ndim
real(kind=wp), intent(in) :: Umin, Umax, R(ndim), PotR(ndim), Temp
character(len=80), intent(in) :: Title
integer(kind=iwp) :: iPrint, i, J, iad, ind, ist, ist1, ist2, ndim1, nv, nv1, nv2, nv11, nv12, nv22, nwork
logical(kind=iwp) :: WarnV, WarnR
real(kind=wp) :: Joule, beta, del, E, O, Obsr, ri, Smax, W, WarnRq, WarnVq, Wfirst, Wlast
integer(kind=iwp), allocatable :: nv1w(:), nv2w(:)
real(kind=wp), allocatable :: Vib(:), X(:), SumO(:), SumW(:), Wmin(:), Wmax(:), Ener(:), S(:), Obs(:), Sw(:)

write(u6,*)
call CollapseOutput(1,'Matrix elements of observable: '//Title)
WarnV = .false.
WarnR = .false.
WarnVq = Zero
WarnRq = Zero
Joule = 1.0e3_wp*auTokJ
beta = One/(kBoltzmann*Temp)
!write(u6,'(a,e12.5)') 'Bolzmann ',kBoltzmann
!write(u6,'(a,e12.5)') 'au to J  ',Joule
!write(u6,'(a,e12.5)') 'Temp     ',Temp
!write(u6,'(a,e12.5)') 'beta     ',beta
ndim1 = ndim+1
nwork = ne*ndim1
call mma_allocate(Vib,nwork,label='Vib')
call mma_allocate(X,ndim,label='X')
call mma_allocate(SumO,ne,label='SumO')
call mma_allocate(SumW,ne,label='SumW')
call mma_allocate(Wmin,ne,label='Wmin')
call mma_allocate(Wmax,ne,label='Wmax')
call mma_allocate(Ener,ne,label='Ener')
call mma_allocate(nv1w,ne*(ne+1)/2,label='nv1w')
call mma_allocate(nv2w,ne*(ne+1)/2,label='nv2w')
call mma_allocate(S,ne*(ne+1)/2,label='S')
call mma_allocate(Sw,ne*(ne+1)/2,label='Sw')
call mma_allocate(Obs,ne*(ne+1)/2,label='Obs')

SumO(:) = Zero
SumW(:) = Zero
Wmin(:) = Zero
Wmax(:) = Zero

del = (Umax-Umin)/(ndim-1)
iPrint = 1

! loop over rotational quantum numbers
do J=J1A,J2A

  ! Read vibrational functions for this J-value. store in vib.
  ist = 1
  iad = iad12(J-J1A+1)
  do nv=1,ne
    call DDafile(Vibwvs,2,Vib(ist),ndim1,iad)
    !write(u6,6668) nv,(Vib(i+ist-1),i=1,ndim+1)
    ist = ist+ndim1
  end do
  call DDafile(Vibwvs,2,Ener,ne,iad)

  ! Compute overlap matrix S
  ist1 = -ndim1
  nv12 = 0
  do nv1=1,ne
    ist1 = ist1+ndim1
    ist2 = -ndim1
    do nv2=1,nv1
      ist2 = ist2+ndim1
      nv12 = nv12+1

      ! Set up scalar product and integrate
      do i=1,ndim
        X(i) = Vib(i+ist1)*Vib(i+ist2)*R(i)**2
      end do
      call Simpsn(X,del,ndim,S(nv12))
    end do
  end do

  ! Check overlap matrix for non-orthogonality
  nv12 = 0
  Smax = Zero
  do nv1=1,ne
    do nv2=1,nv1
      nv12 = nv12+1
      if (nv1 /= nv2) then
        if (abs(S(nv12)) > abs(Smax)) Smax = S(nv12)
      end if
    end do
  end do
  if (abs(Smax) > 1.0e-4_wp) write(u6,1200) Smax

  ! Print overlap matrix
  if (iPrint >= 2) then
    write(u6,'(1x,A)') 'Overlap matrix for vibrational wave functions'
    nv12 = 0
    do nv1=1,ne
      do nv2=1,nv1
        nv12 = nv12+1
        nv1w(nv12) = nv1
        nv2w(nv12) = nv2
        Sw(nv12) = S(nv12)
      end do
    end do
    write(u6,'(6(3X,2I3,F12.6))') (nv1w(i),nv2w(i),Sw(i),i=1,nv12)
  end if

  ! set up integrand array for this Observable
  ist1 = -ndim
  nv12 = 0
  nv11 = 0
  do nv1=1,ne
    nv11 = nv11+nv1
    nv22 = 0
    ist1 = ist1+ndim
    ist2 = -ndim
    do nv2=1,nv1
      nv22 = nv22+nv2
      ist2 = ist2+ndim
      nv12 = nv12+1
      do i=1,ndim
        ri = R(i)
        X(i) = Vib(i+ist1)*Vib(i+ist2)*potR(i)*ri**2
      end do
      call Simpsn(x,del,ndim,Obsr)
      Obs(nv12) = Obsr/sqrt(S(nv11)*S(nv22))
    end do
  end do

  ! write matrix elements
  write(u6,'(">>>>",1x,A,I3)') 'matrix elements over vibrational wave functions (atomic units) for rotational quantum number',J
  nv12 = 0
  do nv1=1,ne
    do nv2=1,nv1
      nv12 = nv12+1
      nv1w(nv12) = nv1
      nv2w(nv12) = nv2
      Sw(nv12) = Obs(nv12)
    end do
  end do
  write(u6,'(6(3X,2I3,F12.6))') (nv1w(i),nv2w(i),Sw(i),i=1,nv12)
  Wfirst = One
  Wlast = Zero
  do nv=1,ne
    ind = nv*(nv+1)/2
    O = Obs(ind)
    E = Ener(nv)*Joule
    W = exp(-beta*E)
    SumW(nv) = SumW(nv)+W
    SumO(nv) = SumO(nv)+W*O
    if (J == J1a) then
      Wmax(nv) = W
      Wmin(nv) = W
    else
      Wmax(nv) = max(W,Wmax(nv))
      Wmin(nv) = min(W,Wmin(nv))
    end if
    if (nv == 1) Wfirst = W
    if (nv == ne) WLast = W
    !write(u6,'(a,3g15.6)') '... E,W,O ',E,W,O
  end do
  !write(u6,'(a,f12.6)') 'Wlast/Wfirst',Wlast/Wfirst
  if (Wlast/Wfirst > 1.0e-3_wp) then
    WarnVq = max(WarnVq,Wlast/Wfirst)
    WarnV = .true.
  end if

! End of loop over rotational quantum number
end do

do nv=1,ne
  !write(u6,'(a,i3)') 'Temp terms for vibration qn ',nv
  !write(u6,'(a,es15.6)') 'SumW(nv) ',SumW(nv)
  !write(u6,'(a,es15.6)') 'SumO(nv) ',SumO(nv)
  !write(u6,'(a,es15.6)') 'Wmax(nv) ',Wmax(nv)
  !write(u6,'(a,es15.6)') 'Wmin(nv) ',Wmin(nv)
  !write(u6,*)
  if (Wmin(nv)/Wmax(nv) > 1.0e-3_wp) then
    WarnRq = max(WarnRq,Wmin(nv)/Wmax(nv))
    WarnR = .true.
  end if
end do
O = Zero
W = Zero
do nv=1,ne
  O = O+SumO(nv)
  W = W+SumW(nv)
end do
O = O/W
write(u6,*)
write(u6,'(a,es15.6,a,f9.3,a)') 'Temperature averaged observable:',O,'  at',Temp,'K'
if (WarnV) then
  write(u6,*)
  write(u6,'(a)') '***'
  write(u6,'(a)') '*** Warning, temperature weighting not converged with respect to vibrational quantum numbers'
  write(u6,'(a)') '***'
  write(u6,'(a,es10.3,a)') '*** Quotient',WarnVq,' should be small'
  write(u6,'(a)') '***'
  write(u6,*)
end if
if (WarnR) then
  write(u6,*)
  write(u6,'(a)') '***'
  write(u6,'(a)') '*** Warning, temperature weighting not converged with respect to rotational quantum numbers'
  write(u6,'(a)') '***'
  write(u6,'(a,es10.3,a)') '*** Quotient',WarnRq,' should be small'
  write(u6,'(a)') '***'
  write(u6,*)
end if
call CollapseOutput(0,'Matrix elements of observable: '//Title)

call mma_deallocate(Vib)
call mma_deallocate(X)
call mma_deallocate(SumO)
call mma_deallocate(SumW)
call mma_deallocate(Wmin)
call mma_deallocate(Wmax)
call mma_deallocate(Ener)
call mma_deallocate(nv1w)
call mma_deallocate(nv2w)
call mma_deallocate(S)
call mma_deallocate(Sw)
call mma_deallocate(OBs)

! End of calculation for this Observable.

return

!6668 format(/1x,'Vib-values',i3/(1x,6f12.6))
1200 format(/1x,'***** Warning: Non-orthogonality between vibrational',1x,'wave functions.' &
            /13x,'Largest overlap matrix element is',f14.6)

end subroutine Vibmat
