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

subroutine MAGN_NO_MF(EXCH,N,X,Y,Z,H,W,dM,sM,nT,T,sopt,WZ,ZB,S,M,DBG)
! this Subroutine computes the magnetisation for zJ==0.0
!
! definition of the variables:
!   EXCH -- total number of exchange states, Integer, input
!      N -- size of the Zeeman matrix, Integer, input, NM <= EXCH !
!  X,Y,Z -- projections of the magnetic field, specIfying the orientation of the applied
!           magnetic field, Real(kind=8) ::, input;  rule: ( X**2 + Y**2 + Z**2 = 1);
!      H -- strength of the magnetic field in Tesla, Real(kind=8) ::, input;
!      W -- energies of the exchange states; Real(kind=8) :: array (EXCH);
!     dM -- matrix of the magnetic moment, Complex(kind=8) :: (3,EXCH,EXCH) array, input;
!     sM -- matrix of the     spin moment, Complex(kind=8) :: (3,EXCH,EXCH) array, input;
!     nT -- number of temperature points for which magnetisation is computed, input;
!      T -- temperature values(in Kelvin) for which magnetisation is computed, input;
!   sopt -- logical parameter. If sopt=.true. Then spin magnetisation is computed.
!                              If sopt=.false. Then spin part is skipped.
!
!     WZ -- Zeeman energies, true values (not shifted to 0), in cm-1, Real(kind=8) :: (N) array, output;
!     ZB -- statistical Boltzmann distribution, for each temperature, Real(kind=8) :: (nT) array, output;
!      S -- spin magnetisation, Real(kind=8) :: (3,nT) array, output;
!      M -- magnetisation, Real(kind=8) :: (3,nT) array, output;
! m_paranoid --  logical parameter.
!            If m_paranoid = .true.  Then the average spin is computed for each temperature point exactly
!            If m_paranoid = .false. Then  the average spin is computed only for the lowest temperature point
!---------
!  temporary (local) variables:
!    MZ -- matrix of the magnetic moment, Complex(kind=8) :: (3,EXCH,EXCH) array
!    SZ -- matrix of the     spin moment, Complex(kind=8) :: (3,EXCH,EXCH) array
!    WM -- array containing the Zeeman eigenstates and, If N<EXCH the exchange eigenstates
!          for the states higher in energy than N, Real(kind=8) ::, (EXCH) array;
!    ZM -- Zeeman eigenvectors, (N,N) Complex(kind=8) :: array,
!   i,j -- labers of the states;
!     l -- labels the cartesian component of the momentum (convention: x=1, y=2, z=3)
!    iT -- labes the temperature points;

implicit none
#include "stdalloc.fh"
integer, parameter :: wp = kind(0.d0)
integer, intent(in) :: EXCH, N, nT
real(kind=8), intent(in) :: X, Y, Z, H
real(kind=8), intent(in) :: W(EXCH), T(nT)
real(kind=8), intent(out) :: ZB(nT), WZ(N)
real(kind=8), intent(out) :: S(3,nT), M(3,nT)
complex(kind=8), intent(in) :: dM(3,EXCH,EXCH)
complex(kind=8), intent(in) :: sM(3,EXCH,EXCH)
logical, intent(in) :: sopt
! local variables:
integer :: i, l, iT
real(kind=8) :: zJ
real(kind=8), allocatable :: WM(:), ST(:), RWORK(:) ! WM(EXCH), ST(3)
complex(kind=8), allocatable :: HZEE(:), WORK(:), W_c(:)
complex(kind=8), allocatable :: ZM(:,:) ! ZM(N,N)
complex(kind=8), allocatable :: SZ(:,:,:), MZ(:,:,:) ! SZ(3,EXCH,EXCH), MZ(3,EXCH,EXCH)
logical :: DBG

! a few checks, before proceeding:
do iT=1,nT
  if (T(iT) == 0.0_wp) return
end do
if (H == 0.0_wp) return
if (N > EXCH) return
zJ = 0.0_wp ! it must be defined for ZEEM
! initialization:
call mma_allocate(WM,exch,'WM')
call mma_allocate(ZM,exch,exch,'ZM')
call mma_allocate(ST,3,'ST')
call mma_allocate(SZ,3,exch,exch,'SZ')
call mma_allocate(MZ,3,exch,exch,'MZ')

! temporary arrays used in ZEEM_SA:
call mma_allocate(RWORK,(3*N-2),'ZEEM_RWORK')
call mma_allocate(HZEE,(N*(N+1)/2),'ZEEM_HZEE')
call mma_allocate(WORK,(2*N-1),'ZEEM_WORK')
call mma_allocate(W_c,N,'ZEEM_W_c')

! zero everything:
call dcopy_(3*N-2,[0.0_wp],0,RWORK,1)
call zcopy_(N*(N+1)/2,[(0.0_wp,0.0_wp)],0,HZEE,1)
call zcopy_(2*N-1,[(0.0_wp,0.0_wp)],0,WORK,1)
call zcopy_(N,[(0.0_wp,0.0_wp)],0,W_c,1)

call dcopy_(3,[0.0_wp],0,ST,1)
call dcopy_(N,[0.0_wp],0,WZ,1)
call dcopy_(exch,[0.0_wp],0,WM,1)
call zcopy_(exch*exch,[(0.0_wp,0.0_wp)],0,ZM,1)
call zcopy_(3*exch*exch,[(0.0_wp,0.0_wp)],0,SZ,1)
call zcopy_(3*exch*exch,[(0.0_wp,0.0_wp)],0,MZ,1)
! start calculations:
if (DBG) then
  write(6,*) 'Enter ZEEM::'
  write(6,*) 'Input data:   N = ',N
  write(6,*) 'Input data:   H = ',H
  write(6,*) 'Input data:   X = ',X
  write(6,*) 'Input data:   Y = ',Y
  write(6,*) 'Input data:   Z = ',Z
  write(6,*) 'Input data:  zJ = ',zJ
  write(6,*) 'Input data: W() = ',W(1:N)
  write(6,*) 'Input data: ST()= ',ST(1:3)
  call prmom('INput data dM:',dM,N)
  call prmom('INput data sM:',sM,N)
end if

! Build and diagonalize the Zeeman Hamiltonian
! most important output are: WM (energies) and ZM (eigenvectors)
call ZEEM_SA(N,H,X,Y,Z,W(1:N),dM(1:3,1:N,1:N),sM(1:3,1:N,1:N),ST,zJ,WM(1:N),ZM,DBG,RWORK,HZEE,WORK,W_c)
if (DBG) write(6,*) 'Exit ZEEM::'

call DCOPY_(N,WM(1:N),1,WZ(1:N),1)
if (N /= EXCH) then
  do i=N+1,EXCH
    WM(i) = W(i)
  end do
end if

! transform the momenta
call UTMU(EXCH,N,ZM,sM,SZ)
call UTMU(EXCH,N,ZM,dM,MZ)

! compute magnetization at different temperatures:
call dcopy_(nT,[0.0_wp],0,ZB,1)
call dcopy_(3*nT,[0.0_wp],0,M,1)
call dcopy_(3*nT,[0.0_wp],0,S,1)

if (N == EXCH) then
  do iT=1,nT
    do l=1,3
      if (sopt) call calcmagn1(EXCH,WM,SZ(l,:,:),T(iT),S(l,iT),ZB(iT))
      call calcmagn1(EXCH,WM,MZ(l,:,:),T(iT),M(l,iT),ZB(iT))
    end do
  end do
else
  do iT=1,nT
    do l=1,3
      if (sopt) call calcmagn2(EXCH,N,WM,T(iT),H,SZ,X,Y,Z,l,S(l,iT),ZB(iT))
      call calcmagn2(EXCH,N,WM,T(iT),H,MZ,X,Y,Z,l,M(l,iT),ZB(iT))
    end do
  end do !iT
end if

! deallocate temporary arrays:
call mma_deallocate(RWORK)
call mma_deallocate(HZEE)
call mma_deallocate(WORK)
call mma_deallocate(W_c)

call mma_deallocate(WM)
call mma_deallocate(ZM)
call mma_deallocate(ST)
call mma_deallocate(SZ)
call mma_deallocate(MZ)

return

end subroutine MAGN_NO_MF
