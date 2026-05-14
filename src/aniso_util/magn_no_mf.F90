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
!  X,Y,Z -- projections of the magnetic field, specifying the orientation of the applied
!           magnetic field, Real(kind=wp) ::, input;  rule: ( X**2 + Y**2 + Z**2 = 1);
!      H -- strength of the magnetic field in tesla, Real(kind=wp) ::, input;
!      W -- energies of the exchange states; Real(kind=wp) :: array (EXCH);
!     dM -- matrix of the magnetic moment, Complex(kind=wp) :: (3,EXCH,EXCH) array, input;
!     sM -- matrix of the     spin moment, Complex(kind=wp) :: (3,EXCH,EXCH) array, input;
!     nT -- number of temperature points for which magnetisation is computed, input;
!      T -- temperature values(in kelvin) for which magnetisation is computed, input;
!   sopt -- logical parameter. If sopt=.true. Then spin magnetisation is computed.
!                              If sopt=.false. Then spin part is skipped.
!
!     WZ -- Zeeman energies, true values (not shifted to 0), in cm-1, Real(kind=wp) :: (N) array, output;
!     ZB -- statistical Boltzmann distribution, for each temperature, Real(kind=wp) :: (nT) array, output;
!      S -- spin magnetisation, Real(kind=wp) :: (3,nT) array, output;
!      M -- magnetisation, Real(kind=wp) :: (3,nT) array, output;
! m_paranoid --  logical parameter.
!            If m_paranoid = .true.  Then the average spin is computed for each temperature point exactly
!            If m_paranoid = .false. Then  the average spin is computed only for the lowest temperature point
!---------
!  temporary (local) variables:
!    MZ -- matrix of the magnetic moment, Complex(kind=wp) :: (3,EXCH,EXCH) array
!    SZ -- matrix of the     spin moment, Complex(kind=wp) :: (3,EXCH,EXCH) array
!    WM -- array containing the Zeeman eigenstates and, If N<EXCH the exchange eigenstates
!          for the states higher in energy than N, Real(kind=wp) ::, (EXCH) array;
!    ZM -- Zeeman eigenvectors, (N,N) Complex(kind=wp) :: array,
!   i,j -- labels of the states;
!     l -- labels the cartesian component of the momentum (convention: x=1, y=2, z=3)
!    iT -- labels the temperature points;

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: EXCH, N, nT
real(kind=wp), intent(in) :: X, Y, Z, H, W(EXCH), T(nT)
complex(kind=wp), intent(in) :: dM(3,EXCH,EXCH), sM(3,EXCH,EXCH)
logical(kind=iwp), intent(in) :: sopt
real(kind=wp), intent(out) :: WZ(N), ZB(nT), S(3,nT), M(3,nT)
logical(kind=iwp), intent(in) :: DBG
integer(kind=iwp) :: iT, l
real(kind=wp) :: ST(3), zJ
real(kind=wp), allocatable :: RWORK(:), WM(:)
complex(kind=wp), allocatable :: dM_TMP(:,:,:), HZEE(:), MZ(:,:,:), sM_TMP(:,:,:), SZ(:,:,:), TMP(:,:), W_c(:), WORK(:), ZM(:,:)

WZ(:) = Zero
ZB(:) = Zero
S(:,:) = Zero
M(:,:) = Zero

! a few checks, before proceeding:
do iT=1,nT
  if (T(iT) == Zero) return
end do
if (H == Zero) return
if (N > EXCH) return
zJ = Zero ! it must be defined for ZEEM
! initialization:
call mma_allocate(WM,exch,'WM')
call mma_allocate(ZM,exch,exch,'ZM')
call mma_allocate(SZ,3,exch,exch,'SZ')
call mma_allocate(MZ,3,exch,exch,'MZ')

! temporary arrays used in ZEEM_SA:
call mma_allocate(RWORK,3*N-2,'ZEEM_RWORK')
call mma_allocate(HZEE,nTri_Elem(N),'ZEEM_HZEE')
call mma_allocate(WORK,2*N-1,'ZEEM_WORK')
call mma_allocate(W_c,N,'ZEEM_W_c')

! zero everything:
RWORK(:) = Zero
HZEE(:) = cZero
WORK(:) = cZero
W_c(:) = cZero

ST(:) = Zero
WM(:) = Zero
ZM(:,:) = cZero
SZ(:,:,:) = cZero
MZ(:,:,:) = cZero
! start calculations:
if (DBG) then
  write(u6,*) 'Enter ZEEM::'
  write(u6,*) 'Input data:   N = ',N
  write(u6,*) 'Input data:   H = ',H
  write(u6,*) 'Input data:   X = ',X
  write(u6,*) 'Input data:   Y = ',Y
  write(u6,*) 'Input data:   Z = ',Z
  write(u6,*) 'Input data:  zJ = ',zJ
  write(u6,*) 'Input data: W() = ',W(1:N)
  write(u6,*) 'Input data: ST()= ',ST(1:3)
  call prmom('Input data dM:',dM,N)
  call prmom('Input data sM:',sM,N)
end if

! Build and diagonalize the Zeeman Hamiltonian
! most important output are: WM (energies) and ZM (eigenvectors)
if (N == EXCH) then
  call ZEEM_SA(N,H,X,Y,Z,W,dM,sM,ST,zJ,WM,ZM,DBG,RWORK,HZEE,WORK,W_c)
else
  call mma_allocate(dM_TMP,3,N,N,label='dM_TMP')
  call mma_allocate(sM_TMP,3,N,N,label='sM_TMP')
  dM_TMP(:,:,:) = dM(:,1:N,1:N)
  sM_TMP(:,:,:) = sM(:,1:N,1:N)
  call ZEEM_SA(N,H,X,Y,Z,W(1:N),dM_TMP,sM_TMP,ST,zJ,WM(1:N),ZM,DBG,RWORK,HZEE,WORK,W_c)
  call mma_deallocate(dM_TMP)
  call mma_deallocate(sM_TMP)
end if
if (DBG) write(u6,*) 'Exit ZEEM::'

WZ(:) = WM(1:N)
WM(N+1:) = W(N+1:)

! transform the momenta
call UTMU(EXCH,N,ZM,sM,SZ)
call UTMU(EXCH,N,ZM,dM,MZ)

! compute magnetization at different temperatures:

if (N == EXCH) then
  call mma_allocate(TMP,EXCH,EXCH,label='TMP')
  do iT=1,nT
    do l=1,3
      if (sopt) then
        TMP(:,:) = SZ(l,:,:)
        call calcmagn1(EXCH,WM,TMP,T(iT),S(l,iT),ZB(iT))
      end if
      TMP(:,:) = MZ(l,:,:)
      call calcmagn1(EXCH,WM,TMP,T(iT),M(l,iT),ZB(iT))
    end do
  end do
  call mma_deallocate(TMP)
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
call mma_deallocate(SZ)
call mma_deallocate(MZ)

return

end subroutine MAGN_NO_MF
