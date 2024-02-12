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

subroutine MAGN_ZJ_PAR(EXCH,N,X,Y,Z,H,W,zJ,dM,sM,nT,T,sopt,WZ,ZB,S,M,thrs,m_paranoid,dbg)
! this Subroutine computes the magnetisation for zJ /= 0.0, paranoid accuracy
!
! definition of the variables:
!   EXCH -- total number of exchange states, Integer, input
!      N -- size of the Zeeman matrix, Integer, input, NM <= EXCH !
!  X,Y,Z -- projections of the magnetic field, specIfying the orientation of the applied
!           magnetic field, Real(kind=wp) ::, input;  rule: ( X**2 + Y**2 + Z**2 = 1);
!      H -- strength of the magnetic field in tesla, Real(kind=wp) ::, input;
!      W -- energies of the exchange states; Real(kind=wp) :: array (EXCH);
!     zJ -- parameter of intermolecular interaction, Real(kind=wp) ::, input;
!   THRS -- threshold for convergence of the spin magnetisation. Real(kind=wp) ::, input;
!           Of any importance only when (zJ /= 0.0_wp), otherwise unused. Default value 1.D-8.
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
!  SCHK -- variable used for checking the convergence; SCHK = ABS(Sx) + ABS(Sy) + ABS(Sz), Real(kind=wp) ::
!   i,j -- labers of the states;
!     l -- labels the cartesian component of the momentum (convention: x=1, y=2, z=3)
!    iT -- labes the temperature points;
!    ST -- value of the average spin of neighboring sites, Real(kind=wp) :: (3) array;

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: EXCH, N, nT
real(kind=wp), intent(in) :: X, Y, Z, H, zJ, W(EXCH), T(nT), thrs
complex(kind=wp), intent(in) :: dM(3,EXCH,EXCH), sM(3,EXCH,EXCH)
logical(kind=iwp), intent(in) :: sopt, m_paranoid, dbg
real(kind=wp), intent(out) :: WZ(N), ZB(nT), S(3,nT), M(3,nT)
integer(kind=iwp) :: i, iT, l
real(kind=wp) :: ST(3), STsave(3)
real(kind=wp), allocatable :: RWORK(:), WM(:)
complex(kind=wp), allocatable :: HZEE(:), MZ(:,:,:), SZ(:,:,:), W_c(:), WORK(:), ZM(:,:)

! a few checks, before proceeding:
do iT=1,nT
  if (T(iT) == Zero) return
end do
if (H == Zero) return
if (N > EXCH) return

! allocate memory:
call mma_allocate(WM,exch,'WM')
call mma_allocate(ZM,exch,exch,'ZM')
call mma_allocate(SZ,3,exch,exch,'SZ')
call mma_allocate(MZ,3,exch,exch,'MZ')

! temporary arrays used in ZEEM_SA:
call mma_allocate(RWORK,3*N-2,'ZEEM_RWORK')
call mma_allocate(HZEE,N*(N+1)/2,'ZEEM_HZEE')
call mma_allocate(WORK,2*N-1,'ZEEM_WORK')
call mma_allocate(W_c,N,'ZEEM_W_c')

! zero everything:
WM(:) = Zero
ZM(:,:) = cZero
SZ(:,:,:) = cZero
MZ(:,:,:) = cZero

RWORK(:) = Zero
HZEE(:) = cZero
WORK(:) = cZero
W_c(:) = cZero

! start calculations:
! code for the case (zJ /= 0):
ZB(:) = Zero
M(:,:) = Zero
S(:,:) = Zero

do iT=1,nT
  ! determine first the average spin of neighboring
  ! molecules for each temperature point ST(1:3)
  if (m_paranoid) then
    ST(:) = Zero
    call mean_field(EXCH,N,H,X,Y,Z,zJ,T(iT),W,thrs,DM,SM,ST,dbg)
    STsave(:) = ST(:)
  else
    ! i.e. when m_paranoid=.false.
    if (iT == 1) then
      ST(:) = Zero
      call mean_field(EXCH,N,H,X,Y,Z,zJ,T(iT),W,thrs,DM,SM,ST,dbg)
      STsave(:) = ST(:)
    else
      ! use the last saved value of the ST:
      ST(:) = STsave(:)
    end if
  end if
  !---------------------------------------------------------------------
  ! here  we have the value of the averaged spin for this temperature
  ! proceed with the computation of magnetism for this temperature
  if (DBG) write(u6,'(A,3ES13.5)') 'Average spin finished. ST on entrance to last ZEEM:',(ST(i),i=1,3)
  WM(:) = Zero
  ZM(:,:) = cZero

  call ZEEM_SA(N,H,X,Y,Z,W(1:N),dM(1:3,1:N,1:N),sM(1:3,1:N,1:N),ST,zJ,WM(1:N),ZM(1:N,1:N),DBG,RWORK,HZEE,WORK,W_c)

  ! move WM energies to WZ:
  WZ(:) = WM(:)
  ! /// calculation of matrix elements of spin momentum in the basis of Zeeman states
  WM(N+1:EXCH) = W(N+1:EXCH)

  ! transform the momenta
  SZ(:,:,:) = cZero
  MZ(:,:,:) = cZero
  call UTMU(EXCH,N,ZM(1:N,1:N),SM,SZ)
  call UTMU(EXCH,N,ZM(1:N,1:N),DM,MZ)

  ! calculation of magnetizations at different temperatures:
  if (N == EXCH) then
    do l=1,3
      if (sopt) call calcmagn1(EXCH,WM,SZ(l,1:EXCH,1:EXCH),T(iT),S(l,iT),ZB(iT))
      call calcmagn1(EXCH,WM,MZ(l,1:EXCH,1:EXCH),T(iT),M(l,iT),ZB(iT))
    end do
  else
    do l=1,3
      if (sopt) call calcmagn2(EXCH,N,WM,T(iT),H,SZ,X,Y,Z,l,S(l,iT),ZB(iT))
      call calcmagn2(EXCH,N,WM,T(iT),H,MZ,X,Y,Z,L,M(l,iT),ZB(iT))
    end do
  end if

end do ! iT

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

end subroutine MAGN_ZJ_PAR
