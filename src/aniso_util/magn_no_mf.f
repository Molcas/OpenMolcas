************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine MAGN_NO_MF( EXCH,N,X,Y,Z,H,W,dM,sM,nT,T,sopt,WZ,ZB,
     &                       S,M, DBG )
c this Subroutine computes the magnetisation for zJ==0.0
c
c definition of the variables:
c   EXCH -- total number of exchange states, Integer, input
c      N -- size of the Zeeman matrix, Integer, input, NM .LE. EXCH !
c  X,Y,Z -- projections of the magnetic field, specIfying the orientation of the applied
c           magnetic field, Real(kind=wp) ::, input;  rule: ( X**2 + Y**2 + Z**2 = 1);
c      H -- strength of the magnetic field in Tesla, Real(kind=wp) ::, input;
c      W -- energies of the exchange states; Real(kind=wp) :: array (EXCH);
c     dM -- matrix of the magnetic moment, Complex(kind=wp) :: (3,EXCH,EXCH) array, input;
c     sM -- matrix of the     spin moment, Complex(kind=wp) :: (3,EXCH,EXCH) array, input;
c     nT -- number of temperature points for which magnetisation is computed, input;
c      T -- temperature values(in Kelvin) for which magnetisation is computed, input;
c   sopt -- logical parameter. If sopt=.true. Then spin magnetisation is computed.
c                              If sopt=.false. Then spin part is skipped.
c
c     WZ -- Zeeman energies, true values (not shifted to 0), in cm-1, Real(kind=wp) :: (N) array, output;
c     ZB -- statistical Boltzmann distribution, for each temperature, Real(kind=wp) :: (nT) array, output;
c      S -- spin magnetisation, Real(kind=wp) :: (3,nT) array, output;
c      M -- magnetisation, Real(kind=wp) :: (3,nT) array, output;
c m_paranoid --  logical parameter.
c            If m_paranoid = .true.  Then the average spin is computed for each temperature point exactly
c            If m_paranoid = .false. Then  the average spin is computed only for the lowest temperature point
c---------
c  temporary (local) variables:
c    MZ -- matrix of the magnetic moment, Complex(kind=wp) :: (3,EXCH,EXCH) array
c    SZ -- matrix of the     spin moment, Complex(kind=wp) :: (3,EXCH,EXCH) array
c    WM -- array containing the Zeeman eigenstates and, If N<EXCH the exchange eigenstates
c          for the states higher in energy than N, Real(kind=wp) ::, (EXCH) array;
c    ZM -- Zeeman eigenvectors, (N,N) Complex(kind=wp) :: array,
c   i,j -- labers of the states;
c     l -- labels the cartesian component of the momentum (convention: x=1, y=2, z=3)
c    iT -- labes the temperature points;
      Implicit None
#include "stdalloc.fh"
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)          :: EXCH, N, nT
      Real(kind=wp), intent(in)    :: X, Y, Z, H
      Real(kind=wp), intent(in)    :: W(EXCH), T(nT)
      Real(kind=wp), intent(out)   :: ZB(nT), WZ(N)
      Real(kind=wp), intent(out)   :: S(3,nT), M(3,nT)
      Complex(kind=wp), intent(in) :: dM(3,EXCH,EXCH)
      Complex(kind=wp), intent(in) :: sM(3,EXCH,EXCH)
      Logical, intent(in)          :: sopt
c local variables:
      Integer          :: i, l, iT
      Real(kind=wp)    :: kB, mB, zJ
      Real(kind=wp), allocatable :: WM(:), ST(:), RWORK(:)
!                                   WM(EXCH), ST(3)
      Complex(kind=wp), allocatable :: HZEE(:), WORK(:), W_c(:)
      Complex(kind=wp), allocatable :: ZM(:,:) !ZM(N,N)
      Complex(kind=wp), allocatable :: SZ(:,:,:), MZ(:,:,:)
!                                      SZ(3,EXCH,EXCH), MZ(3,EXCH,EXCH)
      Logical          :: DBG
      Call qEnter('MAGN_NO_MF')
      kB=0.6950356000_wp   ! Boltzmann constant,  in cm^-1*K-1
      mB=0.4668643740_wp   ! Bohr magneton,       in cm-1*T-1

c a few checks, before proceeding:
      Do iT=1,nT
        If (T(iT)==0.0_wp) Return
      End Do
      If(H==0.0_wp) Return
      If(N > EXCH) Return
      zJ = 0.0_wp ! it must be defined for ZEEM
c initialization:
      Call mma_allocate(WM,exch,'WM')
      Call mma_allocate(ZM,exch,exch,'ZM')
      Call mma_allocate(ST,3,'ST')
      Call mma_allocate(SZ,3,exch,exch,'SZ')
      Call mma_allocate(MZ,3,exch,exch,'MZ')

      ! temporary arrays used in ZEEM_SA:
      Call mma_allocate(RWORK,(3*N-2),'ZEEM_RWORK')
      Call mma_allocate(HZEE,(N*(N+1)/2),'ZEEM_HZEE')
      Call mma_allocate(WORK,(2*N-1),'ZEEM_WORK')
      Call mma_allocate(W_c,N,'ZEEM_W_c')

      ! zero everything:
      Call dcopy_(3*N-2,0.0_wp,0,RWORK,1)
      Call zcopy_(N*(N+1)/2,(0.0_wp,0.0_wp),0,HZEE,1)
      Call zcopy_(2*N-1,(0.0_wp,0.0_wp),0,WORK,1)
      Call zcopy_(N,(0.0_wp,0.0_wp),0,W_c,1)

      Call dcopy_(   3,0.0_wp,0,ST,1)
      Call dcopy_(   N,0.0_wp,0,WZ,1)
      Call dcopy_(exch,0.0_wp,0,WM,1)
      Call zcopy_(  exch*exch,(0.0_wp,0.0_wp),0,ZM,1)
      Call zcopy_(3*exch*exch,(0.0_wp,0.0_wp),0,SZ,1)
      Call zcopy_(3*exch*exch,(0.0_wp,0.0_wp),0,MZ,1)
c start calculations:
      If (DBG) Write(6,*) 'Enter ZEEM::'
      If (DBG) Write(6,*) 'Input data:   N = ', N
      If (DBG) Write(6,*) 'Input data:   H = ', H
      If (DBG) Write(6,*) 'Input data:   X = ', X
      If (DBG) Write(6,*) 'Input data:   Y = ', Y
      If (DBG) Write(6,*) 'Input data:   Z = ', Z
      If (DBG) Write(6,*) 'Input data:  zJ = ', zJ
      If (DBG) Write(6,*) 'Input data: W() = ', W(1:N)
      If (DBG) Write(6,*) 'Input data: ST()= ', ST(1:3)
      If (DBG) Call prmom('INput data dM:',dM,N)
      If (DBG) Call prmom('INput data sM:',sM,N)
      If (DBG) Call xFlush(6)

      ! Build and diagonalize the Zeeman Hamiltonian
      ! most important output are: WM (energies) and ZM (eigenvectors)
      Call ZEEM_SA( N, H, X,Y,Z, W(1:N), dM(1:3,1:N,1:N),
     &             sM(1:3,1:N,1:N), ST, zJ, WM(1:N), ZM,
     &             DBG, RWORK, HZEE, WORK, W_c )
      If (DBG) Write(6,*) 'Exit ZEEM::'
      If (DBG) Call xFlush(6)


      Call DCOPY_(N, WM(1:N), 1, WZ(1:N), 1)
      If(N.ne.EXCH) Then
        Do i=N+1,EXCH
          WM(i)=W(i)
        End Do
      End If


      ! transform the momenta
      Call UTMU( EXCH, N, ZM, sM, SZ )
      Call UTMU( EXCH, N, ZM, dM, MZ )


      ! compute magnetization at different temperatures:
      Call dcopy_(  nT,0.0_wp,0,ZB,1)
      Call dcopy_(3*nT,0.0_wp,0,M,1)
      Call dcopy_(3*nT,0.0_wp,0,S,1)

      If (N==EXCH) Then
        Do iT=1,nT
          Do l=1,3
            If (sopt) Then
              Call calcmagn1( EXCH, WM, SZ(l,:,:), T(iT), S(l,iT),
     &                        ZB(iT) )
            End If
              Call calcmagn1( EXCH, WM, MZ(l,:,:), T(iT), M(l,iT),
     &                        ZB(iT) )
          End Do
        End Do
      Else
        Do iT=1,nT
          Do l=1,3
            If(sopt) Then
              Call calcmagn2( EXCH, N, WM, T(iT), H, SZ, X, Y, Z, l,
     &                        S(l,iT), ZB(iT) )
            End If
              Call calcmagn2( EXCH, N, WM, T(iT), H, MZ, X, Y, Z, l,
     &                        M(l,iT), ZB(iT) )
          End Do
        End Do !iT
      End If

      ! deallocate temporary arrays:
      Call mma_deallocate(RWORK)
      Call mma_deallocate(HZEE)
      Call mma_deallocate(WORK)
      Call mma_deallocate(W_c)

      Call mma_deallocate(WM)
      Call mma_deallocate(ZM)
      Call mma_deallocate(ST)
      Call mma_deallocate(SZ)
      Call mma_deallocate(MZ)

      Call qExit('MAGN_NO_MF')
      Return
      End
