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
      Subroutine MAGN_ZJ_PAR( EXCH,N,X,Y,Z,H,W,zJ, dM,sM,nT,T,sopt,
     &                        WZ, ZB, S, M, thrs, m_paranoid, dbg )
c this Subroutine computes the magnetisation for zJ .ne. 0.0, paranoid accuracy
c
c definition of the variables:
c   EXCH -- total number of exchange states, Integer, input
c      N -- size of the Zeeman matrix, Integer, input, NM .LE. EXCH !
c  X,Y,Z -- projections of the magnetic field, specIfying the orientation of the applied
c           magnetic field, Real(kind=8) ::, input;  rule: ( X**2 + Y**2 + Z**2 = 1);
c      H -- strength of the magnetic field in Tesla, Real(kind=8) ::, input;
c      W -- energies of the exchange states; Real(kind=8) :: array (EXCH);
c     zJ -- parameter of intermolecular interaction, Real(kind=8) ::, input;
c   THRS -- threshold for convergence of the spin magnetisation. Real(kind=8) ::, input;
c           Of any importance only when (zJ.ne.0.0_wp), otherwise unused. Default value 1.D-8.
c     dM -- matrix of the magnetic moment, Complex(kind=8) :: (3,EXCH,EXCH) array, input;
c     sM -- matrix of the     spin moment, Complex(kind=8) :: (3,EXCH,EXCH) array, input;
c     nT -- number of temperature points for which magnetisation is computed, input;
c      T -- temperature values(in Kelvin) for which magnetisation is computed, input;
c   sopt -- logical parameter. If sopt=.true. Then spin magnetisation is computed.
c                              If sopt=.false. Then spin part is skipped.
c
c     WZ -- Zeeman energies, true values (not shifted to 0), in cm-1, Real(kind=8) :: (N) array, output;
c     ZB -- statistical Boltzmann distribution, for each temperature, Real(kind=8) :: (nT) array, output;
c      S -- spin magnetisation, Real(kind=8) :: (3,nT) array, output;
c      M -- magnetisation, Real(kind=8) :: (3,nT) array, output;
c m_paranoid --  logical parameter.
c            If m_paranoid = .true.  Then the average spin is computed for each temperature point exactly
c            If m_paranoid = .false. Then  the average spin is computed only for the lowest temperature point
c---------
c  temporary (local) variables:
c    MZ -- matrix of the magnetic moment, Complex(kind=8) :: (3,EXCH,EXCH) array
c    SZ -- matrix of the     spin moment, Complex(kind=8) :: (3,EXCH,EXCH) array
c    WM -- array containing the Zeeman eigenstates and, If N<EXCH the exchange eigenstates
c          for the states higher in energy than N, Real(kind=8) ::, (EXCH) array;
c    ZM -- Zeeman eigenvectors, (N,N) Complex(kind=8) :: array,
c  SCHK -- variable used for checking the convergence; SCHK = ABS(Sx) + ABS(Sy) + ABS(Sz), Real(kind=8) ::
c   i,j -- labers of the states;
c     l -- labels the cartesian component of the momentum (convention: x=1, y=2, z=3)
c    iT -- labes the temperature points;
c mxIter--defines the maximal number of iterations for determination of the average spin
c    ST -- value of the average spin of neighboring sites, Real(kind=8) :: (3) array;
      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)          :: EXCH, N, nT
      Real(kind=8), intent(in)    :: X, Y, Z, H, zJ
      Real(kind=8), intent(in)    :: W(EXCH), T(nT)
      Real(kind=8), intent(in)    :: thrs
      Complex(kind=8), intent(in) :: dM(3,EXCH,EXCH)
      Complex(kind=8), intent(in) :: sM(3,EXCH,EXCH)
      Logical, intent(in)          :: sopt
      Logical, intent(in)          :: m_paranoid
      Logical, intent(in)          :: dbg

      Real(kind=8), intent(out)   :: ZB(nT)
      Real(kind=8), intent(out)   :: WZ(N)
      Real(kind=8), intent(out)   :: S(3,nT)
      Real(kind=8), intent(out)   :: M(3,nT)

c local variables:
#include "stdalloc.fh"
      Integer, parameter :: mxIter=100
! defines the maximal number of iterations for determination of the average spin
!      parameter (mxIter=100) ! it is of any importance only If ( zJ.ne.0 )
      Integer          :: i, l, iT
      Real(kind=8)    :: ST(3), STsave(3)
      Real(kind=8), allocatable :: WM(:)
      Real(kind=8), allocatable :: RWORK(:)
      Complex(kind=8), allocatable :: HZEE(:), WORK(:), W_c(:)
      Complex(kind=8), allocatable :: ZM(:,:)
      Complex(kind=8), allocatable :: SZ(:,:,:)
      Complex(kind=8), allocatable :: MZ(:,:,:)
!                                      SZ(3,EXCH,EXCH), MZ(3,EXCH,EXCH)

c a few checks, before proceeding:
      Do iT=1,nT
        If (T(iT)==0.0_wp) Return
      End Do
      If(H==0.0_wp) Return
      If(N > EXCH) Return


! allocate memory:
      Call mma_allocate(WM,exch,'WM')
      Call dcopy_(exch,[0.0_wp],0,WM,1)

      Call mma_allocate(ZM,exch,exch,'ZM')
      Call zcopy_(exch*exch,[(0.0_wp,0.0_wp)],0,ZM,1)

      Call mma_allocate(SZ,3,exch,exch,'SZ')
      Call zcopy_(3*exch*exch,[(0.0_wp,0.0_wp)],0,SZ,1)

      Call mma_allocate(MZ,3,exch,exch,'MZ')
      Call zcopy_(3*exch*exch,[(0.0_wp,0.0_wp)],0,MZ,1)

      ! temporary arrays used in ZEEM_SA:
      Call mma_allocate(RWORK,(3*N-2),'ZEEM_RWORK')
      Call mma_allocate(HZEE,(N*(N+1)/2),'ZEEM_HZEE')
      Call mma_allocate(WORK,(2*N-1),'ZEEM_WORK')
      Call mma_allocate(W_c,N,'ZEEM_W_c')

      ! zero everything:
      Call dcopy_(3*N-2,[0.0_wp],0,RWORK,1)
      Call zcopy_(N*(N+1)/2,[(0.0_wp,0.0_wp)],0,HZEE,1)
      Call zcopy_(2*N-1,[(0.0_wp,0.0_wp)],0,WORK,1)
      Call zcopy_(N,[(0.0_wp,0.0_wp)],0,W_c,1)

c start calculations:
      ! code for the case (zJ .ne. 0 ):
      Call dcopy_(  nT,[0.0_wp],0,ZB,1)
      Call dcopy_(3*nT,[0.0_wp],0,M,1)
      Call dcopy_(3*nT,[0.0_wp],0,S,1)


      Do iT=1,nT
         ! determine first the average spin of neighboring
         ! molecules for each temperature point ST(1:3)
         If(m_paranoid) Then
            Call dcopy_(3,[0.0_wp],0,ST,1)
            Call dcopy_(3,[0.0_wp],0,STsave,1)
            Call mean_field( EXCH, N, H, X,Y,Z, zJ, T(iT), W, thrs,
     &                       DM, SM, ST, dbg )
            Call dcopy_(3,ST,1,STsave,1)
         Else
            ! i.e. when m_paranoid=.false.
            If (iT == 1) Then
               Call dcopy_(3,[0.0_wp],0,ST,1)
               Call dcopy_(3,[0.0_wp],0,STsave,1)
               Call mean_field( EXCH, N, H, X,Y,Z, zJ, T(iT), W, thrs,
     &                          DM, SM, ST, dbg )
               Call dcopy_(3,ST,1,STsave,1)
            Else
               ! use the last saved value of the ST:
               Call dcopy_(3,STsave,1,ST,1)
            End If
         End If
!        --------------------------------------------------------------
!        here  we have the value of the averaged spin for this temperature
!        proceed with the computation of magnetism for this temperature
         If(DBG) Write(6,'(A,3E13.5)')
     &                 'Average spin finished. ST on entrance to '//
     &                 'last ZEEM:',(ST(i),i=1,3)
         Call dcopy_(exch,[0.0_wp],0,WM,1)
         Call zcopy_(exch*exch,[(0.0_wp,0.0_wp)],0,ZM,1)


         Call ZEEM_SA( N, H, X,Y,Z, W(1:N),
     &              dM(1:3,1:N,1:N), sM(1:3,1:N,1:N),
     &              ST, zJ, WM(1:N), ZM(1:N,1:N),
     &              DBG, RWORK, HZEE, WORK, W_c )


         ! move WM energies to WZ:
         Call dcopy_(N, [0.0_wp], 0, WZ, 1)
         Call dcopy_(N,       WM, 1, WZ, 1)
c /// calculation of matrix elements of spin momentum in the basis of Zeeman states
        If(N<EXCH) Then
           Do i=N+1,EXCH
              WM(i)=W(i)
           End Do
        End If


        ! transform the momenta
        Call zcopy_(3*exch*exch, [(0.0_wp,0.0_wp)], 0, SZ, 1)
        Call zcopy_(3*exch*exch, [(0.0_wp,0.0_wp)], 0, MZ, 1)
        Call UTMU( EXCH, N, ZM(1:N,1:N), SM, SZ )
        Call UTMU( EXCH, N, ZM(1:N,1:N), DM, MZ )



        ! calculation of magnetizations at different temperatures:
         If (N==EXCH) Then
            Do l=1,3
               If (sopt) Then
                  Call calcmagn1( EXCH, WM, SZ(l,1:EXCH,1:EXCH),
     &                            T(iT), S(l,iT), ZB(iT) )
               End If
                 Call calcmagn1( EXCH, WM, MZ(l,1:EXCH,1:EXCH),
     &                           T(iT), M(l,iT), ZB(iT))
            End Do
         Else
            Do l=1,3
               If(sopt) Then
                  Call calcmagn2( EXCH, N, WM, T(iT), H,
     &                            SZ, X, Y, Z, l,
     &                            S(l,iT), ZB(iT) )
               End If
                  Call calcmagn2( EXCH, N, WM, T(iT), H,
     &                            MZ, X,Y,Z, L,
     &                            M(l,iT), ZB(iT) )
            End Do
         End If

      End Do ! iT

      ! deallocate temporary arrays:
      Call mma_deallocate(RWORK)
      Call mma_deallocate(HZEE)
      Call mma_deallocate(WORK)
      Call mma_deallocate(W_c)

      Call mma_deallocate(WM)
      Call mma_deallocate(ZM)
      Call mma_deallocate(SZ)
      Call mma_deallocate(MZ)
      Return
      End
