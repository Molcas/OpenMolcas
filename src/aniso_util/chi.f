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
      Subroutine chi( M1, M2, E, N, T, Z, X)
!     computes the chi  of one site
! definition of the variables:
!     N -- number of states included in the calculation of chi, Integer, input
!    M1 -- moment 1, Complex*16, (3,N,N) array, input
!    M2 -- moment 2, Complex*16, (3,N,N) array, input
!     E -- energy of the N states, Real(kind=wp) ::, (N) array, input
!     T -- temperature at which the chi is computed, Real(kind=wp) ::, input
!     Z -- statistical sum according to Boltzmann distribution law, Real(kind=wp) ::, output
!     X -- susceptibility tensor, Real(kind=wp) ::, (3,3) array, output
!--------
!  temporary (local) variables:
! iS,jS -- denote states over which the chi is computed
!    pB   -- partial Boltzmann population of a given state, Real(kind=wp) ::
!    dE -- energy diference E(i)-E(j)

      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)          :: N
      Real(kind=wp), intent(in)    :: E(N), T
      Complex(kind=wp), intent(in) :: M1(3,N,N), M2(3,N,N)
      Real(kind=wp), intent(out)   :: Z, X(3,3)
! local variables
      Integer                      :: i,j,iS,jS
      Real(kind=wp)                :: pB, dE, c2(3,3), R, F
      Real(kind=wp)                :: boltz_k

      Call qEnter('CHI')
      boltz_k=0.6950356_wp !   in cm^-1*k-1

      pB=0.0_wp
      dE=0.0_wp
      Z =0.0_wp
      Call dcopy_(3*3,0.0_wp,0,X,1)

      Do iS=1,N
      ! first loop over all states
         Call dcopy_(3*3,0.0_wp,0,c2,1)
         ! pB = statistical sum for state iS at temperature T
         pB = exp( -E(iS)/ boltz_k / T )
         ! accumulate the total statistical sum Z
         Z =  Z + pB
         Do jS=1,N
         ! second loop over all states
            dE=E(iS)-E(jS)
            ! set the multiplication factor:
            If(abs(dE) < 0.001_wp) Then
               F= 1.0_wp
            Else
               F=-2.0_wp*boltz_k * T / dE
            End If
            ! accumulate the contributions to the X tensor:
            Do i=1,3
               Do j=1,3
                  R=dble( M1(i,iS,jS)* conjg(M2(j,iS,jS)) )
                  c2(i,j) = c2(i,j) + F*R
               End Do ! j
            End Do ! i
         End Do ! jS

         ! add the (iS) contribution to the X tensor, according to its
         ! Boltzmann distribution:
         Call daxpy_(3*3,pB,c2,1,X,1)
      End Do ! iS

      ! scale the total tensor by the total statistical sum Z:
      Call dscal_(3*3, 1.0_wp/Z, X, 1)
      Call qExit('CHI')
      Return
      End
