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
      Subroutine calcmagn2(N,NM,W,T,H,M,dX,dY,dZ,L, MT,Z)
!  this subroutine sums the contribution to molar magnetization from
!  different states:
!    -- from the states obtained in exact Zeeman diagonalization (NM) and
!    -- from the states higher in energy via second order perturbation
!----------------------------------------------------------------------
!  Input data:
!     N  -  total number of states
!    NM  -  number of states enterring the exact Zeeman diagonalization
!          (NM <= N)
!   W(N) -  energy of the input states
!
!
      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)          :: N, NM, L
      Real(kind=8), intent(in)    :: W(N), T, dX, dY, dZ, H
      Real(kind=8), intent(out)   :: MT, Z
      Complex(kind=8), intent(in) :: M(3,N,N)

      Integer                    :: i, j
      Real(kind=8)              :: pB, dltw, S2, S1, mB, kB
      Call qEnter('calcmagn2')
c /// constants
      mB=0.4668643740_wp                  ! * in cm-1*T-1
      kB=0.69503560_wp                    !   in cm^-1*K-1
      pB=0.0_wp
      Z =0.0_wp
      MT=0.0_wp
      DLTW=0.0_wp

      Do I=1,N
        pB=EXP(-(W(I)-W(1))/kB/T)
        Z=Z+pB
        If(I.LE.NM) Then
c  case when I <= NM
          S2=0.0_wp
          S2=DBLE(M(L,I,I))
          Do J = NM+1, N
            DLTW=W(I)-W(J)
            S1=0.0_wp
            S1=DBLE( M(L,I,J)*CONJG(M(1,I,J)) )*dX
     &        +DBLE( M(L,I,J)*CONJG(M(2,I,J)) )*dY
     &        +DBLE( M(L,I,J)*CONJG(M(3,I,J)) )*dZ
            S2=S2 - 2.0_wp * mB * H * S1 / DLTW
          End Do ! J
        Else  !  I
c  case when I > NM
          Do J=1,N
            DLTW=W(I)-W(J)
            S1=0.0_wp
            S2=0.0_wp
            S1=DBLE( M(L,I,J)*CONJG(M(1,I,J)) )*dX
     &        +DBLE( M(L,I,J)*CONJG(M(2,I,J)) )*dY
     &        +DBLE( M(L,I,J)*CONJG(M(3,I,J)) )*dZ

            If(abs(DLTW).lt.1.D-3) Then
              S2=S2 + 1.0_wp * mB * H * S1 / kB / T
            Else
              S2=S2 - 2.0_wp * mB * H * S1 / DLTW
            End If
          End Do ! J
        End If  !  I
        MT=MT+pB*S2
      End Do !I

      MT=MT/Z

      Call qExit('calcmagn2')
      Return
      End
