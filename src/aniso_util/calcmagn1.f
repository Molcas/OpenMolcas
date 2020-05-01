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
      Subroutine calcmagn1(N,E,M,T,MT,Z)
!----------------------------------------------------------------------
!   Input data:
!     N =  number of the states to be considered, scalar integer, input
!     E =  Zeeman energy of these states, array(N), real, input
!     M =  momentum of these states ( diagonal), complex array (N,N), input
!          ( only diagonal elements are used )
!     T =  scalar real number denoting temperature in K, real, input
!   Output data:
!     MT=  computed momentum at temperature T, scalar real, output
!     Z =  Boltzmann statistical sum, real scalar
!----------------------------------------------------------------------
      Implicit None
      Integer, parameter          :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)         :: N
      Real(kind=8),intent(in)    :: E(N), T
      Complex(kind=8),intent(in) :: M(N,N)

      Real(kind=8),intent(out)   :: Z, MT
!-- local variables:
      Integer                 :: im,mp1,i
      Real(kind=8)           :: kB
!----------------------------------------------------------------------
      Call qEnter('calcmagn1')
      kB=0.6950356000_wp !   in cm^-1*K-1
      Z =0.0_wp
      MT=0.0_wp
      im=0
      mp1=0
      im=MOD(N,11)

      If (im.ne.0) Then
        Do i = 1,im !N   ! im
          Z = Z + EXP( -(E(i)-E(1))/kB/T )
          MT=MT + EXP( -(E(i)-E(1))/kB/T ) * DBLE( M(I,I) )
        End Do
        If ( N < 11 ) Then
          MT=MT/Z
        Return
        End If
      End If

      mp1 = im + 1
      Do i = mp1, N, 11
        Z = Z + EXP( -(E(i   )-E(1))/kB/T )
     &        + EXP( -(E(i+ 1)-E(1))/kB/T )
     &        + EXP( -(E(i+ 2)-E(1))/kB/T )
     &        + EXP( -(E(i+ 3)-E(1))/kB/T )
     &        + EXP( -(E(i+ 4)-E(1))/kB/T )
     &        + EXP( -(E(i+ 5)-E(1))/kB/T )
     &        + EXP( -(E(i+ 6)-E(1))/kB/T )
     &        + EXP( -(E(i+ 7)-E(1))/kB/T )
     &        + EXP( -(E(i+ 8)-E(1))/kB/T )
     &        + EXP( -(E(i+ 9)-E(1))/kB/T )
     &        + EXP( -(E(i+10)-E(1))/kB/T )

        MT=MT + EXP( -(E(i   )-E(1))/kB/T ) * DBLE( M(i   ,i   ) )
     &        + EXP( -(E(i+ 1)-E(1))/kB/T ) * DBLE( M(i+ 1,i+ 1) )
     &        + EXP( -(E(i+ 2)-E(1))/kB/T ) * DBLE( M(i+ 2,i+ 2) )
     &        + EXP( -(E(i+ 3)-E(1))/kB/T ) * DBLE( M(i+ 3,i+ 3) )
     &        + EXP( -(E(i+ 4)-E(1))/kB/T ) * DBLE( M(i+ 4,i+ 4) )
     &        + EXP( -(E(i+ 5)-E(1))/kB/T ) * DBLE( M(i+ 5,i+ 5) )
     &        + EXP( -(E(i+ 6)-E(1))/kB/T ) * DBLE( M(i+ 6,i+ 6) )
     &        + EXP( -(E(i+ 7)-E(1))/kB/T ) * DBLE( M(i+ 7,i+ 7) )
     &        + EXP( -(E(i+ 8)-E(1))/kB/T ) * DBLE( M(i+ 8,i+ 8) )
     &        + EXP( -(E(i+ 9)-E(1))/kB/T ) * DBLE( M(i+ 9,i+ 9) )
     &        + EXP( -(E(i+10)-E(1))/kB/T ) * DBLE( M(i+10,i+10) )
      End Do

      MT=MT/Z

      Call qExit('calcmagn1')
      Return
      End







