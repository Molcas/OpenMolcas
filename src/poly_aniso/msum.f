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
      Subroutine Msum( N, Mex,Zex, ML,ZL, MR,ZR, iopt, M, Z )
!   computes the total M, provided all input values are provided
!   according to the desired partition scheme  (iopt = 1 or 2 )
!   all M vectors must be given in the general coordinate system
!
! definition of the variables:
!     N  -- total number of magnetic sites, Integer, input
!    Mex -- magnetisation vector arising from the exchange states only, Real(kind=8) ::, (3) array, input
!    Zex -- statistical sum according to Boltzmann distribution law of the exchange states, Real(kind=8) ::, input
!    ML  -- magnetisation vector arising from LOCAL states (all), Real(kind=8) ::, (N,3) array, input
!    ZL  -- statistical sum according to Boltzmann distribution law arising from LOCAL states (all), Real(kind=8) ::, (N) array, input
!    MR  -- magnetisation vector arising from LOCAL states (exchange only), Real(kind=8) ::, (N,3) array, input
!    ZR  -- statistical sum according to Boltzmann distribution law arising from LOCAL states (exchange only), Real(kind=8) ::, (N) array, input
!   iopt -- option allowing to choose the desired formula, (Integer, input):
!           iopt=1  =>  formula for weak exchange limit ( new derivation)
!           iopt=2  =>  formula for strong exchange limit ( simple sumation of moments),size consistent;
!           iopt=3  =>  formula for strong exchange limit ( not to be used...)
!    M   -- total magnwtisation, Real(kind=8) ::, (3) array, output
!    Z   -- total statistical sum according to Boltzmann distribution, Real(kind=8) ::, output
!---------
!  temporary (local) variables:
!
      Implicit None
      Integer, parameter        :: wp=kind(0.d0)
      Integer       :: N, iopt
      Real(kind=8) :: Mex(3), Zex
      Real(kind=8) :: ML(N,3), ZL(N)
      Real(kind=8) :: MR(N,3), ZR(N)
      Real(kind=8) :: M(3), Z
! local variables
      Integer       :: i, ic
      Real(kind=8) :: ZLT, ZRT, MLT(3), MRT(3)
      Logical       :: DBG
      DBG=.false.
      ZLT=1.0_wp
      ZRT=1.0_wp
      Z=0.0_wp
      MLT(:)=0.0_wp
      MRT(:)=0.0_wp
      M(:)=0.0_wp

      If( iopt .eq. 1 ) Then
! compute the total ZT
! my formula (simple):
        Do i=1,N
          ZLT=ZLT*ZL(i)
          ZRT=ZRT*ZR(i)
        End Do
        Z=Zex+ZLT-ZRT
        Do ic=1,3
          Do i=1,N
            MLT(ic) = MLT(ic) + ML(i,ic)
            MRT(ic) = MRT(ic) + MR(i,ic)
          End Do
          M(ic) = Mex(ic) + MLT(ic) - MRT(ic)
        End Do
        If (DBG) Then
          Write(6,'(A,3F10.6)') 'Contribution from exchange states',    &
     &                           (Mex(ic),ic=1,3)
          Write(6,'(A,3F10.6)') 'Contribution from excited states ',    &
     &                           (MLT(ic) - MRT(ic),ic=1,3)
        End If

      Else If( iopt .eq. 2 ) Then
! "thesis formula:"
        Do i=1,N
          ZLT=ZLT*ZL(i)
          ZRT=ZRT*ZR(i)
        End Do
        Z = Zex+ZLT-ZRT
        Do ic=1,3
          Do i=1,N
            MLT(ic) = MLT(ic) + ML(i,ic)*ZLT
            MRT(ic) = MRT(ic) + MR(i,ic)*ZRT
          End Do
          M(ic) = (Mex(ic)*Zex + MLT(ic) - MRT(ic))/Z
        End Do

      Else

        Write(6,'(A)') 'chi_sum: IOPT parameter out of range'
        Write(6,'(A,i8)') 'IOPT = ', IOPT
        Return

      End If

      Return
      End
