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
      Subroutine chi_sum( N, Xex,Zex, XL,ZL, XR,ZR, iopt, X, Z )
c   computes the total CHI, provided all input values are provided
c   according to the desired partition scheme  (iopt = 1 or 2 )
c   all X tensors must be given in the general coordinate system
c
c definition of the variables:
c     N  -- total number of magnetic sites, Integer, input
c    Xex -- susceptibility tensor arising form the exchange states only, Real(kind=wp) :: , (3,3) array, input
c    Zex -- statistical sum according to Boltzmann distribution law of the exchange states, Real(kind=wp) :: , input
c    XL  -- susceptibility tensor arising from LOCAL states (all), Real(kind=wp) :: , (N,3,3) array, input
c    ZL  -- statistical sum according to Boltzmann distribution law arising from LOCAL states (all), Real(kind=wp) :: , (N) array, input
c    XR  -- susceptibility tensor arising from LOCAL states (exchange only), Real(kind=wp) :: , (N,3,3) array, input
c    ZR  -- statistical sum according to Boltzmann distribution law arising from LOCAL states (exchange only), Real(kind=wp) :: , (N) array, input
c   iopt -- option allowing to choose the desired formula, (Integer, input):
c           iopt=1  =>  formula for weak exchange limit ( new derivation)
c           iopt=2  =>  formula for strong exchange limit
c           iopt=3  =>  formula for strong exchange limit
c    X   -- total susceptibility, Real(kind=wp) :: , (3,3) array, output
c    Z   -- total statistical sum according to Boltzmann distribution, Real(kind=wp) :: , output
c---------
c  temporary (local) variables:
c
      Implicit None
      Integer, parameter        :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)       :: N, iopt
      Real(kind=wp), intent(in) :: Xex(3,3), Zex
      Real(kind=wp), intent(in) :: XL(N,3,3), ZL(N)
      Real(kind=wp), intent(in) :: XR(N,3,3), ZR(N)
      Real(kind=wp), intent(out):: X(3,3), Z
c local variables
      Integer       :: i, ic, jc
      Real(kind=wp) :: ZLT, ZRT, XLT(3,3), XRT(3,3)

      ZLT=1.0_wp
      ZRT=1.0_wp
      Z=0.0_wp
      Call dcopy_(3*3,0.0_wp,0,XRT,1)
      Call dcopy_(3*3,0.0_wp,0,XLT,1)
      Call dcopy_(3*3,0.0_wp,0,X,1)
      ! compute the total ZT
      If ( iopt == 1 ) Then
        ! my formula (simple):
        Do i=1,N
          ZLT=ZLT*ZL(i)
          ZRT=ZRT*ZR(i)
        End Do
        Z = Zex+ZLT-ZRT
        Do ic=1,3
          Do jc=1,3
            Do i=1,N
              XLT(ic,jc) = XLT(ic,jc) + XL(i,ic,jc)
              XRT(ic,jc) = XRT(ic,jc) + XR(i,ic,jc)
            End Do
            X(ic,jc) = Xex(ic,jc) + XLT(ic,jc) - XRT(ic,jc)
          End Do
        End Do

      Else If ( iopt == 2 ) Then
        ! "thesis formula:"
        Do i=1,N
          ZLT=ZLT*ZL(i)
          ZRT=ZRT*ZR(i)
        End Do
        Z = Zex+ZLT-ZRT
        Do ic=1,3
          Do jc=1,3
            Do i=1,N
              XLT(ic,jc) = XLT(ic,jc) + XL(i,ic,jc)*ZLT
              XRT(ic,jc) = XRT(ic,jc) + XR(i,ic,jc)*ZRT
            End Do
            X(ic,jc) = (Xex(ic,jc)*Zex + XLT(ic,jc) - XRT(ic,jc))/Z
          End Do
        End Do

      Else If ( iopt == 3 ) Then
        ! "weird formula as implemented in some version of the
        ! code, e.g. in the 2013 version:"
        Do i=1,N
          ZLT=ZLT*ZL(i)
          ZRT=ZRT*ZR(i)
        End Do
        ZLT = Zex*ZLT/ZRT
        ZRT = Zex
        Z   = ZLT - ZRT + Zex
        Do ic=1,3
          Do jc=1,3
            Do i=1,N
              XLT(ic,jc) = XLT(ic,jc) + XL(i,ic,jc)
              XRT(ic,jc) = XRT(ic,jc) + XR(i,ic,jc)
            End Do
                X(ic,jc) = ( XLT(ic,jc)*ZLT
     &                     - XRT(ic,jc)*ZRT
     &                     + Xex(ic,jc)*Zex ) / Z
          End Do
        End Do

      Else

        Write(6,'(A)') 'chi_sum: IOPT parameter out of range'
        Write(6,'(A,i8)') 'IOPT = ', IOPT
        Return

      End If

      Return
      End
