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
      Subroutine Set_knm(knm)

!----------------------------------------!
!  knm( rank, projection )
!  are proportionality coefficients between
!    ITO ( JCP 137, 064112 (2012); doi: 10.1063/1.4739763 ) and
!    Extended Stevens Operators (ESO)
!
! were obtained by direct comparison the matrix elements
! of the ITO and ESO for each rank n and projection m
!----------------------------------------!

      Implicit None

      Integer, parameter :: wp=SELECTED_REAL_KIND(p=15,r=307)

      Integer       :: i, j
      Real(kind=wp) :: knm(12,0:12)

      Do i=1,12
         Do j=0,12
            knm(i,j)=0.0_wp
         End Do
      End Do
!----------------------------------------!
      knm( 1, 0) = 1.0_wp
      knm( 1, 1) = sqrt( 0.500_wp )
!----------------------------------------!
      knm( 2, 0) = 1.0_wp
      knm( 2, 1) = sqrt( 6.000_wp )
      knm( 2, 2) = sqrt( 1.500_wp )
!----------------------------------------!
      knm( 3, 0) = 1.0_wp
      knm( 3, 1) = sqrt( 0.750_wp )
      knm( 3, 2) = sqrt( 7.500_wp )
      knm( 3, 3) = sqrt( 1.250_wp )
!----------------------------------------!
      knm( 4, 0) = 1.0_wp
      knm( 4, 1) = sqrt(  20.000_wp )
      knm( 4, 2) = sqrt(  10.000_wp )
      knm( 4, 3) = sqrt( 140.000_wp )
      knm( 4, 4) = sqrt(  17.500_wp )
!----------------------------------------!
      knm( 5, 0) = 1.0_wp
      knm( 5, 1) = sqrt(   7.500_wp )
      knm( 5, 2) = sqrt( 210.000_wp )
      knm( 5, 3) = sqrt(   8.750_wp )
      knm( 5, 4) = sqrt( 157.500_wp )
      knm( 5, 5) = sqrt(  15.750_wp )
!----------------------------------------!
      knm( 6, 0) = 1.0_wp
      knm( 6, 1) = sqrt(  42.000_wp )
      knm( 6, 2) = sqrt(  26.250_wp )
      knm( 6, 3) = sqrt( 105.000_wp )
      knm( 6, 4) = sqrt(  31.500_wp )
      knm( 6, 5) = sqrt( 693.000_wp )
      knm( 6, 6) = sqrt(  57.750_wp )
!----------------------------------------!
      knm( 7, 0) = 1.0_wp
      knm( 7, 1) = sqrt(   0.875_wp )
      knm( 7, 2) = sqrt(   5.250_wp )
      knm( 7, 3) = sqrt(   2.625_wp )
      knm( 7, 4) = sqrt( 115.500_wp )
      knm( 7, 5) = sqrt(  28.875_wp )
      knm( 7, 6) = sqrt( 750.750_wp )
      knm( 7, 7) = sqrt( 107.250_wp )
!----------------------------------------!
      knm( 8, 0) = 1.0_wp
      knm( 8, 1) = sqrt(    72.000_wp )
      knm( 8, 2) = sqrt(  1260.000_wp )
      knm( 8, 3) = sqrt(  9240.000_wp )
      knm( 8, 4) = sqrt(  1386.000_wp )
      knm( 8, 5) = sqrt( 72072.000_wp )
      knm( 8, 6) = sqrt(  1716.000_wp )
      knm( 8, 7) = sqrt( 51480.000_wp )
      knm( 8, 8) = sqrt(  3217.500_wp )
!----------------------------------------!
      knm( 9, 0) = 1.0_wp
      knm( 9, 1) = sqrt(    22.500_wp )
      knm( 9, 2) = sqrt(  1980.000_wp )
      knm( 9, 3) = sqrt(  1155.000_wp )
      knm( 9, 4) = sqrt( 90090.000_wp )
      knm( 9, 5) = sqrt(  1287.000_wp )
      knm( 9, 6) = sqrt(  8580.000_wp )
      knm( 9, 7) = sqrt(  1608.750_wp )
      knm( 9, 8) = sqrt( 54697.500_wp )
      knm( 9, 9) = sqrt(  3038.750_wp )
!----------------------------------------!
      knm(10, 0) = 1.0_wp
      knm(10, 1) = sqrt(    110.000_wp )
      knm(10, 2) = sqrt(     82.500_wp )
      knm(10, 3) = sqrt(   8580.000_wp )
      knm(10, 4) = sqrt(   4290.000_wp )
      knm(10, 5) = sqrt(   1716.000_wp )
      knm(10, 6) = sqrt(    536.250_wp )
      knm(10, 7) = sqrt(  36465.000_wp )
      knm(10, 8) = sqrt(   6077.500_wp )
      knm(10, 9) = sqrt( 230945.000_wp )
      knm(10,10) = sqrt(  11547.250_wp )
!----------------------------------------!
      knm(11, 0) = 1.0_wp
      knm(11, 1) = sqrt(      8.250_wp )
      knm(11, 2) = sqrt(   1072.500_wp )
      knm(11, 3) = sqrt(   3753.750_wp )
      knm(11, 4) = sqrt(  18018.000_wp )
      knm(11, 5) = sqrt(    160.875_wp )
      knm(11, 6) = sqrt(   1823.250_wp )
      knm(11, 7) = sqrt(   4558.125_wp )
      knm(11, 8) = sqrt( 346417.500_wp )
      knm(11, 9) = sqrt(   5773.625_wp )
      knm(11,10) = sqrt( 242492.250_wp )
      knm(11,11) = sqrt(  11022.375_wp )
!----------------------------------------!
      knm(12, 0) = 1.0_wp
      knm(12, 1) = sqrt(     156.000_wp )
      knm(12, 2) = sqrt(    6006.000_wp )
      knm(12, 3) = sqrt(    4004.000_wp )
      knm(12, 4) = sqrt(    2252.250_wp )
      knm(12, 5) = sqrt(  306306.000_wp )
      knm(12, 6) = sqrt(    2431.000_wp )
      knm(12, 7) = sqrt(  277134.000_wp )
      knm(12, 8) = sqrt(   69283.500_wp )
      knm(12, 9) = sqrt(  646646.000_wp )
      knm(12,10) = sqrt(   88179.000_wp )
      knm(12,11) = sqrt( 4056234.000_wp )
      knm(12,12) = sqrt(  169009.750_wp )
!----------------------------------------!

      Return
      End subroutine set_knm
