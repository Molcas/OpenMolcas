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

module CovRad_Data

use Definitions, only: wp

real*8 :: CovRadT_(0:92) = [ &
          0.00_wp, &
          3.08_wp,0.93_wp, &
          1.23_wp,0.90_wp,0.82_wp,3.59_wp,3.31_wp,3.08_wp,0.72_wp,0.71_wp, &
          1.54_wp,0.36_wp,1.18_wp,1.11_wp,1.06_wp,1.02_wp,0.99_wp,0.98_wp, &
          2.03_wp,0.74_wp, &
          1.44_wp,1.32_wp,1.22_wp,1.18_wp,1.17_wp, &
          1.17_wp,1.16_wp,1.15_wp,1.17_wp,1.25_wp, &
          1.26_wp,1.22_wp,1.20_wp,1.16_wp,1.14_wp,1.89_wp, &
          2.16_wp,0.91_wp, &
          1.62_wp,1.45_wp,1.34_wp,1.30_wp,1.27_wp, &
          1.25_wp,1.25_wp,1.28_wp,1.34_wp,1.41_wp, &
          1.44_wp,1.41_wp,1.40_wp,1.36_wp,1.33_wp,1.31_wp, &
          2.35_wp,0.98_wp,1.25_wp, &
          1.65_wp,1.65_wp,1.64_wp,1.63_wp,1.62_wp,1.85_wp,1.61_wp, &
          1.59_wp,1.59_wp,1.58_wp,1.57_wp,1.56_wp,1.70_wp,1.56_wp, &
          1.44_wp,1.34_wp,1.30_wp,0.28_wp, &
          1.26_wp,1.27_wp,1.30_wp,1.34_wp,1.49_wp, &
          1.48_wp,1.47_wp,1.46_wp,1.53_wp,1.47_wp,1.3_wp, &
          2.50_wp,2.05_wp,1.50_wp, &
          1.65_wp,1.50_wp,1.42_wp]
real*8 :: CovRad_(0:86) = [ &
          0.000_wp, &
          0.643_wp,0.643_wp,2.457_wp,1.909_wp,1.587_wp,1.436_wp,1.209_wp, & ! 1-7
          1.096_wp,1.020_wp,0.945_wp,2.986_wp,2.646_wp,2.400_wp,2.192_wp, & ! 8-14
          2.060_wp,1.890_wp,1.795_wp,1.701_wp,3.836_wp,3.288_wp,2.721_wp, & ! 15-21
          2.494_wp,2.305_wp,2.230_wp,2.211_wp,2.211_wp,2.192_wp,2.173_wp, & ! 22-28
          2.211_wp,2.362_wp,2.381_wp,2.305_wp,2.268_wp,2.192_wp,2.154_wp, & ! 29-35
          2.116_wp,4.082_wp,3.609_wp,3.061_wp,2.740_wp,2.532_wp,2.457_wp, & ! 36-42
          2.400_wp,2.362_wp,2.362_wp,2.419_wp,2.532_wp,2.797_wp,2.721_wp, & ! 43-49
          2.665_wp,2.646_wp,2.570_wp,2.513_wp,2.476_wp,4.441_wp,3.742_wp, & ! 50-56
          3.194_wp,3.118_wp,3.118_wp,3.099_wp,3.080_wp,3.061_wp,3.496_wp, & ! 57-63
          3.042_wp,3.005_wp,3.005_wp,2.986_wp,2.967_wp,2.948_wp,2.948_wp, & ! 64-70
          2.948_wp,2.721_wp,2.532_wp,2.457_wp,2.419_wp,2.381_wp,2.400_wp, & ! 71-77
          2.457_wp,2.532_wp,2.816_wp,2.797_wp,2.778_wp,2.759_wp,2.759_wp, & ! 78-84
          2.740_wp,2.710_wp]                                                ! 85-

end module CovRad_Data
