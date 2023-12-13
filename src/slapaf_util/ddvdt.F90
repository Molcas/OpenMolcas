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

module ddvdt

use Definitions, only: wp

implicit none
private

#define _OLD_CODE_
                            ! Parameters for HMF, updated 2007-11-08
real(kind=wp), parameter :: aAV(3,3) = reshape([1.0000_wp,0.3949_wp,0.3949_wp, &
                                                0.3949_wp,0.2800_wp,0.1200_wp, &
                                                0.3949_wp,0.1200_wp,0.0600_wp],[3,3]), &
                            rAV(3,3) = reshape([1.3500_wp,2.1000_wp,2.5300_wp, &
                                                2.1000_wp,2.8700_wp,3.8000_wp, &
                                                2.5300_wp,3.8000_wp,4.5000_wp],[3,3]), &
                            ! Original parameters from the 1999 paper
                            !aAV(3,3) = reshape([1.0000_wp,0.3949_wp,0.3949_wp, &
                            !                    0.3949_wp,0.2800_wp,0.2800_wp, &
                            !                    0.3949_wp,0.2800_wp,0.2800_wp],[3,3]) &
                            !rAV(3,3) = reshape([1.3500_wp,2.1000_wp,2.5300_wp, &
                            !                    2.1000_wp,2.8700_wp,3.4000_wp, &
                            !                    2.5300_wp,3.4000_wp,3.4000_wp],[3,3]) &
                            f_Const_Min = 1.0e-3_wp, & ! 1.0e-2_wp, Zero
#                           ifdef _OLD_CODE_
                            ! HMF augmented by weak forces, MGS
                            alpha_vdW = 5.0_wp, &
                            r_ref_vdW(3,3) = reshape([0.0_wp,3.6_wp,3.6_wp, &
                                                      3.6_wp,5.3_wp,5.3_wp, &
                                                      3.6_wp,5.3_wp,5.3_wp],[3,3]), &
#                           else
                            ! approximate with sums of average vdW radius
                            alpha_vdW = 10.0_wp, &
                            r_ref_vdW(3,3) = reshape([0.0_wp,5.5_wp,5.5_wp, &
                                                      5.5_wp,6.3_wp,6.3_wp, &
                                                      5.5_wp,6.3_wp,6.3_wp],[3,3]), &
#                           endif
                            ! a la Schlegel
                            A_Str = 1.734_wp, &
                            B_Str(6) = [-0.244_wp,0.352_wp,1.085_wp,0.660_wp,1.522_wp,2.068_wp], &
                            A_StrH(2) = [0.3601_wp,1.944_wp], &
                            ! Parameters from 1999 paper (rkr, rkr_vdw)
                            rkr = 0.45_wp, &
                            rkr_vdW = 0.05_wp, &
                            ! Parameter from 1999 paper (rkf)
                            !rkf = 0.15_wp, &
                            !rkf = 0.2_wp, &
                            rkf = 0.1_wp, &
                            A_Bend(2) = [0.16_wp,0.25_wp], &
                            ! Parameter from 1999 paper (rkt)
                            !rkt = 0.005_wp, &
                            rkt = 0.0025_wp, &
                            A_Trsn(2) = [0.0023_wp,0.0700_wp], &
                            !rko = 0.16_wp, &
                            rko = 0.0016_wp, &
                            Trans_Const = 0.5_wp, &
                            Rot_Const = 0.5_wp

public :: A_Bend, A_Str, A_StrH, A_Trsn, aAV, alpha_vdW, B_Str, f_Const_Min, r_ref_vdW, rAV, rko, rkf, rkr, rkr_vdW, rkt, &
          Rot_Const, Trans_Const

end module ddvdt
