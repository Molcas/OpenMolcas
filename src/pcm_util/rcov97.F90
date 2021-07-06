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

function RCov97(IA,IB)
! This function returns an estimated covalent bond distance (in Ang)
! between atoms of atomic numbers IA and IB. Setting IB to 0 returns
! the covalent radius of IA. Parameters for atoms heavier than At are
! taken from the UFF force field (A.K.Rappe',C.J.Casewit,K.S.Colwell,
! W.A.Goddard III and W.M.Skiff, J.Am.Chem.Soc. 114,10024 (1992)).

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: RCov97
integer(kind=iwp), intent(in) :: IA, IB
real(kind=wp), parameter :: Rii(0:104) = [ &
    Zero  , &
  0.354_wp,0.849_wp, &
  1.336_wp,1.010_wp,0.838_wp,0.757_wp,0.700_wp,0.658_wp,0.668_wp,0.920_wp, &
  1.539_wp,1.421_wp,1.244_wp,1.117_wp,1.101_wp,1.064_wp,1.044_wp,1.032_wp, &
  1.953_wp,1.761_wp, &
           1.513_wp,1.412_wp,1.402_wp,1.345_wp,1.382_wp,1.270_wp,1.241_wp,1.164_wp,1.302_wp,1.193_wp, &
                    1.260_wp,1.197_wp,1.211_wp,1.190_wp,1.192_wp,1.147_wp, &
  2.260_wp,2.052_wp, &
           1.698_wp,1.564_wp,1.473_wp,1.467_wp,1.322_wp,1.478_wp,1.332_wp,1.338_wp,1.386_wp,1.403_wp, &
                    1.459_wp,1.398_wp,1.407_wp,1.386_wp,1.382_wp,1.267_wp, &
  2.570_wp,2.277_wp, &
           1.943_wp,1.841_wp,1.823_wp,1.816_wp,1.801_wp,1.780_wp,1.771_wp, &
           1.735_wp,1.732_wp,1.710_wp,1.696_wp,1.673_wp,1.660_wp,1.637_wp, &
           1.671_wp,1.611_wp,1.511_wp,1.392_wp,1.372_wp,1.372_wp,1.371_wp,1.364_wp,1.262_wp,1.340_wp, &
                    1.518_wp,1.459_wp,1.512_wp,1.500_wp,1.545_wp,1.420_wp, &
  2.880_wp,2.512_wp, &
           1.983_wp,1.721_wp,1.711_wp,1.684_wp,1.666_wp,1.657_wp,1.660_wp, &
           1.801_wp,1.761_wp,1.750_wp,1.724_wp,1.712_wp,1.689_wp,1.679_wp, &
           1.698_wp,1.850_wp &
]

RCov97 = Rii(min(max(IA,lbound(Rii,1)),ubound(Rii,1)))+Rii(min(max(IB,lbound(Rii,1)),ubound(Rii,1)))

return

end function RCov97
