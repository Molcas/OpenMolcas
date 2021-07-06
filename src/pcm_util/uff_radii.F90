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

function UFF_radii(IA)

use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: UFF_radii
integer(kind=iwp), intent(in) :: IA
! Warning: the following are atomic diameters.
real(kind=wp), parameter :: Rii(0:104) = [ &
    Zero  , &
  !  H        He
  2.886_wp,2.362_wp, &
  ! Li        Be         B         C         N         O         F        Ne
  2.451_wp,2.745_wp,4.083_wp,3.851_wp,3.660_wp,3.500_wp,3.364_wp,3.243_wp, &
  ! Na        Mg        Al        Si         P         S        Cl        Ar
  2.983_wp,3.021_wp,4.499_wp,4.295_wp,4.147_wp,4.035_wp,3.947_wp,3.868_wp, &
  !  K        Ca
  3.812_wp,3.399_wp, &
           ! Sc        Ti         V        Cr        Mn        Fe        Co        Ni        Cu        Zn
           3.295_wp,3.175_wp,3.144_wp,3.023_wp,2.961_wp,2.912_wp,2.872_wp,2.834_wp,3.495_wp,2.763_wp, &
                    ! Ga        Ge        As        Se        Br        Kr
                    4.383_wp,4.280_wp,4.230_wp,4.205_wp,4.189_wp,4.141_wp, &
  ! Rb        Sr (+2)
  4.114_wp,3.641_wp, &
           !  Y  (+3)  Zr (+4)   Nb (+5)   Mo (+6)   Tc (+5)   Ru (+2)   Rh (+3)   Pd (+2)   Ag (+1)   Cd (+2)
           3.345_wp,3.124_wp,3.165_wp,3.052_wp,2.998_wp,2.963_wp,2.929_wp,2.899_wp,3.148_wp,2.848_wp, &
                    !       In            Sn           Sb              Te          I       Xe
                    4.463_wp,4.392_wp,4.420_wp,4.470_wp,4.50_wp,4.404_wp, &
  ! Cs        Ba (+2)
  4.517_wp,3.703_wp, &
           ! La (+3)   Ce (+3)   Pr (+3)   Nd (+3)   Pm (+3)   Sm (+3)   Eu (+3)
           ! Gd (+3)   Tb (+3)   Dy (+3)   Ho (+3)   Er (+3)   Tm (+3)   Yb (+3)
           3.522_wp,3.556_wp,3.606_wp,3.575_wp,3.547_wp,3.520_wp,3.493_wp, &
           3.368_wp,3.451_wp,3.428_wp,3.409_wp,3.391_wp,3.374_wp,3.355_wp, &
           ! Lu (+3)   Hf (+4)   Ta (+5)    W(+4,+6) Re(+5,+7) Os (+6)   Ir (+3)   Pt        Au        Hg
           3.640_wp,3.141_wp,3.170_wp,3.069_wp,2.954_wp,3.120_wp,2.840_wp,2.754_wp,3.293_wp,2.705_wp, &
                    ! Tl        Pb        Bi (+3)   Po (+2)   At        Rn (+4)
                    4.347_wp,4.297_wp,4.370_wp,4.709_wp,4.750_wp,4.765_wp, &
  ! Fr        Ra (+2)
  4.900_wp,3.677_wp, &
           ! Ac (+3)   Th (+4)   Pa (+4)    U (+4)   Np (+4)   Pu (+4)   Am (+4)
           ! Cm (+3)   Bk (+3)   Cf (+3)   Es (+3)   Fm (+3)   Md (+3)   No (+3)
           3.478_wp,3.396_wp,3.424_wp,3.395_wp,3.424_wp,3.424_wp,3.381_wp, &
           3.326_wp,3.339_wp,3.313_wp,3.299_wp,3.286_wp,3.274_wp,3.248_wp, &
           ! Lw(+3)
           3.236_wp,3.500_wp &
]

UFF_radii = Rii(min(max(IA,lbound(Rii,1)),ubound(Rii,1)))*Half

return

end function UFF_radii
