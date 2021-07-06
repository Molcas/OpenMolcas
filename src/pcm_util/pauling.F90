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

function Pauling(N)
! Pauling radius for atomic number N (UFF when not defined)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Pauling
integer(kind=iwp), intent(in) :: N
real(kind=wp), parameter :: R(110) = [ &
  !  H      He
  1.20_wp,1.20_wp, &
  ! Li      Be       B       C       N       O       F      Ne
  1.37_wp,1.45_wp,1.45_wp,1.50_wp,1.50_wp,1.40_wp,1.35_wp,1.30_wp, &
  ! Na      Mg      Al      Si       P       S      Cl      Ar
  1.57_wp,1.36_wp,1.24_wp,1.17_wp,1.90_wp,1.85_wp,1.80_wp,1.88_wp, &
  !  K      Ca
  2.75_wp,  Zero , &
          ! Sc      Ti       V      Cr      Mn      Fe      Co      Ni      Cu      Zn
            Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,1.63_wp,1.40_wp,1.39_wp, &
                  ! Ga      Ge      As      Se      Br      Kr
                  1.87_wp,1.86_wp,2.00_wp,2.00_wp,1.95_wp,2.02_wp, &
  ! Rb      Sr
    Zero ,  Zero , &
          !  Y      Zr      Nb      Mo      Tc      Ru      Rh      Pd      Ag      Cd
            Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,1.63_wp,1.72_wp,1.58_wp, &
                  ! In      Sn      Sb      Te       I      Xe
                  1.93_wp,2.17_wp,2.20_wp,2.20_wp,2.15_wp,2.16_wp, &
  ! Cs(n.a.) Ba(n.a.)
    Zero ,  Zero , &
          ! La      Ce      Pr      Nd      Pm      Sm      Eu      Gd      Tb      Dy      Ho      Er      Tm      Yb
            Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero , &
          ! Lu      Hf      Ta       W      Re      Os      Ir      Pt      Au      Hg
            Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,1.72_wp,1.66_wp,1.55_wp, &
                  ! Tl      Pb      Bi      Po      At      Rn
                  1.96_wp,1.02_wp,  Zero ,  Zero ,  Zero ,  Zero , &
  ! Fr(n.a.) Ra(n.a.)
    Zero ,  Zero , &
          ! Ac      Th      Pa       U      Np      Pu      Am      Cm      Bk      Cf      Es      Fm      Md
            Zero ,  Zero ,  Zero ,1.86_wp,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero , &
          !
            Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero &
]
real(kind=wp), external :: UFF_radii

! If the Pauling radius is not defined, resort to the UFF value
if (R(N) == Zero) then
  Pauling = UFF_radii(N)
else
  Pauling = R(N)
end if

return

end function Pauling
