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

function RList(IA)
! Assigns Caillet-Claverie's atomic parameters for dispersion and repulsion.

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: RList
integer(kind=iwp), intent(in) :: IA
real(kind=wp), parameter :: RW(110) = [ &
  !  H      He
  1.20_wp,1.28_wp, &
  ! Li      Be       B       C       N       O       F      Ne
    Zero ,  Zero ,  Zero ,1.70_wp,1.60_wp,1.50_wp,1.45_wp,1.38_wp, &
  ! Na      Mg      Al      Si       P       S      Cl      Ar
  1.20_wp,  Zero ,  Zero ,  Zero ,1.85_wp,1.80_wp,1.76_wp,1.66_wp, &
  !  K      Ca
  1.46_wp,  Zero , &
          ! Sc      Ti       V      Cr      Mn      Fe      Co      Ni      Cu      Zn
            Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero , &
                  ! Ga      Ge      As      Se      Br      Kr
                    Zero ,  Zero ,  Zero ,  Zero ,1.85_wp,1.76_wp, &
  ! Rb      Sr
    Zero ,  Zero , &
          !  Y      Zr      Nb      Mo      Tc      Ru      Rh      Pd      Ag      Cd
            Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero , &
                  ! In      Sn      Sb      Te       I      Xe
                    Zero ,  Zero ,  Zero ,  Zero ,1.96_wp,1.85_wp, &
  ! Cs      Ba
    Zero ,  Zero , &
          ! La      Ce      Pr      Nd      Pm      Sm      Eu      Gd      Tb      Dy      Ho      Er      Tm      Yb
            Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero , &
          ! Lu      Hf      Ta       W      Re      Os      Ir      Pt      Au      Hg
            Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero , &
                  ! Tl      Pb      Bi      Po      At      Rn
                    Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero , &
  ! Fr      Ra
    Zero ,  Zero , &
          ! Ac      Th      Pa       U      Np      Pu      Am      Cm      Bk      Cf      Es      Fm      Md
            Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero , &
          !
            Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero ,  Zero &
]

if ((IA >= lbound(RW,1)) .and. (IA <= ubound(RW,1))) then
  RList = RW(IA)
else
  RList = Zero
  write(u6,'(a)') 'IA out of range in RList.'
  call Abend()
end if

return

end function RList
