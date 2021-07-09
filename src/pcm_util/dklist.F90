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

function DKList(NumAt)
! Assigns Caillet-Claverie's atomic parameters for dispersion and repulsion.

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: DkList
integer(kind=iwp), intent(in) :: NumAt
real(kind=wp), parameter :: K(110) = [ &
  !  H       He
  1.000_wp,0.594_wp, &
  ! Li       Be        B        C        N       O         F       Ne
    Zero  ,  Zero  ,  Zero  ,1.000_wp,1.180_wp,1.360_wp,1.500_wp,1.209_wp, &
  ! Na       Mg       Al       Si        P       S        Cl       Ar
  1.400_wp,  Zero  ,  Zero  ,  Zero  ,2.100_wp,2.400_wp,2.100_wp,2.127_wp, &
  !  K       Ca
  2.900_wp,  Zero  , &
           ! Sc       Ti        V       Cr       Mn       Fe       Co       Ni       Cu       Zn
             Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  , &
                    ! Ga       Ge       As       Se       Br       Kr
                      Zero  ,  Zero  ,  Zero  ,  Zero  ,2.400_wp,2.566_wp, &
  ! Rb       Sr
    Zero  ,  Zero  , &
           !  Y       Zr       Nb       Mo       Tc       Ru       Rh       Pd       Ag       Cd
             Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  , &
                    ! In       Sn       Sb       Te        I       Xe
                      Zero  ,  Zero  ,  Zero  ,  Zero  ,3.200_wp,  Zero  , &
  ! Cs       Ba
    Zero  ,  Zero  , &
           ! La       Ce       Pr       Nd       Pm       Sm       Eu
           ! Gd       Tb       Dy       Ho       Er       Tm       Yb
             Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  , &
             Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  , &
           ! Lu       Hf       Ta        W       Re       Os       Ir       Pt       Au       Hg
             Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  , &
                    ! Tl       Pb       Bi       Po       At       Rn
                      Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  , &
  ! Fr       Ra
    Zero  ,  Zero  , &
           ! Ac       Th       Pa        U       Np       Pu       Am
           ! Cm       Bk       Cf       Es       Fm       Md
             Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  , &
             Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  , &
           !
             Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero  ,  Zero &
]

DKList = K(NumAt)

return

end function DKList
