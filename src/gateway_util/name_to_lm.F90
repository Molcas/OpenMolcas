!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2017, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  Name_to_lm
!
!> @brief
!>   Get \f$ l \f$ and \f$ m \f$ numbers from a basis function name.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Given a basis function name as a character string, extract the \f$ l
!> \f$ and \f$ m \f$ quantum numbers.
!>
!> Spherical harmonics functions are expected with the format `nnLmms`,
!> where `nn` is the shell number, `L` is a letter denoting angular
!> momentum, `mm` is the absolute value of the \f$ m \f$ number, and `s`
!> is the sign of \f$ m \f$.
!>
!> Cartesian functions are expected with the format `Lxxyyzz`, where `L`
!> is a letter denoting angular momentum, and `xx`, `yy`, `zz` are the
!> powers of \f$ x \f$, \f$ y \f$ and \f$ z \f$ (\f$ m_x,m_y,m_z \f$).
!> in this case, the numbers returned are \f$ -l \f$ and
!> \f$ ((m_y+m_z)^2+m_z-m_y)/2-m_x \f$.
!>
!> @param[in]  BName Basis function name
!> @param[out] l     \f$ l \f$ number (\f$ -l \f$ for Cartesians)
!> @param[out] m     \f$ m \f$ number (see details for Cartesians)
!***********************************************************************

subroutine Name_to_lm(BName,l,m)

use define_af, only: AngTp
use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: BName
integer(kind=iwp), intent(out) :: l, m
character :: Letter
integer(kind=iwp) :: i, lx, ly, lz

Letter = BName(3:3)
call LoCase(Letter)
l = 0
m = 0
if (Letter == 's') then
  ! Default is s
  return
else if (Letter == 'p') then
  ! p usually appear as px, py, pz, except when they are contaminants
  l = 1
  if (BName(4:4) /= '0') then
    Letter = BName(4:4)
    call LoCase(Letter)
    if (Letter == 'x') then
      m = 1
    else if (Letter == 'y') then
      m = -1
    else if (Letter == 'z') then
      m = 0
    end if
    return
  end if
end if
! Parse the label for other cases
l = -1
! Find if there is an angular label
do i=sum(lbound(AngTp)),sum(ubound(AngTp))
  if (Letter == AngTp(i)) then
    l = i
    exit
  end if
end do
if (l >= 0) then
  ! If a label is found it is a spherical shell, just read m
  read(BName(4:5),*) m
  if (BName(6:6) == '-') m = -m
else
  ! If no label, this is a Cartesian shell, return -l and some convention for m.
  ! We use m=T(ly+lz)-(lx+ly), where T(n) is the nth triangular number: n*(n+1)/2).
  ! From here, ly+lz can be recovered as the (int) triangular root of m+l: (sqrt(8*(m+l)+1)-1)/2,
  ! lz is m+l-T(ly+lz) and lx is l-(ly+lz). This has the property that all
  ! possible combinations of lx,ly,lz are encoded in -l plus a number from -l to l*(l+1)/2,
  ! in descending order with priority lx>ly>lz: (3,0,0), (2,1,0), (2,0,1), (1,2,0),
  ! (1,1,1), (1,0,2), (0,3,0), (0,2,1), (0,1,2), (0,0,3)
  read(BName(2:3),*) lx
  read(BName(4:5),*) ly
  read(BName(6:7),*) lz
  l = -lx-ly-lz
  m = (ly+lz)*(ly+lz+1)/2-(lx+ly)
end if

return

end subroutine Name_to_lm
