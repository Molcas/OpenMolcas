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
! Copyright (C) 2000, Markus P. Fuelscher                              *
!***********************************************************************

subroutine PrtTim()
!***********************************************************************
!                                                                      *
!     print out timings for the various sections of the program        *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 2000                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use splitcas_data, only: DoSplitCAS
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "timers.fh"
integer(kind=iwp), parameter :: MAX_TIMERS = 40
real(kind=wp) :: F(MAX_TIMERS), T(MAX_TIMERS)
integer(kind=iwp) :: i

T(:) = Zero
F(:) = Zero

T(MAX_TIMERS) = Ebel_3
T(15) = Ebel_3-Ebel_2
T(5) = Ebel_2-Ebel_1
T(1) = Ebel_1
T(4) = Eterna_3-Eterna_2
T(3) = Eterna_2-Eterna_1
T(2) = T(1)-T(3)-T(4)
T(6) = Fortis_3
T(7) = Candino_3
T(8) = Piaget_3
T(9) = Zenith_3
T(10) = Tissot_3
T(11) = Omega_3
T(12) = Rolex_3
T(13) = Rado_3
T(14) = Gucci_3
T(16) = Oris_2
T(17) = Movado_2
T(18) = T(15)-T(16)-T(17)
T(19) = Alfex_3
T(20) = WTC_3
T(21) = Longines_3
T(22) = C_Dress_3
T(23) = W_Dress_3
T(24) = C_get_Cm3
T(25) = W_get_Cm3
T(26) = TSIGMA(1)
T(27) = TSIGMA(2)
T(28) = TSIGMA(3)
T(29) = TSIGMA(4)
T(30) = TSIGMA(5)
T(31) = TSIGMA(6)
T(32) = TDENSI(1)
T(33) = TDENSI(2)
T(34) = TDENSI(3)

do i=1,MAX_TIMERS-1
  if (1000.0_wp*T(i) > One) then
    F(i) = T(i)/T(MAX_TIMERS)
  else
    F(i) = Zero
  end if
end do

write(u6,*)
write(u6,100) 'Timings'
write(u6,100) '-------'
write(u6,*)
write(u6,100) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
write(u6,101) ' ',' ','        time','    fraction'
write(u6,100) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
write(u6,102) '1) Input section',':',T(1),F(1)
write(u6,102) '   - Input processing',':',T(2),F(2)
write(u6,102) '   - Create GUGA tables',':',T(3),F(3)
write(u6,102) '   - Create determinant tables',':',T(4),F(4)
write(u6,102) '2) Wave function optimization',':',T(5),F(5)
write(u6,102) '   - transformation section',':',T(6),F(6)
write(u6,102) '     . AO=>MO integral transformation',':',T(7),F(7)
write(u6,102) '     . Fock-matrix generation',':',T(8),F(8)
write(u6,102) '   - CI optimization',':',T(9),F(9)
write(u6,102) '     . construct Hdiag',':',T(10),F(10)

if (.not. DoSplitCAS) then ! (GLMJ)
  write(u6,102) '     . construct Hsel',':',T(11),F(11)
  write(u6,102) '     . Davidson diagonalization',':',T(19),F(19)

  write(u6,102) '       .. sigma vector generation',':',T(12),F(12)
  write(u6,102) '          |-> aa/bb 1-electron   ',':',T(26),F(26)
  write(u6,102) '          |-> aa/bb 2-electron   ',':',T(27),F(27)
  write(u6,102) '          \-> alpha-beta         ',':',T(28),F(28)
  write(u6,102) '              |-> C prefetch     ',':',T(29),F(29)
  write(u6,102) '              |-> matrix multiply',':',T(30),F(30)
  write(u6,102) '              \-> S scatter      ',':',T(31),F(31)

  write(u6,102) '       .. HCSCE',':',T(21),F(21)
  write(u6,102) '       .. page_in/page_out',':',T(20),F(20)
else
  write(u6,102) '     . U_AA diagonalization',':',T(22),F(22)
  write(u6,102) '     . compute Cm coeff',':',T(24),F(24)
end if

write(u6,102) '     . density matrix generation',':',T(13),F(13)
write(u6,102) '          |-> aa/bb 1-electron  ',':',T(32),F(32)
write(u6,102) '          |-> aa/bb 2-electron  ',':',T(33),F(33)
write(u6,102) '          \-> alpha-beta        ',':',T(34),F(34)

write(u6,102) '   - orbital optimization',':',T(14),F(14)
write(u6,102) '3) Output section',':',T(15),F(15)
write(u6,102) '   - Create/update the file RELAX',':',T(16),F(16)
write(u6,102) '   - Create/update the file RUNFILE',':',T(17),F(17)
write(u6,102) '   - Create/update the file JOBIPH',':',T(18),F(18)
write(u6,*)
write(u6,100) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
write(u6,102) '   Total',':',T(MAX_TIMERS),F(MAX_TIMERS)
write(u6,100) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
write(u6,*)

100 format(2X,A)
101 format(2X,A,T44,A,A,A)
102 format(2X,A,T44,A,2F12.2)

end subroutine PrtTim
