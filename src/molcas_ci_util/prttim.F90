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

use timers, only: C_Dress, C_get_Cm, TimeAoMo, TimeCIOpt, TimeDavid, TimeDens, TimeFock, TimeHCSCE, TimeHDiag, TimeHSel, &
                  TimeInput, TimeOrb, TimeOutput, TimePage, TimeRelax, TimeSigma, TimeTotal, TimeTrans, TimeWfn, W_Dress, W_get_Cm
use lucia_data, only: TDENSI, TSIGMA
use splitcas_data, only: DoSplitCAS
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), parameter :: MAX_TIMERS = 35
real(kind=wp) :: F(MAX_TIMERS), T(MAX_TIMERS)

T(:) = Zero
F(:) = Zero
F(MAX_TIMERS) = One

T(MAX_TIMERS) = TimeTotal
T(1) = TimeInput
!T(2) = Zero
!T(3) = Zero
T(4) = T(1)-T(2)-T(3)
T(5) = TimeWfn
T(6) = TimeTrans
T(7) = TimeAoMo
T(8) = TimeFock
T(9) = TimeCIOpt
T(10) = TimeHDiag
T(11) = TimeHSel
T(12) = TimeSigma
T(13) = TimeDens
T(14) = TimeOrb
T(15) = TimeOutput
T(16) = TimeRelax
!T(17) = Zero
T(18) = T(15)-T(16)-T(17)
T(19) = TimeDavid
T(20) = TimePage
T(21) = TimeHCSCE
T(22) = C_Dress
T(23) = W_Dress
T(24) = C_get_Cm
T(25) = W_get_Cm
T(26:31) = TSIGMA(:)
T(32:34) = TDENSI(:)

if (T(MAX_TIMERS) > 1.0e-6_wp) F(:) = T(:)/T(MAX_TIMERS)

write(u6,*)
write(u6,100) 'Timings'
write(u6,100) '-------'
write(u6,*)
write(u6,100) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
write(u6,101) ' ',' ','        time','    fraction'
write(u6,100) '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
write(u6,102) '1) Input section',':',T(1),F(1)
write(u6,102) '   - Input processing',':',T(4),F(4)
!write(u6,102) '   - Create GUGA tables',':',T(2),F(2)
!write(u6,102) '   - Create determinant tables',':',T(3),F(3)
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
!write(u6,102) '   - Create/update the file RUNFILE',':',T(17),F(17)
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
