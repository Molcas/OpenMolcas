************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2000, Markus P. Fuelscher                              *
************************************************************************
      Subroutine PrtTim_m

************************************************************************
*                                                                      *
*     print out timings for the various sections of the program        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 2000                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************

      Implicit Real*8 (A-H,O-Z)

#include "timers.fh"
#include "output_ras.fh"
#include "splitcas.fh"
      Parameter (ROUTINE='XXXXXXXX')
      Integer, Parameter :: MAX_TIMERS = 40

      Dimension T(MAX_TIMERS),F(MAX_TIMERS)
*
      T = 0.0d0
      F = 0.0d0
*
      T(MAX_TIMERS) = Ebel_3
      T(15) = Ebel_3 - Ebel_2
      T( 5) = Ebel_2 - Ebel_1
      T( 1) = Ebel_1
      T( 4) = Eterna_3 - Eterna_2
      T( 3) = Eterna_2 - Eterna_1
      T( 2) = T(1) - T(3) - T(4)
      T( 6) = Fortis_3
      T( 7) = Candino_3
      T( 8) = Piaget_3
      T( 9) = Zenith_3
      T(10) = Tissot_3
      T(11) = Omega_3
      T(12) = Rolex_3
      T(13) = Rado_3
      T(14) = Gucci_3
      T(16) = Oris_2
      T(17) = Movado_2
      T(18) = T(15) - T(16) - T(17)
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
*
      Do i = 1,(MAX_TIMERS - 1)
        If ( 1000.0d0*T(i).gt.1.0d0 ) then
          F(i) = T(i)/T(MAX_TIMERS)
        Else
          F(i) = 0.0d0
        End If
      End Do
*
      Write(LF,*)
      Write(LF,'(2X,A)') 'Timings'
      Write(LF,'(2X,A)') '-------'
      Write(LF,*)
      Write(LF,'(2X,A)')
     &      '- - - - - - - - - - - - - - - - -'//
     &      ' - - - - - - - - - - - - - - - - -'
      Write(LF,'(2X,A,T44,A,A,A)')
     &      ' ',' ','        time','    fraction'
      Write(LF,'(2X,A)')
     &      '- - - - - - - - - - - - - - - - -'//
     &      ' - - - - - - - - - - - - - - - - -'
      Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '1) Input section',':',T(1),F(1)
      Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '   - Input processing',':',T(2),F(2)
      Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '   - Create GUGA tables',':',T(3),F(3)
      Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '   - Create determinant tables',':',T(4),F(4)
      Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '2) Wave function optimization',':',T(5),F(5)
      Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '   - transformation section',':',T(6),F(6)
      Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '     . AO=>MO integral transformation',':',T(7),F(7)
      Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '     . Fock-matrix generation',':',T(8),F(8)
      Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '   - CI optimization',':',T(9),F(9)
      Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '     . construct Hdiag',':',T(10),F(10)

      If (.not.DoSplitCAS) then ! (GLMJ)
        Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '     . construct Hsel',':',T(11),F(11)
        Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '     . Davidson diagonalization',':',T(19),F(19)

        Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '       .. sigma vector generation',':',T(12),F(12)
        Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '          |-> aa/bb 1-electron   ',':',T(26),F(26)
        Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '          |-> aa/bb 2-electron   ',':',T(27),F(27)
        Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '          \-> alpha-beta         ',':',T(28),F(28)
        Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '              |-> C prefetch     ',':',T(29),F(29)
        Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '              |-> matrix multiply',':',T(30),F(30)
        Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '              \-> S scatter      ',':',T(31),F(31)

        Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '       .. HCSCE',':',T(21),F(21)
        Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '       .. page_in/page_out',':',T(20),F(20)
      else
        Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '     . U_AA diagonalization',':',T(22),F(22)
        Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '     . compute Cm coeff',':',T(24),F(24)
      End if

      Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '     . density matrix generation',':',T(13),F(13)
      Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '          |-> aa/bb 1-electron  ',':',T(32),F(32)
      Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '          |-> aa/bb 2-electron  ',':',T(33),F(33)
      Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '          \-> alpha-beta        ',':',T(34),F(34)

      Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '   - orbital optimization',':',T(14),F(14)
      Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '3) Output section',':',T(15),F(15)
      Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '   - Create/update the file RELAX',':',T(16),F(16)
      Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '   - Create/update the file RUNFILE',':',T(17),F(17)
      Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '   - Create/update the file JOBIPH',':',T(18),F(18)
      Write(LF,*)
      Write(LF,'(2X,A)')
     &      '- - - - - - - - - - - - - - - - -'//
     &      ' - - - - - - - - - - - - - - - - -'
      Write(LF,'(2X,A,T44,A,2F12.2)')
     &      '   Total',':',T(MAX_TIMERS),F(MAX_TIMERS)
      Write(LF,'(2X,A)')
     &      '- - - - - - - - - - - - - - - - -'//
     &      ' - - - - - - - - - - - - - - - - -'
      Write(LF,*)
*
      Return
      End
