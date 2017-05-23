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
* Copyright (C) 2017, Felix Plasser                                    *
************************************************************************

      Subroutine wfa(ireturn)
      USE iso_c_binding
      Implicit None
#include "Molcas.fh"
#include "standard_iounits.fh"
#include "warnings.fh"

      Integer ireturn

#ifdef _WFA_
      Integer iLU, ist, ien
      Character*180 Get_Ln
      External Get_Ln
      Character Line*180
      Character Inp*900
#endif

      ireturn=_RC_NOT_AVAILABLE_

#ifdef _WFA_
!   Read the input file and parse as a string to libwfa
      ist = 1
      Call SpoolInp(iLU)
      Do While (Line(1:4).ne.'END ')
          Line=(Get_ln(iLU))

          ien = ist + len(trim(Line))
          Inp(ist:ien) = trim(Line)
          ist = ien + 1

          Call Normal(Line)
      End Do

      write(6,*) 'Starting wavefunction analysis ...'
      call wfa_driver(ireturn,Inp)
#else
      write(6,*) 'WFA module not installed!'
#endif

      End Subroutine wfa
