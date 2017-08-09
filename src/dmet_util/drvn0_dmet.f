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
* Copyright (C) 1991,1995, Roland Lindh                                *
************************************************************************
      SubRoutine DrvN0_DMET(ireturn)
************************************************************************
*                                                                      *
* Object: to compute the nuclear contibutions to the nuclear potential *
*                                                                      *
* Called from: Input                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             October '91                                              *
*             Modified for various other contributions May 95', RL     *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
      Real*8 A(3), B(3), RB(3)
      Integer iDCRR(0:7), jCoSet(8,8), iStb(0:7), jStb(0:7)
      Logical EQ, NoLoop
*
*     Statement function for Cartesian index
*
*

      iRout = 33
      iPrint = nPrint(iRout)
      Call qEnter('DrvN0')
      r12_Min=0.0D0
*
*
*VB for H2      PotNuc = 0.71510434d0
      PotNuc = 0.000d0

*        PotNuc=PotNuc+PXX
*

         If (Show) Then
            Write (6,*)
            Write (6,*)
            Write (6,'(11X,A,F16.8,A)')
     &     ' Total Nuclear Potential Energy        ',PotNuc,' au'
            Write (6,*)
         End If
*
      Call Put_dScalar('PotNuc',PotNuc)
*
      Call qExit('DrvN0')
      Return
      End
