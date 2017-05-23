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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      SubRoutine Exp_1(Vector,n1,n2,Array,Fact)
************************************************************************
*                                                                      *
* Object: expand an array.                                             *
*                                                                      *
* Called from: Rys2Dg                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             October '91                                              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Real*8 Vector(n1,n2), Array(n1)
*
*     Call qEnter('Exp_1')
*
      Do 10 i2 = 1, n2
         Do 20 i1 = 1, n1
            Vector(i1,i2) = Array(i1)*Fact
 20      Continue
 10   Continue
*
*     Call qExit('Exp_1')
      Return
      End
