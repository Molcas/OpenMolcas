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
* Copyright (C) Per-Olof Widmark                                       *
************************************************************************
************************************************************************
*                                                                      *
*======================================================================*
*                                                                      *
* Author: Per-Olof Widmark                                             *
*         IBM Sweden                                                   *
*                                                                      *
************************************************************************
      Subroutine zlcase(line)
      Implicit Integer (a-z)
      Character*(*) line
      Do 100 k=1,len(line)
         Index=iChar(line(k:k))
*        If(97.le.Index .and. Index.le.122) line(k:k)=Char(Index-32)
         If(65.le.Index .and. Index.le.90) line(k:k)=Char(Index+32)
100   Continue
      Return
      End
