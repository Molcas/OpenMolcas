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
* This routine prints the trailer part to the postscript figure of     *
* occupation numbers.                                                  *
*                                                                      *
*======================================================================*
*                                                                      *
* Author: Per-Olof Widmark                                             *
*         IBM Sweden                                                   *
*                                                                      *
************************************************************************
      Subroutine FigCls(lu)
      Write(lu,'(a)') '%---'
      Write(lu,'(a)') 'grestore'
      Write(lu,'(a)') 'showpage'
      Write(lu,'(a)') '%.AutoFeed'
      Return
      End
