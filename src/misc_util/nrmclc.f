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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
************************************************************************
      SubRoutine NrmClc(Vec,lth,SubNam,MatNam)
************************************************************************
*                                                                      *
*     purpose: compute and print out norms for debuging purposes       *
*                                                                      *
*     input:                                                           *
*       Vec     : vector whose norm is going to be computed (lth)      *
*       SumNam  : name of a subroutine it is called from               *
*       MatNam  : name of a matrix                                     *
*       DeBug   : T - print out norm                                   *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
*
      Implicit Real*8 (a-h,o-z)
*
      Real*8 Vec(lth)
      Character SubNam*(*), MatNam*(*)
*
      R = DDot_(lth,Vec,1,Vec,1)
      Q = DDot_(lth,[1.0D0],0,Vec,1)
      S = 0.0D0
      Do i = 1, lth
         S = S + Vec(i)*DBLE(i)
      End Do
      Write(6,'(5A,3E24.16)')' Norm of ',MatNam,' in ',SubNam,' = ',
     &                        R,Q,S
*
      Return
      End
