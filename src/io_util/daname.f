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
* Copyright (C) 1991, Per-Olof Widmark                                 *
*               1993,1996,1997, Markus P. Fuelscher                    *
*               1996, Luis Serrano-Andres                              *
************************************************************************
      Subroutine DaName(Lu,String)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Open unit Lu for direct access I/O and link the data stream to   *
*     the logical file name Name.                                      *
*                                                                      *
*     calling arguments:                                               *
*     Lu      : integer, input                                         *
*               logical unit number                                    *
*     Name    : character string, input                                *
*               logical file name                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, IBM Sweden, 1991                                   *
*     M.P. Fuelscher, University of Lund, Sweden, 1993, 1996, 1997     *
*     L. Serrano-Andres, University of Lund, Sweden, 1996              *
*                                                                      *
************************************************************************

      Character*(*) String
      Integer Lu
      Logical mf,wa
      mf=.false.
      wa=.false.
      Call Daname_Main(Lu, String, mf, wa)

      Return
      End
