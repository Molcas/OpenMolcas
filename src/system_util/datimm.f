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
* Copyright (C) 1993, Markus P. Fuelscher                              *
************************************************************************
      Subroutine Datimm
************************************************************************
*                                                                      *
*     Print the current date and time.                                 *
*                                                                      *
*     Note:                                                            *
*     The VS-FORTRAN subroutines, Datim and Datimx, are replaced       *
*     by a similar procedure written in C.                             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Character*72 Line
*----------------------------------------------------------------------*
*     start                                                            *
*----------------------------------------------------------------------*
      Line='                                                      '

      Call Datimx(Line)
      Write(6,'(5A)')
     &  '* Started ',Line(1:10),Line(20:24),
     &  ' at ',Line(12:19)
*----------------------------------------------------------------------*
*     exit                                                             *
*----------------------------------------------------------------------*
      Return
      End
