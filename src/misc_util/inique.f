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
      Subroutine IniQue
************************************************************************
*                                                                      *
*     Initialize the 'pseudo dynamic' memory controller                *
*                                                                      *
*     calling arguments: none                                          *
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
      Implicit Integer (A-Z)

#include "SysCtl.fh"
*----------------------------------------------------------------------*
*     Initialize the Common / MemCtl / the first time it is referenced *
*----------------------------------------------------------------------*
      Do 10 i=1,lQueCtl
         QueCtl(i)=0
10    End Do
      QueCtl(ipStat)   = ON
      QueCtl(ipTrace)  = OFF
      QueCtl(ipQuery)  = OFF
      QueCtl(ipSysOut) = 6
      QueCtl(iplTbl)   = 0
      QueCtl(iplStk)   = 0
*----------------------------------------------------------------------*
*     exit                                                             *
*----------------------------------------------------------------------*
      Return
      End
