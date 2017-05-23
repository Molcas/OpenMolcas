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
      Subroutine RelEny(ERelMV,ERelDC,D,VMs,Dar,lth)
************************************************************************
*                                                                      *
*     purpose: Compute relativistic corrections to the energy          *
*                                                                      *
*     input:                                                           *
*       D       : density matrix of length lth                         *
*       VMs     : mass velocity integrals of length lth                *
*       Dar     : Darvin integrals of length lth                       *
*                                                                      *
*     output:                                                          *
*       ERelMV  : mass velocity correction                             *
*       ERelDC  : Darvin correction                                    *
*                                                                      *
*     called from: PrFin                                               *
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
      Dimension D(lth),VMs(lth),Dar(lth)
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*
      ERelMV = DDot_(lth,D,1,VMs,1)
      ERelDC = DDot_(lth,D,1,Dar,1)
*
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
