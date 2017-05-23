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
* Copyright (C) 1993,1999, Markus P. Fuelscher                         *
*               1993, Per Ake Malmqvist                                *
************************************************************************
      Subroutine UPKI4(nData,nByte,InBuf,OutBuf)
************************************************************************
*                                                                      *
*     purpose: unpack Integers                                         *
*                                                                      *
*    Calling parameters:                                               *
*    nData : number of unpacked integrals                              *
*    nByte : length of pack array in bytes                             *
*    InBuf : contains on input pack integers                           *
*    OutBuf: contains unpacked integers on output                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher and P. A. Malmqvist                               *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history:                                                         *
*     replace various older versions of the algorithm                  *
*     M.P. Fuelscher, Univeristy of Lund, Sweden, 1999                 *
*                                                                      *
************************************************************************
*
      Character*1 InBuf(*)
      Integer OutBuf(nData)
*
*----------------------------------------------------------------------*
*
      iOpt = 1
      Call iunzip(iOpt,nData,nByte,InBuf,OutBuf)
*
*----------------------------------------------------------------------*
*
      Return
      End
