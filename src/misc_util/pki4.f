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
      Subroutine PKI4(nData,nByte,InBuf,OutBuf)
************************************************************************
*                                                                      *
*     purpose: pack Integers                                           *
*                                                                      *
*    Calling parameters:                                               *
*    nData : number of unpacked integrals                              *
*    nByte : length of pack array in bytes                             *
*    InBuf : contains on input unpack integers                         *
*    OutBuf: contains packed integers on output                        *
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
*     M.P. Fuelscher, University of Lund, Sweden, 1999                 *
*                                                                      *
************************************************************************
*
      Integer InBuf(nData)
      Integer OutBuf(*)
      Call cPKI4(OutBuf)
*
*     This is to allow type punning without an explicit interface
      Contains
      Subroutine cPKI4(OutBuf)
      Use Iso_C_Binding
      Integer, Target :: OutBuf(*)
      Character, Pointer :: cOutBuf(:)
*
*----------------------------------------------------------------------*
*
      iOpt = 1
      Call C_F_Pointer(C_Loc(OutBuf(1)),cOutBuf,[1])
      Call izip(iOpt,nData,nByte,InBuf,cOutBuf)
      Nullify(cOutBuf)
*
*----------------------------------------------------------------------*
*
      Return
      End Subroutine cPKI4
*
      End
