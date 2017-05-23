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
* Copyright (C) 1996, Markus P. Fuelscher                              *
*               2002, Roland Lindh                                     *
*               2005, Jesper Wisborg Krogh                             *
************************************************************************
      Subroutine SpoolInp(LuSpool)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Copy all input to a local scratch file                           *
*                                                                      *
*     calling arguments:                                               *
*     LuSpool   : integer                                              *
*               logical unit number of the temporary file              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*     Modified by:  R. Lindh to allow for use of old input file.       *
*                   Tokyo, Japan, June 2002.                           *
*                   J.W. Krogh to make auto.plx create the old input   *
*                   file. Lund, Sweden, October 2005.                  *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Integer (A-Z)
      External Get_ProgName
#include "standard_iounits.fh"
      Character*100 ProgName, Get_ProgName
      Character*128 FileName
      logical exist
*
      Call QEnter('SpoolInp')
*
*     Get the name of the module
*
      ProgName=Get_ProgName()
      Call Upcase(ProgName)
      Call LeftAd(ProgName)
*
      iEnd = 1
 99   If (ProgName(iEnd:iEnd).ne.' ') Then
         iEnd=iEnd+1
         Go To 99
      End If
      iEnd = Min(iEnd-1,5)
      FileName=progname(1:iend)//'INP'

*
*     If Spool is true, then the input in StdIn is just used.
*     Else we'll use the latest input. This is created by
*     auto.plx
*
      LuSpool = 17
      If (Spool) Then
         LuSpool = LuRd
      Else
         call f_inquire('LASTEN',exist) ! customized Last_Energy input
         if(exist) then
           LuSpool = IsFreeUnit(LuSpool)
           call Molcas_Open(LuSpool,'LASTEN')
         else
           call f_inquire(Filename,exist)
           if(exist) then
             LuSpool = IsFreeUnit(LuSpool)
             call Molcas_Open(LuSpool,Filename)
           endif
         endif
      End If
*
       Call QExit('SpoolInp')
      Return
*                                                                      *
************************************************************************
*                                                                      *
*

      End
*
      Subroutine Close_LuSpool(LuSpool)
      Implicit Integer (A-Z)
#include "standard_iounits.fh"
*
      If (.Not. Spool) Then
         Close(LuSpool)
      End If
*
      Return
      End
