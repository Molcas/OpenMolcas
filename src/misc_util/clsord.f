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
      Subroutine ClsOrd(rc,option)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Close the two-electron integral file.                            *
*                                                                      *
*     input:                                                           *
*     option : Switch to set options                                   *
*              (not used at present)                                   *
*                                                                      *
*     output:                                                          *
*     rc     : Return code.                                            *
*              A value of 0 (zero) is returned upon successful         *
*              completion of the request. A nonzero value indi-        *
*              cates an error.                                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. P. Fuelscher                                                  *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Integer (A-Z)
*
#include "TwoRc.fh"
#include "TwoFlags.fh"
#include "Molcas.fh"
#include "TwoDat.fh"
*----------------------------------------------------------------------*
*     Start procedure:                                                 *
*----------------------------------------------------------------------*
*     Call qEnter('ClsOrd')
      rc=rc0000
*----------------------------------------------------------------------*
*     Check the file status                                            *
*----------------------------------------------------------------------*
      If( AuxTwo(isStat).ne.1 ) Then
         rc=rcCL01
      Call SysAbendMsg('ClsOrd',
     *  'The ORDINT file has not been opened',' ')
      End If
*----------------------------------------------------------------------*
*     Reset error code,open flag and unit number. Close file.          *
*----------------------------------------------------------------------*
      LuOrd=AuxTwo(isUnit)
      iDisk=0
      Call iDaFile(LuOrd,1,TocTwo,lTocTwo,iDisk)
      Call DaClos(LuOrd)
      AuxTwo(isUnit) = iNoNum
      AuxTwo(isStat) = iNoNum
      AuxTwo(isDaDa) = iNoNum
*
      If (RAMD) Then
         Call GetMem('RAMD','Free','Real',ip_RAMD,RAMD_size)
         RAMD=.False.
      End If
*----------------------------------------------------------------------*
*     Terminate procedure                                              *
*----------------------------------------------------------------------*
*     Call qExit('ClsOrd')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real(option)
      End
