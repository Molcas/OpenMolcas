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
      Subroutine ClsOne(rc,Option)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Close the one-electron integral file.                            *
*                                                                      *
*     input:                                                           *
*     Option : Switch to set options                                   *
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
*     M. P. Fuelscher, University of Lund, Sweden, 1993                *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
*
      Implicit Integer (A-Z)
*
#include "OneRc.fh"
#include "OneFlags.fh"

#include "OneDat.fh"
*---------------------------------------------------------------------*
*     Start procedure                                                 *
*---------------------------------------------------------------------*
*     Call qEnter('ClsOne')
      rc=rc0000
      LuOne=AuxOne(pLu)
*----------------------------------------------------------------------*
*     Check the file status                                            *
*----------------------------------------------------------------------*
      If( AuxOne(pOpen).ne.1 ) Then
         rc=rcCL01
      Call SysAbendMsg('ClsOne',
     *  'The ONEINT file has not been opened',' ')
      End If
      AuxOne(pOpen)=0
*----------------------------------------------------------------------*
*     Dump the TOC upon request                                        *
*----------------------------------------------------------------------*
      If ( iAnd(Option,1024).ne.0 ) Call DmpOne
*----------------------------------------------------------------------*
*     Reset error code,open flag and unit number. Close file.          *
*----------------------------------------------------------------------*
      Call DaClos(LuOne)
      Call iCopy(lAux,[NaN],0,AuxOne,1)
      Call iCopy(lToc,[NaN],0,TocOne,1)
*----------------------------------------------------------------------*
*     Terminate procedure                                              *
*----------------------------------------------------------------------*
*     Call qExit('ClsOne')
      Return
      End
