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
      Subroutine ClsFls_SCF
************************************************************************
*                                                                      *
*     purpose: Close files after SCF calculations                      *
*                                                                      *
*     called from: SCF                                                 *
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
#include "file.fh"

#include "mxdm.fh"
#include "infscf.fh"
#include "scfwfn.fh"
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*
*---  close two-electron integral file --------------------------------*
      If (.Not.DSCF .And. .Not.DoCholesky) Then
         iRc=-1
         iOpt=0
         Call ClsOrd(iRc,iOpt)
         If (iRc.ne.0) Then
            Write (6,*) 'ClsFls: Error closing ORDINT'
            Call Abend()
         End If
      End If
*
*---  close DNSMAT, dVxcdR, TWOHAM and GRADIENT -----------------------*
      Call DaClos(LuDSt)
      Call DaClos(LuOSt)
      Call DaClos(LuTSt)
      Call DaClos(LuGrd)
*
*---  close 2nd order updatefiles       -------------------------------*
      Call DaClos(LuDGd)
      Call DaClos(Lux)
      Call DaClos(LuDel)
      Call DaClos(Luy)
*
#ifdef _HDF5_
      call mh5_close_file(wfn_fileid)
#endif

      End
