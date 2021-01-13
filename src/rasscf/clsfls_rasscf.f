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
      Subroutine ClsFls_RASSCF
************************************************************************
*                                                                      *
*     Close files.                                                     *
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
#ifdef _HDF5_
      use mh5, only: mh5_close_file
#endif
      Implicit Real*8 (A-H,O-Z)
      Logical DoCholesky
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='ClsFls  ')

#include "davctl.fh"
#include "qnctl.fh"
#include "raswfn.fh"
*----------------------------------------------------------------------*
*     Start                                                            *
*-------------------------------------- -------------------------------*
C Local print level (if any)
*---  close the JOBOLD file -------------------------------------------*
      If(JOBOLD.gt.0.and.JOBOLD.ne.JOBIPH) Then
        Call DaClos(JOBOLD)
        JOBOLD=-1
      Else If (JOBOLD.gt.0) Then
        JOBOLD=-1
      End If
*---  close the JOBIPH file -------------------------------------------*
      If(JOBIPH.gt.0) Then
        Call DaClos(JOBIPH)
        JOBIPH=-1
      End If
#ifdef _HDF5_
      call mh5_close_file(wfn_fileid)
#endif
*---  close the ORDINT file -------------------------------------------*
      CALL DecideOnCholesky(DoCholesky)
       If (.not.DoCholesky) then
         iRc=-1
         iOpt=0
         Call ClsOrd(iRc,iOpt)
         If ( iRc.ne.0 ) Then
           Call WarningMessage(1,'Failed to close the ORDINT file.')
         End If
       End If
*---  close the file carrying the transformed two-electron integrals --*
      Call DaClos(LUINTM)
*---  close the DAVID file carrying temporary CI and sigma vectros ----*
      Call DaClos(LUDAVID)
*---  open the file carrying the hessian update vectors ---------------*
      Call DaClos(LuQune)
*---  close the file for storage of informations on CI-iterations
      Close(ITERFILE)
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
      Return
      End
