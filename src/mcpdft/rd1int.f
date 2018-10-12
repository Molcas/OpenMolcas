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
      Subroutine Rd1Int_m
************************************************************************
*                                                                      *
*     Read header and matrices from the one-electron integral file     *
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
      Implicit Real*8 (A-H,O-Z)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
      Call qEnter('Rd1Int')
*---  read file header  -----------------------------------------------*
      Call Get_cArray('Seward Title',Header,144)
*---  read number of symm. species ------------------------------------*
      Call Get_iScalar('nSym',nSym)
*---  read number of basis functions per symmetry species -------------*
      Call Get_iArray('nBas',nBas,nSym)
*---  read nuclear potential ------------------------------------------*
* Get POTNUC from the runfile, where it was set by seward.
* ( Do not trust reading it from JOBIPH).
      Call Get_dScalar('potNuc',PotNuc)
*---  read basis function labels --------------------------------------*
      nBas_tot=0
      Do iSym = 1, nSym
         nBas_tot=nBas_tot+nBas(iSym)
      End Do
      Call Get_cArray('Unique Basis Names',Name,(LENIN8)*nBas_tot)
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
      Call qExit('Rd1Int')
      Return
      End
