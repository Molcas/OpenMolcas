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
* Copyright (C) 1991, Markus P. Fuelscher                              *
************************************************************************
      Subroutine InpCtl_Motra(ipOvlp,ipHOne,ipKine,ipCMO)

************************************************************************
*                                                                      *
*     Purpose:                                                         *
*     Read all information required                                    *
*                                                                      *
***** M.P. Fuelscher, University of Lund, Sweden, 1991 *****************
*
      Implicit Real*8 (A-H,O-Z)


#include "motra_global.fh"
#include "files_motra.fh"
#include "trafo_motra.fh"
#include "WrkSpc.fh"
*
      Call qEnter('InpCtl')
*----------------------------------------------------------------------*
*     Read the content of the one electron integral file               *
*----------------------------------------------------------------------*
      Call Rd1Int_Motra(ipOvlp,ipHOne,ipKine)
*----------------------------------------------------------------------*
*     Read auxiliary  input                                            *
*----------------------------------------------------------------------*
      Call RdInp_Motra
*----------------------------------------------------------------------*
*     Read Reaction field and add to one-electron integrals            *
*----------------------------------------------------------------------*
      If ( iRFpert.eq.1 ) Call RdRfld(ipHOne)
*----------------------------------------------------------------------*
*     Read the MO coefficients and occupations                         *
*----------------------------------------------------------------------*
      Call GetMem('CMO','Allo','Real',ipCMO,nTot2)
      Call RdCmo_motra(Work(ipCMO),Work(ipOvlp))
*----------------------------------------------------------------------*
*     Delete orbitals with occupations samller than a given value      *
*----------------------------------------------------------------------*
      If ( iAutoCut.eq.1 ) Call AutoCut
*----------------------------------------------------------------------*
*     Print the input and orbital definitions                          *
*----------------------------------------------------------------------*
      If (iPrint.GE.0) Call PrInp(Work(ipCMO))
*----------------------------------------------------------------------*
*     Normal termination                                               *
*----------------------------------------------------------------------*
      Call qExit('InpCtl')
*
      Return
      End
