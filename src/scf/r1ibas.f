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
      Subroutine R1IBas
************************************************************************
*                                                                      *
*     purpose: Read basis set informations.                            *
*                                                                      *
*     called from: ReadIn                                              *
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

#include "mxdm.fh"
#include "infscf.fh"

#include "WrkSpc.fh"
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
      Call qEnter('R1IBas')
*
*---- read file header
      Call Get_cArray('Seward Title',Header,144)
*---- read number of symm. species
      Call Get_iScalar('nSym',nSym)
*---- read number of basis functions per symmetry species
      Call Get_iArray('nBas',nBas,nSym)
*---- read basis function labels
      nBas_tot=0
      Do iSym = 1, nSym
         nBas_tot=nBas_tot+nBas(iSym)
      End Do
      Call Get_cArray('Unique Basis Names',Name,(LENIN8)*nBas_tot)
*---- read number of atoms
      Call Get_iScalar('Unique atoms',nAtoms)
*---- read nuclear potential
*     Call get_dScalar('PotNuc',PotNuc)
      Call Peek_dScalar('PotNuc',PotNuc)
*
*---- Compute lengths of matrices
      lthBas = 0
      Do iSym = 1, nSym
         lthBas = lthBas + nBas(iSym)
      End Do
*
*---- Define atom and type
      Do i = 1, lthBas
         Atom(i) = Name(i)(1:LENIN)
         Type(i) = Name(i)(LENIN1:LENIN8)
      End Do
*
      Call qExit('R1IBas')
      Return
      End
