************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Rd1Int_FFPT
*
************************************************************************
*                                                                      *
*     Objective: Read the header of the one-electron integral file     *
*                Extract also symmetry and basis set information.      *
*                In addition read the overlap, the nuclear attraction  *
*                and kinetic integrals.                                *
*                                                                      *
************************************************************************
*
      Implicit Real*8 ( A-H,O-Z )
*
#include "input.fh"
#include "WrkSpc.fh"
*
*
*----------------------------------------------------------------------*
*                                                                      *
*     Start procedure                                                  *
*     Read one-electron integral file header etc.                      *
*                                                                      *
*----------------------------------------------------------------------*
*
      iOpt=0
      iComp=1
      iSyLbl=0
*
      Call Get_cArray('Seward Title',Header,144)
      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)
      Call Get_iScalar('Unique atoms',nAtoms)
      Call Get_dArray('Unique Coordinates',Coor,3*nAtoms)
*
      Return
      End
