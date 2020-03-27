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
************************************************************************
*                                                                      *
*  Integer Function NbfShl    returns # of bf for shell,symmetry       *
*                                                                      *
************************************************************************
c----------------------------------------------------------------------
      Integer Function nbfshl(iSkal,irp)
      use iSD_data
c----------------------------------------------------------------------
      Implicit Real*8 (A-H,O-Z)
*
*  returns number of basis functions for given shell and symmetry
*
#include "itmax.fh"
#include "info.fh"
*
      nbfshl=0
      iShell = iSD(11,iSkal)
      iCmp   = iSD( 2,iSkal)
*     loop over components of shell...
      Do i=1, iCmp
        If (IAND(IrrCmp(IndS(iShell)+i),2**irp).ne.0) nbfshl =
     &           nbfshl + iSD(3,iSkal)
      End Do
      return
      End
