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
*  Subroutine CoordShl        returns coordinates of a given shell     *
*                                                                      *
************************************************************************
      SubRoutine CoordShl(xi,yi,zi,iskal)
c----------------------------------------------------------------------
      use iSD_data
      Implicit Real*8 (A-H,O-Z)
#include "shinf.fh"
#include "nsd.fh"
#include "itmax.fh"
#include "setup.fh"
#include "WrkSpc.fh"
*
************************************************************************
*                                                                      *
* returns coordinates of a shell                                       *
*                                                                      *
************************************************************************
*
      jpcn=iSD(8,iSkal)
      xi=work(jpcn)
      yi=work(jpcn+1)
      zi=work(jpcn+2)
      return
      end
