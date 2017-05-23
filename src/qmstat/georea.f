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
*----------------------------------------------------------------------*
*
*----------------------------------------------------------------------*
      Subroutine GeoRea(nskipp,quantum)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"

      Logical quantum
*-----------------------------------------------------------------------*
* Enter.                                                                *
*-----------------------------------------------------------------------*
      Call Qenter('GeoRea')
*-----------------------------------------------------------------------*
* Read!                                                                 *
*-----------------------------------------------------------------------*
      iDisk=0
      Ind=nCent*nPart
      IndMa=nPol*nPart
      If(nSkipp.ne.0.and.iPrint.ge.4) then      !If we are to skip
        Write(6,*)' Reading from configuration ',nskipp,'.' !something.
      Endif
      Do 11, j=1,nSkipp+1
        If(j.ne.1.and.iRead.ne.9) then
          Call dDaFile(9,2,Etot,1,iDisk)
          Call dDaFile(9,2,Ract,1,iDisk)
          Call dDaFile(9,2,GamOld,1,iDisk)
          Call dDaFile(9,2,Gam,1,iDisk)
          Call dDaFile(9,2,ESub,1,iDisk)
        Endif
        If(iRead.eq.9) then  !If this is a sampfile we do not care about
                          !the induced dipoles, so we just read them to
                          !get rid of them.
c          Do 12, i=1+nPol,IndMa
c
c12        Continue
        Endif
11    Continue
*-----------------------------------------------------------------------*
* Exit.                                                                 *
*-----------------------------------------------------------------------*
      Call Qexit('GeoRea')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_logical(quantum)
      End
