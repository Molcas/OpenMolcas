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
      Dimension Dum(1)
*-----------------------------------------------------------------------*
* Enter.                                                                *
*-----------------------------------------------------------------------*
*-----------------------------------------------------------------------*
* Read!                                                                 *
*-----------------------------------------------------------------------*
      iDisk=0
      If(nSkipp.ne.0.and.iPrint.ge.4) then      !If we are to skip
        Write(6,*)' Reading from configuration ',nskipp,'.' !something.
      Endif
      Do 11, j=1,nSkipp+1
        If(j.ne.1.and.iRead.ne.9) then
          Call dDaFile(9,2,Dum,1,iDisk) !Etot
          Call dDaFile(9,2,Dum,1,iDisk) !Ract
          Call dDaFile(9,2,Dum,1,iDisk) !GamOld
          Call dDaFile(9,2,Dum,1,iDisk) !Gam
          Call dDaFile(9,2,Dum,1,iDisk) !ESub
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
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_logical(quantum)
      End
