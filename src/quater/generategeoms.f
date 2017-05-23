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
      subroutine GenerateGeoms(Q)
      implicit none
#include "geoms.fh"
#include "debug.fh"
#include "WrkSpc.fh"
#include "options.fh"
      Real*8 Q(0:3)
      Real*8 Vtrans(3)
      integer igeom,iLU
      character*6 fname
      logical isfreeunit

      Do igeom=3,ngeoms+2
        call dcopy_(nat(igeom)*3,Work(ipgeo(2)),1,Work(ipgeo(igeom)),1)
      End Do
      if (rotate) Call RotateGeoms(Q)
      if (translate) then
        Call SetVectTrans(nat(1),Work(ipgeo(1)), XYZ1,
     &                            nat(2),Work(ipgeo(2)), XYZ2,Vtrans)
        Call TranslateGeoms(Vtrans)
       End If

      Do iLU=12,99
        if (isfreeunit(iLU)) goto 999
      End Do
      iLU=-1
      Call SysWarnMsg("GenerateGeoms","No free unit number found",
     &                "Geometries not written to files")
999   continue
      Do igeom=1,ngeoms+2
        if (igeom.lt.100) write(fname,'(a4,i2)') "GEOM",igeom
        if (igeom.lt.10) write(fname,'(a5,i1)') "GEOM0",igeom
        call Molcas_Open (iLU,fname)
        call PrintGeom(iLU,nat(igeom),title(igeom),
     &           Work(ipgeo(igeom)),geolbl(1,igeom))
        close (iLU)
      End Do

      end
