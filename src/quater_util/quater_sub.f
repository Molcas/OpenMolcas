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
* Copyright (C) Yannick Carissan                                       *
************************************************************************
      subroutine quater_sub(nAtoms,G1,G2,ireturn)
************************************************************
*
*   <DOC>
*     <Name>quater</Name>
*     <Syntax>quater(ireturn)</Syntax>
*     <Arguments>
*       \Argument{ireturn}{return code}{Integer}{out}
*     </Arguments>
*     <Purpose>Driver for quater</Purpose>
*     <Dependencies>quater util and util</Dependencies>
*     <Author>Y. Carissan</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects>none</Side_Effects>
*     <Description>
*        Driver for quater
*     </Description>
*    </DOC>
*
************************************************************
      Implicit none
#include "WrkSpc.fh"
#include "debug.fh"
#include "options.fh"
#include "geoms.fh"
      Real*8 U1(3),U2(3),V1(3),V2(3)
      Real*8 V1best(3),V2best(3)
      Real*8 Q(0:3),Vtrans(3)
      Integer ireturn,nAtoms,XYZ(3)
      Real*8  G1(3*nAtoms),G2(3*nAtoms)

      Data XYZ/1,2,3/
#ifdef _DEBUG_
      debug=.true.
#else
      debug=.false.
#endif
      rotate=.true.
      translate=.true.
      ngeoms=1
      title(1)="ORIGINAL"
      title(2)="GEOHYPER"
      title(3)="bALIGNED"

      call quaterinit()

      nat(1)=nAtoms
      nat(2)=nAtoms
      nat(3)=nAtoms

      Call GetMem("GEO0","Allo","Real",ipgeo(1),3*nAtoms)
      Call GetMem("GEO1","Allo","Real",ipgeo(2),3*nAtoms)
      Call GetMem("GEOF","Allo","Real",ipgeo(3),3*nAtoms)

      call dcopy_(3*nAtoms,G1,1,Work(ipgeo(1)),1)
      call dcopy_(3*nAtoms,G2,1,Work(ipgeo(2)),1)

      Call SetVect(nAtoms,Work(ipgeo(1)),XYZ,V1,V2)
      Call SetVect(nAtoms,Work(ipgeo(2)),XYZ,U1,U2)


      if (debug) then
        Write(6,*) 'Reference axis'
        Call RecPrt("U1",' ',U1,3,1)
        Call RecPrt("U2",' ',U2,3,1)
        Write(6,*) 'New axis'
        Call RecPrt("V1",' ',V1,3,1)
        Call RecPrt("V2",' ',V2,3,1)
      end if

      call QuaterSolve(U1,U2,V1,V2,Q)

      if (debug) then
        Write(6,*) 'Normalized Reference axis'
        Call RecPrt("U1",' ',U1,3,1)
        Call RecPrt("U2",' ',U2,3,1)
        Write(6,*) 'Normalized New axis'
        Call RecPrt("V1",' ',V1,3,1)
        Call RecPrt("V2",' ',V2,3,1)
        call QuaterRotation(Q,U1,V1best)
        call QuaterRotation(Q,U2,V2best)
        Call RecPrt("Best V1",' ',V1best,3,1)
        Call RecPrt("Best V2",' ',V2best,3,1)
      end if

      call dcopy_(3*nAtoms,Work(ipgeo(2)),1,Work(ipgeo(3)),1)

      Call RotateGeoms(Q)
      Call SetVectTrans(nat(1),Work(ipgeo(1)),XYZ,
     &                  nat(3),Work(ipgeo(3)),XYZ,
     &                  Vtrans)
      Call TranslateGeoms(Vtrans)

      if (debug) then
         Call PrintGeom(-1,nat(3),title(3),
     &           Work(ipgeo(3)),geolbl(1,3))
      end if

      call dcopy_(3*nAtoms,Work(ipgeo(3)),1,G2,1)

      Call GetMem("GEO0","Free","Real",ipgeo(1),3*nAtoms)
      Call GetMem("GEO1","Free","Real",ipgeo(2),3*nAtoms)
      Call GetMem("GEOF","Free","Real",ipgeo(3),3*nAtoms)


      ireturn=0

      return
      end
