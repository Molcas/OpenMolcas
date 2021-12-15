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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      SubRoutine MltNuc(CoOP,Chrg,Coor,nAtm,rNucMm,ir,nComp)
************************************************************************
*                                                                      *
* Object: to compute the multipole moments for the nuclei.             *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
      Real*8 Chrg(nAtm), Coor(3,nAtm), rNucMm((ir+1)*(ir+2)/2), CoOp(3)
*
      iRout = 124
      iPrint = nPrint(iRout)
      If (iPrint.ge.99) Then
         Call RecPrt(' In MltNuc:Coor',' ',Coor,3,nAtm)
         Call RecPrt(' In MltNuc:Chrg',' ',Chrg,nAtm,1)
         Call RecPrt(' In MltNuc:CoOp',' ',CoOp,1,3)
      End If
*
*     Compute the nuclear contribution to the multipole moments
*
      ip = 0
      Do 71 ix = ir, 0, -1
         Do 72 iy = ir-ix, 0, -1
            ip = ip + 1
            iz = ir-ix-iy
            temp = Zero
*           Write (*,*) ' ix,iy,iz=',ix,iy,iz
            Do 73 iAtom = 1, nAtm
               If (ix.eq.0) Then
                  CCoMx=One
               Else
                  CCoMx=(Coor(1,iAtom)-CoOp(1))**ix
               End If
               If (iy.eq.0) Then
                  CCoMy=One
               Else
                  CCoMy=(Coor(2,iAtom)-CoOp(2))**iy
               End If
               If (iz.eq.0) Then
                  CCoMz=One
               Else
                  CCoMz=(Coor(3,iAtom)-CoOp(3))**iz
               End If
*              Write (*,*) CCoMx, CCoMy, CCoMz, temp
               temp = temp + Chrg(iAtom) * CCoMx * CCoMy * CCoMz
 73         Continue
            rNucMm(ip) = temp
 72      Continue
 71   Continue
*
      If (iPrint.ge.99) Call RecPrt(' Nuclear Multipole Moments',
     &                              ' ',rNucMm,ip,1)
      Return
      If (.False.) Call Unused_integer(nComp)
      End
