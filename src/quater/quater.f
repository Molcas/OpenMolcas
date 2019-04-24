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
*  quater
*
*> @brief
*>   Driver for quater
*> @author Y. Carissan
*>
*> @details
*> Driver for quater.
*>
*> @param[out] ireturn return code
************************************************************************
      subroutine quater(ireturn)
      Implicit none
#include "WrkSpc.fh"
#include "debug.fh"
#include "options.fh"
      Real*8 U1(3),U2(3),V1(3),V2(3)
      Real*8 V1best(3),V2best(3)
      Real*8 Q(0:3)
      Integer ireturn

      debug=.false.
      call quaterinit()

      Call RdInput_Quater(U1,U2,V1,V2)

      if (debug) then
        Write(6,*) 'Reference axis'
        Call RecPrt("U1",' ',U1,3,1)
        Call RecPrt("U2",' ',U2,3,1)
        Write(6,*) 'New axis'
        Call RecPrt("V1",' ',V1,3,1)
        Call RecPrt("V2",' ',V2,3,1)
      end if

      call QuaterSolve(U1,U2,V1,V2,Q)

c for test
      Call Add_Info('Quaternion',Q,4,8)
c for test

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

      Call GenerateGeoms(Q)

      Call QuaterReport(V1best,V2best,V1,V2)

      Call QuaterFinish()

      ireturn=0

      return
      end
