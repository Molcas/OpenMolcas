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
      Subroutine SGClose(SGS)
      use stdalloc, only: mma_deallocate
      use Struct, only: SGStruct
      Implicit None
      Type (SGStruct) SGS
      Integer nLev, nVert, lMAW, lLTV
C Unpack structure iSGStruct:
      nLev   =SGS%nLev
      nVert  =SGS%nVert
      lMAW   =SGS%lMAW
      lLTV   =SGS%lLTV
      Call mma_deallocate(SGS%ISm)
      Call mma_deallocate(SGS%DRT)
      Call mma_deallocate(SGS%Down)
      Call mma_deallocate(SGS%Up)
      Call GetMem('MAW','Free','Inte',lMAW,4*nVert)
      Call GetMem('LTV','Free','Inte',lLTV,nLev+2)
      end Subroutine SGClose
