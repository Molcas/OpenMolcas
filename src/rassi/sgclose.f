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

      Call mma_deallocate(SGS%ISm)
      Call mma_deallocate(SGS%DRT)
      Call mma_deallocate(SGS%Down)
      Call mma_deallocate(SGS%Up)
      Call mma_deallocate(SGS%MAW)
      Call mma_deallocate(SGS%LTV)

      end Subroutine SGClose
