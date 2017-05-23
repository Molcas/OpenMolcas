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
* Copyright (C) Eugeniusz Bednarz                                      *
************************************************************************
* <DOC>
*  <NAME>AllocateC</NAME>
*  <Syntax>AllocateC(label,ip,nipM,nip)</Syntax>
*  <Arguments>
*    \Argument{label}{The label of allocated memory}{Character*(*)}{In}
*    \Argument{ip}{The index to proper element of the character vector}{Integer}{Out}
*    \Argument{nipM}{The length of the character}{Integer}{In}
*    \Argument{nip}{The dimension of character vector}{Integer}{In}
*  </Arguments>
*  <Purpose>Pseudo allocation of memory for a string of length nipM of a vector of dimension nip.</Purpose>
*  <Dependencies></Dependencies>
*  <Author>Eugeniusz Bednarz</Author>
*  <Modified_by></Modified_by>
*  <Side_Effects></Side_Effects>
*  <Description>
*    Pseudo allocation of memory for a string of length nipM of a vector of dimension nip.
*    Example:
*      A call
*         Call AllocateC('Label',ip,m,n)
*      will allocate a character vector of n elements, where any element of this vector has a length m.
*      It corresponds to the classical declaration: Character*(m) A(n), where A is a variable name.
*      Moreover, a memory is allocated as multiplicity of real number. The reason of such allocation
*      is the get\_carray subroutine which always returns a characters of length of multiplicity of 8 bytes[real*8].
*  </Description>
* </DOC>
      Subroutine AllocateC(label,ip,nipM,nip)
      Implicit None
#include "WrkSpc.fh"
#include "SysDef.fh"
      !
      Character*(*)
     +     label
      Integer
     +     ip,
     +     nipM,
     +     nip
      !
      Integer
     +     i,
     +     j

      ! fitt character size to the multiplicity of the real number
      i=nipM*nip
      j=mod(i,RtoB)
      If(j.Gt.0) i=i+j

      ! Allocate memory
      Call GetMem(label,'Allo','Char',ip,i)

      Return
      End

* <DOC>
*  <NAME>DeAllocateC</NAME>
*  <Syntax>DeAllocateC(label,ip,nipM,nip)</Syntax>
*  <Arguments>
*    \Argument{label}{The label of allocated memory}{Character*(*)}{In}
*    \Argument{ip}{The index to proper element of the character vector}{Integer}{Out}
*    \Argument{nipM}{The length of the character}{Integer}{In}
*    \Argument{nip}{The dimension of character vector}{Integer}{In}
*  </Arguments>
*  <Purpose>Pseudo deallocation of memory for a string of length nipM of a vector of dimension nip.</Purpose>
*  <Dependencies></Dependencies>
*  <Author>Eugeniusz Bednarz</Author>
*  <Modified_by></Modified_by>
*  <Side_Effects></Side_Effects>
*  <Description>
*    Pseudo deallocation of memory for a string of length nipM of a vector of dimension nip.
*    Example:
*      A call
*         Call DeAllocateC('Label',ip,m,n)
*      will deallocate a character vector of n elements, where any element of this vector has a length m.
*  </Description>
* </DOC>
      Subroutine DeAllocateC(label,ip,nipM,nip)
      Implicit None
#include "WrkSpc.fh"
#include "SysDef.fh"
      !
      Character*(*)
     +     label
      Integer
     +     ip,
     +     nipM,
     +     nip
      !
      Integer
     +     i

      i=nipM*nip
      ! DeAllocate memory
      Call GetMem(label,'Free','Char',ip,i)

      Return
      End
