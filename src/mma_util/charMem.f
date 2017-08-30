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
*  AllocateC
*
*> @brief
*>   Pseudo allocation of memory for a string of length \p nipM of a vector of dimension \p nip
*> @author Eugeniusz Bednarz
*>
*> @details
*> Pseudo allocation of memory for a string of length \p nipM of a vector of dimension \p nip.
*>
*> Example: A call
*> \code
*> Call AllocateC('Label',ip,m,n)
*> \endcode
*> will allocate a character vector of \c n elements, where any element of this vector has a length \c m.
*> It corresponds to the classical declaration: ``Character*(m) A(n)``, where \c A is a variable name.
*> Moreover, a memory is allocated as multiplicity of real number. The reason of such allocation
*> is the ::get_carray subroutine which always returns a characters of length of multiplicity of 8 bytes (``real*8``).
*>
*> @param[in]  label The label of allocated memory
*> @param[out] ip    The index to proper element of the character vector
*> @param[in]  nipM  The length of the character
*> @param[in]  nip   The dimension of character vector
************************************************************************
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

************************************************************************
*  DeAllocateC
*
*> @brief
*>   Pseudo deallocation of memory for a string of length \p nipM of a vector of dimension \p nip
*> @author Eugeniusz Bednarz
*>
*> @details
*> Pseudo deallocation of memory for a string of length \p nipM of a vector of dimension \p nip.
*>
*> Example: A call
*> \code
*> Call DeAllocateC('Label',ip,m,n)
*> \endcode
*> will deallocate a character vector of \p n elements, where any element of this vector has a length \p m.
*>
*> @param[in]  label The label of allocated memory
*> @param[out] ip    The index to proper element of the character vector
*> @param[in]  nipM  The length of the character
*> @param[in]  nip   The dimension of character vector
************************************************************************
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
