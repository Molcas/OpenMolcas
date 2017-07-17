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
*  Init_a_Chunk
*
*> @brief
*>   Initialize a common chunk
*>
*> @details
*> Initialize a common chunk.
*>
*> @param[in] ip The base index for a memory block
*> @param[in] n  The number of elements for a given block
************************************************************************
      Subroutine Init_a_Chunk(ip,n)
      Implicit Real*8 (a-h,o-z)
      Common /chunk/ ip_base,n_tot
*
      ip_base= ip
      n_tot  = n
*
      Return
      End

************************************************************************
*  Get_a_Chunk
*
*> @brief
*>   Calculate the index, \p ip, for a given data type, \p Type_, according to the base index stored at common ``/chunk/``
*>
*> @details
*> The \p Type_ is a string of any size. It is not case sensitive, and only the four first letters matter.
*> The \p n is the number of elements for a given data type.
*>
*> @note
*> The ::Get_a_Chunk subroutine depends from the informations
*> stored on the common chunk.
*>
*> @param[in]  Label A string without meaning
*> @param[in]  Type_ ``REAL`` / ``INTE``
*> @param[out] ip    The index for a memory block
*> @param[in]  n     The number of elements for a given data type
************************************************************************
      Subroutine Get_a_Chunk(Label,Type_,ip,n)
      Implicit Real*8 (a-h,o-z)
#include "SysDef.fh"
#include "WrkSpc.fh"
      Character*(*) Type_, Label
      Character*8   Type
      Common /chunk/ ip_base,n_tot
*
      Type=Type_
      Call Upcase(Type)
*
      If (Type.eq.'REAL') Then
         ip=ip_base + n_tot
         n_tot = n_tot + n
      Else If (Type.eq.'INTE') Then
         ip=ip_of_iWork(Work(ip_base + n_tot))
         n_tot = n_tot + (n-1)/RtoI + 1
      Else
         Write (6,*) 'Get_a_chunk: invalid type!'
         Write (6,'(2A)') 'Type=',Type
         Call QTrace()
         Call Abend()
      End If
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_character(Label)
      End

************************************************************************
*  nChunk
*
*> @brief
*>   Return the number of elements stored on the common chunk
*>
*> @details
*> ::nChunk returns the number of elements stored on the common chunk.
*>
*> @param[out] n A number of elements for a given block
************************************************************************
      Subroutine nChunk(n)
      Implicit Real*8 (a-h,o-z)
      Common /chunk/ ip_base,n_tot
*
      n = n_tot
*
      Return
      End
