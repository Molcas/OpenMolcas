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
************************************************************
*
*   <DOC>
*     <Name>Init\_a\_Chunk</Name>
*     <Syntax>Call Init\_a\_Chunk(ip,n)</Syntax>
*     <Arguments>
*       \Argument{ip}{The base index for a memory block}{Integer}{in}
*       \Argument{n}{The number of elements for a given block}{Integer}{in}
*     </Arguments>
*     <Purpose>
* To initialize a common chunk.
*     </Purpose>
*     <Dependencies></Dependencies>
*     <Author></Author>
*     <Modified_by></Modified_by>
*     <Side_Effects>
*     </Side_Effects>
*     <Description>
* To initialize a common chunk.
*     </Description>
*    </DOC>
*
************************************************************

      Subroutine Init_a_Chunk(ip,n)
      Implicit Real*8 (a-h,o-z)
      Common /chunk/ ip_base,n_tot
*
      ip_base= ip
      n_tot  = n
*
      Return
      End



************************************************************
*
*   <DOC>
*     <Name>Get\_a\_Chunk</Name>
*     <Syntax>Call Get\_a\_Chunk(Label,Type\_,ip,n)</Syntax>
*     <Arguments>
*       \Argument{Label}{A string without meaning}{Character*(*)}{in}
*       \Argument{Type\_}{REAL $|$ INTE}{Character*(*)}{in}
*       \Argument{ip}{The index for a memory block}{Integer}{out}
*       \Argument{n}{The number of elements for a given data type}{Integer}{in}
*     </Arguments>
*     <Purpose>
* To calculate the index, ip, for a given data type, Type\_, according to
* the base index stored at common /chunk/.
*     </Purpose>
*     <Dependencies>
*     The Get\_a\_Chunk subroutine depends from the informations
*     stored on the common chunk.
*     </Dependencies>
*     <Author></Author>
*     <Modified_by></Modified_by>
*     <Side_Effects>
*     </Side_Effects>
*     <Description>
* The Type\_ is a string of any size. It is not case sensitive, and only the four first letters matter.
* The n is the number of elements for a given data type.
*     </Description>
*    </DOC>
*
************************************************************

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

************************************************************
*
*   <DOC>
*     <Name>nChunk</Name>
*     <Syntax>Call nChunk(n)</Syntax>
*     <Arguments>
*       \Argument{n}{A number of elements for a given bloc}{Integer}{out}
*     </Arguments>
*     <Purpose>
* To return the number of elements stored on the common chunk.
*     </Purpose>
*     <Dependencies></Dependencies>
*     <Author></Author>
*     <Modified_by></Modified_by>
*     <Side_Effects>
*     </Side_Effects>
*     <Description>
* nChunk returns the number of elements stored on the common chunk.
*     </Description>
*    </DOC>
*
************************************************************

      Subroutine nChunk(n)
      Implicit Real*8 (a-h,o-z)
      Common /chunk/ ip_base,n_tot
*
      n = n_tot
*
      Return
      End
