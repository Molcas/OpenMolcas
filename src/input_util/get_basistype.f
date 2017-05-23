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
* Copyright (C) Valera Veryazov                                        *
************************************************************************
        Logical function get_BasisType(BasisType)
************************************************************
*
*   <DOC>
*     <Name>get\_BasisType</Name>
*     <Syntax>get\_BasisType(BasisType)</Syntax>
*     <Arguments>
*       \Argument{BasisType}{Basis set type}{character}{in}
*     </Arguments>
*     <Purpose>Logical function to check basis set</Purpose>
*     <Dependencies></Dependencies>
*     <Author>V. Veryazov</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*       Function returns .true. if the basis set in the current
*       calculation has specific type.
*       Only 3 first characters are used, so
*        get\_BasisType('segmented') is the same as get\_BasisType('SEG').
*
*       The list of available basis set types is available at
*       src/Include/basistype.fh and basis\_library/basistype.tbl
*     </Description>
*    </DOC>
*
************************************************************
        character*(*) BasisType
        Character*3 temp, TypeCon, TypeAll,TypeRel
#include "basistype.fh"
        Logical Found
        integer BasisTypes(4)
*
        nData=0
        Found=.False.
        Call Qpg_iArray('BasType',Found,nData)
        if (.not.Found) then
          get_BasisType=.false.
          return
        endif
        call get_bastype(BasisTypes,nData)
        i=BasisTypes(1)
        if(i.le.0) then
          TypeCon='UNK'
        else
          TypeCon=BasTypeCon((i-1)*4+1:(i-1)*4+3)
        endif
        i=BasisTypes(2)
        if(i.le.0) then
          TypeAll='UNK'
        else
          TypeAll=BasTypeAll((i-1)*4+1:(i-1)*4+3)
        endif
        i=BasisTypes(3)
        if(i.le.0) then
          TypeRel='UNK'
        else
          TypeRel=BasTypeRel((i-1)*4+1:(i-1)*4+3)
        endif

        i=len(BasisType)
        if(i.gt.3) i=3
        temp='___'
        temp(1:i)=BasisType(1:i)
        call UpCase(temp)
        get_BasisType=.false.
        if(temp.eq.TypeCon.or.temp.eq.TypeAll.or.temp.eq.TypeRel)
     &  get_BasisType=.true.

        return
        end
