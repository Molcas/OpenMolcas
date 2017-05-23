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
       Function isFreeUnit(iseed)

************************************************************
*
*   <DOC>
*     <Name>isFreeUnit</Name>
*     <Syntax>isFreeUnit(init)</Syntax>
*     <Arguments>
*       \Argument{init}{guess for unit number}{integer}{inout}
*     </Arguments>
*     <Purpose>Find free unit number</Purpose>
*     <Dependencies></Dependencies>
*     <Author>V. Veryazov</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*       Find unused unit number, starting from initial value
*     </Description>
*    </DOC>
*
************************************************************
#include "fio.fh"
c      check free chanal, starting from init
       Logical Opened
*
CVV since more and more developers' calling isfreeunit with constant...
       init=iseed
       if(init.lt.1.or.init.gt.300) then
        write(6,*) '*** Possible bug in opening file'
        write(6,*) '*** isFreeUnit resets the unit number'
        init=12
       endif
       isFreeUnit=-init
       kan=min(init,Mxfile-1)
       kan0=kan
       Go to 2
 1     Continue
       If (kan.eq.MxFile+1) kan=10
       If (kan.eq.kan0) Then
          Call fastio('STATUS')
          Write (6,*) ' isFreeUnit: no available unit!'
          Call QTrace()
          Call Abend()
       End If
 2     Continue
*
*----- Check for Dafile
*
       If (kan.gt.1.and.kan.le.MxFile.and.isOpen(kan).eq.1) Then
          kan=kan+1
          Go To 1
       End If
*
*----- Check for Fortran I/O
*
       Inquire(unit=kan,opened=Opened)
       If (Opened) Then
          kan=kan+1
          Go To 1
       End If
       isFreeUnit=kan
       Return
       End
