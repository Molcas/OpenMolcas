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
* Copyright (C) 2004,2005, Alexander Wolf                              *
*               2004,2005, Markus Reiher                               *
************************************************************************
      subroutine t_create (sused,scrleng,scrchar,tnumber,tcounter,
     *                     ttimes,tmult,t,tscrleng,tscrchar,ttimestot)
c
c******************************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 12.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c******************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer sused,scrleng(maxsnumber),tnumber,tcounter(maxsnumber),
     *        ttimes(maxsnumber),tmult(maxsnumber),tscrleng(maxsnumber),
     *        ttimestot
c
      character*(4) t(maxsnumber)
      character*(9) scrchar(maxsnumber)
      character*(11) tscrchar(maxsnumber)
c
      integer j
      character*(3) dkh_int2char
c
c
      do 10 j=1,sused
        tcounter(j)=1
        ttimes(j)=0
        tmult(j)=0
        t(j)='    '
        tscrchar(j)='           '
        tscrleng(j)=0
  10  continue
      tnumber=0
      ttimestot=0
      do 20 j=1,sused
        t(j)(1:4)='T'//dkh_int2char(j)
        tscrleng(j)=scrleng(j)+2
        tscrchar(j)(1:tscrleng(j))='P'//scrchar(j)(1:scrleng(j))//'P'
 20   continue
c
      if (dbgflg.ge.1) then
        write (dbgunit,1001)
1001    format (/2X,'Creation and initialization of scratch arrays ',
     *          '"t" successful.')
      endif
c
      if (dbgflg.ge.3) then
         write (dbgunit,1156)sused
 1156    format (/2X,'The following ',I3,' Txxx-matrices have been set',
     *           ' up:',//2X,'Txxx',2X,'tscrchar',4X,'tmult',2X,
     *           'tscrleng',/)
         do 6298 j=1,sused
           write (dbgunit,8877) t(j),tscrchar(j),tmult(j),tscrleng(j)
8877       format (2X,A4,2X,A11,2X,I3,5X,I3)
6298     continue
      endif

      return
      end
