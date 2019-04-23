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
      Subroutine Diff_Aux1(nEPotPoints,ipEPCo,nB,OneFile)
      Implicit Real*8 (a-h,o-z)

#include "WrkSpc.fh"
#include "warnings.fh"

      Character*10 Label,OneFile
      Dimension nInt(1)

*
*-- Open One-electron file.
*
      irc=-1
      Lu_One=49
      Lu_One=IsFreeUnit(Lu_One)
      Call OpnOne(irc,0,OneFile,Lu_One)
      If(irc.ne.0) then
        Write(6,*)
        Write(6,*)'ERROR! Could not open one-electron integral file.'
        Call Quit(_RC_IO_ERROR_READ_)
      Endif

*
*-- Loop over all EF0, terminate when return-code is non-zero.
*
      nEPotPoints=0
      maxCen=99999
      Call GetMem('Temporary','Allo','Real',iTmp,maxCen*3)
      Call GetMem('Idiot','Allo','Real',idiot,nB*(nB+1)/2+4)
      Do i=1,maxCen
        Write(Label,'(A3,I5)')'EF0',i
        irc=-1
        iopt=1
        iSmLbl=0
        Call iRdOne(irc,iopt,label,1,nInt,iSmLbl)
        If(irc.ne.0) Go To 9901
        irc=-1
        iopt=0
        iSmLbl=0
        Call RdOne(irc,iopt,label,1,Work(idiot),iSmLbl)
        Work(iTmp+(i-1)*3+0)=Work(idiot+nInt(1)+0)
        Work(iTmp+(i-1)*3+1)=Work(idiot+nInt(1)+1)
        Work(iTmp+(i-1)*3+2)=Work(idiot+nInt(1)+2)
        nEPotPoints=nEPotPoints+1
      Enddo
9901  Continue

*
*-- Put the coordinates and nuclear part in nice and tight vectors.
*
      Call GetMem('PotPointCoord','Allo','Real',ipEPCo,3*nEPotPoints)
      call dcopy_(3*nEPotPoints,Work(iTmp),1,Work(ipEPCo),1)

*
*-- Deallocate.
*
      Call GetMem('Temporary','Free','Real',iTmp,maxCen*3)
      Call GetMem('Idiot','Free','Real',idiot,nB*(nB+1)/2+4)


      Return
      End
