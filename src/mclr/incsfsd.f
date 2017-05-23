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
      Subroutine InCSFSD(iState,State_sym,GUGA)
      Implicit Real*8 (a-h,o-z)
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "Files_mclr.fh"
#include "spinfo_mclr.fh"
#include "WrkSpc.fh"
#include "csfbas_mclr.fh"
      Logical GUGA
      Integer State_sym
*
#include "csfsd.fh"

*     Place pointer
*
      iSym=iEor(State_Sym-1,iState-1)+1
*
      If (isym.eq.1.and.i1.eq.1) Return
      If (isym.eq.iAnders) Return
*
      iAdr=2
      If (iSym.eq.1) iAdr=1
      iad=0
      Do i=1,iState-1
         Call iDafile(LUCSF2SD,0,idum,lldet,iad)
         Call iDafile(LUCSF2SD,0,idum,lconf,iad)
      End Do
*
      If (iSym.ne.1) Then
         If (iAnders.eq.-9)  Then
             Call Getmem('KICTS2','ALLO','INTEGER',Kicts(2),lldet)
             Call Getmem('KICONF2','ALLO','INTEGER',kiconf(2),lConf)
             iAllo=1
         End If
         iAnders=isym
      End If
      If (iSym.eq.1) Then
          If (i1.eq.-9) Then
           Call Getmem('KICTS1','ALLO','INTEGER',Kicts(1),lldet)
           If (.true.)
     &      Call Getmem('KICONF1','ALLO','INTEGER',kiconf(1),lConf)
           i1=1
          End If
      End If

!      open(unit=1422,file="det.index") ! yma
!      do i=1,lldet
!        write(1422,*)iWork(kicts(iAdr)+i-1)
!      end do
!      close(1422)

! calculated from zoo.f, the GUGA number for determinent
      Call iDafile(LUCSF2SD,2,iWork(kicts(iAdr)),lldet,iad)
      Call iDafile(LUCSF2SD,2,iWork(kiconf(iAdr)),lconf,iad)

      Return

c Avoid unused argument warnings
      If (.False.) Call Unused_logical(GUGA)
      End
