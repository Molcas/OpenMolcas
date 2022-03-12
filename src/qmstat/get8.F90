!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!
! ALL SUBROUTINE HERE HAVE THE PURPOSE TO PUT OR GET NUMBERS ON/FROM
! STARTFILES AND SAMPFILES.
!
      Subroutine Get8(Ract,Etot)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "files_qmstat.fh"
#include "WrkSpc.fh"

      Character Head*200

      iDisk=0
      Call DaName(iLuStIn,StFilIn)
      Call WrRdSim(iLuStIn,2,iDisk,iTcSim,64,Etot,Ract,nPart,Gamold     &
     &,GaOld,Esub)
      iDisk=iTcSim(1)
!
!---- In this loop we read the coordinates. The construction of Cordst
!     makes this loop necessary. Maybe we should consider going to
!     dynamic allocation.
!
      Do 1020, i=1,3
        Call GetMem('CTemp','Allo','Real',iCT,nPart*nCent)
        Call dDaFile(iLuStIn,2,Work(iCT),nPart*nCent,iDisk)
        Do 1021, j=1,nCent*nPart
          Cordst(j,i)=Work(iCT+j-1)
1021    Continue
        Call GetMem('CTemp','Free','Real',iCT,nPart*nCent)
        iDisk=iTcSim(i+1)
1020  Continue
      Call DaClos(iLuStIn)
!
!---- If requested, print initial coordinates.
!
      If(iPrint.ge.10) then
        Write(Head,*)'Coordinates read from startfile.'
        Call Cooout(Head,Cordst,nPart,nCent)
      Endif

      Return
      End
