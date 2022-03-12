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
      Subroutine Put8(Ract,Etot,Gamma,Gam,Esav)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "files_qmstat.fh"
#include "WrkSpc.fh"

      Character Head*200

      Call DaName(iLuStUt,StFilUt)  !Here follows a general output
      iDisk=0                   !to the startfile.
      Call WrRdSim(iLuStUt,1,iDisk,iTcSim,64,Etot,Ract,nPart,Gamma      &
     &            ,Gam,Esav)
      iTcSim(1)=iDisk
      Do 1010,i=1,3  !In this loop the coordinates are put on file.
                     !The loop is needed due to how Cordst is
                     !statically allocated.
        Call GetMem('CTemp','Allo','Real',iCT,nPart*nCent)
        Do 1011,j=1,nPart*nCent
          Work(iCT+j-1)=Cordst(j,i)
1011    Continue
        Call dDaFile(iLuStUt,1,Work(iCT),nPart*nCent,iDisk)
        iTcSim(1+i)=iDisk
        Call GetMem('CTemp','Free','Real',iCT,nPart*nCent)
1010  Continue
      iDisk=0
      Call WrRdSim(iLuStUt,1,iDisk,iTcSim,64,Etot,Ract,nPart,Gamma      &
     &            ,Gam,Esav)
      Call DaClos(iLuStUt)
      If(iPrint.ge.10) then  !Print the stored configuration.
        Write(Head,*)' Coordinates put on the startfile solvent configu'&
     &//'ration.'
        Call Cooout(Head,Cordst,nPart,nCent)
      Endif

      Return
      End
