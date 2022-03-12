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
      Subroutine Put9(Etot,Ract,iDT,iHowMSamp,Gamma,Gam,Esav,iDisk)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "files_qmstat.fh"
#include "WrkSpc.fh"

      Dimension iDT(3)
      Character Head*200

      iHowMSamp=iHowMSamp+1
      iDiskOld=iDisk
      Call WrRdSim(iLuSaUt,1,iDisk,iTcSim,64,Etot,Ract,nPart            &
     &            ,Gamma,Gam,Esav)  !A header
      iTcSim(1)=iDisk
      Do 1024, i=1,3
        Call GetMem('CTemp','Allo','Real',iCT,nPart*nCent)
        Do 1025, j=1,nCent*nPart
          Work(iCT+j-1)=Cordst(j,i)
1025    Continue
        Call dDaFile(iLuSaUt,1,Work(iCT),nPart*nCent,iDisk)
                 !The solvent coordinates.
        Call GetMem('CTemp','Free','Real',iCT,nPart*nCent)
        iTcSim(i+1)=iDisk
1024  Continue
!      Do 777, i=1,3
!        Call dDaFile(iLuSaUt,1,-Work(iDT(i)),nPart*nPol,iDisk)
!777   Continue                             !Induced dipoles.
!      iTcSim(5)=iDisk
      iDiskHead=iDiskOld
      Call WrRdSim(iLuSaUt,1,iDiskHead,iTcSim,64,Etot,Ract,nPart        &
     &            ,Gamma,Gam,Esav) !Put header again, but
                         !now with a
                    !meaningful iTcSim vector that contains
             !the table of content which simplifies reading

      If(iPrint.ge.15) then
        Write(Head,*)' Coordinates put on sampfile.'
        Call Cooout(Head,Cordst,nPart,nCent)
      Endif

      Return
! Avoid unused argument warnings
      If (.False.) Call Unused_integer_array(iDT)
      End
