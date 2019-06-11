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
*
* ALL SUBROUTINE HERE HAVE THE PURPOSE TO PUT OR GET NUMBERS ON/FROM
* STARTFILES AND SAMPFILES.
*
      Subroutine Get8(Ract,Etot)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "files_qmstat.fh"
#include "WrkSpc.fh"

      Character Head*200

      iDisk=0
      Call DaName(iLuStIn,StFilIn)
      Call WrRdSim(iLuStIn,2,iDisk,iTcSim,64,Etot,Ract,nPart,Gamold
     &,GaOld,Esub)
      iDisk=iTcSim(1)
*
*---- In this loop we read the coordinates. The construction of Cordst
*     makes this loop necessary. Maybe we should consider going to
*     dynamic allocation.
*
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
*
*---- If requested, print initial coordinates.
*
      If(iPrint.ge.10) then
        Write(Head,*)'Coordinates read from startfile.'
        Call Cooout(Head,Cordst,nPart,nCent)
      Endif

      Return
      End


      Subroutine Get9(Ract,Coord,info_atom,iQ_Atoms,iDiskSa)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "files_qmstat.fh"
#include "WrkSpc.fh"

      Dimension Coord(MxAt*3)
      Dimension info_atom(MxAt)
      Character Head*200

      Call WrRdSim(iLuSaIn,2,iDiskSa,iTcSim,64,Etot,Ract,nPart,Gamold
     &,GaOld,Esub)
      iDiskSa=iTcSim(1)
      Do 1022, i=1,3
        Call GetMem('CTemp','Allo','Real',iCT,nPart*nCent)
        Call dDaFile(iLuSaIn,2,Work(iCT),nPart*nCent,iDiskSa)
        Do 1023, j=1,nCent*nPart
          Cordst(j,i)=Work(iCT+j-1)
1023    Continue
        Call GetMem('CTemp','Free','Real',iCT,nPart*nCent)
        iDiskSa=iTcSim(i+1)
1022  Continue
*
*---- We dummy-read the induced dipoles from the sampfile.
*
*      Call GetMem('Dummy','Allo','Real',iDum,nPart*nPol)
*      Do 1777, i=1,3
*        Call dDaFile(iLuSaIn,2,Work(iDum),nPol*nPart,iDiskSa)
*1777  Continue
*      Call GetMem('Dummy','Free','Real',iDum,nPart*nPol)
*
*---- And now we place the QM-molecule in proper place and set some
*     numbers to zero or one so we only collect configurations from
*     the sampfile.
*
      Call PlaceIt9(Coord,Cordst,info_atom,iQ_Atoms)
      delX=0
      delFi=0
      delR=0
      nMacro=1
      nMicro=1
*
*---- Some printing if requested.
*
      If(iPrint.ge.15) then
        Write(Head,*)'Coordinates after substitution in configuration r'
     &//'ead from sampfile.'
        Call Cooout(Head,Cordst,nPart,nCent)
      Endif

      Return
      End


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
      Call WrRdSim(iLuSaUt,1,iDisk,iTcSim,64,Etot,Ract,nPart
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
*      Do 777, i=1,3
*        Call dDaFile(iLuSaUt,1,-Work(iDT(i)),nPart*nPol,iDisk)
*777   Continue                             !Induced dipoles.
*      iTcSim(5)=iDisk
      iDiskHead=iDiskOld
      Call WrRdSim(iLuSaUt,1,iDiskHead,iTcSim,64,Etot,Ract,nPart
     &            ,Gamma,Gam,Esav) !Put header again, but
                         !now with a
                    !meaningful iTcSim vector that contains
             !the table of content which simplifies reading

      If(iPrint.ge.15) then
        Write(Head,*)' Coordinates put on sampfile.'
        Call Cooout(Head,Cordst,nPart,nCent)
      Endif

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer_array(iDT)
      End


      Subroutine Put8(Ract,Etot,Gamma,Gam,Esav)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "files_qmstat.fh"
#include "WrkSpc.fh"

      Character Head*200

      Call DaName(iLuStUt,StFilUt)  !Here follows a general output
      iDisk=0                   !to the startfile.
      Call WrRdSim(iLuStUt,1,iDisk,iTcSim,64,Etot,Ract,nPart,Gamma
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
      Call WrRdSim(iLuStUt,1,iDisk,iTcSim,64,Etot,Ract,nPart,Gamma
     &            ,Gam,Esav)
      Call DaClos(iLuStUt)
      If(iPrint.ge.10) then  !Print the stored configuration.
        Write(Head,*)' Coordinates put on the startfile solvent configu'
     &//'ration.'
        Call Cooout(Head,Cordst,nPart,nCent)
      Endif

      Return
      End
