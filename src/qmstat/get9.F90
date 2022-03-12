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
      Subroutine Get9(Ract,Coord,info_atom,iQ_Atoms,iDiskSa)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "files_qmstat.fh"
#include "WrkSpc.fh"

      Dimension Coord(MxAt*3)
      Dimension info_atom(MxAt)
      Character Head*200

      Call WrRdSim(iLuSaIn,2,iDiskSa,iTcSim,64,Etot,Ract,nPart,Gamold   &
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
!
!---- We dummy-read the induced dipoles from the sampfile.
!
!      Call GetMem('Dummy','Allo','Real',iDum,nPart*nPol)
!      Do 1777, i=1,3
!        Call dDaFile(iLuSaIn,2,Work(iDum),nPol*nPart,iDiskSa)
!1777  Continue
!      Call GetMem('Dummy','Free','Real',iDum,nPart*nPol)
!
!---- And now we place the QM-molecule in proper place and set some
!     numbers to zero or one so we only collect configurations from
!     the sampfile.
!
      Call PlaceIt9(Coord,Cordst,info_atom,iQ_Atoms)
      delX=0
      delFi=0
      delR=0
      nMacro=1
      nMicro=1
!
!---- Some printing if requested.
!
      If(iPrint.ge.15) then
        Write(Head,*)'Coordinates after substitution in configuration r'&
     &//'ead from sampfile.'
        Call Cooout(Head,Cordst,nPart,nCent)
      Endif

      Return
      End
