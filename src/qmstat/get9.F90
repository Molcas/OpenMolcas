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

! GET NUMBERS FROM SAMPFILE.
subroutine Get9(Ract,Coord,info_atom,iQ_Atoms,iDiskSa)

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "qminp.fh"
#include "files_qmstat.fh"
#include "WrkSpc.fh"

dimension Coord(MxAt*3)
dimension info_atom(MxAt)
character Head*200

call WrRdSim(iLuSaIn,2,iDiskSa,iTcSim,64,Etot,Ract,nPart,Gamold,GaOld,Esub)
iDiskSa = iTcSim(1)
do i=1,3
  call GetMem('CTemp','Allo','Real',iCT,nPart*nCent)
  call dDaFile(iLuSaIn,2,Work(iCT),nPart*nCent,iDiskSa)
  do j=1,nCent*nPart
    Cordst(j,i) = Work(iCT+j-1)
  end do
  call GetMem('CTemp','Free','Real',iCT,nPart*nCent)
  iDiskSa = iTcSim(i+1)
end do

! We dummy-read the induced dipoles from the sampfile.

!call GetMem('Dummy','Allo','Real',iDum,nPart*nPol)
!do i=1,3
!  call dDaFile(iLuSaIn,2,Work(iDum),nPol*nPart,iDiskSa)
!end do
!call GetMem('Dummy','Free','Real',iDum,nPart*nPol)

! And now we place the QM-molecule in proper place and set some
! numbers to zero or one so we only collect configurations from
! the sampfile.

call PlaceIt9(Coord,Cordst,info_atom,iQ_Atoms)
delX = 0
delFi = 0
delR = 0
nMacro = 1
nMicro = 1

! Some printing if requested.

if (iPrint >= 15) then
  write(Head,*) 'Coordinates after substitution in configuration read from sampfile.'
  call Cooout(Head,Cordst,nPart,nCent)
end if

return

end subroutine Get9
