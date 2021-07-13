!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

subroutine Term_David(ICICH,iter,lRoots,nConf,Vector,JOBIPH,LuDavid,iDisk)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Terminate the Davidson diagonalization                           *
!                                                                      *
!     calling arguments:                                               *
!     ICICH   : integer                                                *
!               switch enabling root selection                         *
!     JOBIPH  : integer                                                *
!               logical unit number of the JOBIPH file                 *
!     iDisk   : integer                                                *
!               disk address of the first CI vector on JOBIPH          *
!     iter    : integer                                                *
!               iteration count of the final result                    *
!     nConf   : integer                                                *
!               length of the CI vector in the CSF basis               *
!     Vector  : array of real*8                                        *
!               temporary vector of length nConf                       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

implicit integer(A-Z)
real*8 Vector(nConf)
#include "rasdim.fh"
#include "davctl.fh"
#include "WrkSpc.fh"

! check input arguments
if (nConf < 0) then
  write(6,*) 'Term_David: nConf less than 0'
  write(6,*) 'nConf = ',nConf
  call Abend()
end if
if (iter < 0) then
  write(6,*) 'Term_David: iter less than 0'
  write(6,*) 'iter = ',iter
  call Abend()
end if
if (iter > mxCiIt) then
  write(6,*) 'Term_David: iter greater than mxCiIt'
  write(6,*) 'iter, mxCiIt = ',iter,mxCiIt
  call Abend()
end if

! Restore the final CI vectors and save them for further use.
! If the root selectioning option has been enabled calculate
! also the overlap elemtents with the test vectors
if (ICICH == 1) then
  call GetMem('CIovlp1','Allo','Real',lOvlp1,lRoots*lRoots)
  call dCopy_(lRoots*lRoots,[0.0d0],0,Work(lOvlp1),(1))
  call GetMem('CIovlp2','Allo','Real',lOvlp2,lRoots*lRoots)
  call dCopy_(lRoots*lRoots,[0.0d0],0,Work(lOvlp2),(1))
end if
do iRoot=1,lRoots
  call Load_tmp_CI_vec(iRoot,nConf,Vector,LuDavid)
  call DDaFile(JOBIPH,1,Vector,nConf,iDisk)
  if (ICICH == 1) then
    call CIovlp(iRoot,Work(lOvlp1),Work(lOvlp2),Vector)
  end if
end do

! If the root selectioning option has been enabled
! make a new choice of the current roots
if (ICICH == 1) then
  call CIselect(Work(lOvlp1),Work(lOvlp2))
  call GetMem('CIovlp2','Free','Real',lOvlp2,lRoots*lRoots)
  call GetMem('CIovlp1','Free','Real',lOvlp1,lRoots*lRoots)
end if

! deallocate memory which was used as records of the RAM disk
if (save_mode /= on_disk) then
  do iRecNo=1,MxMemStk
    iMem = memory_address(iRecNo)
    call GetMem(' ','Free','Real',iMem,nConf)
  end do
end if

return

end subroutine Term_David
