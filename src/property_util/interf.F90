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

subroutine Interf(i_root,Ene,isuseene,iscasvb)
!***********************************************************************
!                                                                      *
!     Object: Driver toward MOLDEN interface                           *
!                                                                      *
!***********************************************************************

implicit real*8(a-h,o-z)
#include "rasdim.fh"
#include "general.fh"
#include "casvb.fh"
#include "WrkSpc.fh"
#include "output_ras.fh"
character*10 Filename
character*80 Note
dimension Ene(*)
dimension iDum(7,8)

!                                                                      *
!***********************************************************************
!                                                                      *
! Compute memory requirements and allocate memory

nB = 0
nB2 = 0
do iS=1,nSym
  nB = nB+nBas(iS)
  nB2 = nB2+nBas(iS)**2
end do

call GetMem('OCCA','Allo','Real',ipOccA,nB)
call GetMem('OCCB','Allo','Real',ipOccB,nB)
call GetMem('ENERGY','Allo','Real',ipEA,2*nB)
ipEB = ipEA+nB
call GetMem('CMOA','Allo','Real',ipCA,nB**2)
call GetMem('CMOB','Allo','Real',ipCB,nB**2)
call GetMem('AdCMOA','Allo','Real',mAdCMOA,nB2)
call GetMem('AdCMOB','Allo','Real',mAdCMOB,nB2)
!                                                                      *
!***********************************************************************
!                                                                      *
! For the moment: Orbital energies just zero
if (isuseene /= 0) then
  do i=1,nB
    Work(ipEA+i-1) = Ene(i)
    Work(ipEB+i-1) = Ene(i)
  end do
else
  call FZero(Work(ipEA),2*nB)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Get the coeff. of sym adapted basis functions (ipCA, ipCB) and
! the spin orbital occupations (ipOccA, ipOccB)

call Dens_IF(i_root,Work(ipCA),Work(ipCB),Work(ipOccA),Work(ipOccB))
call Dens_IF_SCF(Work(ipCA),Work(mAdCMOA),'B')
call Dens_IF_SCF(Work(ipCB),Work(mAdCMOB),'B')
!                                                                      *
!***********************************************************************
!                                                                      *
! Write out info on a temporary vector file.

Note = 'Temporary orbital file for the MOLDEN interface.'
LuTmp = 50
LuTmp = IsFreeUnit(LuTmp)
iUHF = IFVB
if (i_root /= 0) iUHF = 1
call WrVec_('TMPORB',LuTmp,'COE',iUHF,nSym,nBas,nBas,Work(mAdCMOA),Work(mAdCMOB),Work(ipOccA),Work(ipOccB),Work(ipEA),Work(ipEB),&
            iDum,Note,0)
!                                                                      *
!***********************************************************************
!                                                                      *
call GetMem('AdCMOB','Free','Real',mAdCMOB,nB2)
call GetMem('AdCMOA','Free','Real',mAdCMOA,nB2)
call GetMem('CMOA','Free','Real',ipCA,nB**2)
call GetMem('CMOB','Free','Real',ipCB,nB**2)
call GetMem('ENERGY','Free','Real',ipEA,nB)
call GetMem('OCCA','Free','Real',ipOccA,nB)
call GetMem('OCCB','Free','Real',ipOccB,nB)
!                                                                      *
!***********************************************************************
!                                                                      *
if (i_root /= 0) then
  if (i_root <= 9) then
    write(filename,'(A7,I1)') 'MD_CAS.',i_root
  else if (i_root <= 99) then
    write(filename,'(A7,I2)') 'MD_CAS.',i_root
  else if (i_root <= 999) then
    write(filename,'(A7,I3)') 'MD_CAS.',i_root
  else
    filename = 'MD_CAS.x'
  end if
else
  filename = 'MD_CAS'
end if
if (iscasvb == 1) filename = 'MD_VB'
!                                                                      *
!***********************************************************************
!                                                                      *
! Call the generic MOLDEN interface

call Molden_Interface(iUHF,'TMPORB',filename)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Interf
