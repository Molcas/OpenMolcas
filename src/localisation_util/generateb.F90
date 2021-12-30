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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine GenerateB(CMO,nBas,nOrb2Loc,ipLbl_AO,ipLbl,nComp,Debug)
! Author: T.B. Pedersen
!
! Purpose: generate the dipole matrices for Boys localisation, i.e.
!          transform from AO to MO basis.

implicit real*8(a-h,o-z)
real*8 CMO(*)
integer ipLbl_AO(nComp), ipLbl(nComp)
logical Debug
#include "WrkSpc.fh"

if ((nBas < 1) .or. (nOrb2Loc < 1)) return

lDbar = nBas*nOrb2Loc
call GetMem('Dbar','Allo','Real',ipDbar,lDbar)
do iComp=1,nComp
  call DGEMM_('N','N',nBas,nOrb2Loc,nBas,1.0d0,Work(ipLbl_AO(iComp)),nBas,CMO,nBas,0.0d0,Work(ipDbar),nBas)
  call DGEMM_('T','N',nOrb2Loc,nOrb2Loc,nBas,1.0d0,CMO,nBas,Work(ipDbar),nBas,0.0d0,Work(ipLbl(iComp)),nOrb2Loc)
end do
call GetMem('Dbar','Free','Real',ipDbar,lDbar)

if (Debug) then
  write(6,*)
  write(6,*) 'In GenerateB'
  write(6,*) '------------'
  write(6,*) '[Assuming doubly occupied orbitals]'
  do iComp=1,nComp
    ip0 = ipLbl(iComp)-1
    Cmp = 0.0d0
    do iMO=1,nOrb2Loc
      Cmp = Cmp+Work(ip0+nOrb2Loc*(iMO-1)+iMO)
    end do
    Cmp = 2.0d0*Cmp
    write(6,'(A,I5,1X,F15.8)') 'Component, Exp. Val.:',iComp,Cmp
    do j=1,nOrb2Loc-1
      do i=j+1,nOrb2Loc
        kij = ip0+nOrb2Loc*(j-1)+i
        kji = ip0+nOrb2Loc*(i-1)+j
        Tst = Work(kij)-Work(kji)
        if (abs(Tst) > 1.0d-14) then
          write(6,*) 'GenerateB: broken symmetry!'
          write(6,*) '  Component: ',iComp
          write(6,*) '  i and j  : ',i,j
          write(6,*) '  Dij      : ',Work(kij)
          write(6,*) '  Dji      : ',Work(kji)
          write(6,*) '  Diff.    : ',Tst
          call SysAbendMsg('GenerateB','Broken symmetry!',' ')
        end if
      end do
    end do
  end do
end if

end subroutine GenerateB
