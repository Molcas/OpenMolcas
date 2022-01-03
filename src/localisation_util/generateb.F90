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

use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: CMO(*)
integer(kind=iwp) :: nBas, nOrb2Loc, nComp, ipLbl_AO(nComp), ipLbl(nComp)
logical(kind=iwp) :: Debug
#include "WrkSpc.fh"
integer(kind=iwp) :: i, iComp, iMO, ip0, ipDbar, j, kij, kji, lDbar
real(kind=wp) :: Cmp, Tst

if ((nBas < 1) .or. (nOrb2Loc < 1)) return

lDbar = nBas*nOrb2Loc
call GetMem('Dbar','Allo','Real',ipDbar,lDbar)
do iComp=1,nComp
  call DGEMM_('N','N',nBas,nOrb2Loc,nBas,One,Work(ipLbl_AO(iComp)),nBas,CMO,nBas,Zero,Work(ipDbar),nBas)
  call DGEMM_('T','N',nOrb2Loc,nOrb2Loc,nBas,One,CMO,nBas,Work(ipDbar),nBas,Zero,Work(ipLbl(iComp)),nOrb2Loc)
end do
call GetMem('Dbar','Free','Real',ipDbar,lDbar)

if (Debug) then
  write(u6,*)
  write(u6,*) 'In GenerateB'
  write(u6,*) '------------'
  write(u6,*) '[Assuming doubly occupied orbitals]'
  do iComp=1,nComp
    ip0 = ipLbl(iComp)-1
    Cmp = Zero
    do iMO=1,nOrb2Loc
      Cmp = Cmp+Work(ip0+nOrb2Loc*(iMO-1)+iMO)
    end do
    Cmp = Two*Cmp
    write(u6,'(A,I5,1X,F15.8)') 'Component, Exp. Val.:',iComp,Cmp
    do j=1,nOrb2Loc-1
      do i=j+1,nOrb2Loc
        kij = ip0+nOrb2Loc*(j-1)+i
        kji = ip0+nOrb2Loc*(i-1)+j
        Tst = Work(kij)-Work(kji)
        if (abs(Tst) > 1.0e-14_wp) then
          write(u6,*) 'GenerateB: broken symmetry!'
          write(u6,*) '  Component: ',iComp
          write(u6,*) '  i and j  : ',i,j
          write(u6,*) '  Dij      : ',Work(kij)
          write(u6,*) '  Dji      : ',Work(kji)
          write(u6,*) '  Diff.    : ',Tst
          call SysAbendMsg('GenerateB','Broken symmetry!',' ')
        end if
      end do
    end do
  end do
end if

end subroutine GenerateB
