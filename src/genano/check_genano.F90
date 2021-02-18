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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine perform some consistency checks.                        *
!                                                                      *
!======================================================================*
!                                                                      *
! Author: Per-Olof Widmark                                             *
!         IBM Sweden                                                   *
!                                                                      *
!***********************************************************************

subroutine Check_genano()

use Genano_globals, only: MxLqn, nSym, nBas, nPrim, pDsym, Center, BasName, symlab
use Constants, only: Zero
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: iBas, ind, iSym, l, icnt(0:MxLqn)

icnt(:) = 0
ind = 0
do iSym=1,nSym
  !write(u6,'(a,i1,a)') ' <<< Symmetry ',iSym,' >>>'
  do iBas=1,nBas(iSym)
    ind = ind+1
    !write(u6,*) BasName(1,ind),BasName(2,ind)
    do l=0,MxLqn
      if (BasName(ind)(1:len(Center)) == Center) then
        if (BasName(ind)(len(Center)+1:) == symlab(l*(l+1)+1)) then
          icnt(l) = icnt(l)+1
        end if
      end if
    end do
  end do
end do
do l=0,MxLqn
  if (icnt(l) /= nPrim(l)) then
    write(u6,*) 'Number of primitives do not match!'
    write(u6,'(1x,a,2i5)') symlab(l*(l+1)+1),nPrim(l),icnt(l)
    call Quit_OnUserError()
  end if
end do
pDsym(:) = Zero

return

end subroutine Check_genano
