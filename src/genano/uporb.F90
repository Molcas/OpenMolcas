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
!======================================================================*
!                                                                      *
! Author: Per-Olof Widmark                                             *
!         IBM Sweden                                                   *
!                                                                      *
!***********************************************************************

subroutine UpOrb(n,o,w,Orb,Lab)

use Genano_globals, only: MxLqn, iSymBk, pDsym, Center, symlab
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: o, w, Orb(n)
character(len=*), intent(in) :: Lab(n)
integer(kind=iwp) :: i, iBas, iBlk, ind, jBas, indx((MxLqn+1)**2), jndx
real(kind=wp) :: add

indx(:) = 0
do iBas=1,n
  !write(u6,'(a,i5)') '     iBas:',iBas
  if (Lab(iBas)(1:len(Center)) == Center) then
    iBlk = 0
    do i=1,(MxLqn+1)**2
      if (symlab(i) == Lab(iBas)(len(Center)+1:)) iBlk = i
    end do
    !write(u6,'(a,i5)') '       iBlk:',iBlk
    if (iBlk == 0) then
      write(u6,*) 'Unknown basis function: ',Lab(iBas)(len(Center)+1:)
      call Quit_OnUserError()
    end if
    indx(iBlk) = indx(iBlk)+1
    !write(u6,'(a,i5)') '     indx:',indx(iBlk)
    jndx = 0
    do jBas=1,iBas
      !write(u6,'(a,i5)') '       jBas:',jBas
      if (Lab(jBas)(1:len(Center)) == Center) then
        if (Lab(iBas)(len(Center)+1:) == Lab(jBas)(len(Center)+1:)) then
          jndx = jndx+1
          !write(u6,'(a,i5)') '     jndx:',jndx
          ind = jndx+indx(iBlk)*(indx(iBlk)-1)/2+iSymBk(iBlk)-1
          !add = w*o*Cmo(iBas)*Cmo(jBas)
          add = w*o*Orb(iBas)*Orb(jBas)
          pDsym(ind) = pDsym(ind)+add
          !write(u6,'(a,f12.6,a,i5)') 'add',add,' to',ind
        end if
      end if
    end do
  end if
end do

return

end subroutine UpOrb
