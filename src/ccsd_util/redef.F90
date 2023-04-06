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

subroutine redef()
! this routine redefines ideffab vector using information about
! idle time on each nodes selected for sumoverab process

use ccsd_global, only: idab, ideffab, ididle, idtmab, nprocab
use Para_Info, only: nProcs
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, ii
real(kind=wp) :: eff, tabtot, tdisp, tdisptot, tdole, tidletot, tmin, tminab

!0 escape, if nprocab=1, there is nothing to redistribute

if (nprocab == 1) return

!1 distribute idtmab and ididle to all nodes

!tmp do i=1,nProcs
!tmp   call MPI_BCAST(idtmab(i),1,MPI_DOUBLE_PRECISION,(i-1),MPI_COMM_WORLD,rc)
!tmp   call MPI_BCAST(ididle(i),1,MPI_DOUBLE_PRECISION,(i-1),MPI_COMM_WORLD,rc)
!tmp end do
! prepis pomocou GA, bcast ekvivalent nevieme
call gadgop(idtmab(1),nProcs,'+')
call gadgop(ididle(1),nProcs,'+')

!2 def real idle time for all nodes

tmin = ididle(1)
tminab = Zero

do i=2,nProcs
  if (tmin > ididle(i)) tmin = ididle(i)
end do

ididle(1:nProcs) = ididle(1:nProcs)-tmin

do i=1,nProcs
  if (tminab < idtmab(i)) tminab = idtmab(i)
end do

!3 calc total time, used in prev. iteration for sumoverab process
!  + total idle time in prev. iter. (on nodes dedicated to ab process)

tabtot = Zero
tidletot = Zero
do ii=1,nprocab
  i = idab(ii)+1
  if (ideffab(ii) > Zero) then
    tabtot = tabtot+idtmab(i)
    tidletot = tidletot+ididle(i)
    if (tminab > idtmab(i)) tminab = idtmab(i)
  else
    tidletot = tidletot+ididle(i)
  end if
end do

tdole = tidletot/nprocab
!? if (tdole > tminab) tdole = tminab

!4 calc new redistribution

tdisptot = Zero
do ii=1,nprocab
  i = idab(ii)+1
  tdisp = idtmab(i)+ididle(i)-tdole
  if (tdisp < Zero) tdisp = Zero
  if (ideffab(ii) == Zero) then
    eff = One
  else
    eff = ideffab(ii)/(idtmab(i)/tabtot)
  end if
  !Stare tdisp = tdisp*eff
  tdisptot = tdisptot+tdisp
end do
write(u6,*) 'Tab   ',tabtot
write(u6,*) 'Tidle ',tidletot
write(u6,*) 'Tdisp ',tdisptot
write(u6,*) 'Tddole',tdole
write(u6,*) 'Tminab',tminab

do ii=1,nprocab
  i = idab(ii)+1
  tdisp = idtmab(i)+ididle(i)-tdole
  if (tdisp < Zero) tdisp = Zero
  if (ideffab(ii) == Zero) then
    eff = One
  else
    eff = ideffab(ii)/(idtmab(i)/tabtot)
  end if
  write(u6,*) ii,idtmab(i),ideffab(ii)
  write(u6,*) eff,tdisp

  !Stare tdisp = tdisp*eff
  ideffab(ii) = tdisp/tdisptot

  if (ideffab(ii) <= 0.02_wp) ideffab(ii) = Zero

end do

!5 renormalization of ideffab

tabtot = sum(ideffab(1:nprocab))
ideffab(1:nprocab) = ideffab(1:nprocab)/tabtot
do ii=1,nprocab
  write(u6,*) ii,ideffab(ii)
end do

! what are these numbers?
! what's the point of the code above?
ideffab(1) = 0.116904633172297_wp
ideffab(2) = 0.129270185505803_wp
ideffab(3) = 0.140060191767431_wp
ideffab(4) = 0.120813846670906_wp
ideffab(5) = 0.086763031703814_wp
ideffab(6) = 0.173676115414579_wp
ideffab(7) = 0.232511995765169_wp

return

end subroutine redef
