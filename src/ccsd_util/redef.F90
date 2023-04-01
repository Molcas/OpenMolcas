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

use Para_Info, only: nProcs
implicit none
#include "parallel.fh"
! help variables
integer i, ii
!LD integer i,ii,rc
real*8 tmin, eff, tabtot, tdisp, tdisptot, tidletot
real*8 tminab, tdole

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
tminab = 0.0d0

do i=2,nProcs
  if (tmin > ididle(i)) then
    tmin = ididle(i)
  end if
end do

do i=1,nProcs
  ididle(i) = ididle(i)-tmin
  if (tminab < idtmab(i)) then
    tminab = idtmab(i)
  end if
end do

!3 calc total time, used in prev. iteration for sumoverab process
!  + total idle time in prev. iter. (on nodes dedicated to ab process)

tabtot = 0.0d0
tidletot = 0.0d0
do ii=1,nprocab
  i = idab(ii)+1
  if (ideffab(ii) > 0.0d0) then
    tabtot = tabtot+idtmab(i)
    tidletot = tidletot+ididle(i)
    if (tminab > idtmab(i)) then
      tminab = idtmab(i)
    end if
  else
    tidletot = tidletot+ididle(i)
  end if
end do

tdole = tidletot/nprocab
!? if (tdole > tminab) tdole = tminab

!4 calc new redistribution

tdisptot = 0.0d0
do ii=1,nprocab
  i = idab(ii)+1
  tdisp = idtmab(i)+ididle(i)-tdole
  if (tdisp < 0.0d0) tdisp = 0.0d0
  if (ideffab(ii) == 0.0d0) then
    eff = 1.0d0
  else
    eff = ideffab(ii)/(idtmab(i)/tabtot)
  end if
  !Stare tdisp = tdisp*eff
  tdisptot = tdisptot+tdisp
end do
write(6,*) 'Tab   ',tabtot
write(6,*) 'Tidle ',tidletot
write(6,*) 'Tdisp ',tdisptot
write(6,*) 'Tddole',tdole
write(6,*) 'Tminab',tminab

do ii=1,nprocab
  i = idab(ii)+1
  tdisp = idtmab(i)+ididle(i)-tdole
  if (tdisp < 0.0d0) tdisp = 0.0d0
  if (ideffab(ii) == 0.0d0) then
    eff = 1.0d0
  else
    eff = ideffab(ii)/(idtmab(i)/tabtot)
  end if
  write(6,*) ii,idtmab(i),ideffab(ii)
  write(6,*) eff,tdisp

  !Stare tdisp = tdisp*eff
  ideffab(ii) = tdisp/tdisptot

  if (ideffab(ii) <= 0.02) then
    ideffab(ii) = 0.0d0
  end if

end do

!5 renormalization of ideffab

tabtot = 0.0d0
do ii=1,nprocab
  tabtot = tabtot+ideffab(ii)
end do
do ii=1,nprocab
  ideffab(ii) = ideffab(ii)/tabtot
  write(6,*) ii,ideffab(ii)
end do

ideffab(1) = 0.116904633172297
ideffab(2) = 0.129270185505803
ideffab(3) = 0.140060191767431
ideffab(4) = 0.120813846670906
ideffab(5) = 0.086763031703814
ideffab(6) = 0.173676115414579
ideffab(7) = 0.232511995765169

return

end subroutine redef
