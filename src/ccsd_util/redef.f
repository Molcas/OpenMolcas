************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
        subroutine redef
c
c        this routine redefine ideffab vector using infromation about
c        idle time on each nodes selected for sumoverab process
c
        use Para_Info, only: nProcs
        implicit none
#include "parallel.fh"
c
c        help variables
c
        integer i,ii
CLD        integer i,ii,rc
        REAL*8 tmin,eff,tabtot,tdisp,tdisptot,tidletot
        REAL*8 tminab,tdole
c
c0        escape, if nprocab=1, there is nothing to redistribute
c
        if (nprocab.eq.1) return
c
c1        distribute idtmab and ididle to all nodes
c
ctmp        do i=1,nProcs
ctmp    call MPI_BCAST (idtmab(i),1,
ctmp c  MPI_DOUBLE_PRECISION,(i-1),MPI_COMM_WORLD,rc)
ctmp    call MPI_BCAST (ididle(i),1,
ctmp c  MPI_DOUBLE_PRECISION,(i-1),MPI_COMM_WORLD,rc)
ctmp    end do
c        prepis pomocou GA, bcast ekvivalent nevieme
        call gadgop (idtmab(1),nProcs,'+')
        call gadgop (ididle(1),nProcs,'+')
c
c2        def real idle time for all nodes
c
        tmin=ididle(1)
        tminab=0.0d0
c
        do i=2,nProcs
        if (tmin.gt.ididle(i)) then
          tmin=ididle(i)
        end if
        end do
c
        do i=1,nProcs
          ididle(i)=ididle(i)-tmin
          if (tminab.lt.idtmab(i)) then
            tminab=idtmab(i)
          end if
        end do
c
c3        calc total time, used in prev. iteration for sumoverab process
c       + total idle time in prev. iter. (on nodes dedicated to ab
c         process)
c
        tabtot=0.0d0
        tidletot=0.0d0
        do ii=1,nprocab
        i=idab(ii)+1
          if (ideffab(ii).gt.0.0d0) then
            tabtot=tabtot+idtmab(i)
            tidletot=tidletot+ididle(i)
            if (tminab.gt.idtmab(i)) then
              tminab=idtmab(i)
            end if
          else
            tidletot=tidletot+ididle(i)
          end if
        end do
c
        tdole=tidletot/nprocab
c?      if (tdole.gt.tminab) tdole=tminab
c
c4        calc new redistribution
c
        tdisptot=0.0d0
        do ii=1,nprocab
        i=idab(ii)+1
          tdisp=idtmab(i)+ididle(i)-tdole
          if (tdisp.lt.0.0d0) tdisp=0.0d0
          if (ideffab(ii).eq.0.0d0) then
            eff=1.0d0
          else
            eff=ideffab(ii)/(idtmab(i)/tabtot)
          end if
cStare    tdisp=tdisp*eff
          tdisptot=tdisptot+tdisp
        end do
        write (6,*) 'Tab   ',tabtot
        write (6,*) 'Tidle ',tidletot
        write (6,*) 'Tdisp ',tdisptot
        write (6,*) 'Tddole',tdole
        write (6,*) 'Tminab',tminab
c
        do ii=1,nprocab
        i=idab(ii)+1
          tdisp=idtmab(i)+ididle(i)-tdole
          if (tdisp.lt.0.0d0) tdisp=0.0d0
          if (ideffab(ii).eq.0.0d0) then
            eff=1.0d0
          else
            eff=ideffab(ii)/(idtmab(i)/tabtot)
          end if
          write (6,*) ii,idtmab(i),ideffab(ii)
          write (6,*) eff,tdisp
c
cStare    tdisp=tdisp*eff
          ideffab(ii)=tdisp/tdisptot
c
           if (ideffab(ii).le.0.02) then
            ideffab(ii)=0.0d0
           end if
c
        end do
c
c5        renormalization of ideffab
c
        tabtot=0.0d0
        do ii=1,nprocab
             tabtot=tabtot+ideffab(ii)
        end do
        do ii=1,nprocab
             ideffab(ii)=ideffab(ii)/tabtot
          write (6,*) ii,ideffab(ii)
        end do
c

        ideffab(1)=0.116904633172297
        ideffab(2)=0.129270185505803
        ideffab(3)=0.140060191767431
        ideffab(4)=0.120813846670906
        ideffab(5)=0.086763031703814
        ideffab(6)=0.173676115414579
        ideffab(7)=0.232511995765169
c
        return
        end
