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
c
c        this package contains following files
c        (mainly for parallel stuff)
c
c        set0
c        joinamplitudes
c        reajalovy
c        distnodes
c       sumabdistt
c       redef
c
c        ----------------------------------------------
c
        subroutine set0 (wrk,wrksize,
     c                   mapd,mapi)
c
c        this routine vanish given mediate
c
        implicit none
c
#include "wrk.fh"
#include "paralell.fh"
c
       integer mapd(0:512,1:6)
       integer mapi(1:8,1:8,1:8)
c
c        help variables
c
       integer poss0,lenght,ii
c
c1        def poss0, legth
        poss0=mapd(1,1)
        ii=mapd(0,5)
        lenght=mapd(ii,1)+mapd(ii,2)-poss0
c
c2        set appropriate wrk = 0
        call mv0zero (lenght,lenght,wrk(poss0))
c
        return
c Avoid unused argument warnings
      if (.false.) call Unused_integer_array(mapi)
        end
c
c        ----------------------------------------------
c
        subroutine joinamplitudes (wrk,wrksize)
c
c        this routine join contributions to amplitudes (Tn) from all
c        processors. Sum will be placed back to Tn on all machines
c
c        N.B. FREE ZONE in purpose of this routine is the space of
c       free working files V1-V4, .....
c        (free zone is not used in GA)
c        Since T2 amplitudes are no more then oovv, at most V1-V3 space
c        will be demaged. (acutally V1 and V2 only)
c
        implicit none
c
#include "wrk.fh"
#include "ccsd1.fh"
#include "ccsd2.fh"
#include "paralell.fh"
c
c        help variables
c
        integer ii,lenght
CLD        integer ii,lenght,rc,i
c
        if (nProcs.eq.1) return
c
c1        join t13,t14
c
c1.1    calc overall length of t13 and t14 (they must be one after the other)
        ii=mapdt14(0,5)
        lenght=mapdt14(ii,1)+mapdt14(ii,2)-posst130
c
c1.2        vanish required part in free zone
c        call mv0zero (lenght,lenght,wrk(possv10))
c
c1.3        allreduce t13 and t14 into free zone
c        call MPI_ALLREDUCE (wrk(posst130),wrk(possv10),lenght,
c    c  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,rc)
c
c1.4        put joined t13,t14 back to t13,t14 place from free zone
c        do i=0,lenght-1
c        wrk(posst130+i)=wrk(possv10+i)
c        end do
c
c1.om   allreduge t13,t14 together
        call gadgop (wrk(posst130),lenght,'+')
c
c
c2        join t21,t22,t23
c
c2.1    calc overall length of t21-t23 (they must be one after the other)
        ii=mapdt23(0,5)
        lenght=mapdt23(ii,1)+mapdt23(ii,2)-posst210
c
c2.2        vanish required part in free zone
c        call mv0zero (lenght,lenght,wrk(possv10))
c
c2.3        allreduce t13 and t14 into free zone
c        call MPI_ALLREDUCE (wrk(posst210),wrk(possv10),lenght,
c    c  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,rc)
c
c2.4        put joined t21-t23 back to t21-t23 place from free zone
c        do i=0,lenght-1
c        wrk(posst210+i)=wrk(possv10+i)
c        end do
c
c2.om   allreduge t21-t23 together
        call gadgop (wrk(posst210),lenght,'+')
c
        return
        end
c
c        ----------------------------------------------
c
       subroutine reajalovy (lun,lenght,vector)
c
c     this routine read blank card
c     with number lun form the given possition and update pointers
c
c     lun    - Logical unit number of file, where mediate is stored (Input)
c     lenght - # of R8 numbers to be read  (Input)
c     vector - space, where numbers are stored after reading  (Output)

c
#include "filemgr.fh"
#include "ccsd1.fh"

#include "SysDef.fh"
c
       integer lun,lenght
       real*8 vector(1:1)
c
       if (iokey.eq.1) then
c      Fortran IO
       read (lun)
c
       else
c      MOLCAS IO
       call ddafile (lun,0,vector,lenght,daddr(lun))
       end if
c
       return
       end
c
c     ----------------------------
c
        subroutine distnodes
c
c        this routine distribute nodes to different parts
c
        implicit none
#include "paralell.fh"
c
        integer i
        REAL*8 efftot
c
ctmp ta zatial takto
        if (nProcs.eq.1) then
c
cI        def nodes for sumoverab
        nprocab=1
        idab(1)=0
        ideffab(1)=1.0
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=0
        idbaab=0
        idbbaa=0
        idbbbb=0
        idaabb=0
        idabba=0
c
cIII        def node for finale
        idfin=0
c
c
        else if (nProcs.eq.2) then
c
cI        def nodes for sumoverab
        nprocab=1
        idab(1)=0
        ideffab(1)=1.0
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=1
        idbaab=1
        idbbaa=1
        idbbbb=1
        idaabb=1
        idabba=1
c
cIII        def node for finale
        idfin=1
c
c
        else if (nProcs.eq.3) then
c
cI        def nodes for sumoverab
        nprocab=1
        idab(1)=0
        ideffab(1)=1.0
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=1
        idbaab=1
        idbbaa=1
        idbbbb=2
        idaabb=2
        idabba=2
c
cIII        def node for finale
        idfin=1
c
c
        else if (nProcs.eq.4) then
c
cI        def nodes for sumoverab
        nprocab=4
        idab(1)=0
        idab(2)=1
        idab(3)=2
        idab(4)=3
        ideffab(1)=0.25
        ideffab(2)=0.25
        ideffab(3)=0.25
        ideffab(4)=0.25
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=0
        idbaab=1
        idbbaa=1
        idbbbb=2
        idaabb=2
        idabba=3
c
cIII        def node for finale
        idfin=3
c
c
        else if (nProcs.eq.5) then
c
cI        def nodes for sumoverab
        nprocab=1
        idab(1)=0
        ideffab(1)=1.0
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=1
        idbaab=1
        idbbaa=2
        idbbbb=3
        idaabb=3
        idabba=4
c
cIII        def node for finale
        idfin=2
c
c
        else if (nProcs.eq.6) then
c
cI        def nodes for sumoverab
        nprocab=6
        idab(1)=0
        idab(2)=1
        idab(3)=2
        idab(4)=3
        idab(5)=4
        idab(6)=5
        ideffab(1)=1.0
        ideffab(2)=1.0
        ideffab(3)=1.0
        ideffab(4)=1.0
        ideffab(5)=1.0
        ideffab(6)=1.0
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=0
        idbaab=1
        idbbaa=2
        idbbbb=3
        idaabb=4
        idabba=5
c
cIII        def node for finale
        idfin=3
c
        else if (nProcs.eq.10) then
c
cI        def nodes for sumoverab
        nprocab=4
          idab(1)=0
          idab(2)=1
          idab(3)=2
          idab(4)=3
          ideffab(1)=1.0
          ideffab(2)=1.0
          ideffab(3)=1.0
          ideffab(4)=1.0
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=4
        idbaab=5
        idbbaa=6
        idbbbb=7
        idaabb=8
        idabba=9
c
cIII        def node for finale
        idfin=5
c
        else
c
cI        def nodes for sumoverab
        nprocab=nProcs
        do i=1,nprocab
          idab(i)=i-1
          ideffab(i)=1.0
        end do
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=0
        idbaab=1
        idbbaa=2
        idbbbb=3
        idaabb=4
        idabba=5
c
cIII        def node for finale
        idfin=6

        end if
c
         return
c
ctmp         koniec tmp riesenia
c
c
c
        if (nProcs.eq.1) then
c
cI        def nodes for sumoverab
        nprocab=1
        idab(1)=0
        ideffab(1)=1.0
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=0
        idbaab=0
        idbbaa=0
        idbbbb=0
        idaabb=0
        idabba=0
c
cIII        def node for finale
        idfin=0
c
c
        else if (nProcs.eq.2) then
c
cI        def nodes for sumoverab
        nprocab=2
        idab(1)=0
        idab(2)=1
        ideffab(1)=0.5
        ideffab(2)=0.5
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=0
        idbaab=0
        idbbaa=0
        idbbbb=1
        idaabb=1
        idabba=1
c
cIII        def node for finale
        idfin=0
c
c
        else if (nProcs.eq.3) then
c
cI        def nodes for sumoverab
        nprocab=3
        idab(1)=0
        idab(2)=1
        idab(3)=2
        ideffab(1)=0.333
        ideffab(2)=0.333
        ideffab(3)=0.333
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=1
        idbaab=1
        idbbaa=1
        idbbbb=2
        idaabb=2
        idabba=2
c
cIII        def node for finale
        idfin=0
c
c
        else if (nProcs.eq.4) then
c
cI        def nodes for sumoverab
        nprocab=2
        idab(1)=2
        idab(2)=3
        ideffab(1)=1.0
        ideffab(2)=1.0
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=0
        idbaab=0
        idbbaa=0
        idbbbb=1
        idaabb=1
        idabba=1
c
cIII        def node for finale
        idfin=0
c
c
        else
c
cI        def nodes for sumoverab
c
        nprocab=nProcs-2
        idab(1)=1
        idab(2)=2
        idab(3)=4
        idab(4)=5
        idab(5)=6
        idab(6)=7
        do i=1,nprocab
        ideffab(i)=1.0d0/nprocab
        end do
        ideffab(1)=ideffab(1)/2
        ideffab(2)=ideffab(2)/2
c
cII        def nodes for sumoverb and intmezzo
        idaaaa=1
        idbaab=0
        idbbaa=0
        idbbbb=2
        idaabb=3
        idabba=3
c
cIII        def node for finale
        idfin=1
c
        end if
c
c        renormalize ideffab
c
        efftot=0.0d0
        do i=1,nprocab
        efftot=efftot+ideffab(i)
        end do
c
        do i=1,nprocab
        ideffab(i)=ideffab(i)/efftot
        end do

c
        return
        end
c
c     ----------------------------
c
        subroutine redef
c
c        this routine redefine ideffab vector using infromation about
c        idle time on each nodes selected for sumoverab process
c
        implicit none
#include "paralell.fh"
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
c
c     ----------------------------
c
        subroutine sumabdistt (n,idtot)
c
c        this routine distribute work for n records among
c        nprocab processors with the frequency, corresponding
c        to ideffab values
c
c        n - # of records to be distributed (I)
c        idtot - distribution vector (O)
c               (idtot(i) -  # of records to be realized by i-th node)
c
        implicit none
#include "paralell.fh"
c
        integer n
        integer idtot(1)
c
c        help parameters
c
        integer i,ntot,max,imax
CLD        integer i,j,ntot,max,imax
        real*8 sum
c
c
c1        distribute recordsc according to eff. coeficients
c
        sum=0.0d0
        do i=1,nprocab
        sum=sum+ideffab(i)
        end do
c
        do i=1,nprocab
        idtot(i)=int(((ideffab(i)*n)/sum)+0.5d0)
        end do
c
c2        do corrections, if roundoff errors caused some diferences
c
1        ntot=0
        do i=1,nprocab
        ntot=ntot+idtot(i)
        end do

        if (ntot.gt.n) then
c        ubrat treba (z najvacsieho dielu)
          max=idtot(1)
          imax=1
          do i=1,nprocab
          if (max.lt.idtot(i)) then
          max=idtot(i)
          imax=i
          end if
          end do
          idtot(imax)=idtot(imax)-1
          goto 1
        else if (ntot.lt.n) then
c        pridat treba (k najvacsiemu dielu)
          max=idtot(1)
          imax=1
          do i=1,nprocab
          if (max.lt.idtot(i)) then
          max=idtot(i)
          imax=i
          end if
          end do
          idtot(imax)=idtot(imax)+1
          goto 1
        end if
c
        return
        end
c
c     ----------------------------
c
