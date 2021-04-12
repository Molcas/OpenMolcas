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
        subroutine o2v4ctl (wrk,wrksize,NvGrp,NvSGrp,LunAux)
c
c!      drajver o2v4 procesu
c
        use Para_Info, only: nProcs
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
#include "wrk.fh"
#include "parcc.fh"
#include "chcc_casy.fh"
c
        integer NvGrp,NvSGrp,LunAux
c
c        help variables
        integer NaGrp,NbeGrp,NaSGrp,NbeSgrp
        integer mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe
c##
        integer aGrp,bGrp,proc,i,j
        integer nJobs,addJobs,actJobs
c
c
c1      Inicializacia premennych (Predbezna)
c
c@c
c        return
c@c
        NaGrp=NvGrp
        NbeGrp=NvGrp
        NaSGrp=NvSGrp
        nbeSGrp=NvSGrp
c
c
c2      define all groups and ssungroup parameters
c
        call DefParo2v4 (NaGrp,NbeGrp,NaSGrp,NbeSgrp,
     c                   mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe)
        if (printkey.ge.10) then
        write (6,*) NaGrp,NbeGrp,NaSGrp,NbeSgrp,
     c                   mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe
        end if
cmp
        Call CWTime(TCpu,TWall)
        if (printkey.gt.1) then
        write (6,*)
        write (6,'(A,f18.1)') ' Cpu last call [s] = ',
     & TCpu-TCpu_l
        write (6,'(A,f18.1)') 'Wall last call [s] = ',
     & TWall-TWall_l
        write (6,*)
        write (6,'(A,f18.1)') 'Total Cpu  [s] = ',
     & TCpu
        write (6,'(A,f18.1)') 'Total Wall [s] = ',
     & TWall-TWall0
        write (6,'(A,f18.2)') 'TCpu/TWall [%] = ',
     & 100.0d0*TCpu/(TWall-TWall0)
        write (6,*)
        end if
        TCpu_l=TCpu
        TWall_l=TWall
cmp
c
c3        distribute work among nodes (def ABID)
c
        do proc=0,(nProcs-1)
        do aGrp=1,NaGrp
        do bGrp=1,NaGrp
          ABID(proc,aGrp,bGrp)=0
        end do
        end do
        end do
c
c
        if ((nProcs.eq.1).or.(NaGrp.eq.1)) then
c3.1        nP = <1>
c        all jobs to one node
c
          do aGrp=1,NaGrp
          do bGrp=1,NaGrp
            ABID(0,aGrp,bGrp)=1
          end do
          end do
c
c
        else
c3.2        nP = <2, inf)
c        N.B. cele zatial trochu odflaknute, lebo sa neberie do uvahy,
c        ze pre half joby je lepsie ked nadvazuju na plne joby s
c        rovnakymi indexami
c
c3.2.1        Full Jobs
c
          i=(NaGrp*(NaGrp-1))/2
          nJobs=int(i/nProcs)
          addJobs=mod(i,nProcs)
c
c          first nodes: 0-addJobs
c             will have int(N'(N'-1)/2 /nProc)+1 Full Jobs
c          rest nodes: addJobs-nProcs-1
c               will have int(N'(N'-1)/2 /nProc)
c
          proc=0
          actJobs=nJobs
          do aGrp=2,NaGrp
          do bGrp=1,aGrp-1
            ABID(proc,aGrp,bGrp)=1
            actJobs=actJobs-1
            if (actJobs.eq.-1) then
              proc=proc+1
              actJobs=nJobs
            else if (actJobs.eq.0) then
              if (addJobs.gt.0) then
                addJobs=addJobs-1
              else
                proc=proc+1
                actJobs=nJobs
              end if
            end if
          end do
          end do
c
c
c3.2.2        Half Jobs
c
c        distribution of Half Jobs
c          - first distribution - addjobs-nProcs-1
c          - second distribution - addjobs-nProcs-1
c          - other distributions - 0 - nProcs-1
c
          addJobs=mod(i,nProcs)
          proc=addJobs
          j=0
c
          do aGrp=1,NaGrp
            ABID(proc,aGrp,aGrp)=1
            proc=proc+1
            if (proc.eq.nProcs) then
              if (j.eq.0) then
c              first distribution - addjobs-nProcs-1
                proc=addJobs
                j=1
              else if (j.eq.1) then
c              second distribution - addjobs-nProcs-1
                proc=addJobs
                j=2
              else
c              other distributions - 0 - nProcs-1
                proc=0
              end if
            end if
          end do
c
        end if
c
c
c        Printing ABID
        do proc=0,nProcs-1
        if (printkey.ge.10) then
        write (6,*)   ' For myRank = ',proc
        end if
          do aGrp=1,NaGrp
          do bGrp=1,aGrp
          if (ABID(proc,aGrp,bGrp).eq.1) then
          if (printkey.ge.10) then
          write (6,*) '    aGrp,bGrp ',aGrp,bGrp
          end if
          end if
          end do
          end do
        end do
c
c
c
c4      A ideme na to
c
        call o2v4 (wrk(1),wrksize,
     c             NaGrp,NbeGrp,NaSGrp,NbeSgrp,
     c             mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe,
     c             LunAux)

c
        return
        end
