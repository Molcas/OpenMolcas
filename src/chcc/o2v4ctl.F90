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

subroutine o2v4ctl(wrk,wrksize,NvGrp,NvSGrp,LunAux)
! drajver o2v4 procesu

use Index_Functions, only: nTri_Elem
use chcc_global, only: ABID, printkey, TCpu, TCpu_l, TWall, TWall_l, TWall0
use Para_Info, only: nProcs
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: wrksize, NvGrp, NvSGrp, LunAux
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp) :: actJobs, addJobs, aGrp, bGrp, i, j, mdGrpa, mdGrpbe, mdSGrpa, mdSGrpbe, NaGrp, NaSGrp, NbeGrp, NbeSgrp, &
                     nJobs, proc

!1 Inicializacia premennych (Predbezna)

!@c
! return
!@c
NaGrp = NvGrp
NbeGrp = NvGrp
NaSGrp = NvSGrp
nbeSGrp = NvSGrp

!2 define all groups and ssungroup parameters

call DefParo2v4(NaGrp,NbeGrp,NaSGrp,NbeSgrp,mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe)
if (printkey >= 10) write(u6,*) NaGrp,NbeGrp,NaSGrp,NbeSgrp,mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe
!mp
call CWTime(TCpu,TWall)
if (printkey > 1) then
  write(u6,*)
  write(u6,'(A,f18.1)') ' Cpu last call [s] = ',TCpu-TCpu_l
  write(u6,'(A,f18.1)') 'Wall last call [s] = ',TWall-TWall_l
  write(u6,*)
  write(u6,'(A,f18.1)') 'Total Cpu  [s] = ',TCpu
  write(u6,'(A,f18.1)') 'Total Wall [s] = ',TWall-TWall0
  write(u6,'(A,f18.2)') 'TCpu/TWall [%] = ',100.0d0*TCpu/(TWall-TWall0)
  write(u6,*)
end if
TCpu_l = TCpu
TWall_l = TWall
!mp

!3 distribute work among nodes (def ABID)

ABID(0:nProcs-1,1:NaGrp,1:NaGrp) = 0

if ((nProcs == 1) .or. (NaGrp == 1)) then
  !3.1 nP = <1>
  ! all jobs to one node

  ABID(0,1:NaGrp,1:NaGrp) = 1

else
  !3.2 nP = <2, inf)
  ! N.B. cele zatial trochu odflaknute, lebo sa neberie do uvahy,
  ! ze pre half joby je lepsie ked nadvazuju na plne joby s
  ! rovnakymi indexami

  !3.2.1 Full Jobs

  i = nTri_Elem(NaGrp-1)
  nJobs = int(i/nProcs)
  addJobs = mod(i,nProcs)

  ! first nodes: 0-addJobs
  !    will have int(N'(N'-1)/2 /nProc)+1 Full Jobs
  ! rest nodes: addJobs-nProcs-1
  !      will have int(N'(N'-1)/2 /nProc)

  proc = 0
  actJobs = nJobs
  do aGrp=2,NaGrp
    do bGrp=1,aGrp-1
      ABID(proc,aGrp,bGrp) = 1
      actJobs = actJobs-1
      if (actJobs == -1) then
        proc = proc+1
        actJobs = nJobs
      else if (actJobs == 0) then
        if (addJobs > 0) then
          addJobs = addJobs-1
        else
          proc = proc+1
          actJobs = nJobs
        end if
      end if
    end do
  end do

  !3.2.2 Half Jobs
  !
  ! distribution of Half Jobs
  !   - first distribution - addjobs-nProcs-1
  !   - second distribution - addjobs-nProcs-1
  !   - other distributions - 0 - nProcs-1

  addJobs = mod(i,nProcs)
  proc = addJobs
  j = 0

  do aGrp=1,NaGrp
    ABID(proc,aGrp,aGrp) = 1
    proc = proc+1
    if (proc == nProcs) then
      if (j == 0) then
        ! first distribution - addjobs-nProcs-1
        proc = addJobs
        j = 1
      else if (j == 1) then
        ! second distribution - addjobs-nProcs-1
        proc = addJobs
        j = 2
      else
        ! other distributions - 0 - nProcs-1
        proc = 0
      end if
    end if
  end do

end if

! Printing ABID
if (printkey >= 10) then
  do proc=0,nProcs-1
    write(u6,*) ' For myRank = ',proc
    do aGrp=1,NaGrp
      do bGrp=1,aGrp
        if (ABID(proc,aGrp,bGrp) == 1) write(u6,*) '    aGrp,bGrp ',aGrp,bGrp
      end do
    end do
  end do
end if

!4 A ideme na to

call o2v4(wrk,wrksize,NaGrp,NbeGrp,NaSGrp,NbeSgrp,mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe,LunAux)

return

end subroutine o2v4ctl
