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

subroutine o3v3ctl(wrk,wrksize,NvGrp,LunAux)
! docasny drajver o3v3 procesov

use chcc_global, only: BeAID, BetaID, printkey
use Para_Info, only: nProcs
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: wrksize, NvGrp, LunAux
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp) :: actJobs, actNode, addJobs, aGrp, beGrp, maxdim, nJobs, proc

!1 Def parameters for o3v3 processes

call DefParo3v3(NvGrp,maxdim)

!2 Distribute work among nodes (def BetaID, BeAID)

!2.1 vanish BetaID, BeAID
BetaID(0:nProcs-1,1:NvGrp) = 0
BeAID(0:nProcs-1,1:NvGrp,1:NvGrp) = 0

if (nProcs == 1) then
  !2.2.1 single node

  BetaID(0,1:NvGrp) = 1
  BeAID(0,1:NvGrp,1:NvGrp) = 1

else
  !2.2.2 multi node
  ! (N.B. trochu odflaknute, dalo by sa este zohladnit ze tie nody,
  ! ktore maju o jeden job navyse nebudu tie, kde je BetaID=1,
  ! resp. nebudu tie, kde sa realizuje X0.1 prispevok a pod)

  nJobs = int((NvGrp*NvGrp)/nProcs)
  addJobs = mod((NvGrp*NvGrp),nProcs)

  actNode = 0
  actJobs = nJobs

  do beGrp=1,NvGrp
    BetaID(actNode,beGrp) = 1
    if (printkey >= 10) write(u6,*) 'BetaID',actnode,beGrp
    do aGrp=1,NvGrp
      BeAID(actNode,beGrp,aGrp) = 1
      if (printkey >= 10) write(u6,*) 'BeAID',actnode,beGrp,aGrp
      actJobs = actJobs-1
      if (actJobs == -1) then
        actNode = actNode+1
        actJobs = nJobs
      else if (actJobs == 0) then
        if (addJobs > 0) then
          addJobs = addJobs-1
        else
          actNode = actNode+1
          actJobs = nJobs
        end if
      end if
    end do
  end do

end if
!@@
if (printkey >= 10) then
  do proc=0,nProcs-1
    do beGrp=1,NvGrp
      write(u6,99) proc,beGrp,(BeAID(proc,beGrp,aGrp),aGrp=1,NvGrp)
    end do
  end do
end if
!@@

!3 A ideme na to

call o3v3jk(wrk,wrksize,NvGrp,maxdim,LunAux)
if (printkey > 1) write(u6,*) ' o3v3jk done'

call o3v3chol(wrk,wrksize,NvGrp,maxdim,LunAux)
if (printkey > 1) write(u6,*) ' o3v3chol done'

call o3v3t2(wrk,wrksize,NvGrp,maxdim,LunAux)
if (printkey > 1) write(u6,*) ' o3v3t2 done'

return

99 format(1x,i3,1x,i2,5x,24(i1,1x))

end subroutine o3v3ctl
