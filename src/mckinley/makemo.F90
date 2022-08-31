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
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************

subroutine MakeMO(AOInt,Temp,nTemp,n_Int,iCmp,iCmpa,ibasi,jbasj,kbask,lbasl,nGr,Indx,moip,naco,nop,indgrd,ishll, &
                  ishell,rmoin,nmoin,iuvwx,iaost,buffer,ianga)
! this is the driver for the two index transformation
! it is not very efficent, but on the other hand it
! usually to take more than a few percent of the total
! CPU time, if someone notice something else I will
! rewrite it, in the mean time, dont worry.

use Basis_Info, only: Shells
use Symmetry_Info, only: nIrrep
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nTemp, n_Int, icmp(4), iCmpa(4), ibasi, jbasj, kbask, lbasl, nGr, Indx(3,4), moip(0:7), naco, &
                                 nop(4), indgrd(3,4,0:nirrep-1), ishll(4), ishell(4), nmoin, iuvwx(4), iAOST(4), ianga(4)
real(kind=wp), intent(in) :: AOInt(n_Int)
real(kind=wp), intent(out) :: Temp(nTemp)
real(kind=wp), intent(inout) :: rmoin(nmoin), buffer(*)
integer(kind=iwp) :: ibas(4), iCar, iCent, iCnt, iGr, ii, iIrrep, iMax, ip, ip0, ip1, ip2, ip5, ipc, ipck, ipcl, ipFin, mSum, &
                     nabcd, nCk, nCl, nij, nijkl, nkl, nScrtch, ntot
logical(kind=iwp) :: lc, pert(0:7)

iMax = 0
mSum = 0
nabcd = iBasi*jBasj*kBask*lBasl
nijkl = icmp(1)*icmp(2)*icmp(3)*icmp(4)
ntot = nabcd*nijkl
iBas(1) = iBasi
iBas(2) = jBasj
iBas(3) = kBask
iBas(4) = lBasl
do ii=1,4
  imax = max(iBas(ii)*iCmp(ii),imax)
  mSum = mSum+iBas(ii)*iCmp(ii)
end do
imax = max(iMax,nAco)

ip = 1
ip0 = ip
ip = ip+nGr*ntot
ip1 = ip
nScrtch = imax**4
ip = ip+nScrtch
ip2 = ip
ip = ip+nScrtch
ip = ip+nScrtch ! some unused memory here?
ip5 = ip
ip = ip+iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4)*iBas(1)*iBas(2)*iBas(3)*iBas(4)
if (ip-1 > nTemp) then
  write(u6,*) 'MakeMO: ip-1 > nTemp'
  write(u6,*) 'ip,nTemp=',ip,nTemp
  call Abend()
end if
!ip = 2
!Temp(ip-1) = Zero
!ip0 = ip
!ip = ip+nGr*ntot+1
!Temp(ip-1) = Zero
!ip1 = ip
!nScrtch = imax**4+1
!ip = ip+nScrtch
!Temp(ip-1) = Zero
!ip2 = ip
!ip = ip+nScrtch
!Temp(ip-1) = Zero
!ip3 = ip
!ip = ip+nScrtch
!Temp(ip-1) = Zero
!ip5 = ip
!ip = ip+iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4)*iBas(1)*iBas(2)*iBas(3)*iBas(4)+1
!Temp(ip-1) = Zero

ipc = 1
!ipD = ipc
!nD = nACO**4
!ipC = ipC+nd
!ipci = ipc
!nCi = iBas(1)*iCmp(1)*nACO
!ipcj = ipc
!nCj = iBas(2)*iCmp(2)*nACO
ipck = ipc
nCk = iBas(3)*iCmp(3)*nACO
ipc = ipc+nCk
ipcl = ipc
nCl = iBas(4)*iCmp(4)*nAcO
ipc = ipc+nCl
nij = iCmp(1)*iBas(1)*iBas(2)*iCmp(2)
nkl = iCmp(3)*iBas(3)*iBas(4)*iCmp(4)
if (ipc-1 /= nMoIn) then
  write(u6,*) 'MakeMO: ipc-1 /= nMoIn'
  write(u6,*) 'ipc,nMoIn=',ipc,nMoIn
  call Abend()
end if

call Sort_mck(AOInt,Temp(ip0),iBas(1),iBas(2),iBas(3),iBas(4),iCmp(1),iCmp(2),iCmp(3),iCmp(4),iBas(1),iBas(2),iBas(3),iBas(4), &
              iCmpa(1),iCmpa(2),iCmpa(3),iCmpa(4),nGr,nop,ianga,ishll)

do iCent=1,4
  lc = .false.
  do iCar=1,3
    pert(0:nIrrep-1) = .false.
    lC = .false.
    do iIrrep=0,nIrrep-1
      if (IndGrd(iCar,iCent,iIrrep) /= 0) pert(iIrrep) = .true.
      if (IndGrd(iCar,iCent,iIrrep) /= 0) lC = .true.
    end do
    if (lc) then
      if (Indx(iCar,iCent) > 0) then

        iGr = Indx(icar,icent)
        call MOAcc(Temp(ip0+(iGr-1)*ntot),Temp(ip1),Temp(ip2),nScrtch,ishell,rmoin(ipCk),nCk,rmoin(ipCl),nCl,Moip,nACO,pert,nOp, &
                   ibas,icmpa,iCar,icent,indgrd,real(iuvwx(iCent),kind=wp)/real(nIrrep,kind=wp),iaost,buffer,nij,nkl, &
                   Shells(ishll(1))%nBasis,Shells(ishll(2))%nBasis,icmpa(1),icmpa(2))

      else if (Indx(iCar,iCent) < 0) then
        Temp(ip5:ip5+ntot-1) = Zero
        do iCnt=1,4
          iGr = Indx(iCar,iCnt)
          if (iGr > 0) then
            ipFin = (iGr-1)*ntot+ip0
            Temp(ip5:ip5+ntot-1) = Temp(ip5:ip5+ntot-1)-Temp(ipFin:ipFin+ntot-1)
          end if
        end do
        call MOAcc(Temp(ip5),Temp(ip1),Temp(ip2),nScrtch,ishell,rmoin(ipCk),nCk,rmoin(ipCl),nCl,moip,nACO,pert,nOp,ibas,icmpa, &
                   iCar,icent,indgrd,real(iuvwx(iCent),kind=wp)/real(nIrrep,kind=wp),iaost,buffer,nij,nkl,Shells(ishll(1))%nBasis, &
                   Shells(ishll(2))%nBasis,icmpa(1),icmpa(2))

      end if
    end if
  end do
end do

return

end subroutine MakeMO
