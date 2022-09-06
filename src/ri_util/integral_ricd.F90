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

subroutine Integral_RICD(iCmp,iShell,MapOrg,iBas,jBas,kBas,lBas,kOp,Shijij,IJeqKL,iAO,iAOst,ijkl,AOInt,SOInt,nSOint,iSOSym,nSkal, &
                         nSOs,TInt,nTInt,iTOffs,nSym)

implicit real*8(A-H,O-Z)
real*8 AOInt(*), SOInt(*), TInt(nTInt)
integer iCmp(4), iShell(4), iAO(4), iAOst(4), kOp(4), iSOSym(2,nSOs), iTOffs(0:7,0:7,0:7), MapOrg(4)
logical Shijij, IJeqKL

if (nSym == 1) then
  call PLF_RICD(AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),iShell,iAO,iAOst,Shijij .and. IJeqKL,iBas,jBas,kBas,lBas,kOp,TInt, &
                iTOffs(1,0,0),iTOffs(2,0,0),iTOffs(0,0,0),iTOffs(3,0,0),iTOffs(4,0,0))
else
  write(6,*) 'Integral_RICD: fatal error!'
  call Abend()
end if

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer_array(MapOrg)
  call Unused_real_array(SOInt)
  call Unused_integer(nSOint)
  call Unused_integer_array(iSOSym)
  call Unused_integer(nSkal)
  call Unused_integer(nSym)
end if

end subroutine Integral_RICD
