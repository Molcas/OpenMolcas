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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************

subroutine CalcAXk2(AXk,D1,D2,PUVX,NPUVX,IndPUVX,Off_Act,Off_Orb)

use MCLR_Data, only: ipMat, nDens, nNA
use input_mclr, only: nAsh, nIsh, nOrb, nSym, ntAsh, ntBas
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: AXk(nDens)
integer(kind=iwp), intent(in) :: NPUVX, IndPUVX(ntBas,ntAsh,ntAsh,ntAsh), Off_Act(nSym), Off_Orb(nSym)
real(kind=wp), intent(in) :: D1(ntAsh**2), D2(ntAsh**2), PUVX(NPUVX)
integer(kind=iwp) :: ip, iq, iSym, it, loc1, loc2, p, q, t, u, v, x
real(kind=wp) :: tempa
real(kind=wp), allocatable :: Opu(:)

call mma_allocate(Opu,ntBas*ntAsh)
do p=1,ntBas
  do u=1,ntAsh
    tempa = Zero
    do v=1,ntAsh
      do x=1,ntAsh
        if (IndPUVX(p,u,v,x) /= 0) tempa = tempa+PUVX(IndPUVX(p,u,v,x))*D2((v-1)*nnA+x)
      end do
    end do
    Opu((p-1)*ntAsh+u) = tempa
  end do
end do
do iSym=1,nSym
  do iq=1,nAsh(iSym)
    q = iq+Off_Act(ISym)
    do ip=1,nIsh(iSym)  ! p is inactive
      p = ip+Off_Orb(ISym)
      tempa = Zero
      do it=1,nAsh(ISym)
        t = it+Off_Act(ISym)
        tempa = tempa+(D1((t-1)*nnA+q)+D1((q-1)*nnA+t))*Opu((p-1)*nnA+t)
      end do
      !write(u6,*) 'tempa after sum over t',tempa
      loc1 = ipMat(iSym,iSym)+(iq-1)*nOrb(iSym)+ip-1+nOrb(iSym)*NIsh(iSym)
      loc2 = ipMat(iSym,iSym)+(ip-1)*nOrb(iSym)+iq-1+NIsh(iSym)
      AXK(loc1) = AXK(loc1)+tempa
      AXK(loc2) = AXK(loc2)-tempa
    end do
    do ip=nIsh(iSym)+nAsh(iSym)+1,nOrb(iSym)  ! p is virtual
      p = ip+Off_Orb(ISym)
      tempa = Zero
      do it=1,nAsh(ISym)
        t = it+Off_Act(ISym)
        tempa = tempa+(D1((t-1)*nnA+q)+D1((q-1)*nnA+t))*Opu((p-1)*nnA+t)
      end do
      loc1 = ipMat(iSym,iSym)+(iq-1)*nOrb(iSym)+ip-1+nOrb(iSym)*NIsh(iSym)
      loc2 = ipMat(iSym,iSym)+(ip-1)*nOrb(iSym)+iq-1+NIsh(iSym)
      AXK(loc1) = AXK(loc1)+tempa
      AXK(loc2) = AXK(loc2)-tempa
    end do
  end do
end do
call mma_deallocate(Opu)

end subroutine CalcAXk2
