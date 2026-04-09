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
! Copyright (C) 1989, Per Ake Malmqvist                                *
!***********************************************************************
!****************************************************************
!  PROGRAM RASSI        PER-AAKE MALMQVIST
!  SUBROUTINE GETH1     IBM-3090 RELEASE 89 01 30
!  READ THE ONE-ELECTRON HAMILTONIAN MATRIX ELEMENTS AND RETURN
!  IT AS HONEAO IN SYMMETRY-BLOCKED SQUARED FORMAT.
!  Also reads and adds reaction field contribution.
!  Also ERFNUC, reaction field contribution to nuclear repulsion.
!****************************************************************

subroutine GETH1_RASSI(HONEAO)

use OneDat, only: sNoNuc, sNoOri
use Cntrl, only: ERFNuc, RFPert
use Symmetry_Info, only: nIrrep
use rassi_data, only: NBASF, NBSQ, NBTRI
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: HONEAO(NBSQ)
integer(kind=iwp) :: iBuf, ICMP, IOPT, IP, IPQ, IQ, IQP, IRC, ISTQ, iSyLab, ISYM, NB
logical(kind=iwp) :: Found
real(kind=wp), allocatable :: H1(:), Tmp(:)
character(len=*), parameter :: OneLbl = 'OneHam'

call mma_allocate(H1,NBTRI,Label='H1')
iRc = -1
iOpt = ibset(ibset(0,sNoOri),sNoNuc)
iCmp = 1
iSyLab = 1
call RdOne(iRc,iOpt,OneLbl,iCmp,H1,iSyLab)
if (IRC /= 0) then
  write(u6,*)
  write(u6,*) '      *** ERROR IN SUBROUTINE  GETH1 ***'
  write(u6,*) '   BARE NUCLEI HAMILTONIAN IS NOT AVAILABLE'
  write(u6,*)
  call ABEND()
end if
ERFNuc = Zero
if (RFpert) then
  call f_Inquire('RUNOLD',Found)
  if (Found) call NameRun('RUNOLD')
  call mma_allocate(Tmp,nBtri,Label='Tmp')
  call Get_dScalar('RF Self Energy',ERFNuc)
  call Get_dArray('Reaction field',Tmp,nBtri)
  if (Found) call NameRun('#Pop')
  H1(:) = H1(:)+Tmp(:)
  call mma_deallocate(Tmp)
end if
IBUF = 1
ISTQ = 0
do ISYM=1,nIrrep
  NB = NBASF(ISYM)
  if (NB == 0) cycle
  do IP=1,NB
    do IQ=1,IP
      IPQ = NB*(IP-1)+IQ+ISTQ
      IQP = NB*(IQ-1)+IP+ISTQ
      HONEAO(IPQ) = H1(IBUF)
      HONEAO(IQP) = H1(IBUF)
      IBUF = IBUF+1
    end do
  end do
  ISTQ = ISTQ+NB**2
end do
call mma_deallocate(H1)

end subroutine GETH1_RASSI
