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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine prtopt2_cvb(iopt1,ioptim,italter,noptim,iorts,ifxorb,ifxstr,idelstr)

use Index_Functions, only: nTri_Elem
use casvb_global, only: icrit, ifinish, imethod, ipr, isaddle, kbasis, lfxvb, lzrvb, mxiter, nfxorb, nfxvb, norb, nort, nvb, &
                        nzrvb, projcas, projsym, spinb
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iopt1, ioptim, italter, noptim, iorts(2,nTri_Elem(norb-1)), ifxorb(norb), ifxstr(nvb), idelstr(nvb)
integer(kind=iwp) :: i, ifx, io
character(len=9) :: sbformat
character(len=3) :: ayn
integer(kind=iwp), allocatable :: tmp(:)
character(len=8), parameter :: methkw(12) = ['Fletcher','    TRIM','Trustopt','Davidson','   Steep','  Vb2cas',' AugHess', &
                                             'AugHess2','   Check',' dFletch','    None','Super-CI']

if (ifinish == 0) then
  if ((ipr(3) >= 1) .or. ((ipr(3) == 0) .and. ((iopt1 == 0) .or. (italter == 1)))) then
    if (noptim == 1) then
      write(u6,'(/,a)') ' -- Starting optimization ------------------'
    else
      write(u6,'(/,a,i3,a)') ' -- Starting optimization - step',ioptim,' --------'
    end if
  end if
  if ((ipr(3) >= 1) .and. ((iopt1 == 0) .or. (italter == 1))) then
    if (icrit == 1) then
      write(u6,'(/,a)') ' Overlap-based optimization (Svb).'
    else if (icrit == 2) then
      write(u6,'(/,a)') ' Energy-based optimization (Evb).'
    end if
    write(u6,'(/,a,11x,a)') ' Optimization algorithm:',methkw(imethod)
    write(u6,'(a,i13)') ' Maximum number of iterations:',mxiter
    ayn = ' No'
    if (projcas) ayn = 'Yes'
    if (projcas) write(u6,'(a,31x,a)') ' Casproj:',ayn
    ayn = ' No'
    if (projsym) ayn = 'Yes'
    if (projsym) write(u6,'(a,31x,a)') ' Symproj:',ayn
    ayn = ' No'
    sbformat = '(a,19x,a)'
    write(sbformat(4:5),'(i2)') 31-len_trim(spinb(kbasis))
    write(u6,sbformat) ' Spin basis:',trim(spinb(kbasis))
    if (isaddle > 0) write(u6,'(/,a,i9)') ' Saddle-point optimization, order:',isaddle
    if (nort > 0) then
      write(u6,'(/,i4,a)') nort,' orthogonalization pairs defined :'
      write(u6,6100) (io,iorts(1,io),iorts(2,io),io=1,nort)
    end if
    if (nfxorb == norb) then
      write(u6,'(/,a)') ' All orbitals will be frozen.'
    else if (nfxorb > 0) then
      write(u6,'(/,a)') ' Following orbitals will be frozen :'
      call mma_allocate(tmp,nfxorb,label='tmp')
      ifx = 0
      do i=1,norb
        if ((ifxorb(i) >= 0) .and. (ifxorb(i) <= norb)) then
          ifx = ifx+1
          tmp(ifx) = i
        end if
      end do
      nfxorb = ifx
      write(u6,'(14i3)') tmp(1:nfxorb)
      call mma_deallocate(tmp)
    end if
    if ((nfxvb > 0) .and. (lfxvb == 0)) then
      write(u6,'(/,a)') ' Following structures will be frozen :'
      write(u6,'(14i3)') ifxstr(1:nfxvb)
    else if ((nfxvb == 0) .and. (lfxvb == 1)) then
      write(u6,'(/,a)') ' All structures will be frozen.'
    else if ((nfxvb > 0) .and. (lfxvb == 1)) then
      write(u6,'(/,a)') ' Following structures coefficients will be optimized :'
      write(u6,'(14i3)') ifxstr(1:nfxvb)
    end if
    if ((nzrvb > 0) .and. (lzrvb == 0)) then
      write(u6,'(/,a)') ' Following structures will be deleted :'
      write(u6,'(14i3)') idelstr(1:nzrvb)
    else if ((nzrvb == 0) .and. (lzrvb == 1)) then
      write(u6,'(/,a)') ' All structures will be deleted.'
    else if ((nzrvb > 0) .and. (lzrvb == 1)) then
      write(u6,'(/,a)') ' Following structures will not be deleted :'
      write(u6,'(14i3)') idelstr(1:nzrvb)
    end if
    write(u6,'(/,a)') ' -------------------------------------------'
  end if
  if ((iopt1 == 0) .or. (italter == 1)) call tuneprint_cvb()
else if (ifinish < 3) then
  if ((ipr(3) >= 1) .or. ((ipr(3) == 0) .and. ((iopt1 == 0) .or. (italter == 1)))) then
    if (noptim == 1) then
      write(u6,'(/,a)') ' -- Wavefunction summary -------------------'
    else
      write(u6,'(/,a,i3,a)') ' -- Wavefunction summary - step',ioptim,' ---------'
    end if
  end if
end if

return
6100 format(3(i4,': ',i2,' -',i2))

end subroutine prtopt2_cvb
