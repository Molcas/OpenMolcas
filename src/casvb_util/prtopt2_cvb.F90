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

implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
#include "spinb_cvb.fh"
#include "WrkSpc.fh"
parameter(nmeth=12)
character*3 ayn
character*8 methkw(nmeth)
character*9 sbformat
dimension iorts(2*norb*(norb-1)/2)
dimension ifxorb(norb), ifxstr(nvb), idelstr(nvb)
save methkw
data methkw/'Fletcher','    TRIM','Trustopt','Davidson','   Steep','  Vb2cas',' AugHess','AugHess2','   Check',' dFletch', &
            '    None','Super-CI'/

if (ifinish == 0) then
  if ((ip(3) >= 1) .or. ((ip(3) == 0) .and. ((iopt1 == 0) .or. (italter == 1)))) then
    if (noptim == 1) then
      write(6,'(/,a)') ' -- Starting optimization ------------------'
    else
      write(6,'(/,a,i3,a)') ' -- Starting optimization - step',ioptim,' --------'
    end if
  end if
  if ((ip(3) >= 1) .and. ((iopt1 == 0) .or. (italter == 1))) then
    if (icrit == 1) then
      write(6,'(/,a)') ' Overlap-based optimization (Svb).'
    else if (icrit == 2) then
      write(6,'(/,a)') ' Energy-based optimization (Evb).'
    end if
    write(6,'(/,a,11x,a)') ' Optimization algorithm:',methkw(imethod)
    write(6,'(a,i13)') ' Maximum number of iterations:',mxiter
    ayn = ' No'
    if (projcas) ayn = 'Yes'
    if (projcas) write(6,'(a,31x,a)') ' Casproj:',ayn
    ayn = ' No'
    if (projsym) ayn = 'Yes'
    if (projsym) write(6,'(a,31x,a)') ' Symproj:',ayn
    ayn = ' No'
    sbformat = '(a,19x,a)'
    write(sbformat(4:5),'(i2)') 31-len_trim_cvb(spinb(kbasis))
    write(6,sbformat) ' Spin basis:',spinb(kbasis)(1:len_trim_cvb(spinb(kbasis)))
    if (isaddle > 0) write(6,'(/,a,i9)') ' Saddle-point optimization, order:',isaddle
    if (nort > 0) then
      write(6,'(/,i4,a)') nort,' orthogonalization pairs defined :'
      write(6,6100) (ior,iorts(1+(ior-1)*2),iorts(2+(ior-1)*2),ior=1,nort)
    end if
    if (nfxorb == norb) then
      write(6,'(/,a)') ' All orbitals will be frozen.'
    else if (nfxorb > 0) then
      write(6,'(/,a)') ' Following orbitals will be frozen :'
      itmp = mstacki_cvb(nfxorb)
      ifx = 0
      do i=1,norb
        if ((ifxorb(i) >= 0) .and. (ifxorb(i) <= norb)) then
          ifx = ifx+1
          iwork(ifx+itmp-1) = i
        end if
      end do
      nfxorb = ifx
      write(6,'(14i3)') (iwork(ii+itmp-1),ii=1,nfxorb)
      call mfreei_cvb(itmp)
    end if
    if ((nfxvb > 0) .and. (lfxvb == 0)) then
      write(6,'(/,a)') ' Following structures will be frozen :'
      write(6,'(14i3)') (ifxstr(ii),ii=1,nfxvb)
    else if ((nfxvb == 0) .and. (lfxvb == 1)) then
      write(6,'(/,a)') ' All structures will be frozen.'
    else if ((nfxvb > 0) .and. (lfxvb == 1)) then
      write(6,'(/,a)') ' Following structures coefficients will be optimized :'
      write(6,'(14i3)') (ifxstr(ii),ii=1,nfxvb)
    end if
    if ((nzrvb > 0) .and. (lzrvb == 0)) then
      write(6,'(/,a)') ' Following structures will be deleted :'
      write(6,'(14i3)') (idelstr(ii),ii=1,nzrvb)
    else if ((nzrvb == 0) .and. (lzrvb == 1)) then
      write(6,'(/,a)') ' All structures will be deleted.'
    else if ((nzrvb > 0) .and. (lzrvb == 1)) then
      write(6,'(/,a)') ' Following structures will not be deleted :'
      write(6,'(14i3)') (idelstr(ii),ii=1,nzrvb)
    end if
    write(6,'(/,a)') ' -------------------------------------------'
  end if
  if ((iopt1 == 0) .or. (italter == 1)) call tuneprint_cvb()
else if (ifinish < 3) then
  if ((ip(3) >= 1) .or. ((ip(3) == 0) .and. ((iopt1 == 0) .or. (italter == 1)))) then
    if (noptim == 1) then
      write(6,'(/,a)') ' -- Wavefunction summary -------------------'
    else
      write(6,'(/,a,i3,a)') ' -- Wavefunction summary - step',ioptim,' ---------'
    end if
  end if
end if

return
6100 format(3(i4,': ',i2,' -',i2))

end subroutine prtopt2_cvb
