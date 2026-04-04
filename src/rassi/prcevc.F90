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

subroutine PRCEVC(NSS,FRAC,SOENE,MAPST,MAPSP,MAPMS,UMATR,UMATI)

use rassi_aux, only: ipglob
use Cntrl, only: NSTATE
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: NSS, MAPST(NSS), MAPSP(NSS), MAPMS(NSS)
real(kind=wp) :: FRAC, SOENE(NSS), UMATR(NSS,NSS), UMATI(NSS,NSS)
integer(kind=iwp) :: I, IMAXSTATE, ISFS, ISS, ISTATE, JEND, JSS, JSTA, nmax(5), NW
real(kind=wp) :: CFFLIM, S, smax(5), SSMAX, SZ, TST, WGTMAX, wmax(5), XMAX
real(kind=wp), allocatable :: sstate(:), weight(:)

call mma_allocate(weight,nss)
call mma_allocate(sstate,nss)

if (IPGLOB >= 3) then
  ! Write out the complex eigenvectors of the spin-orbit states.
  ! Four states at a time are written out.
  ! They are written as complex numbers, four on each line.
  ! Coefficients are written out if any of the four complex coeffs
  ! on the same line have absolute value at least as large as
  ! a certain fraction (FRAC) of the largest such value for the
  ! foursome of states.
  do JSTA=1,NSS,4
    JEND = min(NSS,JSTA+3)
    write(u6,*)
    write(u6,'(1X,A16,F16.8,3(2X,F16.8))') '    Energy (au) ',(SOENE(JSS),JSS=JSTA,JEND)
    write(u6,'(1X,A16,6X,I4,3(14X,I4))') ' SFS  S     Ms  ',(JSS,JSS=JSTA,JEND)
    ! Scan coefficients to pick out the largest:
    WGTMAX = Zero
    do ISS=1,NSS
      do JSS=JSTA,JEND
        WGTMAX = max(WGTMAX,UMATR(ISS,JSS)**2+UMATI(ISS,JSS)**2)
      end do
    end do

    ! Scan coefficients, write if large enough:
    CFFLIM = FRAC*sqrt(WGTMAX)
    do ISS=1,NSS
      ISTATE = MAPST(ISS)
      S = Half*real(MAPSP(ISS)-1,kind=wp)
      SZ = Half*real(MAPMS(ISS),kind=wp)
      TST = Zero
      do JSS=JSTA,JEND
        TST = max(TST,UMATR(ISS,JSS)**2+UMATI(ISS,JSS)**2)
      end do
      if (TST >= CFFLIM**2) &
        write(u6,'(I4,1X,F4.1,1X,F5.1,3X,4(A1,F7.4,A1,F7.4,A1,1x))') ISTATE,S,SZ, &
                                                                     ('(',UMATR(ISS,JSS),',',UMATI(ISS,JSS),')',JSS=JSTA,JEND)
    end do
  end do
else if (IPGLOB >= 2) then
  ! Write out the weights of the five most important spin-free states
  ! For each spin-orbit state (BOR in Krapperup 070226)

  write(u6,*)
  write(u6,*)
  write(u6,*) 'Weights of the five most important spin-orbit-free states for each spin-orbit state.'
  write(u6,*)
  write(u6,*) 'SO State  Total energy (au)           Spin-free states, spin, and weights'
  write(u6,*) '-------------------------------------------------------------------------------------------------------'
  do iss=1,nss
    do isfs=1,nstate
      weight(isfs) = Zero
    end do
    do jss=1,nss
      istate = mapst(jss)
      sstate(istate) = Half*real(mapsp(jss)-1,kind=wp)
      weight(istate) = weight(istate)+umatr(jss,iss)**2+umati(jss,iss)**2
    end do
    ! Sort the weights
    imaxstate = 0
    do
      xmax = Zero
      ssmax = Zero
      nw = 0
      do istate=1,nstate
        if (weight(istate) >= xmax) then
          nw = istate
          xmax = weight(istate)
          ssmax = sstate(istate)
        end if
      end do
      weight(nw) = -One
      imaxstate = imaxstate+1
      nmax(imaxstate) = nw
      wmax(imaxstate) = xmax
      smax(imaxstate) = ssmax
      if (imaxstate >= min(nstate,5)) exit
    end do

    write(u6,'(i5,1x,f16.6,3x,5(i5,f4.1,f8.4))') iss,soene(iss),(nmax(i),smax(i),wmax(i),i=1,(min(nstate,5)))

  end do
  write(u6,*) '-------------------------------------------------------------------------------------------------------'
end if

! Added by Ungur Liviu on 04.11.2009
! Addition of UMATR and UMATI on RunFile

call Put_dArray('UMATR_SINGLE',UMATR,NSS**2)
call Put_dArray('UMATI_SINGLE',UMATI,NSS**2)

call mma_deallocate(weight)
call mma_deallocate(sstate)

end subroutine PRCEVC
