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

subroutine CalcAMt(iOpt,LUQRP,MPLbl,lMax,iSRShll,nProj,iCoShll,rcharge)
!***********************************************************************
!                                                                      *
!...       calculates the non-diagonal spectral representation         *
!          A matrix for an atom. Note that its signs is such           *
!          that the spectral representation must be ADDED to the       *
!          one-electron Hamiltonian.                                   *
!                                                                      *
!***********************************************************************

use Basis_Info, only: Shells
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Four, Pi
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iOpt, LUQRP, lMax, iSRShll, nProj, iCoShll
character(len=20), intent(in) :: MPLbl
real(kind=wp), intent(in) :: rcharge
integer(kind=iwp) :: i, IJ, IJAM0, iprint, iq, iSRSh, J, K, L, LAM, lP1, lpq, maxprim, maxprimt, N, nmat, nnexp, nP, nP1, nP2, &
                     Nrel, nscr
real(kind=wp) :: ADUM, AuxLs, PreFac, ZI, ZJ
integer(kind=iwp), allocatable :: iScratch(:)
real(kind=wp), allocatable :: AUXI(:), COREK(:,:,:), EVN1(:), hcorr(:), OVL(:,:), PVPT(:), RE1R(:), rel(:), Scratch(:), srel(:), &
                              tnrel(:), trel(:), unrel(:), urel(:), VEXTT(:), W1W1(:)
!integer(kind=iwp), parameter :: iExch = 1, iMVPot = 2, iDWPot = 4, iNPPot = 8
integer(kind=iwp), parameter :: iExch = 0, iMVPot = 1, iDWPot = 2, iNPPot = 3
real(kind=wp), external :: OVLMP, Vexch

iprint = 0
call Agin()

! calculate relativistic integrals if needed
lpq = 0
maxprim = 0
do i=1,lmax+1
  nnExp = Shells(iSRShll+i-1)%nExp
  maxprim = max(maxprim,nnExp)
  lpq = lpq+nnExp*(nnExp+1)/2
end do
maxprimt = maxprim*(maxprim+1)/2
Nrel = max(4*lpq,maxprimt)
nscr = 3*maxprimt+5*maxprim*maxprim+5*maxprim
#ifdef _DEBUGPRINT_
write(u6,*) ' basis:',(Shells(iSRShll+i-1)%nExp,i=1,lmax+1)
do i=1,lmax+1
  nnExp = Shells(iSRShll+i-1)%nExp
  write(u6,*) ' number of exponents',nnExp
  write(u6,*) ' exponents, symmetry',i
  do j=1,nnExp
    write(u6,*) Shells(iSRShll+i-1)%Exp(j)
  end do
end do
#endif
call mma_allocate(rel,Nrel,label='rel')
call mma_allocate(srel,maxprimt,label='srel')
call mma_allocate(trel,maxprimt,label='trel')
call mma_allocate(urel,maxprimt,label='urel')
call mma_allocate(unrel,maxprimt,label='unrel')
call mma_allocate(tnrel,maxprimt,label='tnrel')
call mma_allocate(hcorr,maxprimt,label='hcorr')
call mma_allocate(VEXTT,maxprimt,label='VEXTT')
call mma_allocate(PVPT,maxprimt,label='PVPT')
call mma_allocate(EVN1,maxprim**2,label='EVN1')
call mma_allocate(AUXI,maxprim**2,label='AUXI')
call mma_allocate(RE1R,maxprim**2,label='RE1R')
call mma_allocate(W1W1,maxprim**2,label='W1W1')
call mma_allocate(OVL,maxprim,maxprim,label='OVL')
call mma_allocate(COREK,maxprim,maxprim,2,label='COREK')
call mma_allocate(iScratch,maxprimt,label='iScratch')
call mma_allocate(Scratch,nscr,label='Scratch')
COREK(:,:,:) = Zero
do lP1=1,lMax+1
  iSRSh = iSRShll+lP1-1
  nP = Shells(iSRSh)%nExp
  N = lP1
  LAM = lP1
  if (nP <= 0) cycle

  Rel(:) = Zero
  if (btest(iOpt,iMVPot) .and. btest(iOpt,iDWPot)) then
    ! Mass-velocity and/or Darwin potentials
    call Vqr(LUQRP,MPLbl,lP1,Shells(iSRSh)%Exp,nP,rel)
  else if (btest(iOpt,iMVPot) .or. btest(iOpt,iDWPot)) then
    write(u6,*) 'Mass-Velocity and Darwin potentials must be'
    write(u6,*) 'active simultaneosly.'
    call Abend()
  end if

  if (btest(iOpt,iNPPot)) then   ! Zero
    nP1 = 3*nP*(nP+1)/2
    nP2 = nP1+5*nP*nP
    call oeisg(rel,srel,trel,urel,Shells(iSRSh)%Exp,rCharge,lp1,nP,unrel,tnrel,hcorr,iprint,VEXTT,PVPT,EVN1,RE1R,AUXI,W1W1, &
               iScratch,Scratch(1:nP1),Scratch(nP1+1:nP2),Scratch(nP2+1:))
    nmat = nP*(nP+1)/2
    if (iprint >= 10) then
      write(u6,*) ' relativistic integrals'
      write(u6,'(4d19.12)') (hcorr(i),i=1,nmat)
    end if
    Rel(:) = Zero
  end if

  ! Overlap and, if neccesary, exchange.
  IJ = 0
  do I=1,NP
    ZI = Shells(iSRSh)%Exp(I)
    do J=1,I
      IJ = IJ+1
      ZJ = Shells(iSRSh)%Exp(J)
      COREK(I,J,1) = rel(ij)
      if (btest(iOpt,iNPPot)) COREK(I,J,2) = hcorr(ij)
      if (btest(iOpt,iExch)) then
        ! minus exchange potential
        AuxLs = VExch(ZI,N,ZJ,N,LAM,nProj,iCoShll)
        COREK(I,J,1) = COREK(I,J,1)-AuxLs
      end if
      OVL(I,J) = OVLMP(N,ZI,N,ZJ)
      OVL(J,I) = OVL(I,J)
      COREK(J,I,1) = COREK(I,J,1)
      COREK(J,I,2) = COREK(I,J,2)
    end do
  end do

  call MATINV(OVL,rel,NP,0,maxprim)

  PreFac = sqrt((Two/PI)**3)*Four**(lP1-1)

  call mma_Allocate(Shells(ISRSh)%Akl,NP,NP,2,Label='Akl')
  Shells(ISRSh)%nAkl = NP
  do iq=1,2
    do I=1,NP
      ZI = Shells(iSRSh)%Exp(I)
      do L=1,NP
        rel(L) = Zero
        do K=1,NP
          rel(L) = rel(L)+OVL(I,K)*COREK(K,L,iq)
        end do
      end do
      do J=1,I
        ZJ = Shells(iSRSh)%Exp(J)
        ADUM = Zero
        do L=1,NP
          ADUM = ADUM+rel(L)*OVL(L,J)
        end do
        ! in MOLCAS3:X1
        ! multiply by the (radial) normalization constants
        ! of the primitives i and j, so that the spectral
        ! representation coeffients correspond to the
        ! (radially) unnormalized primitives.
        ADUM = ADUM*PreFac*sqrt(sqrt((ZI*ZJ)**(3+2*(lP1-1))))
        Shells(iSRSh)%Akl(I,J,iq) = ADUM
        Shells(iSRSh)%Akl(J,I,iq) = ADUM
      end do
    end do
    IJAM0 = IJAM0+NP**2
  end do

end do

call mma_deallocate(rel)
call mma_deallocate(srel)
call mma_deallocate(trel)
call mma_deallocate(urel)
call mma_deallocate(unrel)
call mma_deallocate(tnrel)
call mma_deallocate(hcorr)
call mma_deallocate(VEXTT)
call mma_deallocate(PVPT)
call mma_deallocate(EVN1)
call mma_deallocate(AUXI)
call mma_deallocate(RE1R)
call mma_deallocate(W1W1)
call mma_deallocate(OVL)
call mma_deallocate(COREK)
call mma_deallocate(iScratch)
call mma_deallocate(Scratch)

return

end subroutine CalcAMt
