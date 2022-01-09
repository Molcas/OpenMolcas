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
!          one-electron hamiltonian.                                   *
!                                                                      *
!     Internal matrices fixed the maximum number of primitives         *
!     per symmetry to 'maxprim'                                        *
!                                                                      *
!***********************************************************************

use Basis_Info

implicit real*8(A-H,O-Z)
external Agin, Ovlmp, Vexch, Vqr
#include "relmp.fh"
#include "stdalloc.fh"
character*20 MPLbl
! working variables (change this)
parameter(maxprim=40)
parameter(mx100=100)
parameter(Mxlpq=(mx100*(mx100+1)/2))
parameter(Nrel=7*Mxlpq+5*mx100*mx100+5*mx100)
real*8 COREK(maxprim,maxprim,2), OVL(maxprim,maxprim), rel(Nrel), srel(maxprim*(maxprim+1)/2), trel(maxprim*(maxprim+1)/2), &
       urel(maxprim*(maxprim+1)/2), tnrel(maxprim*(maxprim+1)/2), unrel(maxprim*(maxprim+1)/2), hcorr(maxprim*(maxprim+1)/2)
real*8 VEXTT(maxprim*(maxprim+1)/2), PVPT(maxprim*(maxprim+1)/2), EVN1(maxprim,maxprim), EVN2(maxprim,maxprim), &
       RE1R(maxprim,maxprim), AUXI(maxprim,maxprim), W1W1(maxprim,maxprim), W1E0W1(maxprim,maxprim)
data iExch/1/,iMVPot/2/,iDWPot/4/,iNPPot/8/

PI = 2d0*acos(0d0)
iprint = 0
call Agin()
call dcopy_(2*MaxPrim**2,[0.0d0],0,Corek,1)

! calculate relativistic integrals if needed
lpq = 0
do i=1,lmax+1
  nnexp = Shells(isrshll+i-1)%nExp
  lpq = lpq+nnexp*(nnexp+1)/2
end do
if (4*lpq > Nrel) then
  write(6,*) ' problem. Nrel and lpq are',nrel,lpq
  write(6,*) ' The dimension of rel must somehow be increased.'
  call Abend()
end if
#ifdef _DEBUGPRINT_
write(6,*) ' basis:',(Shells(isrshll+i-1)%nExp,i=1,lmax+1)
do i=1,lmax+1
  nnExp = Shells(isrshll+i-1)%nExp
  write(6,*) ' number of exponents',nnExp
  write(6,*) ' exponents, symmetry',i
  do j=1,nnExp
    write(6,*) Shells(iSRShll+i-1)%exp(j)
  end do
end do
#endif
do lP1=1,lMax+1
  iSRSh = iSRShll+lP1-1
  nP = Shells(iSRSh)%nExp
  if (np > maxprim) then
    write(6,*) 'CalcAMt: np > maxprim',np,maxprim
    write(6,*) 'Abend: Increase MaxPrim !'
    call Abend()
  end if
  N = lP1
  LAM = lP1
  if (nP <= 0) Go To 1000

  call dcopy_(nRel,[0.0d0],0,Rel,1)
  if ((iand(iOpt,iMVPot) /= 0) .and. (iand(iOpt,iDWPot) /= 0)) then
    ! Mass-velocity and/or Darwin potentials
    call Vqr(LUQRP,MPLbl,lP1,Shells(iSRSh)%Exp,nP,rel)
  else if ((iand(iOpt,iMVPot) /= 0) .or. (iand(iOpt,iDWPot) /= 0)) then
    write(6,*) 'Mass-Velocity and Darwin potentials must be'
    write(6,*) 'active simultaneosly.'
    call Abend()
  end if

  if (iand(iOpt,iNPPot) /= 0) then   ! Zero
    call oeisg(rel,srel,trel,urel,Shells(iSRSh)%Exp,rCharge,mx100,lp1,Shells(iSRSh)%nExp,unrel,tnrel,hcorr,iprint,VEXTT,PVPT,EVN1, &
               EVN2,RE1R,AUXI,W1W1,W1E0W1)
    nmat = (Shells(iSRSh)%nExp*(Shells(iSRSh)%nExp+1))/2
    if (iprint >= 10) then
      write(6,*) ' relativistic integrals'
      write(6,12) (hcorr(i),i=1,nmat)
12    format(4d19.12)
    end if
    call dcopy_(nRel,[0.0d0],0,Rel,1)
  end if

  ! Overlap and, if neccesary, exchange.
  IJ = 0
  do I=1,NP
    ZI = Shells(iSRSh)%exp(I)
    do J=1,I
      IJ = IJ+1
      ZJ = Shells(iSRSh)%exp(J)
      COREK(I,J,1) = rel(ij)
      if (iand(iOpt,iNPPot) /= 0) then
        corek(i,j,2) = hcorr(ij)
      end if
      if (iand(iOpt,iExch) /= 0) then
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

  PreFac = sqrt((2d0/PI)**3)*4d0**(lP1-1)

  call mma_Allocate(Shells(ISRSh)%Akl,NP,NP,2,Label='Akl')
  Shells(ISRSh)%nAkl = NP
  do iq=1,2
    do I=1,NP
      ZI = Shells(iSRSh)%exp(I)
      do L=1,NP
        rel(L) = 0.d0
        do K=1,NP
          rel(L) = rel(L)+OVL(I,K)*COREK(K,L,iq)
        end do
      end do
      do J=1,I
        ZJ = Shells(iSRSh)%exp(J)
        ADUM = 0.d0
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

1000 continue
end do

return

end subroutine CalcAMt
