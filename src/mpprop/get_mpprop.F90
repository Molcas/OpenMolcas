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

subroutine Get_MpProp(nPrim,nBas,nAtoms,nCenters,nMltPl,ip_D_p,ECENTX,ECENTY,ECENTZ,LNearestAtom,LFirstRun,LLumOrb)
! nOcOb,oNum,nOrb,oCof

use MPProp_globals, only: BondMat, Cor, CordMltPl, Frac, iAtBoMltPlAd, iAtBoMltPlAdCopy, iAtMltPlAd, iMltPlAd, iAtPrTab, iQnuc, &
                          Labe, Method, nAtomPBas
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nPrim, nBas, nAtoms, nCenters, nMltPl, ip_D_p
real(kind=wp), intent(in) :: ECENTX(nPrim*(nPrim+1)/2), ECENTY(nPrim*(nPrim+1)/2), ECENTZ(nPrim*(nPrim+1)/2) !, oNum(nOrb), oCof(nBas,nPrim)
logical(kind=iwp), intent(in) :: LNearestAtom, LFirstRun, LLumOrb
#include "WrkSpc.fh"
integer(kind=iwp) :: i, iA, iComp, ii, il, iMltpl, ip, iPBas, iq, iStdout, j, jj, jPBas, k, nA, nB, nl, np, nq
real(kind=wp) :: CorP(3), CorN(3), FracA, FracB, FracN, FracP, Qn, Qp, Qs, R, RA, RB, rnloveril, rnPoveriP, rnqoveriq, Rtot, Rwei, &
                 Smallest, rsum, sum_a, sum_b, xfac, xfac_a, xfac_b, yfac, yfac_a, yfac_b, zfac, zfac_a, zfac_b
integer(kind=iwp), allocatable :: iCompMat(:,:,:)
real(kind=wp), allocatable :: Qexp(:,:) !, DenMat(:,:)

call mma_allocate(iCompMat,[0,nMltPl],[0,nMltPl],[0,nMltPl],label='iCompMat')
call mma_allocate(Qexp,nPrim,nPrim,label='Qexp')
!call mma_allocate(DenMat,nPrim,nPrim,label='DenMat')

!                                                                      *
!***********************************************************************
!                                                                      *
! SOME KIND OF PROLOG

iStdout = u6
write(iStdOut,*)
if (LLumOrb) then
  write(iStdOut,*) ' Using the INPORB file'
  if (LFirstRun) then
    write(iStdOut,*) ' This is first run of Get_MpProp'
    !write(iStdOut,*)
    !if (Method == 'UHF-SCF') then
    !  write(iStdOut,*) ' Number of alpha electrons',nOcOB
    !else
    !  write(iStdOut,*) ' Occupation Number ',nOcOB
    !end if
    !write(iStdOut,*) ' Number of Orbitals' ,nOrb
  else
    write(iStdOut,*) ' Running Get_MpProp for the second time'
    !write(iStdOut,*)
    !write(iStdOut,*) ' Number beta electrons ',nOcOB
    !write(iStdOut,*) ' Number of Orbitals' ,nOrb
  end if
else
  write(iStdOut,*) ' Using the densities from a Molcas calculation'
end if
write(iStdOut,*)
! Set the variable that knows the component
do iMltpl=0,nMltPl
  iComp = 0
  do np=iMltpl,0,-1
    do nq=iMltpl-np,0,-1
      nl = iMltpl-np-nq
      iComp = iComp+1
      iCompMat(np,nq,nl) = iComp
    end do
  end do
end do
! Translate dipoles to the center for the quadrupoles
if (LFirstRun) then
  !do iMltpl=0,1
  ! An error written by me DH
  do i=1,nPrim
    do j=1,i
      do k=1,3
        Work(iWork(iMltPlAd(1)+k-1)+i*(i-1)/2+j-1) = Work(iWork(iMltPlAd(1)+k-1)+i*(i-1)/2+j-1)+ &
                                                     (CordMltPl(k,0)-CordMltPl(k,2))*Work(iWork(iMltPlAd(0))+i*(i-1)/2+j-1)
      end do
    end do
  end do
  do iMltpl=0,1
    CordMltPl(1,iMltPl) = CordMltPl(1,2)
    CordMltPl(2,iMltPl) = CordMltPl(2,2)
    CordMltPl(3,iMltPl) = CordMltPl(3,2)
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! GET THE DENSITY MATRIX
!
!if (LLumOrb) then
!  do K=1,nPrim
!    do L=1,K
!      DenMat(K,L) = Zero
!    end do
!  end do
!  do K=1,nPrim
!    do L=1,K
!      do i=1,nOcOb
!        DenMat(K,L) = DenMat(K,L)+oNum(I)*oCof(I,K)*oCof(I,L)
!      end do
!      DenMat(L,K) = DenMat(K,L)
!    end do
!  end do
!else
!  do K=1,nPrim
!    do L=1,K
!      DenMat(L,K) = Work(ip_D_p+K*(K-1)/2+L-1)
!      if (K /= L) then
!        DenMat(K,L) = DenMat(L,K)
!      end if
!    end do
!  end do
!end if
!                                                                      *
!***********************************************************************
!                                                                      *
! GET THE INTERACTION SITES
do i=1,nPrim
  do j=1,i
    Qexp(i,j) = work(iwork(iMltPlAd(0))+i*(i-1)/2+j-1)*Work(ip_D_p+i*(i-1)/2+j-1)
    Qexp(j,i) = Qexp(i,j)
  end do
end do

iA = -99
do nA=1,nAtoms

  ! NOW SUM UP THE OTHER MULTIPOLE MOMENTS

  do iMltpl=0,nMltPl
    iComp = 0
    do np=iMltpl,0,-1
      do nq=iMltpl-np,0,-1
        nl = iMltpl-np-nq
        iComp = iComp+1
        rsum = Zero
        do ip=0,np
          call NoverP(np,ip,rnPoveriP)
          if (np == ip) then
            xfac = rnpoverip
          else
            xfac = rnPoveriP*(CordMltPl(1,ip)-Cor(1,nA,nA))**(np-ip)
          end if
          do iq=0,nq
            call NoverP(nq,iq,rnqoveriq)
            if (nq == iq) then
              yfac = rnqoveriq
            else
              yfac = rnqoveriq*(CordMltPl(2,iq)-Cor(2,nA,nA))**(nq-iq)
            end if
            do il=0,nl
              call NoverP(nl,il,rnloveril)
              if (nl == il) then
                zfac = rnloveril
              else
                zfac = rnloveril*(CordMltPl(3,il)-Cor(3,nA,nA))**(nl-il)
              end if
              if ((xfac == Zero) .or. (yfac == Zero) .or. (zfac == Zero)) cycle
              do iPBas=1,nAtomPBas(nA)
                i = iAtPrTab(iPBas,nA)
                do jPBas=1,nAtomPBas(nA)
                  !!!! Check i>j
                  if (iAtPrTab(jPBas,nA) > i) then
                    jj = i
                    ii = iAtPrTab(jPBas,nA)
                  else
                    jj = iAtPrTab(jPBas,nA)
                    ii = i
                  end if
                  rsum = rsum+xfac*yfac*zfac*Work(ip_D_p+ii*(ii-1)/2+jj-1)*Work(iWork(iMltPlAd(ip+iq+il)+ &
                        iCompMat(ip,iq,il)-1)+ii*(ii-1)/2+jj-1)
                end do
              end do
            end do
          end do
        end do
        Work(iAtBoMltPlAd(iMltpl)+nCenters*(iComp-1)+nA*(nA+1)/2-1) = Work(iAtBoMltPlAd(iMltpl)+nCenters*(iComp-1)+nA*(nA+1)/2-1)- &
                                                                      rsum ! minus from the negative sign of the electron
        Work(iAtBoMltPlAdCopy(iMltpl)+nCenters*(iComp-1)+nA*(nA+1)/2-1) = Work(iAtBoMltPlAd(iMltpl)+ &
                                                                          nCenters*(iComp-1)+nA*(nA+1)/2-1)
        Work(iAtMltPlAd(iMltPl)+nAtoms*(iComp-1)+nA-1) = Work(iAtMltPlAd(iMltPl)+nAtoms*(iComp-1)+nA-1)-rsum
      end do
    end do
  end do

  ! THE NUCLEAR CHARGE IS ADDED ON

  if (LFirstRun) then
    Work(iAtMltPlAd(0)+nA-1) = Work(iAtMltPlAd(0)+nA-1)+Work(iQnuc+nA-1)
    Work(iAtBoMltPlAd(0)+nA*(nA+1)/2-1) = Work(iAtBoMltPlAd(0)+nA*(nA+1)/2-1)+Work(iQnuc+nA-1)
    Work(iAtBoMltPlAdCopy(0)+nA*(nA+1)/2-1) = Work(iAtBoMltPlAdCopy(0)+nA*(nA+1)/2-1)+Work(iQnuc+nA-1)
  end if

  do nB=1,nA-1

    ! MULTIPOLE BETWEEN NA,NB

    Qp = Zero
    Qn = Zero
    do i=1,3
      CorP(i) = Zero
      CorN(i) = Zero
    end do
    do iPBas=1,nAtomPBas(NA)
      i = iAtPrTab(iPBas,NA)
      do jPBas=1,nAtomPBas(NB)
        j = iAtPrTab(jPBas,NB)
        !!!! Check i>j
        if (iAtPrTab(jPBas,nB) > i) then
          jj = i
          ii = iAtPrTab(jPBas,nB)
        else
          jj = iAtPrTab(jPBas,nB)
          ii = i
        end if
        if (Qexp(i,j) > Zero) then
          Qp = Qp+Qexp(i,j)
          CorP(1) = CorP(1)+QEXP(I,J)*ECENTX(ii*(ii-1)/2+jj)
          CorP(2) = CorP(2)+QEXP(I,J)*ECENTY(ii*(ii-1)/2+jj)
          CorP(3) = CorP(3)+QEXP(I,J)*ECENTZ(ii*(ii-1)/2+jj)
        else
          Qn = Qn+Qexp(i,j)
          CorN(1) = CorN(1)+QEXP(I,J)*ECENTX(ii*(ii-1)/2+J)
          CorN(2) = CorN(2)+QEXP(I,J)*ECENTY(ii*(ii-1)/2+J)
          CorN(3) = CorN(3)+QEXP(I,J)*ECENTZ(ii*(ii-1)/2+J)
        end if
      end do
    end do

    if (Qp /= Zero) then
      CorP(1) = CorP(1)/Qp
      CorP(2) = CorP(2)/Qp
      CorP(3) = CorP(3)/Qp
    end if
    if (Qn /= Zero) then
      CorN(1) = CorN(1)/Qn
      CorN(2) = CorN(2)/Qn
      CorN(3) = CorN(3)/Qn
    end if
    Qs = Qp-Qn
    if (Qs /= Zero) then
      FracP = Qp/Qs
      FracN = One-FracP
    else
      FracP = Zero
      FracN = Zero
    end if
    Cor(1,nA,nB) = CorP(1)*FracP+CorN(1)*FracN
    Cor(2,nA,nB) = CorP(2)*FracP+CorN(2)*FracN
    Cor(3,nA,nB) = CorP(3)*FracP+CorN(3)*FracN

    ! CALCULATE THE WEIGTH OF EACH ATOMIC CENTER BY A SIMPLE GEOMETRIC RATIO

    Rwei = sqrt((Cor(1,nA,nA)-Cor(1,nA,nB))**2+(Cor(2,nA,nA)-Cor(2,nA,nB))**2+(Cor(3,nA,nA)-Cor(3,nA,nB))**2)
    Rtot = sqrt((Cor(1,nA,nA)-Cor(1,nB,nB))**2+(Cor(2,nA,nA)-Cor(2,nB,nB))**2+(Cor(3,nA,nA)-Cor(3,nB,nB))**2)

    ! FRACTION OF MULTIPOLE EXPANSION TO BE RELATED TO THE PAIR ATOMS

    FracB = Rwei/Rtot
    FracA = One-FracB
    Frac(nA,nB) = FracA

    ! Find the closest atom
    if (LNearestAtom .and. (.not. BondMat(nA,nB))) then
      Smallest = huge(Smallest)
      do i=1,nAtoms
        R = sqrt((Cor(1,nA,nB)-Cor(1,i,i))**2+(Cor(2,nA,nB)-Cor(2,i,i))**2+(Cor(3,nA,nB)-Cor(3,i,i))**2)
        if (R < Smallest) then
          Smallest = R
          iA = i
        end if
      end do
      RA = sqrt((Cor(1,nA,nB)-Cor(1,nA,nA))**2+(Cor(2,nA,nB)-Cor(2,nA,nA))**2+(Cor(3,nA,nB)-Cor(3,nA,nA))**2)
      RB = sqrt((Cor(1,nA,nB)-Cor(1,nB,nB))**2+(Cor(2,nA,nB)-Cor(2,nB,nB))**2+(Cor(3,nA,nB)-Cor(3,nB,nB))**2)
      if (((iA /= nA) .and. (iA /= nB)) .and. ((Smallest < RA) .and. (Smallest < RB))) then
        iA = iA
        FracA = One
        FracB = Zero
        write(iStdOut,*)
        write(iStdOut,*) ' Moving bond between the atoms  ',Labe(nA),Labe(nB)
        write(iStdOut,*) ' to the atom                    ',Labe(iA)
        write(iStdOut,*)
      else
        iA = nA
      end if
    else
      iA = nA
    end if

    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! TRANSFORM THE MULTIPOLES TO BONDS AND ATOMS
    !
    ! Do multipoles

    do iMltpl=0,nMltPl
      iComp = 0
      do np=iMltpl,0,-1
        do nq=iMltpl-np,0,-1
          nl = iMltpl-np-nq
          iComp = iComp+1
          rsum = Zero
          do ip=0,np
            call NoverP(np,ip,rnPoveriP)
            if (np == ip) then
              xfac = rnpoverip
            else
              xfac = rnPoveriP*(CordMltPl(1,ip)-Cor(1,nA,nB))**(np-ip)
            end if
            do iq=0,nq
              call NoverP(nq,iq,rnqoveriq)
              if (nq == iq) then
                yfac = rnqoveriq
              else
                yfac = rnqoveriq*(CordMltPl(2,iq)-Cor(2,nA,nB))**(nq-iq)
              end if
              do il=0,nl
                call NoverP(nl,il,rnloveril)
                if (nl == il) then
                  zfac = rnloveril
                else
                  zfac = rnloveril*(CordMltPl(3,il)-Cor(3,nA,nB))**(nl-il)
                end if
                if ((xfac == Zero) .or. (yfac == Zero) .or. (zfac == Zero)) cycle
                do iPBas=1,nAtomPBas(nA)
                  i = iAtPrTab(iPBas,nA)
                  do jPBas=1,nAtomPBas(nB)
                    !!!! Check i>j
                    if (iAtPrTab(jPBas,nB) > i) then
                      jj = i
                      ii = iAtPrTab(jPBas,nB)
                    else
                      jj = iAtPrTab(jPBas,nB)
                      ii = i
                    end if
                    rsum = rsum+xfac*yfac*zfac*Two*Work(ip_D_p+ii*(ii-1)/2+jj-1)*Work(iWork(iMltPlAd(ip+iq+il)+ &
                           iCompMat(ip,iq,il)-1)+ii*(ii-1)/2+jj-1)
                  end do
                end do
              end do
            end do
          end do
          Work(iAtBoMltPlAd(iMltpl)+nCenters*(iComp-1)+nA*(nA-1)/2+nB-1) = Work(iAtBoMltPlAd(iMltpl)+ &
                                                                           nCenters*(iComp-1)+nA*(nA-1)/2+nB-1)-rsum
                                                                           ! minus from the negative sign of the electron
          ! Copy the multipole arrays
          Work(iAtBoMltPlAdCopy(iMltpl)+nCenters*(iComp-1)+nA*(nA-1)/2+nB-1) = Work(iAtBoMltPlAd(iMltpl)+nCenters*(iComp-1)+ &
                                                                               nA*(nA-1)/2+nB-1)
        end do
      end do
    end do

    ! If one should move the bond to the closest atom iA would be equal to
    ! that otherwise iA should be equal to nA

    if (((Method == 'UHF-SCF') .and. (LFirstRun .eqv. .false.)) .or. ((Method /= 'UHF-SCF') .and. (LFirstRun .eqv. .true.))) then
      do iMltpl=0,nMltPl
        iComp = 0
        do np=iMltpl,0,-1
          do nq=iMltpl-np,0,-1
            nl = iMltpl-np-nq
            iComp = iComp+1
            sum_a = Zero
            sum_b = Zero
            do ip=0,np
              call NoverP(np,ip,rnPoveriP)
              if (np == ip) then
                xfac_a = rnPoveriP
                xfac_b = rnPoveriP
              else
                xfac_a = (Cor(1,nA,nB)-Cor(1,iA,iA))**(np-ip)*rnPoveriP
                xfac_b = (Cor(1,nA,nB)-Cor(1,nB,nB))**(np-ip)*rnPoveriP
              end if
              do iq=0,nq
                call NoverP(nq,iq,rnqoveriq)
                if (nq == iq) then
                  yfac_a = rnqoveriq
                  yfac_b = rnqoveriq
                else
                  yfac_a = (Cor(2,nA,nB)-Cor(2,iA,iA))**(nq-iq)*rnqoveriq
                  yfac_b = (Cor(2,nA,nB)-Cor(2,nB,nB))**(nq-iq)*rnqoveriq
                end if
                do il=0,nl
                  call NoverP(nl,il,rnloveril)
                  if (nl == il) then
                    zfac_a = rnloveril
                    zfac_b = rnloveril
                  else
                    zfac_a = (Cor(3,nA,nB)-Cor(3,iA,iA))**(nl-il)*rnloveril
                    zfac_b = (Cor(3,nA,nB)-Cor(3,nB,nB))**(nl-il)*rnloveril
                  end if
                  sum_a = sum_a+xfac_a*yfac_a*zfac_a*FracA*Work(iAtBoMltPlAd(ip+iq+il)+nCenters*(iCompMat(ip,iq,il)-1)+ &
                          nA*(nA-1)/2+nB-1)
                  sum_b = sum_b+xfac_b*yfac_b*zfac_b*FracB*Work(iAtBoMltPlAd(ip+iq+il)+nCenters*(iCompMat(ip,iq,il)-1)+ &
                          nA*(nA-1)/2+nB-1)
                end do
              end do
            end do
            if (BondMat(nA,nB)) then
              ! If bonding
              ! Do atoms
              Work(iAtMltPlAd(iMltpl)+nAtoms*(iComp-1)+iA-1) = Work(iAtMltPlAd(iMltpl)+nAtoms*(iComp-1)+iA-1)+sum_a
              Work(iAtMltPlAd(iMltpl)+nAtoms*(iComp-1)+nB-1) = Work(iAtMltPlAd(iMltpl)+nAtoms*(iComp-1)+nB-1)+sum_b
            else
              ! If not bonding
              ! Do Atoms
              Work(iAtMltPlAd(iMltpl)+nAtoms*(iComp-1)+iA-1) = Work(iAtMltPlAd(iMltpl)+nAtoms*(iComp-1)+iA-1)+sum_a
              Work(iAtMltPlAd(iMltpl)+nAtoms*(iComp-1)+nB-1) = Work(iAtMltPlAd(iMltpl)+nAtoms*(iComp-1)+nB-1)+sum_b
              ! Do bonds
              Work(iAtBoMltPlAd(iMltpl)+nCenters*(iComp-1)+iA*(iA+1)/2-1) = Work(iAtBoMltPlAd(iMltpl)+nCenters*(iComp-1)+ &
                                                                            iA*(iA+1)/2-1)+sum_a
              Work(iAtBoMltPlAd(iMltpl)+nCenters*(iComp-1)+nB*(nB+1)/2-1) = Work(iAtBoMltPlAd(iMltpl)+nCenters*(iComp-1)+ &
                                                                            nB*(nB+1)/2-1)+sum_b
            end if
          end do
        end do
      end do
    end if
  end do
end do

!call mma_deallocate(DenMat)
call mma_deallocate(Qexp)
call mma_deallocate(iCompMat)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nBas)

end subroutine Get_MpProp
