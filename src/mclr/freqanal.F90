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

subroutine Freqanal(nDeg,nrvec,H,converged,ELEC,iel,elout,ldisp,Lu_10)

use Index_Functions, only: nTri_Elem
use input_mclr, only: ChIrr, nDisp, nSRot, nSym, nUserPT, UserP, UserT
use temperatures, only: DefTemp
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Five
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nDeg(*), nrvec(*), iel(3), ldisp(nsym), Lu_10
real(kind=wp), intent(in) :: H(*), elec(*)
logical(kind=iwp), intent(in) :: converged(8)
real(kind=wp), intent(_OUT_) :: elout(*)
integer(kind=iwp) :: i, i1, i3, iCtl, ii, iNeg, ipNx, iSym, iT, ix, j, jpNx, jx, k, kk, ll, lModes, lnm_molpac, nEig, nModes, nx
real(kind=wp) :: Fact, rNorm, Tmp
logical(kind=iwp) :: Do_Molden
real(kind=wp), allocatable :: EVal(:), EVal2(:), EVec(:), EVec2(:,:), Intens(:), NMod(:), RedMas(:), Temp(:), Tmp3(:)
integer(kind=iwp), external :: IsFreeUnit

call mma_allocate(NMod,nDisp**2,Label='NMod')
call mma_allocate(EVec,nDisp**2,Label='EVec')
call mma_allocate(EVec2,2,nDisp**2,Label='EVec2')
call mma_allocate(EVal,nDisp,Label='EVal')
call mma_allocate(Intens,nDisp*2,Label='Intens')
call mma_allocate(RedMas,nDisp,Label='RedMas')
ipNx = 1
nModes = 0
lModes = 0

write(u6,*)
write(u6,*) '     ************************************'
write(u6,*) '     *                                  *'
write(u6,*) '     * Harmonic frequencies in cm-1     *'
write(u6,*) '     * Intensities in km/mole           *'
write(u6,*) '     *                                  *'
write(u6,*) '     * No correction due to curvilinear *'
write(u6,*) '     * representations has been done    *'
write(u6,*) '     *                                  *'
write(u6,*) '     ************************************'
write(u6,*)
i1 = 1
i3 = 1
j = 0
ii = 1
write(Lu_10,'(A)') '*PERTURBATIONS'
write(Lu_10,*) ldisp
write(Lu_10,'(A)') '*BEGIN NORMAL MODES'
write(Lu_10,'(A)') '*NOTICE THAT THEY ARE SYMMETRY ADAPTED'
write(Lu_10,'(A)') '*USING ORTHOGONAL TRANSFORMATIONS'
write(Lu_10,'(A)') '*AND NOT ALASKA TYPE'

! open normal mode file for "normal mode molpac" ! yma
lnm_molpac = 60
lnm_molpac = isFreeUnit(lnm_molpac)
call Molcas_Open(lnm_molpac,'normal_modes_molpac')

Do_Molden = .true.
do iSym=1,nSym
  nX = ldisp(isym)
  if (nX /= 0) then
    write(u6,*)
    write(u6,*) '   Symmetry ',chirr(isym)
    write(u6,*) '  =============='
    write(u6,*)

    if (converged(isym)) then
      call mma_allocate(EVal2,2*nX,Label='EVal2')
      call mma_allocate(Tmp3,nX**2,Label='Tmp3')
      call FREQ(nX,H(i3),nDeg(i1),nrvec(i1),Tmp3,EVec2,EVal2,RedMas,iNeg)
      EVec(:) = EVec2(1,:)
      EVal(i1:i1+nX-1) = EVal2(1:nX)
      call mma_deallocate(EVal2)
      call mma_deallocate(Tmp3)

      iCtl = 0
      ll = 0
      kk = j+1
      do i=1,3
        if (iel(i) == isym) then
          iCtl = 1
          ll = ll+1
          do k=1,nx
            j = j+1
            tmp = Zero
            do it=0,nx-1
              Fact = sqrt(real(nDeg(i1+it),kind=wp))
              tmp = tmp+EVec(1+(k-1)*nx+it)*elec(ii+it)*Fact
            end do
            elout(j) = tmp
          end do
          ii = ii+nx
        end if
      end do
      write(Lu_10,'(A,I1)') '*NORMAL MODES SYMMETRY: ',isym

      ! Save normal modes for later generation of Molden input.

      NMod(ipNx:ipNx+nX**2-1) = EVec(1:nX**2)
      jpNx = ipNx

      do iX=1,nX

        ! Transform from mass-weighted cartesian to cartesian for Molden.

        rNorm = Zero
        do jX=0,nX-1
          Fact = sqrt(real(nDeg(jX+1),kind=wp))
          NMod(ipNx+jX) = NMod(ipNx+jX)/Fact
          rNorm = rNorm+real(nDeg(jX+1),kind=wp)*NMod(ipNx+jX)**2
        end do
        NMod(ipNx:ipNx+nX-1) = NMod(ipNx:ipNx+nX-1)/sqrt(rNorm)

        ipNx = ipNx+nX
        lModes = lModes+nX
      end do
      nModes = nModes+nX
      EVec(1:nX**2) = NMod(jpNx:jpNx+nX**2-1)
      call GF_Print(EVal(i1),EVec,elout(kk),ll,nX,nX,iCtl,Intens(i1),RedMas,Lu_10,i1-1)
    else
      write(u6,*)
      write(u6,*) '     NOT CONVERGED'
      write(u6,*)
      do i=1,3
        if (iel(i) == isym) then
          j = j+1
          ii = ii+nx
          elout(j) = -99999999.0_wp
        end if
      end do
      Do_Molden = .false.
    end if
  end if
  i3 = i3+nTri_Elem(nx)
  i1 = i1+nx
end do
nEig = i1-1

! close the normal mode file
close(lnm_molpac) ! This is for normal_modes_molpac -- yma
!call NM_MOPAC_SNF(nsym,ldisp,natoms) ! f90 not support ....

if (nsym == 1) call Print_Mode_Components(NMod,EVal,nModes,lModes,lDisp)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(Temp,nEig,Label='Temp')
Temp(:) = Eval(1:nEig)

! For verification purpose we skip frequencies close to zero.

do i=1,nEig
  if (abs(Temp(i)) < Five) Temp(i) = Zero
end do
call Add_Info('Harm_Freq',Temp,nEig,1)
call mma_deallocate(Temp)

do i=1,nEig
  if (abs(Intens(i)) < One) Intens(i) = Zero
end do
call Add_Info('IR_Intensities',Intens,nEig,1)
!                                                                      *
!***********************************************************************
!                                                                      *
write(Lu_10,'(A)') '*END NORMAL MODES'

! Calculate thermodynamic properties----------

if ((nUserPT == 0) .and. (nsRot == 0)) then
  UserP = One
  nUserPT = size(DefTemp)
  UserT(1:nUserPT) = DefTemp(1:nUserPT)
  !call ThermoData(EVal,nEig)
end if
call Thermo_Driver(UserT,UserP,nUserPT,nsRot,EVal,nEig,.false.)

! Write stuff on Molden input file

if (Do_Molden) call Freq_Molden(EVal,nModes,NMod,lModes,nSym,Intens,lDisp,RedMas)

call mma_deallocate(NMod)
call mma_deallocate(evec)
call mma_deallocate(evec2)
call mma_deallocate(eval)
call mma_deallocate(intens)
call mma_deallocate(redmas)

end subroutine Freqanal
