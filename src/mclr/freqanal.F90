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

use stdalloc, only: mma_allocate, mma_deallocate
use input_mclr, only: nSym, nDisp, nUserPT, nSRot, UserP, ChIrr, UserT
use temperatures, only: DefTemp

implicit none
integer nDeg(*), nrvec(*)
real*8 H(*)
logical converged(8)
real*8 elec(*)
integer iel(3)
real*8 elout(*)
integer ldisp(nsym), Lu_10
! local variables
logical Do_Molden
real*8, allocatable :: NMod(:), EVec(:), EVec2(:,:), EVal(:), EVal2(:), Intens(:), RedMas(:), Tmp3(:), Temp(:)
integer ipNx, nModes, lModes, i1, i3, j, ii, lnm_molpac, iSym, nx, iCtl, ll, kk, i, k, iT, jpNx, ix, jx, nEig, iNeg
integer, external :: IsFreeUnit
real*8 Tmp, Fact, rNorm

call mma_allocate(NMod,nDisp**2,Label='NMod')
call mma_allocate(EVec,nDisp**2,Label='EVec')
call mma_allocate(EVec2,2,nDisp**2,Label='EVec2')
call mma_allocate(EVal,nDisp,Label='EVal')
call mma_allocate(Intens,nDisp*2,Label='Intens')
call mma_allocate(RedMas,nDisp,Label='RedMas')
ipNx = 1
nModes = 0
lModes = 0

write(6,*)
write(6,*) '     ************************************'
write(6,*) '     *                                  *'
write(6,*) '     * Harmonic frequencies in cm-1     *'
write(6,*) '     * Intensities in km/mole           *'
write(6,*) '     *                                  *'
write(6,*) '     * No correction due to curvilinear *'
write(6,*) '     * representations has been done    *'
write(6,*) '     *                                  *'
write(6,*) '     ************************************'
write(6,*)
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
    write(6,*)
    write(6,*) '   Symmetry ',chirr(isym)
    write(6,*) '  =============='
    write(6,*)

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
            tmp = 0.0d0
            do it=0,nx-1
              Fact = sqrt(dble(nDeg(i1+it)))
              tmp = tmp+EVec(1+(k-1)*nx+it)*elec(ii+it)*Fact
            end do
            elout(j) = tmp
          end do
          ii = ii+nx
        end if
      end do
      write(Lu_10,'(A,I1)') '*NORMAL MODES SYMMETRY: ',isym

      ! Save normal modes for later generation of Molden input.

      call dcopy_(nX**2,EVec,1,NMod(ipNx),1)
      jpNx = ipNx

      do iX=1,nX

        ! Transform from mass-weighted cartesian to cartesian for Molden.

        rNorm = 0.0d0
        do jX=0,nX-1
          Fact = sqrt(dble(nDeg(jX+1)))
          NMod(ipNx+jX) = NMod(ipNx+jX)/Fact
          rNorm = rNorm+dble(nDeg(jX+1))*NMod(ipNx+jX)**2
        end do
        call DScal_(nX,1.0d0/sqrt(rNorm),NMod(ipNx),1)

        ipNx = ipNx+nX
        lModes = lModes+nX
      end do
      nModes = nModes+nX
      call dcopy_(nX**2,NMod(jpNx),1,EVec,1)
      call GF_Print(EVal(i1),EVec,elout(kk),ll,nX,nX,iCtl,Intens(i1),RedMas,Lu_10,i1-1)
    else
      write(6,*)
      write(6,*) '     NOT CONVERGED'
      write(6,*)
      do i=1,3
        if (iel(i) == isym) then
          j = j+1
          ii = ii+nx
          elout(j) = -99999999d0
        end if
      end do
      Do_Molden = .false.
    end if
  end if
  i3 = i3+nx*(nx+1)/2
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
call dcopy_(nEig,Eval,1,Temp,1)

! For verification purpose we skip frequencies close to zero.

do i=1,nEig
  if (abs(Temp(i)) < 5.0d0) Temp(i) = 0.0d0
end do
call Add_Info('Harm_Freq',Temp,nEig,1)
call mma_deallocate(Temp)
!
do i=1,nEig
  if (abs(Intens(i)) < 1.0d0) Intens(i) = 0.0d0
end do
call Add_Info('IR_Intensities',Intens,nEig,1)
!                                                                      *
!***********************************************************************
!                                                                      *
write(Lu_10,'(A)') '*END NORMAL MODES'

! Calculate thermodynamic properties----------

if ((nUserPT == 0) .and. (nsRot == 0)) then
  UserP = 1.0d0
  nUserPT = size(DefTemp)
  do i=1,nUserPT
    UserT(i) = DefTemp(i)
  end do
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
