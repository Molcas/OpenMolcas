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
! Copyright (C) Chad E. Hoyer                                          *
!***********************************************************************
!  DQVDiabat
!
!> @brief
!>   Compute diabats with DQV
!> @author C. E. Hoyer
!>
!> @details
!> This subroutine takes in properties that have been
!> computed by RASSI and uses them to compute diabatic
!> states and thus diabatic energies and couplings.
!> Currently, the user must compute x, y, z, xx,
!> yy, zz, and 1/r, and they must be computed in that
!> order.
!>
!> @param[in] PROP Properties computed in RASSI
!> @param[in] HAM
!***********************************************************************

subroutine DQVDiabat(PROP,HAM)

use Cntrl, only: AlphZ, BetaE, ICOMP, NPROP, NSTATE, PNAME
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half, Quart, Pi
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: PROP(NSTATE,NSTATE,NPROP), HAM(NSTATE,NSTATE)
integer :: I, IPROP, ISTA, JSTA, K, PNUM(7)
real(kind=wp) :: ATerm, BTerm, Chng, CosO, CTerm, RotAngF, RotAngO, SinO, ThrSch, TII, TIJ, TJJ
real(kind=wp), allocatable :: HDIA(:,:), TI(:), TJ(:), TMP(:,:), TROT(:,:), TRQ(:,:)
integer(kind=iwp), parameter :: itMAX = 50
real(kind=wp), parameter :: MTE = 1.0e-8_wp, MTF = 1.0e-14_wp, THRS = 1.0e-8_wp

! Printing some stuff

write(u6,*)
write(u6,*)
write(u6,'(6X,A)') repeat('*',100)
write(u6,'(6X,A,98X,A)') '*','*'
write(u6,'(6X,A,33X,A,34X,A)') '*',' The DQV Diabatization Section ','*'
write(u6,'(6X,A,98X,A)') '*','*'
write(u6,'(6X,A)') repeat('*',100)
write(u6,*)
write(u6,*)

call CollapseOutput(1,'DQV Diabatization section')
write(u6,'(3X,A)') '-------------------------'

! Find the properties we need
PNUM(:) = 0
do IPROP=1,NPROP
  if (PNAME(IPROP) == 'MLTPL  1') then
    if (ICOMP(IPROP) == 1) PNUM(1) = IPROP
    if (ICOMP(IPROP) == 2) PNUM(2) = IPROP
    if (ICOMP(IPROP) == 3) PNUM(3) = IPROP
  end if
  if (PNAME(IPROP) == 'MLTPL  2') then
    if (ICOMP(IPROP) == 1) PNUM(4) = IPROP
    if (ICOMP(IPROP) == 4) PNUM(5) = IPROP
    if (ICOMP(IPROP) == 6) PNUM(6) = IPROP
  end if
  if (PNAME(IPROP) == 'EF0    1') then
    if (ICOMP(IPROP) == 1) PNUM(7) = IPROP
  end if
end do
do IPROP=1,7
  if (PNUM(IPROP) == 0) then
    write(u6,*) 'DQVDiabat: Some required properties are not available'
    call Abend()
  end if
end do

! CEH Now we maximize f_DQV through a Jacobi sweep algorithm

! Initialize rotation matrix to unit matrix (theta=0)
call mma_allocate(TROT,NSTATE,NSTATE,Label='TROT')
call unitmat(TROT,NSTATE)

!Now, we have to form the trace of the quadrupole
!The user has to specify dipole then quadrupole
call mma_allocate(TRQ,NSTATE,NSTATE,Label='TRQ')
TRQ(:,:) = PROP(:,:,PNUM(4))+PROP(:,:,PNUM(5))+PROP(:,:,PNUM(6))

call RecPrt('The TRQ matrix in DQVDiabat','',TRQ,NSTATE,NSTATE)
write(u6,*)

call mma_allocate(TI,NSTATE,Label='TI')
call mma_allocate(TJ,NSTATE,Label='TJ')
do i=1,itMAX
  THRSCH = Zero

  do ISTA=1,NSTATE
    do JSTA=ISTA+1,NSTATE

      ! Equations 9 and 10 of Kleir et al. with \mu written
      ! as a vector in addition to the quadrupole trace

      ATERM = PROP(ISTA,JSTA,PNUM(1))**2+PROP(ISTA,JSTA,PNUM(2))**2+PROP(ISTA,JSTA,PNUM(3))**2+ALPHZ*TRQ(ISTA,JSTA)**2+ &
              BETAE*PROP(ISTA,JSTA,PNUM(7))**2-Quart*(PROP(ISTA,ISTA,PNUM(1))**2+PROP(ISTA,ISTA,PNUM(2))**2+ &
              PROP(ISTA,ISTA,PNUM(3))**2+ALPHZ*TRQ(ISTA,ISTA)**2+BETAE*PROP(ISTA,ISTA,PNUM(7))**2)- &
              Quart*(PROP(JSTA,JSTA,PNUM(1))**2+PROP(JSTA,JSTA,PNUM(2))**2+PROP(JSTA,JSTA,PNUM(3))**2+ALPHZ*TRQ(JSTA,JSTA)**2+ &
              BETAE*PROP(JSTA,JSTA,PNUM(7))**2)+Half*(PROP(ISTA,ISTA,PNUM(1))*PROP(JSTA,JSTA,PNUM(1))+ &
              PROP(ISTA,ISTA,PNUM(2))*PROP(JSTA,JSTA,PNUM(2))+PROP(ISTA,ISTA,PNUM(3))*PROP(JSTA,JSTA,PNUM(3))+ &
              ALPHZ*TRQ(ISTA,ISTA)*TRQ(JSTA,JSTA)+BETAE*PROP(ISTA,ISTA,PNUM(7))*PROP(JSTA,JSTA,PNUM(7)))

      BTERM = PROP(ISTA,JSTA,PNUM(1))*(PROP(ISTA,ISTA,PNUM(1))-PROP(JSTA,JSTA,PNUM(1)))+ &
              PROP(ISTA,JSTA,PNUM(2))*(PROP(ISTA,ISTA,PNUM(2))-PROP(JSTA,JSTA,PNUM(2)))+ &
              PROP(ISTA,JSTA,PNUM(3))*(PROP(ISTA,ISTA,PNUM(3))-PROP(JSTA,JSTA,PNUM(3)))+ &
              ALPHZ*TRQ(ISTA,JSTA)*(TRQ(ISTA,ISTA)-TRQ(JSTA,JSTA))+BETAE*PROP(ISTA,JSTA,PNUM(7))* &
              (PROP(ISTA,ISTA,PNUM(7))-PROP(JSTA,JSTA,PNUM(7)))

      if (abs(BTERM) > MTF) then
        if ((ATERM > MTE) .or. (ATERM < -MTE)) then
          ! This part of the code computes the rotation angle.
          CTERM = -BTERM/ATERM
          ROTANGF = atan(CTERM)

          if ((ROTANGF < Zero) .and. (BTERM > Zero)) ROTANGF = ROTANGF+PI

          if ((ROTANGF > Zero) .and. (BTERM < Zero)) ROTANGF = ROTANGF-PI

          ROTANGO = ROTANGF*Quart

          if ((ROTANGO > PI*Quart) .or. (ROTANGO < -PI*Quart)) write(u6,*) 'Quadrants 2 and 3'

          COSO = cos(ROTANGO)
          SINO = sin(ROTANGO)

        else

          write(u6,*) 'A is pretty close to zero.  Something is fishy.'
          COSO = Zero
          SINO = Zero

        end if

        !Equation 13 of Kleir et al.
        CHNG = ATERM+sqrt(ATERM**2+BTERM**2)

        if (CHNG > THRSCH) THRSCH = CHNG

        if (CHNG > MTE) then

          !Update T rotation matrix
          TI(:) = COSO*TROT(:,ISTA)+SINO*TROT(:,JSTA)
          TJ(:) = COSO*TROT(:,JSTA)-SINO*TROT(:,ISTA)
          TROT(:,ISTA) = TI(:)
          TROT(:,JSTA) = TJ(:)

          ! Rotates transition dipole matrices
          ! These equations can be found in a book on Jacobi sweeps

          do k=1,3
            TII = PROP(ISTA,ISTA,PNUM(k))
            TJJ = PROP(JSTA,JSTA,PNUM(k))
            TIJ = PROP(ISTA,JSTA,PNUM(k))
            TI(:) = COSO*PROP(:,ISTA,PNUM(k))+SINO*PROP(:,JSTA,PNUM(k))
            TJ(:) = COSO*PROP(:,JSTA,PNUM(k))-SINO*PROP(:,ISTA,PNUM(k))
            PROP(ISTA,:,PNUM(k)) = TI(:)
            PROP(:,ISTA,PNUM(k)) = TI(:)
            PROP(JSTA,:,PNUM(k)) = TJ(:)
            PROP(:,JSTA,PNUM(k)) = TJ(:)
            PROP(ISTA,JSTA,PNUM(k)) = (COSO**2-SINO**2)*TIJ+COSO*SINO*(TJJ-TII)
            PROP(JSTA,ISTA,PNUM(k)) = PROP(ISTA,JSTA,PNUM(k))
            PROP(ISTA,ISTA,PNUM(k)) = COSO**2*TII+SINO**2*TJJ+Two*COSO*SINO*TIJ
            PROP(JSTA,JSTA,PNUM(k)) = SINO**2*TII+COSO**2*TJJ-Two*COSO*SINO*TIJ
          end do

          ! Rotates trace of quadrupole matrix

          TII = TRQ(ISTA,ISTA)
          TJJ = TRQ(JSTA,JSTA)
          TIJ = TRQ(ISTA,JSTA)
          TI(:) = COSO*TRQ(:,ISTA)+SINO*TRQ(:,JSTA)
          TJ(:) = COSO*TRQ(:,JSTA)-SINO*TRQ(:,ISTA)
          TRQ(ISTA,:) = TI(:)
          TRQ(:,ISTA) = TI(:)
          TRQ(JSTA,:) = TJ(:)
          TRQ(:,JSTA) = TJ(:)
          TRQ(ISTA,JSTA) = (COSO**2-SINO**2)*TIJ+COSO*SINO*(TJJ-TII)
          TRQ(JSTA,ISTA) = TRQ(ISTA,JSTA)
          TRQ(ISTA,ISTA) = COSO**2*TII+SINO**2*TJJ+Two*COSO*SINO*TIJ
          TRQ(JSTA,JSTA) = SINO**2*TII+COSO**2*TJJ-Two*COSO*SINO*TIJ

          ! Rotates electric potential

          k = 7
          TII = PROP(ISTA,ISTA,PNUM(k))
          TJJ = PROP(JSTA,JSTA,PNUM(k))
          TIJ = PROP(ISTA,JSTA,PNUM(k))
          TI(:) = COSO*PROP(:,ISTA,PNUM(k))+SINO*PROP(:,JSTA,PNUM(k))
          TJ(:) = COSO*PROP(:,JSTA,PNUM(k))-SINO*PROP(:,ISTA,PNUM(k))
          PROP(ISTA,:,PNUM(k)) = TI(:)
          PROP(:,ISTA,PNUM(k)) = TI(:)
          PROP(JSTA,:,PNUM(k)) = TJ(:)
          PROP(:,JSTA,PNUM(k)) = TJ(:)
          PROP(ISTA,JSTA,PNUM(k)) = (COSO**2-SINO**2)*TIJ+COSO*SINO*(TJJ-TII)
          PROP(JSTA,ISTA,PNUM(k)) = PROP(ISTA,JSTA,PNUM(k))
          PROP(ISTA,ISTA,PNUM(k)) = COSO**2*TII+SINO**2*TJJ+Two*COSO*SINO*TIJ
          PROP(JSTA,JSTA,PNUM(k)) = SINO**2*TII+COSO**2*TJJ-Two*COSO*SINO*TIJ

        end if
      end if
    end do
  end do

  if (THRSCH < THRS) then
    write(u6,*) 'Converged in this many iterations:'
    write(u6,*) i
    write(u6,*)
    exit
  end if

end do
if (i > itMAX) write(u6,*) 'Max iterations'
call mma_deallocate(TRQ)
call mma_deallocate(TI)
call mma_deallocate(TJ)

call RecPrt('Diabatic Coefficients','',TROT,NSTATE,NSTATE)

call mma_allocate(TMP,NSTATE,NSTATE,Label='TMP')
TMP(:,:) = TROT(:,:)**2
call RecPrt('Weights of adiabatic states','',TMP,NSTATE,NSTATE)

call DGEMM_('T','N',NSTATE,NSTATE,NSTATE,One,TROT,NSTATE,HAM,NSTATE,Zero,TMP,NSTATE)

call mma_allocate(HDIA,NSTATE,NSTATE,Label='HDIA')
call DGEMM_('N','N',NSTATE,NSTATE,NSTATE,One,TMP,NSTATE,TROT,NSTATE,Zero,HDIA,NSTATE)
call mma_deallocate(TMP)
call mma_deallocate(TROT)

write(u6,*)
write(u6,*) 'The eigenvectors may no longer be in'
write(u6,*) 'energetic order'
write(u6,*)
call RecPrt('Diabatic Hamiltonian','',HDIA,NSTATE,NSTATE)

call CollapseOutput(0,'DQV Diabatization section')
write(u6,*)

! Molcas verify calls
call Add_Info('Adiabat1',HAM(1,1),1,4)
call Add_Info('Adiabat2',HAM(2,2),1,4)
call Add_Info('Adiabat3',HAM(3,3),1,4)
call Add_Info('DQVHam11',HDIA(1,1),1,4)
call Add_Info('DQVHam12',HDIA(1,2),1,4)
call Add_Info('DQVHam13',HDIA(1,3),1,4)
call Add_Info('DQVHam22',HDIA(2,2),1,4)
call Add_Info('DQVHam23',HDIA(2,3),1,4)
call Add_Info('DQVHam33',HDIA(3,3),1,4)
! End Molcas verify calls

call mma_deallocate(HDIA)

end subroutine DQVDiabat
