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
! Copyright (C) 1992,1994, Per Ake Malmqvist                           *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine FMAT_CASPT2(FIFA,nFIFA,FIMO,NFIMO,DREF,NDREF,HONE,nHONE)
! COMPUTE THE INACTIVE FOCK MATRIX IN MO BASIS.
! COMPUTE THE ACTIVE FOCK MATRIX IN MO BASIS.
! DEFINITIONS ARE:
! FIMO(P,Q)= H(P,Q)+ SUM( (2*(PQ,KK)-(PK,QK)) OVER INACTIVE K)
! WHERE H(P,Q) IS THE CORE FOCK MATRIX OBTAINED FROM THE
! TRANSFORMATION PART, I.E., THE CI 1-ELECTRON HAMILTONIAN.
! SIMILARLY,
! FAMO(P,Q)= (1/2)*SUM( (2*(PQ,TU)-(PU,QT))*DREF(T,U) OVER ACTIVE T,U)
! THIS ROUTINE USES DIRECTLY THE TRANSFORMED INTEGRALS TO SET UP
! FIMO AND FAMO.
! CODED 92-12-04 BY MALMQVIST FOR CASPT2, MOLCAS-3 VERSION.

use definitions, only: iwp, wp, u6
use constants, only: Zero, Half, One, Two
use caspt2_global, only: LUINTM
use caspt2_module, only: NSYM, NORB, NISH, NOSH, NAES, NoMx, NoTri
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer(kind=iwp), intent(in) :: nFIFA, NFIMO, NDREF, nHONE
real(kind=wp), intent(inout) :: FIFA(nFIFA), FIMO(NFIMO)
real(kind=wp), intent(in) :: DREF(NDREF), HONE(nHONE)
integer(kind=iwp) IAD2M(3,36*36)
integer(kind=iwp) NDIM2M, IDISK, IFSTA, NBR, NB3, NBNB, ISYS, IS3RS, NIP, NOP, NAESP, ISYQ, IS3PQ, ISADDR, IDISK1, ITABS, ITU, &
                  IUABS, nInts, nBUF
real(kind=wp) DTU
real(kind=wp), allocatable :: BUF(:)

! Inactive and active Fock matrices:
FIMO(:) = Zero
FIFA(:) = Zero ! ued temporarily as FAMO

! notri=Size of an array with symmetry-blocked triangular
! submatrices, using non-frozen, non-deleted MO indices.
! NBUF=Max size of a LUINTM buffer.
NBUF = max(NOMX**2,notri)
call mma_allocate(BUF,NBUF,Label='BUF')

NDIM2M = (NSYM*(NSYM+1))/2
IDISK = 0
call IDAFILE(LUINTM,2,IAD2M,3*36*36,IDISK)

call Do_Loops(1) ! Do Coulomb contributions
call Do_Loops(2) ! Do exchange contributions

call mma_deallocate(BUF)

! FIMO comes from contractions over inactive orbitals, while FAMO from
! contractions over active orbitals and therefore are summed up together
! here into FIFA.

FIMO(1:notri) = FIMO(1:notri)+HONE(1:notri)
FIFA(1:notri) = FIFA(1:notri)+FIMO(1:notri)

contains

subroutine Do_Loops(icase)

  integer(kind=iwp), intent(in) :: iCase
  integer(kind=iwp) ISYR, ISYP, IT, IU

  IFSTA = 1
  do ISYR=1,NSYM
    NBR = NORB(ISYR)
    if (NBR == 0) cycle
    NB3 = (NBR**2+NBR)/2
    NBNB = NBR**2
    ISYS = ISYR
    IS3RS = (ISYR**2+ISYR)/2
    do ISYP=1,NSYM
      NIP = NISH(ISYP)
      NOP = NOSH(ISYP)
      NAESP = NAES(ISYP)
      if (NOP == 0) cycle
      ISYQ = ISYP
      IS3PQ = (ISYP*(ISYP+1))/2
      ISADDR = IS3RS+NDIM2M*(IS3PQ-1)
      IDISK1 = IAD2M(iCase,ISADDR)
      if (IDISK1 == 0) then
        if (iCase == 1) then
          write(u6,*) ' FMAT: COULOMB INTEGRAL BUFFER MISSING!'
        else
          write(u6,*) ' FMAT: EXCH-1  INTEGRAL BUFFER MISSING!'
        end if
        write(u6,'(1X,A,4I4)') 'SYMMETRY BLOCK:',ISYP,ISYQ,ISYR,ISYS
        call ABEND()
      end if
      if (iCase == 1) then
        nInts = NB3
      else
        nInts = NBNB
      end if

      do IT=1,NOP
        do IU=1,IT

          ! Process here contributions to FIMO or FAMO.

          if ((IT <= NIP) .and. (IT == IU)) then
            call DDAFILE(LUINTM,2,BUF,nInts,IDISK1)
            if (iCase == 1) then
              ! ADD 2* BUFFER INTO CORRECT POSITION OF  FIMO.
              call DAXPY_(NB3,Two,BUF,1,FIMO(IFSTA),1)
            else
              call TRIANG(NBR,BUF)
              ! ADD -1* BUFFER INTO CORRECT POSITION OF  FIMO.
              call DAXPY_(NB3,-One,BUF,1,FIMO(IFSTA),1)
            end if
          else if ((IT > NIP) .and. (IT <= NOP) .and. (IU > NIP) .and. (IU <= NOP)) then
            ITABS = NAESP+(IT-NIP)
            IUABS = NAESP+(IU-NIP)
            ITU = (ITABS*(ITABS-1))/2+IUABS
            DTU = DREF(ITU)
            if (IT == IU) DTU = half*DTU
            call DDAFILE(LUINTM,2,BUF,nInts,IDISK1)
            if (iCase == 1) then
              call DAXPY_(NB3,Two*DTU,BUF,1,FIFA(IFSTA),1)
            else
              call TRIANG(NBR,BUF)
              call DAXPY_(NB3,-DTU,BUF,1,FIFA(IFSTA),1)
            end if
          else
            call DDAFILE(LUINTM,0,BUF,nInts,IDISK1)
          end if

        end do
      end do

    end do
    IFSTA = IFSTA+NB3
  end do
end subroutine Do_Loops

end subroutine FMAT_CASPT2
