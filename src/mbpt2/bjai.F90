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

subroutine BJAI(IAD,EPSI,EPSE,E2BJAI,VECL2)

use MBPT2_Global, only: LuIntM
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: IAD(3888)
real(kind=wp), intent(in) :: EPSI(*), EPSE(*)
real(kind=wp), intent(out) :: E2BJAI, VECL2
integer(kind=iwp) :: i, iA, IAD1, IAD13, IAD2, IADA, IADAB, IADB, iB, iI, iJ, ISPQRS, iSymA, iSymB, iSymI, iSymJ, j, LA, LAA, LAB, &
                     LAB1, nExtA, nExtB, nOccI, nOccJ, nOrbA, nOrbB, NRA, NRB, NRI, NRJ
real(kind=wp) :: A, B, CKAL1, CKAL2, dVectorA(255), dVectorB(255), EDIF, EDIFIJ, SKAL1, SKAL2
logical(kind=iwp) :: DoCholesky
real(kind=wp), allocatable :: INT1(:), INT2(:), AIBJ(:), AJBI(:)
logical(kind=iwp), parameter :: Debug = .false.
real(kind=wp), external :: ddot_
#include "corbinf.fh"

SKAL2 = -huge(SKAL2)
IAD13 = 0

call DecideOnCholesky(DoCholesky)
if (Debug) then
  write(u6,*)
  write(u6,'(A,8I3)') '      nOcc:',(nOcc(i),i=1,nSym)
  write(u6,'(A,8I3)') '      nExt:',(nExt(i),i=1,nSym)
  write(u6,'(A,8I3)') '      nOrb:',(nOrb(i),i=1,nSym)
end if

call iDAFILE(LUINTM,2,IAD,3888,IAD13)

VECL2 = One
E2BJAI = Zero
ISPQRS = 0
NRI = 0
do iSymI=1,nSym
  nOccI = nOcc(iSymI)
  NRJ = 0
  do iSymJ=1,iSymI
    nOccJ = nOcc(iSymJ)
    NRA = 0
    do iSymA=1,nSym
      nExtA = nExt(iSymA)
      nOrbA = nOrb(iSymA)
      NRB = 0
      do iSymB=1,iSymA
        nExtB = nExt(iSymB)
        nOrbB = nOrb(iSymB)
        ISPQRS = ISPQRS+1
        if (nOccI*nOccJ*nExtA*nExtB /= 0) then
          if (Mul(iSymI,iSymJ) == Mul(iSymA,iSymB)) then
            IAD1 = IAD(3*ISPQRS-1)
            IAD2 = IAD(3*ISPQRS)

            if (Debug) then
              write(u6,*)
              write(u6,*) ' [BJAI] Integrals <A,B|I,J> : ',iSymA,iSymB,iSymI,iSymJ
            end if
            ! --- iSymI /= iSymJ ---
            if (iSymI /= iSymJ) then
              LAB = nOrbA*nOrbB
              if (DoCholesky) LAB = nExtA*nExtB
              LAB1 = nExtA*nExtB
              call mma_allocate(INT1,LAB,label='INT1')
              call mma_allocate(INT2,LAB,label='INT2')
              call mma_allocate(AIBJ,LAB1,label='AIBJ') ! AB|IJ
              call mma_allocate(AJBI,LAB1,label='AJBI') ! AB|JI
              do iI=1,nOccI
                do iJ=1,nOccJ
                  if (Debug) then
                    write(u6,*)
                    write(u6,*) ' *  i,j = ',iI,iJ
                  end if
                  call dDAFile(LUINTM,2,INT1,LAB,IAD1)
                  call dDAFile(LUINTM,2,INT2,LAB,IAD2)
                  if (Debug) then
                    if (DoCholesky) then
                      call RecPrt('Int1:','(8F10.6)',Int1,nExtB,nExtA)
                      call RecPrt('Int2:','(8F10.6)',Int2,nExtB,nExtA)
                    else
                      call RecPrt('Int1:','(8F10.6)',Int1,nOrbB,nOrbA)
                      call RecPrt('Int2:','(8F10.6)',Int2,nOrbB,nOrbA)
                    end if
                  end if
                  EDIFIJ = EPSI(NRI+iI)+EPSI(NRJ+iJ)
                  IADA = nOcc(iSymA)*nOrbB+nOcc(iSymB)
                  if (DoCholesky) IADA = 0  ! CGG
                  I = 1
                  if (Debug) then
                    write(u6,*)
                  end if
                  do iA=1,nExtA
                    do iB=1,nExtB
                      A = INT1(IADA+iB)
                      B = INT2(IADA+iB)
                      if (Debug) dVectorA(iB) = A
                      if (Debug) dVectorB(iB) = B
                      AIBJ(I) = A+B
                      AJBI(I) = A-B
                      I = I+1
                    end do
                    if (Debug) then
                      write(u6,'(A,I3,6F10.6)') 'A:',iA,(dVectorA(j),j=1,nExtB)
                      write(u6,'(A,I3,6F10.6)') 'B:',iA,(dVectorB(j),j=1,nExtB)
                      write(u6,*) ' -------'
                    end if
                    if (DoCholesky) then
                      IADA = IADA+nExtB
                    else
                      IADA = IADA+nOrbB
                    end if
                  end do
                  I = 1
                  do iA=1,nExtA
                    do iB=1,nExtB
                      EDIF = EPSE(NRA+iA)+EPSE(NRB+iB)-EDIFIJ
                      INT1(I) = AIBJ(I)/EDIF
                      INT2(I) = AJBI(I)/EDIF
                      I = I+1
                    end do
                  end do
                  SKAL1 = DDOT_(LAB1,INT1,1,AIBJ,1)
                  SKAL2 = DDOT_(LAB1,INT2,1,AJBI,1)
                  E2BJAI = E2BJAI-SKAL1-Three*SKAL2
                  CKAL1 = DDOT_(LAB1,INT1,1,INT1,1)
                  CKAL2 = DDOT_(LAB1,INT2,1,INT2,1)
                  VECL2 = VECL2+CKAL1+Three*CKAL2
                end do
              end do
              call mma_deallocate(INT1)
              call mma_deallocate(INT2)
              call mma_deallocate(AIBJ)
              call mma_deallocate(AJBI)
            end if

            ! --- iSymI == iSymJ ---
            if (iSymI == iSymJ) then
              LA = nExtA
              LAB = nOrbA**2
              if (DoCholesky) LAB = nExtA**2
              LAA = (LA*(LA+1))/2
              call mma_allocate(INT1,2*LAB,label='INT1')
              call mma_allocate(AIBJ,2*LAA,label='AIBJ')
              call mma_allocate(AJBI,2*LAA,label='AJBI')
              IADAB = nOcc(iSymA)*(nOrbA+1)
              if (DoCholesky) IADAB = 0
              do iI=1,nOccI
                do iJ=1,iI
                  if (Debug) then
                    write(u6,*)
                    write(u6,*) ' *  i,j = ',iI,iJ
                  end if
                  call dDAFile(LUINTM,2,INT1,LAB,IAD1)
                  if (Debug) then
                    if (DoCholesky) then
                      call RecPrt('Int:','(8F10.6)',Int1,nExtA,nExtA)
                    else
                      call RecPrt('Int:','(8F10.6)',Int1,nOrbA,nOrbA)
                    end if
                  end if
                  EDIFIJ = EPSI(NRI+iI)+EPSI(NRJ+iJ)
                  I = 1
                  IADA = IADAB
                  do iA=1,nExtA
                    IADB = IADAB
                    if (Debug) then
                      write(u6,*) ' iA=',iA,'  IADA:'
                      write(u6,'(8F10.6)') (INT1(IADA+j),j=1,iA)
                    end if
                    do iB=1,iA
                      A = INT1(IADA+iB)
                      B = INT1(IADB+iA)
                      AIBJ(I) = A+B
                      AJBI(I) = A-B
                      I = I+1
                      if (DoCholesky) then
                        IADB = IADB+nExtA
                      else
                        IADB = IADB+nOrbA
                      end if
                    end do
                    if (DoCholesky) then
                      IADA = IADA+nExtA
                    else
                      IADA = IADA+nOrbA
                    end if
                  end do
                  I = 1
                  do iA=0,LA-1
                    do iB=0,iA
                      EDIF = EPSE(NRA+iA+1)+EPSE(NRA+iB+1)-EDIFIJ
                      if (iA == iB) EDIF = Two*EDIF
                      INT1(I) = AIBJ(I)/EDIF
                      INT1(LAA+I) = AJBI(I)/EDIF
                      if ((iA == iB) .and. (iI /= iJ)) then
                        VECL2 = VECL2+Two*INT1(I)**2
                      else if ((iA /= iB) .and. (iI == iJ)) then
                        VECL2 = VECL2+Half*INT1(I)**2
                      else
                        VECL2 = VECL2+INT1(I)**2
                      end if
                      I = I+1
                    end do
                  end do
                  SKAL1 = DDOT_(LAA,INT1,1,AIBJ,1)
                  if (iI /= iJ) then
                    SKAL2 = Three*DDOT_(LAA,INT1(LAA+1),1,AJBI,1)
                  end if
                  if (iI == iJ) E2BJAI = E2BJAI-Half*SKAL1
                  if (iI /= iJ) E2BJAI = E2BJAI-SKAL1-SKAL2
                  if (iI /= iJ) then
                    CKAL2 = DDOT_(LAA,INT1(LAA+1),1,INT1(LAA+1),1)
                    VECL2 = VECL2+Three*CKAL2
                  end if
                end do
              end do
              call mma_deallocate(INT1)
              call mma_deallocate(AIBJ)
              call mma_deallocate(AJBI)
            end if
          end if
        end if
        NRB = NRB+nExtB
      end do
      NRA = NRA+nExtA
    end do
    NRJ = NRJ+nOccJ
  end do
  NRI = NRI+nOccI
end do
VECL2 = sqrt(One/VECL2)

return

end subroutine BJAI
