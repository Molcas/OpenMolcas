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

implicit real*8(A-H,O-Z)

dimension EPSI(*), EPSE(*), IAD(3888)
logical DoCholesky, Debug
data Debug/.false./
dimension dVectorA(255) ! CGG
dimension dVectorB(255) ! CGG
!data Debug /.true./    ! CGG
#include "mxdim.fh"
#include "corbinf.fh"
#include "files_mbpt2.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"

MUL(I,J) = 1+ieor(I-1,J-1)

SKAL2 = -9999999.9
IAD13 = 0
LIADUT = 3888

call DecideOnCholesky(DoCholesky)
if (Debug) then
  write(6,*)
  write(6,'(A,8I3)') '      nOcc:',(nOcc(i),i=1,nSym)
  write(6,'(A,8I3)') '      nExt:',(nExt(i),i=1,nSym)
  write(6,'(A,8I3)') '      nOrb:',(nOrb(i),i=1,nSym)
end if

call iDAFILE(LUINTM,2,IAD,LIADUT,IAD13)

VECL2 = 1.0D+00
E2BJAI = 0.0D+00
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
        if (nOccI*nOccJ*nExtA*nExtB == 0) goto 41
        if (MUL(iSymI,iSymJ) == MUL(iSymA,iSymB)) then
          IAD1 = IAD(3*ISPQRS-1)
          IAD2 = IAD(3*ISPQRS)

          if (Debug) then
            write(6,*)
            write(6,*) ' [BJAI] Integrals <A,B|I,J> : ',iSymA,iSymB,iSymI,iSymJ
          end if
          ! --- iSymI /= iSymJ ---
          if (iSymI /= iSymJ) then
            LAB = nOrbA*nOrbB
            if (DoCholesky) LAB = nExtA*nExtB
            LAB1 = nExtA*nExtB
            call GetMem('INT1','ALLO','REAL',LINT1,LAB)
            call GetMem('INT2','ALLO','REAL',LINT2,LAB)
            call GetMem('AIBJ','ALLO','REAL',LAIBJ,LAB1) ! AB|IJ
            call GetMem('AJBI','ALLO','REAL',LAJBI,LAB1) ! AB|JI
            do iI=1,nOccI
              do iJ=1,nOccJ
                if (Debug) then
                  write(6,*)
                  write(6,*) ' *  i,j = ',iI,iJ
                end if
                call dDAFile(LUINTM,2,Work(LINT1),LAB,IAD1)
                call dDAFile(LUINTM,2,Work(LINT2),LAB,IAD2)
                if (Debug) then
                  if (DoCholesky) then
                    call RecPrt('Int1:','(8F10.6)',Work(lInt1),nExtB,nExtA)
                    call RecPrt('Int2:','(8F10.6)',Work(lInt2),nExtB,nExtA)
                  else
                    call RecPrt('Int1:','(8F10.6)',Work(lInt1),nOrbB,nOrbA)
                    call RecPrt('Int2:','(8F10.6)',Work(lInt2),nOrbB,nOrbA)
                  end if
                end if
                EDIFIJ = EPSI(NRI+iI)+EPSI(NRJ+iJ)
                IADA = nOcc(iSymA)*nOrbB+nOcc(iSymB)
                if (DoCholesky) IADA = 0  ! CGG
                IADB = LINT2+IADA
                IADA = LINT1+IADA
                I = 0
                if (Debug) then
                  write(6,*)
                end if
                do iA=0,nExtA-1
                  do iB=0,nExtB-1
                    A = Work(IADA+iB)
                    B = Work(IADB+iB)
                    if (Debug) dVectorA(iB+1) = A
                    if (Debug) dVectorB(iB+1) = B
                    Work(LAIBJ+I) = A+B
                    Work(LAJBI+I) = A-B
                    I = I+1
                  end do
                  if (Debug) then
                    write(6,'(A,I3,6F10.6)') 'A:',iA+1,(dVectorA(j),j=1,nExtB)
                    write(6,'(A,I3,6F10.6)') 'B:',iA+1,(dVectorB(j),j=1,nExtB)
                    write(6,*) ' -------'
                  end if
                  if (DoCholesky) then
                    IADA = IADA+nExtB
                    IADB = IADB+nExtB
                  else
                    IADA = IADA+nOrbB
                    IADB = IADB+nOrbB
                  end if
                end do
                I = 0
                do iA=1,nExtA
                  do iB=1,nExtB
                    EDIF = EPSE(NRA+iA)+EPSE(NRB+iB)-EDIFIJ
                    Work(LINT1+I) = Work(LAIBJ+I)/EDIF
                    Work(LINT2+I) = Work(LAJBI+I)/EDIF
                    I = I+1
                  end do
                end do
                SKAL1 = DDOT_(LAB1,Work(LINT1),1,Work(LAIBJ),1)
                SKAL2 = DDOT_(LAB1,Work(LINT2),1,Work(LAJBI),1)
                E2BJAI = E2BJAI-SKAL1-3.0D+00*SKAL2
                CKAL1 = DDOT_(LAB1,Work(LINT1),1,Work(LINT1),1)
                CKAL2 = DDOT_(LAB1,Work(LINT2),1,Work(LINT2),1)
                VECL2 = VECL2+CKAL1+3.0D+00*CKAL2
              end do
            end do
            call GetMem('INT1','FREE','REAL',LINT1,LAB)
            call GetMem('INT2','FREE','REAL',LINT2,LAB)
            call GetMem('AIBJ','FREE','REAL',LAIBJ,LAB1)
            call GetMem('AJBI','FREE','REAL',LAJBI,LAB1)
          end if

          ! --- iSymI == iSymJ ---
          if (iSymI == iSymJ) then
            LA = nExtA
            LAB = nOrbA**2
            if (DoCholesky) LAB = nExtA**2
            LAA = (LA*(LA+1))/2
            call GetMem('INT','ALLO','REAL',LINT,2*LAB)
            call GetMem('AIBJ','ALLO','REAL',LAIBJ,2*LAA)
            call GetMem('AJBI','ALLO','REAL',LAJBI,2*LAA)
            IADAB = LINT+nOcc(iSymA)*(nOrbA+1)
            if (DoCholesky) IADAB = LINT
            do iI=1,nOccI
              do iJ=1,iI
                if (Debug) then
                  write(6,*)
                  write(6,*) ' *  i,j = ',iI,iJ
                end if
                call dDAFile(LUINTM,2,Work(LINT),LAB,IAD1)
                if (Debug) then
                  if (DoCholesky) then
                    call RecPrt('Int:','(8F10.6)',Work(lInt),nExtA,nExtA)
                  else
                    call RecPrt('Int:','(8F10.6)',Work(lInt),nOrbA,nOrbA)
                  end if
                end if
                EDIFIJ = EPSI(NRI+iI)+EPSI(NRJ+iJ)
                I = 0
                IADA = IADAB
                do iA=0,nExtA-1
                  IADB = IADAB
                  if (Debug) then
                    write(6,*) ' iA=',iA+1,'  IADA:'
                    write(6,'(8F10.6)') (WORK(IADA+j),j=0,iA)
                  end if
                  do iB=0,iA
                    A = Work(IADA+iB)
                    B = Work(IADB+iA)
                    Work(LAIBJ+I) = A+B
                    Work(LAJBI+I) = A-B
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
                I = 0
                do iA=0,LA-1
                  do iB=0,iA
                    EDIF = EPSE(NRA+iA+1)+EPSE(NRA+iB+1)-EDIFIJ
                    if (iA == iB) EDIF = 2.0D+00*EDIF
                    Work(LINT+I) = Work(LAIBJ+I)/EDIF
                    Work(LINT+LAA+I) = Work(LAJBI+I)/EDIF
                    if ((iA == iB) .and. (iI /= iJ)) then
                      VECL2 = VECL2+2.0D+00*Work(LINT+I)**2
                    else if ((iA /= iB) .and. (iI == iJ)) then
                      VECL2 = VECL2+0.5D+00*Work(LINT+I)**2
                    else
                      VECL2 = VECL2+Work(LINT+I)**2
                    end if
                    I = I+1
                  end do
                end do
                SKAL1 = DDOT_(LAA,Work(LINT),1,Work(LAIBJ),1)
                if (iI /= iJ) then
                  SKAL2 = 3.0D+00*DDOT_(LAA,Work(LINT+LAA),1,Work(LAJBI),1)
                end if
                if (iI == iJ) E2BJAI = E2BJAI-0.5*SKAL1
                if (iI /= iJ) E2BJAI = E2BJAI-SKAL1-SKAL2
                if (iI /= iJ) then
                  CKAL2 = DDOT_(LAA,Work(LINT+LAA),1,Work(LINT+LAA),1)
                  VECL2 = VECL2+3.0D+00*CKAL2
                end if
              end do
            end do
            call GetMem('INT','FREE','REAL',LINT,LAB+LA)
            call GetMem('AIBJ','FREE','REAL',LAIBJ,LAA)
            call GetMem('AJBI','FREE','REAL',LAJBI,LAA)
          end if
        end if
41      NRB = NRB+nExtB
      end do
      NRA = NRA+nExtA
    end do
    NRJ = NRJ+nOccJ
  end do
  NRI = NRI+nOccI
end do
VECL2 = sqrt(1.0D+00/VECL2)

return

end subroutine BJAI
