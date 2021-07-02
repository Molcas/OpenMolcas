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
! Copyright (C) 1998, Jun-ya Hasegawa                                  *
!               2004, Giovanni Ghigo                                   *
!***********************************************************************

subroutine RDINT2(IPRX,DoTCVA)
! SECOND ORDER TWO-ELECTRON TRANFORMATION PROGRAM. TEST SECTION
!
! THIS SUBROUTINE READS AND CHECKS THE RESULT OF THE SECOND ORDER
! TWO-ELECTRON TRANSFORMATION PROGRAM TRA2. IT CAN BE CALLED BY
! TR2CTL IMMEDIATELY AFTER THE CALL TO TRA2
!  IPRX=0 Do not print Coulomb nor Exchange Integrals
!  IPRX=1 Print only Coulomb Integrals
!  IPRX=2 Print all Integrals
!  IPRX=3 Print only Exchange Integrals
!
! Modified for caspt2 by Jun-ya Hasegawa(98/08/18)
! Modified again by G. Ghigo - September-December 2004
!
! A,B are MO indices, counting only non-frozen and non-deleted.
! T,U are occupied MO indices, only non-frozen and non-deleted.
!
! <AB/TU> ARE ALWAYS GENERATED
! EXCHANGE INTEGRALS <AT/BU> ARE GENERATED AS FOLLOWS:
! <AT/BU> IF ISP >= ISR
! <AT/UB> IF ISP > ISS AND ISP /= ISQ
! <TA/BU> IF ISQ > ISR AND ISP /= ISQ
! <TA/UB> IF ISQ >= ISS AND ISP /= ISQ
!
! IAD2M CONTAINS START ADRESS FOR EACH TYPE OF INTEGRALS:
! IAD2M(1,iSymIJAB)   COULOMB INTEGRALS <AB|TU>
! IAD2M(2,iSymIJAB)   EXCHANGE INTEGRALS <AB|TU> FOR SYM T > SYM U
! IAD2M(3,iSymIJAB)   EXCHANGE INTEGRALS <AB|TU> FOR SYM T < SYM U
! THE LAST ADRESS IS ZERO IF SYM T = SYM U

implicit integer(i-n)
implicit real*8(A-H,O-Z)
parameter(MxFro=50)
parameter(MxDel=50)
#include "rasdim.fh"
#include "caspt2.fh"
#include "trafo.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
real*8, allocatable :: Tmp(:)
dimension IAD2M(3*36*36)
logical DoTCVA
logical Found

! GG-Dec04  The following informations must be passed to the Cholesky
! transformation section through RunFile. COMMON blocks could not be
! used due to several conflicts.
call Get_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)
call Get_iArray('nFroPT',nFro,nSym)
call Get_iArray('nDelPT',nDel,nSym)
call Get_iArray('nIsh',nIsh,nSym)
! PAM 2007 call Get_iArray('nAsh',nAsh,nSym)
! Replaced by:
do i=1,nSym
  nAsh(i) = 0
end do
call qpg_iArray('nAsh',Found,nData)
if (Found .and. (nData == nSym)) then
  call Get_iArray('nAsh',nAsh,nSym)
end if
! End of replacement
do i=1,nSym
  nOrb(i) = nBas(i)-nFro(i)-nDel(i)
  nIsh(i) = nIsh(i)-nFro(i)
  nOsh(i) = nIsh(i)+nAsh(i)
  nSsh(i) = nOrb(i)-nOsh(i)
end do

write(6,*)
write(6,*) 'SECOND ORDER TWO-ELECTRON TRANFORMATION PROGRAM. TEST SECTION:'
write(6,*) ' A,B are MO-symmetry indices, counting only non-frozen and non-deleted.'
write(6,*) ' I,J are occupied MO-symmetry indices, only non-frozen and non-deleted.'
write(6,*) ' i,j are occupied MO indices'
!write(6,*) ' <AB/TU> COULOMB INTEGRALS ARE ALWAYS GENERATED.'
!write(6,*) ' <AB|TU> EXCHANGE INTEGRALS ARE GENERATED AS FOLLOWS:'
!write(6,*) '   <AT/BU> IF ISP >= ISR'
!write(6,*) '   <AT/UB> IF ISP > ISS AND ISP /= ISQ'
!write(6,*) '   <TA/BU> IF ISQ > ISR AND ISP /= ISQ'
!write(6,*) '   <TA/UB> IF ISQ >= ISS AND ISP /= ISQ'
!write(6,*) ' IAD2M(1,iSymIJAB)  COULOMB INTEGRALS <AB|TU>'
!write(6,*) ' IAD2M(2,iSymIJAB)  EXCHANGE INTEGRALS <AB|TU> FOR SYM T > SYM U'
!write(6,*) ' IAD2M(3,iSymIJAB)  EXCHANGE INTEGRALS <AB|TU> FOR SYM T < SYM U'
!write(6,*) ' THE LAST ADRESS IS ZERO IF SYM T = SYM U'
write(6,*)
write(6,'(A,8I3)') '        Symmetries :',(i,i=1,nSym)
write(6,*)
write(6,'(A,8I3)') '           Frozen  :',(nFro(i),i=1,nSym)
write(6,'(A,8I3)') '      Inactive (I) :',(nIsh(i),i=1,nSym)
write(6,'(A,8I3)') '        Active (A) :',(nAsh(i),i=1,nSym)
write(6,'(A,8I3)') '     Secondary (S) :',(nSsh(i),i=1,nSym)
write(6,'(A,8I3)') '          Deleted  :',(nDel(i),i=1,nSym)
write(6,*)
write(6,'(A,8I3)') '  Total correlated :',(nOrb(i),i=1,nSym)
call XFlush(6)

! READ ADDRESS RECORD ON UNIT LUINTM

IAD13 = 0
LIADUT = 3888
call iDAFILE(LUINTM,2,IAD2M,LIADUT,IAD13)

! LOOP OVER QUADRUPLES OF SYMMETRIES (iSymI,iSymJ,iSymA,iSymB)

iSymIJAB = 0
LTotInt = 0
LTotEx1 = 0
LTotEx2 = 0
do iSymI=1,NSYM
  nOrbI = nOrb(iSymI)
  nOccI = nOsh(iSymI)
  do iSymJ=1,iSymI
    nOrbJ = nOrb(iSymJ)
    nOccJ = nOsh(iSymJ)
    iSymIJ = MUL(iSymI,iSymJ)
    do iSymA=1,NSYM
      nOrbA = nOrb(iSymA)
      nOccA = nOsh(iSymA)
      iSymIJA = MUL(iSymIJ,iSymA)
      ISR = iSymA
      do iSymB=1,iSymA
        nOrbB = nOrb(iSymB)
        nOccB = nOsh(iSymB)
        iSymIJAB = iSymIJAB+1

        if (iSymIJA /= iSymB) GO TO 101
        if ((nOccI*nOccJ) == 0) GO TO 101
        if (nOrbI*nOrbJ*nOrbA*nOrbB == 0) GO TO 101

        ! FIND ADDRESSES FOR THIS SYMMETRY BLOCK

        IADC = IAD2M(3*iSymIJAB-2)
        IADX1 = IAD2M(3*iSymIJAB-1)
        IADX2 = IAD2M(3*iSymIJAB)
        write(6,1000) iSymA,iSymB,iSymI,iSymJ
1000    format(/1X,'SYMMETRY BLOCK < A B | I J >',4I4)

        if (IADC == 0) then
          !write(6,*)
          write(6,1100)
1100      format(1X,'NO COULOMB INTEGRALS FOR THIS SYMMETRY BLOCK')
        else
          write(6,1200) IADC
1200      format(1X,'ADDRESS FOR COULOMB INTEGRALS',I8)
          IAD13C = IADC
        end if
        if (IADX1 == 0) then
          !write(6,*)
          write(6,1110)
1110      format(1X,'NO EXCHAN1 INTEGRALS FOR THIS SYMMETRY BLOCK')
        else
          write(6,1210) IADX1
1210      format(1X,'ADDRESS FOR EXCHAN1 INTEGRALS',I8)
          IAD131 = IADX1
        end if
        if (IADX2 == 0) then
          !write(6,*)
          write(6,1120)
1120      format(1X,'NO EXCHAN2 INTEGRALS FOR THIS SYMMETRY BLOCK')
        else
          write(6,1220) IADX2
1220      format(1X,'ADDRESS FOR EXCHAN2 INTEGRALS',I8)
          IAD132 = IADX2
        end if
        LREC = nOrbA*nOrbB
        if (iSymA == iSymB) LREC = (nOrbA**2+nOrbA)/2
        LRECX = nOrbA*nOrbB
        if (.not. DoTCVA) LRECX = (nOrbA-nOccA)*(nOrbB-nOccB)
        LInt = 0
        LEx1 = 0
        LEx2 = 0
        do NI=1,nOccI
          NUM = nOccJ
          if (iSymI == iSymJ) NUM = NI
          do NJ=1,NUM

            ! THE LOOP ABOVE OVER T AND U RECOVERS ONE BLOCK OF INTEGRALS <AB|TU>
            ! FOR EACH PAIR T,U. TO PROCESS ONE BLOCK FOR ALL A AND B THE
            ! FOLLOWING LOOP STRUCTURE IS USED:

            if (IADC /= 0) then
              LInt = LInt+LREC
              if ((IPRX > 0) .and. (IPRX < 3)) then
                call mma_allocate(Tmp,LREC,Label='Tmp')
                call dDAFILE(LUINTM,2,Tmp,LREC,IAD13C)
                write(6,1300) NI,NJ,IAD13C-LREC,(Tmp(I),I=1,LREC)
                call mma_deallocate(Tmp)
1300            format(/1X,'<AB|IJ> COULOMB INTEGR.S FOR |ij> PAIR',2I3,'  DiskAdd=',I8/(8F10.6))
              end if
            end if

            ! THE EXCHANGE INTEGRALS OF TYPE 1 ,(AT|BU) ARE PROCESSED AS THE
            ! COULOMB INTEGRALS. IF NST /= NSU THERE ARE ALSO EXCHANGE
            ! INTEGRALS OF TYPE 2, (AU|BT). THE ORDERING IS STILL T,U AND A,B
            ! BUT T IS NOW THE FOURTH INDEX AND U THE SECOND
            ! EXCHANGE INTEGRALS ARE ALWAYS QUADRATIC IN A,B

            if (IADX1 /= 0) then
              LEx1 = LEx1+LRECX
              if (IPRX > 1) then
                call mma_allocate(Tmp,LRECX,Label='Tmp')
                call dDAFILE(LUINTM,2,Tmp,LRECX,IAD131)
                write(6,1310) NI,NJ,IAD131-LRECX,(Tmp(I),I=1,LRECX)
                call mma_deallocate(Tmp)
1310            format(/1X,'EXCHAN1 INTEGRALS FOR |ij> PAIR',2I3,'  DiskAdd=',I8/(8F10.6))
              end if
            end if

            if (IADX2 /= 0) then
              LEx2 = LEx2+LRECX
              if (IPRX > 1) then
                call mma_allocate(Tmp,LRECX,Label='Tmp')
                call dDAFILE(LUINTM,2,Tmp,LRECX,IAD132)
                write(6,1320) NI,NJ,IAD132-LRECX,(Tmp(I),I=1,LRECX)
                call mma_deallocate(Tmp)
1320            format(/1X,'EXCHAN2 INTEGRALS FOR |ij> PAIR',2I3,'  DiskAdd=',I8/(8F10.6))
              end if
            end if
          end do
        end do

        write(6,*)
        write(6,1010) LInt,LEx1,LEx2
1010    format(3X,'LCou=',I8,' , LEx1=',I8,' , LEx2=',I8)
        LTotInt = LTotInt+LInt
        LTotEx1 = LTotEx1+LEx1
        LTotEx2 = LTotEx2+LEx2

        ! ALL INTEGRALS FOR SYMMETRY BLOCK iSymI,iSymJ,iSymA,iSymB ARE READ

101     continue
      end do
    end do
  end do
end do
write(6,*)
write(6,2000) LTotInt,LTotEx1,LTotEx2
2000 format(3X,'LTotCou=',I8,' , LTotEx1=',I8,' , LTotEx2=',I8)
write(6,*) '   LTotTot=',LTotInt+LTotEx1+LTotEx2
write(6,*)

return

end subroutine RDINT2
