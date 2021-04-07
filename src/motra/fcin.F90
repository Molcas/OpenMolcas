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

subroutine FCIN(FLT,nFLT,DLT,FSQ,DSQ,EMY,CMO)

implicit real*8(A-H,O-Z)
#include "motra_global.fh"
#include "trafo_motra.fh"
#include "WrkSpc.fh"
real*8 CMO(*)
dimension DLT(*), FLT(nFLT), DSQ(*), FSQ(*)
logical DoCholesky

! Construct the one-electron density matrix for frozen space

call DONEI(DLT,DSQ,CMO)

! Compute the one-electron energy contribution to EMY

EONE = 0.0d0
do NPQ=1,NTOT1
  EONE = EONE+DLT(NPQ)*FLT(NPQ)
end do
EMY = EONE
if ((IPRINT >= 5) .or. (DEBUG /= 0)) then
  write(6,'(6X,A,E20.10)') 'ONE-ELECTRON CORE ENERGY:',EONE
end if

! Quit here if there are no frozen orbitals

NTFRO = 0
n_Bas = 0
do ISYM=1,NSYM
  NTFRO = NTFRO+NFRO(ISYM)
  n_Bas = max(n_Bas,nBas(iSym))
end do
if (NTFRO == 0) then
  return
end if

! Compute the two-electron contribution to the Fock matrix

call Allocate_Work(ipTemp,nFlt)
call FZero(Work(ipTemp),nFlt)

call DecideOnCholesky(DoCholesky)

if (DoCholesky) then
  call Cho_Fock_MoTra(nSym,nBas,nFro,DLT,DSQ,FLT,nFLT,FSQ,1.0d0)

  ! Print the Fock-matrix

  if ((IPRINT >= 5) .or. (DEBUG /= 0)) then
    write(6,'(6X,A)') 'Fock matrix in AO basis'
    ISTLTT = ipTemp
    do ISYM=1,NSYM
      NB = NBAS(ISYM)
      if (NB > 0) then
        write(6,'(6X,A,I2)') 'symmetry species:',ISYM
        call TRIPRT(' ',' ',WORK(ISTLTT),NB)
        ISTLTT = ISTLTT+NB*(NB+1)/2
      end if
    end do
  end if

  goto 99  ! jump over the conventional ERIs calculation

end if

call GETMEM('FCIN2','ALLO','REAL',LW2,n_Bas**2)
call GETMEM('FCIN1','MAX','REAL',LW1,LBUF)
LBUF = max(LBUF-LBUF/10,0)
call GETMEM('FCIN1','ALLO','REAL',LW1,LBUF)

call FTWOI(DLT,DSQ,Work(ipTemp),nFlt,FSQ,LBUF,WORK(LW1),WORK(LW2))

call GETMEM('FCIN1','FREE','REAL',LW1,LBUF)
call GETMEM('FCIN2','FREE','REAL',LW2,n_Bas**2)

99 continue
call DaXpY_(nFlt,1.0D+0,Work(ipTemp),1,Flt,1)
call Free_Work(ipTemp)

! Add the two-electron contribution to EMY
ETWO = -EONE
do NPQ=1,NTOT1
  ETWO = ETWO+DLT(NPQ)*FLT(NPQ)
end do
EMY = EONE+0.5d0*ETWO
if ((IPRINT >= 5) .or. (DEBUG /= 0)) then
  write(6,'(6X,A,E20.10)') 'TWO-ELECTRON CORE ENERGY:',ETWO
end if

! Exit

return

end subroutine FCIN
