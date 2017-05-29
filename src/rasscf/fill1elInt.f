************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2015, Giovanni Li Manni                                *
************************************************************************
      SUBROUTINE Fill1elInt(Fock,nacpar,emy,DIAF,iter)
************************************************************************
*                                                                      *
* Purpose: generate 1-electron integral entries in FCIDUMP file (NECI).*
*          This should be working with Choleski as well.               *
*                                                                      *
*     Author:    G. Li Manni, Max-Planck Institute January 2015        *
*                                                                      *
************************************************************************
*
#include "fciqmc_global.fh"
#include "files_fciqmc.fh"
#include "fciqmc.fh"
#include "WrkSpc.fh"
*     Declaring explicitely the needed ones... I would like to get read of the IMPLICIT
      Real*8 Fock(nacpar), Fval,emy,emyn,DIAF(*)
      Real*8 CPT,CPE,TIOT,TIOE
      integer nacpar,i,k,iorb,jorb,nactel,ioff,icount,isym,ITER
      integer itotnbas,Dummy, iErr
      logical okay, Dbg
*
      Call qEnter('Fill1elInt')
      Dbg=.false.
*
*     Set time at start of transformation
      CALL SETTIM
      if(Dbg) then
        CALL TRIPRT('Fock matrix in Fill1ElInt',' ',Fock,NAC)
      end if
      Call Get_iScalar('nActel',nActEl)
      if(Dbg) write(6,*) 'nActEl = ', nActEl

      IF(NACTEL.NE.0) THEN
        EMYN= EMY/DBLE(NACTEL)
      ELSE
        EMYN = 0.0d0
      END IF

      call Add_Info('Fock-elements',Fock,nacpar,8)

      do k = 1, nacpar
* We are going to compute i and j from k=ij....
        iorb = ceiling(-0.5d0+sqrt(2.0d0*k))
        jorb = k - (iorb-1)*iorb/2
        Fval = Fock(k)
        if(iorb.eq.jorb) Fval = Fock(k)- EMYN
        if(abs(Fval).ge.1.0d-11)then
          write(LuFCI,'(1X,G20.11,4I5)')  Fval,iorb,jorb,0,0
          if(Dbg)then
             write(6,'(1X,G20.11,4I5)')   Fval,iorb,jorb,0,0
          end if
        endif
      enddo

c Orbital energies section ....
        IF(ITER.eq.1) then
c read orbital energies from INPORB file
         call f_Inquire (FnInpOrb,okay)
         If ( okay ) Then
          itotnbas = 0
          do i = 1, nsym
            itotnbas= itotnbas + nbas(i)
          end do
          call getmem('EORB','Allo','Real',ipEOrb,itotnbas)
          Call RdVec(FnInpOrb,LuInpOrb,'E',nSym,nBas,nBas,
     &          Dummy, Dummy, work(ipEOrb), iDummy,
     &          VecTit, 0, iErr)
         Else
          Write (6,*) 'RdCMO: Error finding MO file'
          Call QTrace()
          Call Abend()
         End If

c write orbital energies into FCIDUMP file
         ioff   = 0
         icount = 0
         do ISYM=1,NSYM
           IF ( NORB(ISYM).GT.0 ) THEN
            do i = 1,norb(isym)
             write(LuFCI,'(1X,G20.11,4I5)') WORK(ipEOrb+ioff+nfro(isym)+
     &                                    i-1),i+icount,0,0,0
             if(Dbg)then
               write(6,'(1X,G20.11,4I5)')   WORK(ipEOrb+ioff+nfro(isym)+
     &                                    i-1),i+icount,0,0,0
             endif
            enddo
           END IF
           ioff   = ioff + nbas(isym)
           icount = icount + norb(isym)
         end do
         call getmem('EORB','Free','Real',ipEOrb,itotnbas)
        ELSE
c Here orbital energies are updated from the latest Fock Matrix
c We assume that the Fock Matrix is diagonal and that diagonal elements represent
c orbital energies
         ioff   = 0
         icount = 0
         do isym=1,nsym
           IF ( NORB(ISYM).GT.0 ) THEN
            do i = 1,norb(isym)
             write(LuFCI,'(1X,G20.11,4I5)') DIAF(ioff+nfro(isym)+
     &                                   i),i+icount,0,0,0
             if(Dbg)then
                write(6,'(1X,G20.11,4I5)') DIAF(ioff+nfro(isym)+
     &                                     i),i+icount,0,0,0
             end if
            enddo
           END IF
           ioff   = ioff + nbas(isym)
           icount = icount + norb(isym)
         end do
        END IF
c write core energy into FCIDMP file

*     write the core energy to dump file ....
      write(LuFCI,'(1X,G20.11,4I5)') EMY,0,0,0,0

      call Add_Info('core energy',EMY,1,8)

      if(Dbg)then
        write(6,*) 'Core energy...'
        write(6,'(1X,G20.11,4I5)') EMY,0,0,0,0
      end if

      CALL TIMING(CPT,CPE,TIOT,TIOE)
      If (iPrint.GE.5) WRITE(6,2200) CPT,TIOT
2200  FORMAT(/6X,' TOTAL CPU TIME(SEC)',F8.2,'TOTAL I/O TIME(SEC)',F8.2)
c
*     exits
      Call qExit('Fill1elInt')
      RETURN
      END
