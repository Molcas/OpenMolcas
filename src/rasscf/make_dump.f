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
* Copyright (C) 2014, Giovanni Li Manni                                *
************************************************************************
      subroutine make_dump(TUVX,nacpr2,LW1,nacpar,EMY,
     & LCMO,DIAF,LuIntA,ITER)
************************************************************************
*                                                                      *
*     Objective: routine to generate FCIDUMP file needed by NECI       *
*                                                                      *
*     Author:    G. Li Manni, Max-Planck Institute December 2014       *
*                                                                      *
************************************************************************
#include "trafo_fciqmc.fh"
#include "fciqmc_global.fh"
#include "WrkSpc.fh"
#include "fciqmc.fh"
      Real*8 LCMO(*),DIAF(*),TUVX(*)
      Logical DoCholesky,int2check,dbg
      integer ITER

*----------------------------------------------------------------------*
*---- Start program and say Hello                                      *
*----------------------------------------------------------------------*
      Call qEnter('make_dump')
      Dbg=.false.
*----------------------------------------------------------------------*
*---- Run through the input section                                    *
*----------------------------------------------------------------------*
      Call init_dump
*----------------------------------------------------------------------*
*---- Check on existence of 1-electron integrals                       *
*----------------------------------------------------------------------*
      if(dbg) write(6,*) 'EMY in make_dump= ', EMY
      call f_Inquire('ONEINT',int2check)
      If (.not.int2check) then
       write(6,*)' One-electron integrals file not found!'
       write(6,*)' One-electron integrals are required by FCIQMC jobs'
       write(6,*)' See you later ;)'
       Call QTrace()
       Call Abend()
      End If
*----------------------------------------------------------------------*
*     Read the content of the one electron integral file               *
*     to collect information about:                                    *
*     Overlap, core Hamiltonian, kinetic integrals                     *
*----------------------------------------------------------------------*
c ----- I need to work this out.... GLM

        Call Get_iScalar('nSym',nSym)
        Call Get_iArray('nBas',nBas,nSym)
        Call Get_iArray('nOrb',nOrb,nSym)
        Call Get_iArray('nFro',nFro,nSym)
        Call Get_iArray('nISh',nISh,nSym)
        Call Get_iArray('nASh',nAsh,nSym)
        Call Get_iArray('nDel',nDel,nSym)
        Call Get_iScalar('nActel',nActEl)
        Call Get_iScalar('LSym',lSym)
        Call get_iscalar('iSpin',iSpin)
        if(dbg) then
          write(6,'(6X,A5,I2)')  'nSym:', nSym
          write(6,'(6X,A5,8I3)') 'nBas:', (nBas(i),i=1,nSym)
          write(6,'(6X,A5,8I3)') 'nOrb:', (nOrb(i),i=1,nSym)
          write(6,'(6X,A5,8I3)') 'nFro:', (nFro(i),i=1,nSym)
          write(6,'(6X,A5,8I3)') 'nIsh:', (nIsh(i),i=1,nSym)
          write(6,'(6X,A5,8I3)') 'nAsh:', (nAsh(i),i=1,nSym)
          write(6,'(6X,A5,8I3)') 'nDel:', (nDel(i),i=1,nSym)
          write(6,'(6X,A7, I3)') 'nActEl:', nActEl
          write(6,'(6X,A7, I3)') 'lSym:', lSym
        end if

        nTot=0
        nTot1=0
        nTot2=0
        n2max=0
        Do iSym=1,nSym
          iBas=nBas(iSym)
          nTot=nTot+iBas
          nTot1=nTot1+iBAs*(iBas+1)/2
          nTot2=nTot2+iBas*iBas
          n2max=Max(n2max,iBas*iBas)
        End Do

        Call Rd1Int_fciqmc(ipOvlp,ipHOne,ipKine)
        if(dbg) then
          write(6,*) 'Overlap integrals in make_dump.f'
          write(6,*) (Work(ipOvlp+i),i=0,nTot1-1)
          write(6,*)
          write(6,*) 'Core integrals in make_dump.f'
          write(6,*) (Work(ipHOne+i),i=0,nTot1-1)
          write(6,*)
          write(6,*) 'Kinetic terms in make_dump.f'
          write(6,*) (Work(ipKine+i),i=0,nTot1-1)
        end if
*----------------------------------------------------------------------*
*---- backup nOrb(isym) to nOrb2(isym)                                 *
*---- backup nFro(isym) to nFro2(isym)                                 *
*---- backup nDel(isym) to nDel2(isym)                                 *
*---- We want to transform only active space... therefore for the MO   *
*---- transformation suite actual Frozen + Inactive would be set to be *
*---- Frozen.... Accordingly, nOrb is set to nAsh only.                *
*---- Easiest way, maybe!?                                             *
*----                                                                  *
*---- Original values will be restored for later uses                  *
*----                                                                  *
*---- nOrbt and nOrbtt computed as needed in tr1ctl_rasscf             *
*----------------------------------------------------------------------*
        nOrbt=0
        nOrbtt=0
        Do iSym=1,nSym
          nOrb2(iSym)=nOrb(iSym)
          nFro2(iSym)=nFro(iSym)
          nDel2(iSym)=nDel(iSym)
          nFro(iSym)=nFro(iSym)+nIsh(iSym)
          nOrb(iSym)=nAsh(iSym)
          nDel(iSym)=nBas(iSym)-nFro(iSym)-nOrb(iSym)
          nOrbt=nOrbt+nOrb(iSym)
          nOrbtt=nOrbtt+nOrb(iSym)*(nOrb(iSym)+1)/2
        end do

        if(dbg) then
          write(6,*) 'nOrb being bucked up...'
          write(6,'(6X,8I3)') (nOrb2(i),i=1,nSym)
          write(6,*) 'nFro being bucked up...'
          write(6,'(6X,8I3)') (nFro2(i),i=1,nSym)
          write(6,*) 'nDel being bucked up...'
          write(6,'(6X,8I3)') (nDel2(i),i=1,nSym)
          write(6,*) 'new nFro...'
          write(6,'(6x,8i3)') (nFro(i),i=1,nSym)
          write(6,*) 'new nOrb...'
          write(6,'(6x,8i3)') (nOrb(i),i=1,nSym)
          write(6,*) 'new nDel...'
          write(6,'(6x,8i3)') (nDel(i),i=1,nSym)
        end if

*----------------------------------------------------------------------*
*     Print MO coefficients and occupations                            *
*----------------------------------------------------------------------*
        If (dbg) then
         Write(6,*)
         Write(6,*) ' CMO in make_dump'
         Write(6,*) ' ---------------------'
         Write(6,*)
         ioff=1
         Do iSym = 1,nSym
          iBas = nBas(iSym)
          if(iBas.ne.0) then
            write(6,*) 'Sym =', iSym
            do i= 1,iBas
              write(6,*)(LCMO(ioff+iBas*(i-1)+j),j=0,iBas-1)
            end do
            iOff = iOff + (iBas*iBas)
          end if
         End Do
        End If

*----------------------------------------------------------------------*
*---- Check on Cholesky approach...                                    *
*----------------------------------------------------------------------*
      Call DecideOnCholesky(DoCholesky)
*----------------------------------------------------------------------*
*---- Check on existence of 2-electron integrals in AO basis           *
*----------------------------------------------------------------------*
      call f_Inquire('ORDINT',int2check)
      If (.not.(int2check.or.DoCholesky)) then
       write(6,*)' Two-electron integrals file not found!'
       write(6,*)' IF Choleski is used than keep going... else stop!'
      End If

c I let the master node do this to avoid a bug with FockTwo later on.
c The bug has to do with reading Buf of integrals in parallel.
c As final effect of the bug the Fock matrix in AO and MO basis (FLT) are wrong
c (differences to the 10th digit).
c      write(6,*) 'myRank at make_dump:', myRank
c      IF(myRank.eq.0) then
*----------------------------------------------------------------------*
*----  Initialize LuFCI for NECI interface
*----------------------------------------------------------------------*
       LuFCI= isFreeUnit(38)
c      call OpnFl('FCIDMP',LuFCI,Exist)
       call molcas_open(LuFCI,'FCIDMP')
c
       write(LuFCI,'(1X,A11,I3,A7,I3,A5,I3,A)') ' &FCI NORB=',norbt,
     &',NELEC=',nActEl,',MS2=',int((ISPIN-1.0d0)),','

       write(LuFCI,'(A,500(I2,","))')'  ORBSYM=',((j,i=1,norb(j)),
     &  j=1,nsym)
       write(LuFCI,'(2X,A5,I1)') 'ISYM=',LSYM -1
       write(LuFCI,'(A)') ' &END'

       LuTwoAO = LuIntA
*----------------------------------------------------------------------*
*      Transform the two-electron integrals                            *
*----------------------------------------------------------------------*
c       Call Tr2Ctl_rasscf(LCMO)
       Call Fill2elInt(TUVX,nacpr2)
*----------------------------------------------------------------------*
*      Transform the one-electron integrals                            *
*----------------------------------------------------------------------*
c       Call Tr1Ctl_rasscf(Work(ipOvlp),Work(ipHOne),Work(ipKine),LCMO
c     &,DIAF,ITER)
       Call Fill1elInt(Work(LW1),nacpar,EMY,DIAF,ITER)
*----------------------------------------------------------------------*
*      Restore Old value of nOrb and nFro                              *
*----------------------------------------------------------------------*
c      END IF
c      IF(myRank.ne.0) then
c        call bcast_2RDM("FCIDMP")
c      END IF

        Do iSym=1,nSym
          nFro(iSym)=nFro2(iSym)
          nOrb(iSym)=nOrb2(iSym)
          nDel(iSym)=nDel2(iSym)
        end do

        if(dbg) then
          write(6,*) 'nOrb being restored...'
          write(6,'(6X,8I3)') (nOrb(i),i=1,nSym)
          write(6,*) 'nFro being restored up...'
          write(6,'(6X,8I3)') (nFro(i),i=1,nSym)
          write(6,*) 'nDel being restored up...'
          write(6,'(6X,8I3)') (nDel(i),i=1,nSym)
        end if
*----------------------------------------------------------------------*
*     Clean up and  termination                                        *
*----------------------------------------------------------------------*
      Call GetMem('Kine','Free','Real',ipKine,nTot1+4)
      Call GetMem('HOne','Free','Real',ipHOne,nTot1+4)
      Call GetMem('Ovlp','Free','Real',ipOvlp,nTot1+4)
*
      close(LuFCI)
      Call FastIO('STATUS')
      ireturn=0
      Call qExit('make_dump')
      return
      End
