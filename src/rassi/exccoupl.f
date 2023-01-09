************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine exccoupl()

      use constants, only: Zero, Half
      use definitions, only: wp, iwp, u6
      use frenkel_global_vars, only: iTyp, jTyp, valst, nestla,
     &                               nestlb, doexch, excl, eNucB
      use stdalloc, only: mma_allocate, mma_deallocate
      implicit real(kind=wp) (A-H,O-Z)
#include "rasdim.fh"
#include "cntrl.fh"
#include "Files.fh"
#include "Morsel.fh"
#include "Struct.fh"
#include "SysDef.fh"
#include "rassi.fh"
#include "prgm.fh"
#include "rasdef.fh"
#include "jobin.fh"
#include "symmul.fh"
      integer(kind=iwp) :: nstat1, nstat2,
     &                     run, lWKX, dimn, a, b
      integer(kind=iwp), external :: isFreeUnit
      character(len=13) :: filnam1, filnam2, filnam3, filnam4,
     &            filnam5, filnam6,  filnam9
      logical :: WK_C_exists, WK_X_exists, states1, states2
      character(len=1) :: labi, labj
      real(kind=wp), allocatable:: rBvA(:), Charge(:),
     &                         WKX1(:),WKX2(:), WK1(:),
     &                  Frenkeltri(:), Frenkelunknwn(:)
      real(kind=wp) :: exCoupl, rlWKX(1)
#ifdef _DEBUGPRINT_RASSI_
      logical :: debug_rassi_code = .true.
#else
      logical :: debug_rassi_code = .false.
#endif

      call StatusLine('RASSI:','Starting Excitonic Coupling Section')
      LuT_ = 10
      labi = 'A'
      labj = 'B'
      if (iTyp == 0) then
        write(u6,*) ' MONA or MONB keyword missing! '
        call Abend()
      else
        jTyp=1+mod(iTyp,2) ! 1 or 2
        if (jTyp .eq. 1) then
           labi='B'
           labj='A'
        end if
      end if
      write(u6,*) ' '
      write(u6,*) ' ************************************************'
      write(u6,*) '  Excitonic couplings section:  '
      write(u6,*) ' '
      write(u6,*) '   ( v:  eA*eB + eA*nucB + eB*nucA + nucA*nucB ) '
      write(u6,*) ' ************************************************'
      write(u6,*) ' '


      call mma_allocate(rBvA,mxroot*(mxroot+1)/2,Label='Wnn_ab')
      rBvA(:) = Zero
!      read in <rhoB|VnucA> from already renamed and saved RunFile called rhoBnucA
      call NameRun('AUXRFIL1')
      call Get_dArray('<rhoB|VnucA>',rBvA,mxroot*(mxroot+1)/2)
      call Get_iScalar('Unique atoms', nAtoms)
      call mma_allocate(Charge,nAtoms,Label='Zcharge')
      call Get_dArray('Effective nuclear Charge',Charge,nAtoms)
      VNN_AB = Zero
      call NameRun('AUXRFIL2')
      call PotNuc_nad(nSym,nAtoms,Charge,VNN_AB)
      if (debug_rassi_code) then
        write(u6,*) 'VNN_AB (nuc-nuc interaction)', VNN_AB
      end if
      call mma_deallocate(Charge)
      call NameRun('#Pop')    ! switch back to old RUNFILE

! read in number of states on each monomer
      LuT_ = 10
      write(filnam5,'(A,I1)') 'states_', iTyp
      inquire(file=filnam5, exist=states1)
      LuT = isFreeUnit(LuT_)
      call molcas_open(LuT, filnam5)
      read(LuT,*) nstat1
      close(LuT)

      write(filnam6,'(A,I1)') 'states_',jTyp
      inquire(file=filnam6, exist=states2)
      LuT = isFreeUnit(LuT_)
      call molcas_open(LuT,filnam6)
      read(LuT,*) nstat2
      close(LuT)

      if (EXCL) then
        dimn = nstat1 * nstat2
        write(u6,*) 'dimn EXCL', nstat1*nstat2
      else
        dimn = 1 + (nstat1-1) + (nstat2-1)
      end if
      ! allocate array with more entries then needed first
      call mma_allocate(Frenkelunknwn,dimn*(dimn+1)/2)
      Frenkelunknwn(:) = zero


      if (DoExch) then
        LuTX4 = isFreeUnit(LuT_)
        write(filnam9,'(A)') 'nRedX'
        call DANAME(LuTX4,filnam9)
        iAddr = 0
        call dDaFile(LuTX4,2,rlWKX,1,iAddr)
        lWKX = int(rlWKX(1), kind=iwp)
        call DACLOS(LuTX4)
      end if

      nijkl = 0

      b = 0
      a = 0
      run = 0
      !loop over I
      do istate=1,nstat1
        !loop over K
        do kstate=1,nstat2
          if ((istate /= 1) .or. (kstate /= 1)) then
            if (EXCL) then
              if ((ALL(nestla /= istate)) .or.
     &           (ALL(nestlb /= kstate))) then
                cycle
              end if
            else
              if ((istate <= valst) .and. (kstate <= valst)) cycle
              if ((istate > valst) .and. (kstate > valst)) cycle
            end if
          end if
          b = b + 1
          a = 0
         !loop over J
         do jstate=1,nstat1
           ! * Start Coulomb term
           if (istate > jstate ) then
           write(filnam1,'(A,I1,A,I3.3,A,I3.3)') 'WK_C',iTyp,'_',
     &            istate,'_',jstate
           else
             write(filnam1,'(A,I1,A,I3.3,A,I3.3)') 'WK_C',iTyp,'_',
     &            jstate,'_',istate
           end if
           inquire(file=filnam1, exist=WK_C_exists)
           LuT1 = isFreeUnit(LuT_)
           call molcas_open(LuT1, filnam1)
           read(LuT1,*) lWK1
           call mma_allocate(WK1,2*lWK1,Label='WK12')
           read(LuT1,*) (WK1(i),i=1,lWK1)
           ! End Coulomb term
           ! Start-Exchange-term
           if (DoExch) then
             if (istate > jstate) then
               write(filnam3,'(A,I1,I2.2,I2.2)') 'X',iTyp,
     &             istate,jstate
             else
               write(filnam3,'(A,I1,I2.2,I2.2)') 'X',iTyp,
     &             jstate,istate
             end if
             inquire(file=filnam3, exist=WK_X_exists)
             LuTX1 = isFreeUnit(LuT_)
             call DANAME(LuTX1,filnam3)
             iAddr = 0
             call mma_allocate(WKX1,lWKX)
             iAddr = 0
             call dDaFile(LuTX1,2,WKX1,lWKX,iAddr)
             call DACLOS(LuTX1)
           end if
          ! End-Exchange-term

           !loop over L
           do lstate=1,nstat2
              if ((jstate /= 1) .or. (lstate /= 1)) then
                if (EXCL) then
                  if ((ALL(nestla /= jstate)) .or.
     &               (ALL(nestlb /= lstate))) then
                    cycle
                  end if
                else
                  if ((jstate <= valst) .and. (lstate <= valst)) cycle
                  if ((jstate > valst) .and. (lstate > valst)) cycle
                end if
              end if
              a = a + 1
              if (a <= b) then
                if (kstate > lstate) then
                  write(filnam2,'(A,I1,A,I3.3,A,I3.3)') 'WK_C',jTyp,
     &                     '_',KSTATE,'_',LSTATE
                else
                  write(filnam2,'(A,I1,A,I3.3,A,I3.3)') 'WK_C',jTyp,
     &                     '_',LSTATE,'_',KSTATE
                end if
               inquire(file=filnam2, exist=WK_C_exists)
               if (WK_C_exists) then
                 AB_nuc = Zero
                 eBnucA = Zero
                 if (istate == jstate) then
                   eBnucA = rBvA(kstate*(kstate-1)/2+lstate)
                 end if
                 eAnucB = Zero
                 if (kstate == lstate) then
                   eAnucB = eNucB(ISTATE*(ISTATE-1)/2+JSTATE)
                   if (istate == jstate) AB_nuc = VNN_AB
                 end if

                LuT2 = isFreeUnit(LuT_)
                call molcas_open(LuT2, filnam2)
                read(LuT2,*) lWK2
                if (lWK2 == LWK1) then
                  read(LuT2,*) (WK1(lWK1+i),i=1,lWK2)
                  eeCoupl = ddot_(lWK1,WK1(1),1,WK1(lWK1+1),1)
                  if (debug_rassi_code) then
                    write(u6,*) 'J = ', eeCoupl
                  endif
      !           Start-Exchange-term
                  if (DoExch) then
                    if (kstate > lstate) then
                      write (filnam4,'(A,I1,I2.2,I2.2)') 'X',jTyp,
     &                       kstate,lstate
                    else
                      write(filnam4,'(A,I1,I2.2,I2.2)')'X',jTyp,
     &                       lstate,kstate
                    end if
                    inquire(file=filnam4, exist=WK_X_exists)

                    LuTX2 = isFreeUnit(LuT_)
                    call DANAME(LuTX2,filnam4)
                    call mma_allocate(WKX2,lWKX)
                    iAddr = 0
                    call dDaFile(LuTX2,2,WKX2,lWKX,iAddr)
                    call DACLOS(LuTX2)

                    exCoupl = -Half*ddot_(lWKX,WKX1(1),1,WKX2(1),1)

                    call mma_deallocate(WKX2)
                    if (debug_rassi_code) then
                      write(u6,*) 'K = ', exCoupl
                    end if
                    eeCoupl = eeCoupl + exCoupl
                  end if
                ! End-Exchange-term
                  write(u6,'(3X,A,I3.3,A,A,A,I3.3,A,A,A,
     &                                  I3.3,A,A,A,I3.3,A,A,A,
     &                                  E18.8)')
     &                         '<(',istate,')',labi,'(',kstate,')',labj,
     &                       '|v|(',jstate,')',labi,
     &                          '(',lstate,')',labj,'> = ', eeCoupl
     &                                                     +eAnucB
     &                                                     +eBnucA
     &                                                     +AB_nuc

                  run = run + 1
                  if (debug_rassi_code) then
                    write(u6,*) 'eeCoupl', eeCoupl, 'eAnucB', eAnucB,
     &                          'eBnucA', eBnucA, 'AB_nuc', AB_nuc
                  end if
                  ! write into the Frenkel basis matrix (triangular)
                  Frenkelunknwn(run) = eeCoupl+eAnucB+eBnucA+AB_nuc
                  nijkl = nijkl+1
                else
                  write(u6,*) ' Size of WK array mismatch! '
                  call abend()
                ! endif lWK2 eq LWK1
                end if
              ! endif WK_C exists
              end if
            ! end IF (a .le. b)
            end if
            close(LuT2)
          ! end L
          end do
          close(LuT1)
          if (DoExch) then
            call mma_deallocate(WKX1)
          end if
          call mma_deallocate(WK1)
         ! end J
         end do
        ! end K
        end do
      ! end I
      end do

      write(u6,*) ' '
      write(u6,*) '  Total Nr of couplings: ', nijkl
      write(u6,*) ' ************************************************'
      write(u6,*) ' '
      !determine dimensions of Hamiltonian
      discrim = 1+4*2*nijkl
      dimn = int((-1+sqrt(discrim))/2, kind=iwp)
      write(u6,'(A,I3.3,A,I3.3)') 'determined Hamiltonian dimensions:',
     &                            dimn,"x",dimn

      call mma_allocate(Frenkeltri,dimn*(dimn+1)/2,label='frenkeltri')
      Frenkeltri(:) = zero
      n = dimn*(dimn+1)/2
      do i=1,n
        Frenkeltri(i) = Frenkelunknwn(i)
      end do

      call frenkelexc(Frenkeltri,dimn,nstat1,nstat2)

      call mma_deallocate(Frenkeltri)
      call mma_deallocate(Frenkelunknwn)
      call mma_deallocate(rBvA)
      if (excl) then
        call mma_deallocate(NESTLA)
        call mma_deallocate(NESTLB)
      end if

      end subroutine exccoupl
