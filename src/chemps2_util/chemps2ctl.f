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
* Copyright (C) 2016, Sebastian Wouters                                *
*               2016, Quan Phung                                       *
************************************************************************
! CheMPS2-Molcas main interface
! Based on Block interface, written by N. Nakatani
! Written by Quan Phung and Sebastian Wouters, Leuven, Aug 2016
! Adapted for Molcas 8.1 by Quan Phung, Leuven, Oct 2016

      Subroutine Chemps2Ctl( LW1, TUVX, IFINAL, IRST )

#ifdef _MOLCAS_MPP_
      Use MPI
#endif

      Implicit Real*8 (A-H,O-Z)

      Dimension LW1(*), TUVX(*)

      Integer iChMolpro(8)
      Character*3 Label
      Integer LINSIZE, NUM_TEI, dtemp, nooctemp, labelpsi4
      Integer conversion(8)
      Integer activesize(8)
      Real*8  chemps2_totale_4d, revdiff, chemps2_conv
      Logical fiedler, mps0
      Integer chemroot, chemps2_info
      character(len=10) :: rootindex
      character(len=100) :: imp1, imp2
      Integer :: iOper(0:7), ihfocc

#ifdef _MOLCAS_MPP_
      Integer*4 IERROR4
      External King, Is_Real_Par
      Logical King
      Logical Is_Real_Par
#endif



#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='CHEMPS2CTL')
      Call qEnter(ROUTINE)

! Quan: FIXME: Do we need this?
* Load symmetry info from RunFile
      iOper = 0
      Call Get_iScalar('NSYM',nIrrep)
      Call Get_iArray('Symmetry operations',iOper,nIrrep)
      Call Get_iScalar('Rotational Symmetry Number',iSigma)

* Get character table to convert MOLPRO symmetry format
      Call MOLPRO_ChTab_BIS(nSym,Label,iChMolpro)

* Convert orbital symmetry into MOLPRO format
      Call Getmem('OrbSym','Allo','Inte',lOrbSym,NAC)
      iOrb=1
      Do iSym=1,nSym
        Do jOrb=1,NASH(iSym)
          iWork(lOrbSym+iOrb-1)=iChMolpro(iSym)
          iOrb=iOrb+1
        End Do
      End Do
      lSymMolpro=iChMolpro(lSym)

      NRDM_ORDER=2
      If (NACTEL.EQ.1) NRDM_ORDER=1


**********************
*  WRITEOUT FCIDUMP  *
**********************

      LINSIZE = ( NAC * ( NAC + 1 ) ) / 2
      NUM_TEI = ( LINSIZE * ( LINSIZE + 1 ) ) / 2
      Call FCIDUMP_OUTPUT( NAC, NACTEL, ISPIN-1,
     &                     lSymMolpro, iWork(lOrbSym),
     &                     0.0d0, LW1, TUVX,
     &                     LINSIZE, NUM_TEI )


      Call Getmem('OrbSym','Free','Inte',lOrbSym,NAC)

**************************
*  WRITEOUT ACTIVE FOCK  *
**************************

!      Write(6,*) "Currently the Fock matrix is printed in fckpt2.f"

*************************
*  WRITEOUT INPUT FILE  *
*************************

#ifdef _MOLCAS_MPP_
      if ( KING() ) then
#endif
      IF (IRST.EQ.0) THEN
! Cleanup chemps2.log.total
        call c_remove("chemps2.log.total")
!        call system('rm -f chemps2.log.total')
! Check if checkpoint files exist
        if (chemps2_restart.EQV..TRUE.) then
           call f_inquire('CHEMNATFIE',fiedler)
           call f_inquire('CHEMNATMPS0',mps0)
!           INQUIRE(FILE='molcas_natorb_fiedler.txt', EXIST=fiedler)
!           INQUIRE(FILE='CheMPS2_natorb_MPS0.h5', EXIST=mps0)
           if (fiedler .and. mps0) then
             write(6,*) 'CHEMPS2> Found checkpoint files for DMRG-SCF'
! Copy CheMPS2_natorb_MPSxxx.h5 to CheMPS2_MPSxxx.h5
             call fcopy('CHEMNATFIE','CHEMFIE',iErr)
!             write(6,*) 'CHEMPS2> DB: 118', iErr
             do chemroot=1,lroots
              write (rootindex, "(I2)") chemroot-1
              imp1="CheMPS2_natorb_MPS"//trim(adjustl(rootindex))//".h5"
              imp2="CheMPS2_MPS"//trim(adjustl(rootindex))//".h5"
              call fcopy(imp1,imp2,iErr)
!              write(6,*) 'CHEMPS2> DB: 124', iErr
             enddo
!              call system("cp molcas_natorb_fiedler.txt
!     &                     molcas_fiedler.txt")
!              do chemroot=1,lroots
!                write (rootindex, "(I2)") chemroot-1
!                andrea=""
!             andrea="cp CheMPS2_natorb_MPS"//trim(adjustl(rootindex))//
!     &              ".h5 CheMPS2_MPS"//trim(adjustl(rootindex))//".h5"
!                call system(andrea)
!              enddo
           else
! Reset chemps2_restart = .false. if not checkpoint files
              write(6,*) 'CHEMPS2> No checkpoint files for DMRG-SCF'
              chemps2_restart=.false.
           endif
        endif
! Check if checkpoint files for 3RDM exist
        if (chemps2_lrestart.EQ.1.) then
           call f_inquire('CHEMCANFIE',fiedler)
           call f_inquire('CHEMCANMPS0',mps0)
!           INQUIRE(FILE='molcas_canorb_fiedler.txt', EXIST=fiedler)
!           INQUIRE(FILE='CheMPS2_canorb_MPS0.h5', EXIST=mps0)
           if (fiedler .and. mps0) then
             write(6,*) 'CHEMPS2> Found checkpoint files for n-RDM'
           else
             write(6,*) 'CHEMPS2> No checkpoint files for n-RDM'
             chemps2_lrestart=0
           endif
        endif

      ENDIF
#ifdef _MOLCAS_MPP_
      endif
#endif

      LUCHEMIN=isFreeUnit(29)
      call molcas_open(LUCHEMIN,'chemps2.input')
      write(LUCHEMIN,*) 'FCIDUMP = FCIDUMP_CHEMPS2'

      call group_psi4number(Label,Labelpsi4)
      write(LUCHEMIN,'(1X,A8,I1)') 'GROUP = ', Labelpsi4
      write(LUCHEMIN,*)

      write(LUCHEMIN,'(1X,A13,I2)') 'EXCITATION = ', lRoots-1
      write(LUCHEMIN,*)

      IF ((ABS(CBLBM)>chemps2_blb .AND. IFINAL.NE.2) .OR.
     &   (IRST.EQ.0 .AND. (chemps2_restart.EQV..FALSE.)) .OR.
     &   (IFINAL.EQ.2 .AND. (Do3RDM.EQV..TRUE.)
     &                .AND. (chemps2_lrestart.EQ.0)) .OR.
     &   (IFINAL.EQ.2 .AND. iOrbTyp.EQ.2
     &                .AND. (chemps2_lrestart.EQ.0))) THEN

            call c_remove("molcas_fiedler.txt")
!            call system('rm -f molcas_fiedler.txt')
!Quan: FIXME: how to remove CheMPS2_MPS0.h5, etc with c_remove
            call systemf("rm -f CheMPS2_MPS*.h5", iErr)
!            write(6,*) 'CHEMPS2> DB: 187', iErr
            write(LUCHEMIN,*) 'MOLCAS_MPS     = TRUE'
            write(6,*) 'CHEMPS2> Start DMRG from scratch'

        write(LUCHEMIN,'(1X,A21)',ADVANCE='NO') 'SWEEP_STATES       = '
        dtemp = 500
        do
          write(LUCHEMIN,('(I7,A2)'),ADVANCE='NO') dtemp, ','
          dtemp = dtemp + min(dtemp,1000)
          if (dtemp .GE. MxDMRG) then
            write(LUCHEMIN,'(I7,A2,I7)') MxDMRG, ',', MxDMRG
            exit
          endif
        enddo

        write(LUCHEMIN,'(1X,A21)',ADVANCE='NO') 'SWEEP_ENERGY_CONV  = '
        dtemp = 500
        do
          if (dtemp .EQ. 500) then
            write(LUCHEMIN,'(E12.5,A3)',ADVANCE='NO') THRE*1000.0, ','
         else
           write(LUCHEMIN,'(E12.5,A3)',ADVANCE='NO') THRE*100.0, ','
         endif

         dtemp = dtemp + min(dtemp,1000)
         if (dtemp .GE. MxDMRG) then
           write(LUCHEMIN,'(E12.5,A3,E12.5)') THRE*5.0, ',', THRE/2.0
           exit
         endif
        enddo

        write(LUCHEMIN,'(1X,A21)',ADVANCE='NO') 'SWEEP_MAX_SWEEPS   = '
        dtemp = 500
        do
          if (dtemp .EQ. 500) then
            write(LUCHEMIN,'(I7,A3)',ADVANCE='NO') max_sweep, ','
          else
            write(LUCHEMIN,'(I7,A3)',ADVANCE='NO') max_sweep/2, ','
          endif

          dtemp = dtemp + min(dtemp,1000)
          if (dtemp .GE. MxDMRG) then
           if (IFINAL.EQ.2) then
             if (Do3RDM .OR. (iOrbTyp.EQ.2)) THEN
              write(LUCHEMIN,'(I7,A3,I7)') max_sweep/2, ','
     &                                   , max_canonical
             else
              write(LUCHEMIN,'(I7,A3,I7)') max_sweep/2, ','
     &                                   , max_sweep*5
             endif
           else
            write(LUCHEMIN,'(I7,A3,I7)') max_sweep/2, ',', max_sweep
           endif
           exit
          endif
        enddo

        write(LUCHEMIN,'(1X,A21)',ADVANCE='NO') 'SWEEP_NOISE_PREFAC = '
        dtemp = 500
        do
         write(LUCHEMIN,'(E12.5,A3)',ADVANCE='NO') chemps2_noise, ','
         dtemp = dtemp + min(dtemp,1000)
         if (dtemp .GE. MxDMRG) then
           write(LUCHEMIN,'(E12.5,A10)') chemps2_noise, ' ,   0.00'
           exit
         endif
        enddo

        write(LUCHEMIN,'(1X,A21)',ADVANCE='NO') 'SWEEP_DVDSON_RTOL  = '
        dtemp = 500
        do
         if (dtemp .EQ. 500) then
           write(LUCHEMIN,'(A9)',ADVANCE='NO') '1.0e-3 ,'
         else
           write(LUCHEMIN,'(A9)',ADVANCE='NO') '1.0e-4 ,'
         endif

         dtemp = dtemp + min(dtemp,1000)
         if (dtemp .GE. MxDMRG) then
           write(LUCHEMIN,'(A9,E12.5)') '1.0e-4 ,', davidson_tol
           exit
         endif
        enddo
        write(LUCHEMIN,*)

      ELSE
! DMRG restart with fixed orbital order
        IF ((ABS(CBLBM)>chemps2_blb/10.0 .AND. IFINAL.NE.2) .OR.
     &   (IRST.EQ.0 .AND. (chemps2_restart.EQV..TRUE.))) THEN
          write(6,*) 'CHEMPS2> Partial restart DMRG ',
     &                  'from previous step'

          write(LUCHEMIN,*) 'MOLCAS_MPS     = TRUE'
          write(LUCHEMIN,'(1X,A21)',ADVANCE='NO')
     &                                      'SWEEP_STATES       = '
          write(LUCHEMIN,'(I7,A2,I7)') MxDMRG, ',', MxDMRG

          write(LUCHEMIN,'(1X,A21)',ADVANCE='NO')
     &                                      'SWEEP_ENERGY_CONV  = '
          write(LUCHEMIN,'(E12.5,A3,E12.5)') THRE*5.0, ',', THRE/2.0

          write(LUCHEMIN,'(1X,A21)',ADVANCE='NO')
     &                                      'SWEEP_MAX_SWEEPS   = '
          if (IFINAL.EQ.2) then
             if (Do3RDM .OR. (iOrbTyp.EQ.2)) THEN
               write(LUCHEMIN,'(I7,A3,I7)') max_sweep/2, ','
     &                                    , max_canonical
             else
               write(LUCHEMIN,'(I7,A3,I7)') max_sweep/2, ','
     &                                    , max_sweep*5
             endif
          else
            write(LUCHEMIN,'(I7,A3,I7)') max_sweep/2, ',', max_sweep
          endif

          write(LUCHEMIN,'(1X,A21)',ADVANCE='NO')
     &                                      'SWEEP_NOISE_PREFAC = '
          write(LUCHEMIN,'(E12.5,A10)') chemps2_noise, ' ,   0.00'

          write(LUCHEMIN,'(1X,A21)',ADVANCE='NO')
     &                                      'SWEEP_DVDSON_RTOL  = '
          write(LUCHEMIN,'(A9,E12.5)') '1.0e-4 ,', davidson_tol
          write(LUCHEMIN,*)
        ELSE
!          write(6,*) 'Full Restart'
          write(6,*) 'CHEMPS2> Fully restart DMRG from previous step'
          write(LUCHEMIN,*) 'MOLCAS_MPS     = TRUE'
          write(LUCHEMIN,'(1X,A21)',ADVANCE='NO')
     &                                      'SWEEP_STATES       = '
          write(LUCHEMIN,'(I7)') MxDMRG

          write(LUCHEMIN,'(1X,A21)',ADVANCE='NO')
     &                                      'SWEEP_ENERGY_CONV  = '
          write(LUCHEMIN,'(E12.5)') THRE/2.0

          write(LUCHEMIN,'(1X,A21)',ADVANCE='NO')
     &                                      'SWEEP_MAX_SWEEPS   = '
          if (IFINAL.EQ.2 .OR. (IFINAL.EQ.1 .AND. iCIonly.EQ.1)) then
            if (IFINAL.EQ.2 .AND. (Do3RDM .OR. (iOrbTyp.EQ.2))) THEN
              write(LUCHEMIN,'(I7)') max_canonical
            else
              write(LUCHEMIN,'(I7)') max_sweep*5
            endif
          else
            write(LUCHEMIN,'(I7)') max_sweep
          endif

          write(LUCHEMIN,'(1X,A21)',ADVANCE='NO')
     &                                      'SWEEP_NOISE_PREFAC = '
          write(LUCHEMIN,'(A16)') '   0.00'

          write(LUCHEMIN,'(1X,A21)',ADVANCE='NO')
     &                                      'SWEEP_DVDSON_RTOL  = '
          write(LUCHEMIN,'(E12.5)') davidson_tol
          write(LUCHEMIN,*)
        ENDIF
      ENDIF


      write(LUCHEMIN,'(A8)',ADVANCE='NO') 'NOCC = '
      do nooctemp=1,NSYM-1
        write(LUCHEMIN,'(A5)',ADVANCE='NO') '0 ,'
      enddo
      write(LUCHEMIN,'(A3)') '0'

      write(LUCHEMIN,'(A8)',ADVANCE='NO') 'NACT = '
      call molpro2psi(Label,conversion)
      do iSym=1,nSym
        activesize( conversion( iChMolpro( iSym ) ) ) = NASH( iSym )
      end do
      do nooctemp=1,NSYM-1
        write(LUCHEMIN,'(I3,A2)',ADVANCE='NO')
     &                        activesize(nooctemp), ' ,'
      enddo
      write(LUCHEMIN,'(I3)')  activesize(NSYM)

      write(LUCHEMIN,'(A8)',ADVANCE='NO') 'NVIR = '
      do nooctemp=1,NSYM-1
        write(LUCHEMIN,'(A5)',ADVANCE='NO') '0 ,'
      enddo
      write(LUCHEMIN,'(A3)') '0'
      write(LUCHEMIN,*)

      write(LUCHEMIN,*) 'MOLCAS_FIEDLER = TRUE'
      write(LUCHEMIN,*) 'MOLCAS_STATE_AVG = TRUE'
      write(LUCHEMIN,*) 'MOLCAS_2RDM    = molcas_2rdm.h5'

      if (sum(hfocc) .NE. 0) then
        write(LUCHEMIN,'(A13)',ADVANCE='NO') 'MOLCAS_OCC ='
        do ihfocc=1,NAC-1
          write(LUCHEMIN,'(I3,A2)', ADVANCE='NO') HFOCC(ihfocc), ', '
        enddo
        write(LUCHEMIN,'(I3)') HFOCC(NAC)
      endif

      If (IFINAL.EQ.2 .AND. Do3RDM .AND. NACTEL.GT.2) Then
         write(6,*)  'CHEMPS2> Running 3-RDM and F.4-RDM'
         write(LUCHEMIN,*) 'MOLCAS_3RDM    = molcas_3rdm.h5'
         write(LUCHEMIN,*) 'MOLCAS_F4RDM   = molcas_f4rdm.h5'
         write(LUCHEMIN,*) 'MOLCAS_FOCK    = FOCK_CHEMPS2'
      endif

      write(LUCHEMIN,*)

      write(LUCHEMIN,*) 'PRINT_CORR = TRUE'
      write(LUCHEMIN,*) 'TMP_FOLDER = ./'

      close(LUCHEMIN)

#ifdef _MOLCAS_MPP_
      write(6,'(1X,A21,I3)') 'CHEMPS2> ITERATION : ', ITER
      if ( KING() ) then
#endif

! Quan: overwrite CheMPS2_xxxorb_MPSX.h5 to CheMPS2_MPSX.h5
         if (
     &   (IFINAL.EQ.2 .AND. Do3RDM
     &                .AND. (chemps2_lrestart>0)) .OR.
     &   (IFINAL.EQ.2 .AND. iOrbTyp.EQ.2
     &                .AND. (chemps2_lrestart>0))) THEN
           if (chemps2_lrestart.EQ.1) then
             write(6,*) 'CHEMPS2> Using user-supplied checkpoint files'
             call fcopy('CHEMCANFIE','CHEMFIE',iErr)
!             write(6,*) 'CHEMPS2> DB: 390', iErr
!             call system("cp molcas_canorb_fiedler.txt
!     &                   molcas_fiedler.txt")
             do chemroot=1,lroots
              write (rootindex, "(I2)") chemroot-1
              imp1="CheMPS2_canorb_MPS"//trim(adjustl(rootindex))//".h5"
              imp2="CheMPS2_MPS"//trim(adjustl(rootindex))//".h5"
              call fcopy(imp1,imp2,iErr)
!              write(6,*) 'CHEMPS2> DB: 398', iErr
!               andrea=""
!              andrea="cp CheMPS2_canorb_MPS"//trim(adjustl(rootindex))//
!     &               ".h5 CheMPS2_MPS"//trim(adjustl(rootindex))//".h5"
!               call system(andrea)
             enddo
           endif

           if (chemps2_lrestart.EQ.2) then
             write(6,*) 'CHEMPS2> Using checkpoint files from',
     &                   ' previous step (not recommended)'
             call fcopy('CHEMNATFIE','CHEMFIE',iErr)
!             write(6,*) 'CHEMPS2> DB: 410', iErr
!             call system("cp molcas_natorb_fiedler.txt
!     &                   molcas_fiedler.txt")
             do chemroot=1,lroots
              write (rootindex, "(I2)") chemroot-1
              imp1="CheMPS2_natorb_MPS"//trim(adjustl(rootindex))//".h5"
              imp2="CheMPS2_MPS"//trim(adjustl(rootindex))//".h5"
              call fcopy(imp1,imp2,iErr)
!              write(6,*) 'CHEMPS2> DB: 418', iErr
!               andrea=""
!              andrea="cp CheMPS2_natorb_MPS"//trim(adjustl(rootindex))//
!     &               ".h5 CheMPS2_MPS"//trim(adjustl(rootindex))//".h5"
!               call system(andrea)
             enddo
           endif
         endif



! Quan: save CANORB before actually calculating
         if (
     &   (IFINAL.EQ.2 .AND. Do3RDM
     &                ) .OR.
     &   (IFINAL.EQ.2 .AND. iOrbTyp.EQ.2
     &                )) THEN
            write(6,*) 'CHEMPS2> Save CANORB'
! Quan: FIXME: Bug!
!            Call OrbFiles(JOBIPH,IPRLEV)
         endif

         call systemf("chemps2 --file=chemps2.input > chemps2.log",iErr)
!         write(6,*) 'CHEMPS2> DB: 441', iErr
! Quan: FIXME: How to save chemps2.log file
         call systemf("cat chemps2.log >> chemps2.log.total",iErr)
!         write(6,*) 'CHEMPS2> DB: 444', iErr

! Quan: save natorb checkpoint file in all iteration
         if (IFINAL<2) then
           if (IFINAL.EQ.1) then
             write(6,*) 'CHEMPS2> Save natorb checkpoint files'
           endif
           call fcopy('CHEMFIE','CHEMNATFIE',iErr)
!           write(6,*) 'CHEMPS2> DB: 452', iErr
             do chemroot=1,lroots
              write (rootindex, "(I2)") chemroot-1
              imp1="CheMPS2_natorb_MPS"//trim(adjustl(rootindex))//".h5"
              imp2="CheMPS2_MPS"//trim(adjustl(rootindex))//".h5"
              call fcopy(imp2,imp1,iErr)
!              write(6,*) 'CHEMPS2> DB: 458', iErr
             enddo


!           call system("cp molcas_fiedler.txt
!     &                     molcas_natorb_fiedler.txt")
!           do chemroot=1,lroots
!             write (rootindex, "(I2)") chemroot-1
!             andrea=""
!             andrea="cp CheMPS2_MPS"//trim(adjustl(rootindex))//
!     &       ".h5 CheMPS2_natorb_MPS"//trim(adjustl(rootindex))//".h5"
!             call system(andrea)
!           enddo
         endif

! Quan: save canorb checkpoint file if possible
         if (
     &   (IFINAL.EQ.2 .AND. Do3RDM
     &                ) .OR.
     &   (IFINAL.EQ.2 .AND. iOrbTyp.EQ.2
     &                )) THEN

           write(6,*) 'CHEMPS2> Save canorb checkpoint files'
           call fcopy('CHEMFIE','CHEMCANFIE',iErr)
!           write(6,*) 'CHEMPS2> DB: 482', iErr
             do chemroot=1,lroots
              write (rootindex, "(I2)") chemroot-1
              imp1="CheMPS2_canorb_MPS"//trim(adjustl(rootindex))//".h5"
              imp2="CheMPS2_MPS"//trim(adjustl(rootindex))//".h5"
              call fcopy(imp2,imp1,iErr)
!              write(6,*) 'CHEMPS2> DB: 488', iErr
             enddo
!           call system("cp molcas_fiedler.txt
!     &                     molcas_canorb_fiedler.txt")
!           do chemroot=1,lroots
!             write (rootindex, "(I2)") chemroot-1
!             andrea=""
!             andrea="cp CheMPS2_MPS"//trim(adjustl(rootindex))//
!     &       ".h5 CheMPS2_canorb_MPS"//trim(adjustl(rootindex))//".h5"
!             call system(andrea)
!           enddo

         endif

! Quan: Cleanup checkpoint files
         if (IFINAL.EQ.2) then
            call c_remove("molcas_fiedler.txt")
!            call system('rm -f molcas_fiedler.txt')
!Quan: FIXME: how to remove CheMPS2_MPS0.h5, etc with c_remove
            call systemf("rm -f CheMPS2_MPS*.h5",iErr)
!            write(6,*) 'CHEMPS2> DB: 508', iErr
         endif

#ifdef _MOLCAS_MPP_
      end if


      if ( Is_Real_Par() ) then
          CALL MPI_Barrier(MPI_COMM_WORLD, IERROR4)
      end if

!Quan: FIXME: softlink all the n-RDM files
      if ( Is_Real_Par().AND.( KING().EQV..false. ) ) then
!         call system("ln -sf ../molcas_2rdm.h5 .")
!         call system("ln -sf ../molcas_3rdm.h5 .")
!         call system("ln -sf ../molcas_f4rdm.h5 .")
!         call system("ln -sf ../chemps2.log .")
        do chemroot=1,lroots
          write(rootindex,"(I2)") chemroot-1
          imp1="ln -sf ../molcas_2rdm.h5.r"//
     &           trim(adjustl(rootindex))//" ."
          call systemf(imp1,iErr)
!          write(6,*) 'CHEMPS2> DB: 529', iErr
          imp1="ln -sf ../molcas_3rdm.h5.r"//
     &           trim(adjustl(rootindex))//" ."
          call systemf(imp1,iErr)
!          write(6,*) 'CHEMPS2> DB: 533', iErr
          imp1="ln -sf ../molcas_f4rdm.h5.r"//
     &           trim(adjustl(rootindex))//" ."
          call systemf(imp1,iErr)
!          write(6,*) 'CHEMPS2> DB: 537', iErr
        enddo
        call systemf("ln -sf ../chemps2.log .",iErr)
!        write(6,*) 'CHEMPS2> DB: 541', iErr
      end if
#endif

!Quan: a very dirty way to extract the total energy
       call systemf(
     & 'grep "***  2-RDM" -B 5 chemps2.log | grep "all instructions" '//
     & '| cut -c 61- > chemps2_totale_4d', iErr)
!       write(6,*) 'CHEMPS2> DB: 549', iErr

!Quan: fix bug E(FCI) != E(CASSCF)
      call systemf(
     & 'grep "Econst" chemps2.log | cut -c 39- > chemps2_totale', iErr)
!      write(6,*) 'CHEMPS2> DB: 554', iErr

      LUTOTE = isFreeUnit(30)
      call molcas_open(LUTOTE,'chemps2_totale')

!Quan: write energy to ENER
      do chemroot=1,lroots
        read(LUTOTE,*) ENER(chemroot,ITER)
      enddo
      close(LUTOTE)

!Quan: check the difference between ener and ener_4d
      LUTOTE = isFreeUnit(30)
      call molcas_open(LUTOTE,'chemps2_totale_4d')
      do chemroot=1,lroots
        read(LUTOTE,*) chemps2_totale_4d
        revdiff = abs(chemps2_totale_4d -
     &                  ENER(chemroot,ITER))/chemps2_totale_4d
!        write(6,*) 'CHEMPS2> DB: revdiff', revdiff
        if (revdiff > 1.0D-9) then
          write(6,*) 'CHEMPS2> large (E(4m) - E(m))/E(4m) = ',
     &      revdiff, 'for root', chemroot, ', consider increasing m!'
        endif
      enddo
      close(LUTOTE)


!Quan: check chemps2 convergence
      write (rootindex, "(I2)") lroots+8
      imp1=""
      imp1="grep ""***  2-RDM"" -B "//trim(adjustl(rootindex))//
     & " chemps2.log | grep ""Energy difference"""//
     & " | cut -c 69- > chemps2_conv"
      call systemf(imp1,iErr)
!      write(6,*) 'CHEMPS2> DB: 587', iErr

      LUCONV = isFreeUnit(30)
      call molcas_open(LUCONV,'chemps2_conv')
      do chemroot=1,lroots
         read(LUCONV,*) chemps2_conv
         write(6,'(1X,A14,I3,A30,E10.2)') 'CHEMPS2> Root ', chemroot,
     &              ' :: DMRG energy convergence : ', chemps2_conv
         if (abs(chemps2_conv) > THRE/2.0) then
            write(6,*) 'CHEMPS2> DMRG not converged, ',
     &                 'consider increasing MXSWeep'
         endif
      enddo
      close(LUCONV)

!Quan: check if CheMPS2 finished without error
      imp1="grep ""Info on DMRG"" chemps2.log | "//
     &       "cut -c 43- > chemps2_info"
      call systemf(imp1,iErr)
!      write(6,*) 'CHEMPS2> DB: 687', iErr

      LUCONV = isFreeUnit(30)
      call molcas_open(LUCONV,'chemps2_info')
      read(LUCONV,*) chemps2_info
      close(LUCONV)
      if ( chemps2_info .NE. 0 ) then
       write(6,*) 'CHEMPS2> CheMPS2 ends abnormally, check calculation'
      endif


      Call qExit(ROUTINE)

      Return
      End
