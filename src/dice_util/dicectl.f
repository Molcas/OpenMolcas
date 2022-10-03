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
* Copyright (C) 2017, Quan Phung                                       *
************************************************************************
* Main control file for DICE. Template from CheMPS2 interface.

      Subroutine DiceCtl( LW1, TUVX, IFINAL, IRST )

      Implicit Real*8 (A-H,O-Z)

      Dimension LW1(*), TUVX(*)

      Integer iChMolpro(8)
      Character*3 Label
      Integer LINSIZE, NUM_TEI, dtemp, nooctemp, labelpsi4
      Integer conversion(8)
      Integer activesize(8)
      Real*8  revdiff, pt2ener
      Logical fiedler, mps0
      Integer chemroot
      character(len=10) :: rootindex
      character(len=3) :: dice_nprocs
      character(len=150) :: imp1, imp2
      Integer :: iOper(0:7)


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
#ifdef _MOLCAS_MPP_
#include "mpif.h"
#endif

! Quan: FIXME: Do we need this?
* Load symmetry info from RunFile
      iOper = 0
      Call Get_iScalar('NSYM',nIrrep)
      Call Get_iArray('Symmetry operations',iOper,nIrrep)
      Call Get_iScalar('Rotational Symmetry Number',iSigma)

* Get character table to convert MOLPRO symmetry format
      Call MOLPRO_ChTab(nSym,Label,iChMolpro)

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
      Call systemf('ln -sf FCIDUMP_CHEMPS2 FCIDUMP',iErr)

      Call Getmem('OrbSym','Free','Inte',lOrbSym,NAC)

*************************
*  WRITEOUT INPUT FILE  *
*************************
#ifdef _MOLCAS_MPP_
      if ( KING() ) then
#endif

! Cleanup output.dat.total
      if (IRST.EQ.0) then
        call c_remove("output.dat.total")
      endif
      write(6,*) 'DICE> INTERATION : ', ITER
      LUTOTE = isFreeUnit(30)
      call molcas_open(LUTOTE,'input.dat')

      write(LUTOTE,'(A4,I4)') 'nocc', NACTEL
      do iref_dice=1,nref_dice
          write(LUTOTE,'(A)') trim(diceocc(iref_dice))
      enddo
      write(LUTOTE,'(A3)') 'end'
      write(LUTOTE,'(A6,I3)') 'nroots', lroots
      write(LUTOTE,*)
      write(LUTOTE,'(A8)') 'schedule'
      write(LUTOTE,'(A1,E12.5)') '0', dice_eps1*1.0d1
      write(LUTOTE,'(A1,E12.5)') '3', dice_eps1*1.0d1
      write(LUTOTE,'(A1,E12.5)') '6', dice_eps1
      write(LUTOTE,'(A3)') 'end'
      write(LUTOTE,'(A7,I6)') 'maxiter', dice_iter
      write(LUTOTE,'(A5)') 'DoRDM'
      write(LUTOTE,'(A8)') 'dE 1.e-8'
      write(LUTOTE,*)
      write(LUTOTE,'(A7,I6)') 'SampleN', dice_sampleN
      write(LUTOTE,'(A8,E12.5)') 'epsilon2', dice_eps2
      write(LUTOTE,'(A18)') 'targetError 8.0e-5'
      if (IRST>0 .or. dice_restart.eqv..true.) then
        write(LUTOTE,'(A11)') 'fullrestart'
      endif

      if (dice_stoc.EQV..False.) then
        write(LUTOTE,'(A13)') 'deterministic'
      else
        write(LUTOTE,'(A13,E12.5)') 'epsilon2Large', dice_eps2*10.0d0
      endif

      close(LUTOTE)
#ifdef _MOLCAS_MPP_
      endif
#endif


******************************
*  RUNNING DICE AND ANALYSE  *
******************************

#ifdef _MOLCAS_MPP_
      if ( KING() ) then
#endif
        call get_environment_variable("MOLCAS_DICE",
     &                                dice_nprocs, status=ierr)
        if (ierr.NE.0) then
          imp2 = "Dice >output.dat 2>dice.err"
        else
          imp2 = "mpirun -np "//trim(adjustl(dice_nprocs))//
     &                " Dice >output.dat 2>dice.err"
        endif

        call systemf(imp2,iErr)
        call systemf("cat output.dat >> output.dat.total",iErr)
        if (iErr.NE.0) then
          write(6,*) 'DICE> DICE crashed'
        endif
#ifdef _MOLCAS_MPP_
      end if

      if ( Is_Real_Par() ) then
          CALL MPI_Barrier(MPI_COMM_WORLD, IERROR4)
      end if

      if ( Is_Real_Par().AND.( KING().EQV..false. ) ) then
        do chemroot=1,lroots
          write(rootindex,"(I2)") chemroot-1
          imp1="ln -sf ../spatialRDM."//
     &           trim(adjustl(rootindex))//"."//
     &           trim(adjustl(rootindex))//".txt ."
          call systemf(imp1,iErr)
        enddo
          call systemf("ln -sf ../output.dat .", iErr)

      end if
#endif

      if (dice_stoc.EQV..False.) then
        imp1 = "grep PTEnergy output.dat | cut -c 10- > dice.energy"
        call systemf(imp1,iErr)
      else
        do chemroot=1,lroots
          write(rootindex,"(I2)") chemroot+2
          imp1 = 'grep -A '//trim(adjustl(rootindex))//
     &    ' "VARIATIONAL CALCULATION RESULT" output.dat
     &         | tail -n 1 | cut -c 6-30 >> dice.energy'
          call systemf(imp1,iErr)
        enddo
        imp2 = "grep +/- output.dat | cut -c 10- > PT2.energy"
        call systemf(imp2,iErr)
      endif

      LUTOTE = isFreeUnit(30)
      call molcas_open(LUTOTE,'dice.energy')

      do chemroot=1,lroots
        read(LUTOTE,*) ENER(chemroot,ITER)
      enddo
      close(LUTOTE)
      call c_remove("dice.energy")

      if (dice_stoc.EQV..True.) then
        call molcas_open(LUTOTE,'PT2.energy')

        do chemroot=1,lroots
          read(LUTOTE,*) PT2ENER
          write(6,*) 'DICE> PT2 Energy: ', PT2ENER
        enddo
        close(LUTOTE)
        call c_remove("PT2.energy")
      endif

      Return
      End
