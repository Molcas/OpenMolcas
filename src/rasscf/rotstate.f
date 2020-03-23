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
* Copyright (C) 2020, Jie J. Bao                                       *
************************************************************************
      Subroutine RotState()
* ****************************************************************
* history:                                                       *
* Jie J. Bao, on Mar. 13, 2020, created this file.               *
* ****************************************************************
#include "rasdim.fh"
#include "rasscf.fh"
#include "splitcas.fh"
#include "general.fh"
#include "gas.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='CICTL   ')
#include "csfbas.fh"
#include "gugx.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "rctfld.fh"
#include "timers.fh"
#include "casvb.fh"
#include "wadr.fh"
#include "rasscf_lucia.fh"
#include "pamint.fh"
#include "input_ras.fh"
#include "stdalloc.fh"



      Integer LHrot,NHrot                ! storing info in H0_Rotate.txt
      Integer LRCIVec,LRCItmp,NRCIVec,LRCIScr ! storing CIVec
      Integer LRState,LRSttmp,NRState    ! storing info in Do_Rotate.txt
      Integer LHScr                      ! calculating rotated H
      Integer rcidisk
      INTEGER LURot,IsFreeUnit
      EXTERNAL IsFreeUnit
      Logical ReadHAM
      INTEGER JRoot,Kroot,IPRLEV

      write(LF,*)
      write(LF,*) ('=',i=1,61)
      write(LF,*)
      write(LF,'(6X,A)')'Do_Rotate.txt is found in scratch directory.'
      write(LF,'(6X,A)')'Following properties are for rotated states.'
      write(LF,*)

      NRState=lRoots**2
      NHRot=NRState
      NRCIVec=lRoots*NConf

      CALL GETMEM('RCIVEC','ALLO','REAL',LRCIVec,NRCIVec)
      CALL GETMEM('RCIScr','ALLO','REAL',LRCIScr,NRCIVec)
      CALL GETMEM('HScr','ALLO','REAL',LHScr,NHRot)
      CALL GETMEM('RState','ALLO','REAL',LRState,NRState)
      CALL GETMEM('HRot','ALLO','REAL',LHRot,NHRot)


      IPRLEV=IPRLOC(3)

*JB   read rotation matrix in Do_Rotate.txt
      LUROT=183
      LUROT=IsFreeUnit(LURot)
      CALL Molcas_Open(LURot,'ROT_VEC')
      LRSttmp=LRState
      Do jRoot = 1,lRoots
       read(LURot,*) (Work(LRSttmp+kRoot-1),kRoot=1,lRoots)
       LRSttmp=LRSttmp+lRoots
      End Do
      close(LURot)
      iF(IPRLEV.GE.DEBUG) Then
        write(LF,*)'rotation matrix'
        LRSttmp=LRState
        Do jRoot = 1,lRoots
         write(LF,*) (Work(LRSttmp+kRoot-1),kRoot=1,lRoots)
         LRSttmp=LRSttmp+lRoots
        End Do
      eND iF
      NHRot=lRoots**2
      ReadHAM=.false.
      CALL f_inquire('ROT_HAM',ReadHAM)
      iF (ReadHAM) then
      write(LF,'(6X,A)')'H0_Rotate.txt is found in scratch directory.'
      write(LF,'(6X,2A)')'Reading rotated Hamiltonian from ',
     & 'H0_Rotate.txt'
        LUROT=IsFreeUnit(LURot)
        CALL Molcas_Open(LURot,'ROT_HAM')
        Do Jroot=1,lroots
          read(LUROT,*) (Work(LHRot+Jroot-1+(Kroot-1)*lroots)
     &                 ,kroot=1,lroots)
        End Do
        Close(LUROT)
        if(IPRLEV.GE.DEBUG) Then
         write(LF,'(6X,A)') 'Rotated H matrix read from scratch'
         write(LF,'(6X,F8.4)') (Work(LHRot+jroot),jroot=0,NHRot-1)
        End if
      eLSE
        write(LF,'(6X,2A)')'H0_Rotate.txt is not found in scratch ',
     &  'directory.'
        write(LF,'(6X,2A)')'Generating Hamiltonian matrix for ',
     &  'rotated states'
        write(LF,'(6X,A)')'and storing the matrix in H0_Rotate.txt'
       CALL DCOPY_(NHRot,[0.0d0],0,WORK(LHRot),1)
       Do I=1,lRoots
         WORK(LHRot+(I-1)*(lRoots+1))=ENER(I,ITER)
       End Do
       Call DGEMM_('n','n',lRoots,lRoots,lRoots,1.0D0,Work(LRState),
     &      lRoots,Work(LHRot),lRoots,0.0D0,Work(LHScr),lRoots)
       Call DGEMM_('n','t',lRoots,lRoots,lRoots,1.0D0,Work(LHScr),
     &      lRoots,Work(LRState),lRoots,0.0D0,Work(LHRot),lRoots)
       LUROT=IsFreeUnit(LURot)
       CALL Molcas_Open(LURot,'ROT_HAM')
       Do Jroot=1,lroots
         write(LUROT,*) (Work(LHRot+Jroot-1+(Kroot-1)*lroots)
     &                ,kroot=1,lroots)
       End Do
       Close(LUROT)
      eND iF
      if(IPRLEV.GE.DEBUG) Then
       write(LF,'(6X,A)') 'Rotated Hamialtonian matrix '
       write(LF,*) (Work(LHRot+jroot),jroot=0,NHRot-1)
      End if

*JB   read CI vector from jobiph
      rcidisk=IADR15(4)
      LRCItmp=LRCIScr
      Do jRoot = 1,lRoots
        Call DDafile(JOBIPH,2,Work(LRCItmp),nConf,rcidisk)
        LRCItmp=LRCItmp+NConf
      End Do
      Call DGEMM_('n','n',NConf,lRoots,lRoots,1.0D0,Work(LRCIScr),
     &     nConf,Work(LRState),lRoots,0.0D0,Work(LRCIVec),nConf)
C      Call DGEMM_('n','t',lRoots,NConf,lRoots,1.0D0,Work(LRState),
C    &       lRoots,Work(LRCIVec),nConf,0.0D0,Work(LRCIScr),lRoots)

C     updating final energies as those for rotated states
      rcidisk=IADR15(4)
      LRCItmp=LRCIVec
      Do I=1,lRoots
        Call DDafile(JOBIPH,1,Work(LRCItmp),nConf,rcidisk)
        ENER(I,ITER)=WORK(LHRot+(I-1)*(lRoots+1))
        LRCItmp=LRCItmp+NConf
      End Do

      IF(IPRLEV.GE.DEBUG) Then
      write(LF,'(2A)')'Printing the coeff of the first CSF',
     &' for each state'
      Do I=1,lRoots
        write(LF,*)WORK(LRCIVec+(I-1)*NConf)
      End Do
      End If

      CALL GETMEM('HScr','FREE','REAL',LHScr,NHRot)
      CALL GETMEM('RCIScr','Free','REAL',LRCIScr,NRCIVec)
      CALL GETMEM('RState','Free','REAL',LRState,NRState)
      CALL GETMEM('RCIVEC','FREE','REAL',LRCIVec,NRCIVec)
      CALL GETMEM('HRot','FREE','REAL',LHRot,NHRot)

      write(LF,*)
      write(LF,*) ('=',i=1,61)

      Return
      End Subroutine
