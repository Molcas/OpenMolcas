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
      INTEGER JRoot,Kroot,IPRLEV
      CHARACTER(Len=18)::MatInfo
      INTEGER ReadStat
      write(LF,*)
      write(LF,*) ('=',i=1,71)
      write(LF,*)
      write(LF,'(11X,A)')'Do_Rotate.txt is found in scratch directory.'
      IF(IXMSP.eq.1) THEN
       write(LF,'(11X,A)')
     & 'Following properties are for XMS intermediate states.'
      ELSE IF(ICMSP.eq.1) THEN
       write(LF,'(11X,A)')
     & 'Following properties are for CMS intermediate states.'
      ELSE
       write(LF,'(11X,A)')
     & 'Following properties are for intermediate states'
       write(LF,'(11X,A)')
     & ' obtained from the user-supplied rotation matrix'
      ENDIF

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
      Read(LURot,*,iostat=ReadStat) MatInfo
      IF(ReadStat.eq.-1) MatInfo='an unknown method'
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
      CALL DCOPY_(NHRot,[0.0d0],0,WORK(LHRot),1)
      Do I=1,lRoots
        WORK(LHRot+(I-1)*(lRoots+1))=ENER(I,ITER)
      End Do
      Call DGEMM_('n','n',lRoots,lRoots,lRoots,1.0D0,Work(LRState),
     &     lRoots,Work(LHRot),lRoots,0.0D0,Work(LHScr),lRoots)
      Call DGEMM_('n','t',lRoots,lRoots,lRoots,1.0D0,Work(LHScr),
     &     lRoots,Work(LRState),lRoots,0.0D0,Work(LHRot),lRoots)
      LUROT=IsFreeUnit(LURot)
      CALL Molcas_Open(LURot,'ROT_HAM')
      Do Jroot=1,lroots
        write(LUROT,*) (Work(LHRot+Jroot-1+(Kroot-1)*lroots)
     &               ,kroot=1,lroots)
      End Do
      write(LURot,*) MatInfo
      Close(LUROT)
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
      Call DGEMM_('n','t',NConf,lRoots,lRoots,1.0D0,Work(LRCIScr),
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
      write(LF,*) ('=',i=1,71)

      Return
      End Subroutine
