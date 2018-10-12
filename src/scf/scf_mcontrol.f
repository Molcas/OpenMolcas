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

      SUBROUTINE Scf_Mcontrol(id_call)

      Implicit Real*8 (a-h,o-z)
#include "mxdm.fh"
#include "infscf.fh"
#include "para_info.fh"
*
      Integer ALGO,NSCREEN
      Logical REORD,DECO,timings
      Real*8  dmpk,dFKmat
*
      Common /CHOSCF / REORD,DECO,dmpk,dFKmat,ALGO,NSCREEN
      COMMON /CHOTIME / timings
*
      Integer id_call
      Character*512 List
      Character*32  Value
*******************************************************

      icount=0

      If (id_call .eq. 1) Then

* --- Label definitions

         write(List,100) 'SCF_started_OK:(-:-):',ALGO,timings,dmpK,
     &                   EThr,DThr,FThr,nIter(1),nScreen


* --- Initialize the control system

         Call MolcasControlInit(List)
         Return

      Else

* --- Read the molcas control file

*1
         Call MolcasControl('Cho_ALGO',Value)
         If (Value(1:4).ne.'    ') Then
            read(Value,*,err=101,end=102) ALGO
            write(6,*)'--- Warning: Cho_ALGOrithm changed by user to ',
     & 'the value ',ALGO
            icount = icount + 1
         EndIf
*2
         Call MolcasControl('Chotime',Value)
         If (Value(1:4).ne.'    ') Then
            read(Value,*,err=101,end=102) timings
            write(6,*)'--- Warning: Cholesky timings visualization ',
     &'changed by user to the value ',timings
            icount = icount + 1
         EndIf
*3
         Call MolcasControl('En_thr',Value)
         If (Value(1:4).ne.'    ') Then
            read(Value,*,err=101,end=102) Ethr
            write(6,*)'--- Warning: SCF Energy threshold changed by ',
     &'user to the value ',Ethr
            icount = icount + 1
         EndIf
*4
         Call MolcasControl('D_thr',Value)
         If (Value(1:4).ne.'    ') Then
            read(Value,*,err=101,end=102) Dthr
            write(6,*)'--- Warning: SCF Density threshold changed by ',
     &'user to the value ',Dthr
            icount = icount + 1
         EndIf
*5
         Call MolcasControl('F_thr',Value)
         If (Value(1:4).ne.'    ') Then
            read(Value,*,err=101,end=102) Fthr
            write(6,*)'--- Warning: SCF Fmat threshold changed by user',
     &' to the value ',Fthr
            icount = icount + 1
         EndIf
*6
         Call MolcasControl('MaxIter',Value)
         If (Value(1:4).ne.'    ') Then
            read(Value,*,err=101,end=102) nIter(1)
            write(6,*)'--- Warning: SCF Max # iterations changed by ',
     &'user to the value ',nIter(1)
            icount = icount + 1
         EndIf
*7
         Call MolcasControl('nScreen',Value)
         If (Value(1:4).ne.'    ') Then
            read(Value,*,err=101,end=102) nScreen
            write(6,*)'--- Warning: Cholesky LK option nSCREEN changed',
     &' by user to the value ',nScreen
            icount = icount + 1
         EndIf
*8
         Call MolcasControl('dmpK',Value)
         If (Value(1:4).ne.'    ') Then
            read(Value,*,err=101,end=102) dmpK
            write(6,*)'--- Warning: Cholesky LK option DMPK changed by',
     &' user to the value ',dmpK
            icount = icount + 1
         EndIf

      EndIf

      icount0=icount

* --- Get the true updated counter in parallel runs
*
      Call gaIgOP_SCAL(icount,'max')

      If (MyRank.eq.0) Then
         If (icount.gt.icount0) Then
           write(6,*)' Steering will NOT be activated this time because'
           write(6,*)' molcas.control file must be changed on node_0 !!'
           Call gaIgOP_SCAL(icount,'min')
         EndIf
      EndIf

      If (icount.gt.0) Then

* --- Trick to broadcast the values across nodes
*
         If (MyRank.ne.0) Then
            ALGO=0
            nIter(1)=0
            nScreen=0
            dmpK=0.0d0
            EThr=0.0d0
            DThr=0.0d0
            FThr=0.0d0
         EndIf
         Call gaIgOP_SCAL(ALGO,'+')
         Call gaIgOP_SCAL(nScreen,'+')
         Call gaIgOP(nIter(1),1,'+')
         Call gadgOP_SCAL(dmpK,'+')
         Call gadgOP_SCAL(EThr,'+')
         Call gadgOP_SCAL(DThr,'+')
         Call gadgOP_SCAL(FThr,'+')

* --- Update label values (note that "timings" is locally updated!)
*
         write(List,100) 'SCF_modified_by_user:',ALGO,timings,dmpK,
     &                   EThr,DThr,FThr,nIter(1),nScreen


* --- Initialize the control system with the new values
*
         Call MolcasControlInit(List)
         Return

      EndIf

      Return

100   FORMAT(A21,
     &',Cho_ALGO=',I2,
     &',Chotime=',L2,
     &',dmpK=',E11.4,
     &',En_thr=',E11.4,
     &',D_thr=',E11.4,
     &',F_thr=',E11.4,
     &',MaxIter=',I4,
     &',nScreen=',I4)

101   write(6,*) 'Scf_Mcontrol: error in data Input. ( icount= ',
     &           icount,' )'
102   write(6,*) 'Scf_Mcontrol: reached end of file. ( icount= ',
     &           icount,' )'

      End
