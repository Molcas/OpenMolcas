* eafmod $ this file belongs to the Molcas repository $
************************************************************************
* Wrapper routines agains the MPI-IO package.                          *
*                                                                      *
* In none mpp installation the AIXIO facility will be used.            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* History:                                                             *
*     MPI-IO: V.P. Vysotskiy, University of Lund, Sweden, 2012-2013    *
*                                                                      *
************************************************************************



      Subroutine EAFOpen(Lu,FName)
      Implicit None
      Character*(*) FName
      Integer Lu
#ifdef _MOLCAS_MPP_
      Character*200 FN
      Integer  iRC,n
      Integer  molcas_eaf_open
      Logical  Is_Real_Par
      External molcas_eaf_open,Is_Real_Par
#endif
      Integer  isfreeunit
      External isfreeunit
#include "fio.fh"
*
      Call QEnter('EAFOpen')
#ifdef _MOLCAS_MPP_
      If(Is_Real_Par()) Then
      FN=FName
      Call Ext_Rank(FN)
      n=LEN(TRIM(FN));
      FN(n+1:n+1)=char(0)

      iRC=molcas_eaf_open(FN,Lu)
      If (iRC.ne.0) Then
         Write (6,*) 'EAFOpen: Abort!'
         Write (6,*) 'iRC=',iRC
         Write (6,*) 'Lu=',Lu
         Write (6,*) 'Fname=',FName
         Write (6,*) 'molcas_eaf_Err_Code =',iRC
         Call Abend
      End If
      Else
#endif
      Lu=isfreeunit(7)
      Call DaName_MF(Lu,FName)
#ifdef _MOLCAS_MPP_
      End If
#endif
      Call QExit('EAFOpen')
      Return
      End
************************************************************************
      Subroutine EAFClose(Lu)
      Implicit None
      Integer Lu
#ifdef _MOLCAS_MPP_
      Integer  iRC
      Integer  molcas_eaf_close
      Logical  Is_Real_Par
      External molcas_eaf_close,Is_Real_Par
*
      If(Is_Real_Par()) Then
      Call QEnter('EAFClose')
      iRC=molcas_eaf_close(Lu)
      If (iRC.ne.0) Then
         Write (6,*) 'EAFClose: Abort!'
         Write (6,*) 'molcas_eaf_Err_Code =',iRC
         Call Abend
      End If
      Else
#endif
      Call QEnter('EAFClose')
      Call DaClos(Lu)
#ifdef _MOLCAS_MPP_
      End If
#endif
      Call QExit('EAFClose')
      Return
      End
************************************************************************
      Subroutine EAFAWrite(Lu,Buf,nBuf,Disk,id)
      Implicit None
      Integer Lu, nBuf, Buf(nBuf), id
      Real*8 Disk
#ifdef _MOLCAS_MPP_
      Integer  iRC
      Integer  molcas_eaf_awrite
      Logical  Is_Real_Par
      External molcas_eaf_awrite,Is_Real_Par
#include "SysDef.fh"
#endif
      Call QEnter('EAFAWrite')
#ifdef _MOLCAS_MPP_
      If(Is_Real_Par()) Then
      iRC=molcas_eaf_awrite(Lu,Disk,Buf,nBuf*ItoB,id)
      If (iRC.ne.0) Then
         Write (6,*) 'EAFAWrite: Abort!'
         Write (6,*) 'molcas_eaf_Err_Code =',iRC
         Call Abend
      End If
      Else
#endif
      id=0
      Call EAFWrite(Lu,Buf,nBuf,Disk)
#ifdef _MOLCAS_MPP_
      End If
#endif

      Call QExit('EAFAWrite')
      Return
      End
************************************************************************
      Subroutine EAFARead(Lu,Buf,nBuf,Disk,id)
      Implicit None
      Integer Lu, nBuf, Buf(nBuf), id
      Real*8 Disk
#ifdef _MOLCAS_MPP_
      Integer  iRC
      Integer  molcas_eaf_aread
      Logical  Is_Real_Par
      External molcas_eaf_aread,Is_Real_Par
#include "SysDef.fh"
#endif
      Call QEnter('EAFARead')
#ifdef _MOLCAS_MPP_
      If(Is_Real_Par()) Then
      iRC=molcas_eaf_aread(Lu,Disk,Buf,nBuf*ItoB,id)
      If (iRC.ne.0) Then
         Write (6,*) 'EAFARead: Abort!'
         Write (6,*) 'molcas_eaf_Err_Code =',iRC
         Call Abend
      End If
      Else
#endif
      id=0
      Call EAFRead(Lu,Buf,nBuf,Disk)
#ifdef _MOLCAS_MPP_
      End If
#endif
      Call QExit('EAFARead')
      Return
      End
************************************************************************
      Subroutine EAFWrite(Lu,Buf,nBuf,Disk)
      Implicit None
      Integer Lu, nBuf, Buf(nBuf)
      Real*8 Disk
      Integer iDisk
#ifdef _MOLCAS_MPP_
      Integer  iRC
      Integer  molcas_eaf_write
      Logical  Is_Real_Par
      External molcas_eaf_write,Is_Real_Par
#include "SysDef.fh"
#endif
      Call QEnter('EAFWrite')
#ifdef _MOLCAS_MPP_
      If(Is_Real_Par()) Then
      iRC=molcas_eaf_write(Lu,Disk,Buf,nBuf*ItoB)
      If (iRC.ne.0) Then
         Write (6,*) 'EAFWrite: Abort!'
         Write (6,*) 'molcas_eaf_Err_Code =',iRC
         Call Abend
      End If
      Else
#endif
      iDisk=Int(Disk)
      Call iDaFile(Lu,1,Buf,nBuf,iDisk)
      Disk=DBLE(iDisk)
#ifdef _MOLCAS_MPP_
      End If
#endif
      Call QExit('EAFWrite')
      Return
      End
************************************************************************
      Subroutine EAFRead(Lu,Buf,nBuf,Disk)
      Implicit None
      Integer Lu, nBuf, Buf(nBuf)
      Real*8 Disk
      Integer iDisk
#ifdef _MOLCAS_MPP_
      Integer  iRC
      Integer  molcas_eaf_read
      Logical  Is_Real_Par
      External molcas_eaf_read,Is_Real_Par
#include "SysDef.fh"
#endif
c.... buf and offset length in bytes
      Call QEnter('EAFRead')
#ifdef _MOLCAS_MPP_
      If(Is_Real_Par()) Then
      iRC=molcas_eaf_read(Lu,Disk,Buf,nBuf*ItoB)
      If (iRC.ne.0) Then
         Write (6,*) 'EAFRead: Abort!'
         Write (6,*) 'molcas_eaf_Err_Code =',iRC
         Call Abend
      End If
      Else
#endif
      iDisk=Int(Disk)
      Call iDaFile(Lu,2,Buf,nBuf,iDisk)
      Disk=DBLE(iDisk)
#ifdef _MOLCAS_MPP_
      End If
#endif
      Call QExit('EAFRead')
      Return
      End
************************************************************************
      Subroutine EAFWait(LU,id)
      Implicit None
      Integer LU, id
#ifdef _MOLCAS_MPP_
      Integer iRC
      Integer  molcas_eaf_wait
      External molcas_eaf_wait
      Logical  Is_Real_Par
      External Is_Real_Par
#endif
*
*
      Call QEnter('EAFWait')
#ifdef _MOLCAS_MPP_
      if(Is_Real_Par()) Then
      iRC=molcas_eaf_wait(LU,id)
      If (iRC.ne.0) Then
         Write (6,*) 'EAFWait: Abort!'
         Write (6,*) 'molcas_eaf_Err_Code =',iRC
         Call Abend
      End If
      End If
#else
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(LU)
         Call Unused_integer(id)
      End If
#endif
*
      Call QExit('EAFWait')
      Return
      End

      Function Disk2Byte(Disk)
      Real*8 Disk2Byte, Disk
#ifdef _MOLCAS_MPP_
      Disk2Byte=Disk
#else
#include "blksize.fh"
      Disk2Byte=Disk*DBLE(min_block_length)
#endif
      Return
      End
