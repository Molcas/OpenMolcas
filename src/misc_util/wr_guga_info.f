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
      Subroutine WR_GUGA(Lu,iOpt,iDisk,
     &                   NFREF,S,N,LN,NSYM,IR1,IR2,IFIRST,INTNUM,
     &                   LSYM,NREF,LN1,NRLN1,MUL,nMUL,NSH,NISH,MxSym,
     &                   JRC,nJRC,JJS,nJJS,NVAL,IOCR,nIOCR)
      Implicit Real*8 (a-h,o-z)
      Integer MUL(nMUL), NSH(MxSym), NISH(MxSym), JRC(nJRC),
     &        JJS(nJJS), NVAL(MxSym), IOCR(nIOCR)
*
      Call s_iDaFile_guga(Lu,iOpt,NFREF,    1,iDisk)
      Call s_dDaFile_guga(Lu,iOpt,S,        1,iDisk)
      Call s_iDaFile_guga(Lu,iOpt,N,        1,iDisk)
      Call s_iDaFile_guga(Lu,iOpt,LN,       1,iDisk)
      Call s_iDaFile_guga(Lu,iOpt,NSYM,     1,iDisk)
      Call s_iDaFile_guga(Lu,iOpt,IR1,      1,iDisk)
      Call s_iDaFile_guga(Lu,iOpt,IR2,      1,iDisk)
      Call s_iDaFile_guga(Lu,iOpt,IFIRST,   1,iDisk)
      Call s_iDaFile_guga(Lu,iOpt,INTNUM,   1,iDisk)
      Call s_iDaFile_guga(Lu,iOpt,LSYM,     1,iDisk)
      Call s_iDaFile_guga(Lu,iOpt,NREF,     1,iDisk)
      Call s_iDaFile_guga(Lu,iOpt,LN1,      1,iDisk)
      Call s_iDaFile_guga(Lu,iOpt,NRLN1,    1,iDisk)
      Call iDaFile(Lu,iOpt,MUL,   nMUL,iDisk)
      Call iDaFile(Lu,iOpt,NSH,  MxSym,iDisk)
      Call iDaFile(Lu,iOpt,NISH, MxSym,iDisk)
      Call iDaFile(Lu,iOpt,JRC,   nJRC,iDisk)
      Call iDaFile(Lu,iOpt,JJS,   nJJS,iDisk)
      Call iDaFile(Lu,iOpt,NVAL, MxSym,iDisk)
      Call iDaFile(Lu,iOpt,IOCR, nIOCR,iDisk)
*
      Return
*
*     This is to allow type punning without an explicit interface
      Contains
      Subroutine s_iDaFile_guga(Lu,iOpt,Buf,lBuf_,iDisk_)
      Use Iso_C_Binding
      Integer Lu, iOpt, lBuf_, iDisk_
      Integer, Target :: Buf
      Integer, Pointer :: pBuf(:)
      Call C_F_Pointer(C_Loc(Buf),pBuf,[1])
      Call iDaFile(Lu,iOpt,pBuf,lBuf_,iDisk_)
      Nullify(pBuf)
      End Subroutine s_iDaFile_guga
      Subroutine s_dDaFile_guga(Lu,iOpt,Buf,lBuf_,iDisk_)
      Use Iso_C_Binding
      Integer Lu, iOpt, lBuf_, iDisk_
      Real*8, Target :: Buf
      Real*8, Pointer :: pBuf(:)
      Call C_F_Pointer(C_Loc(Buf),pBuf,[1])
      Call dDaFile(Lu,iOpt,pBuf,lBuf_,iDisk_)
      Nullify(pBuf)
      End Subroutine s_dDaFile_guga
*
      End
