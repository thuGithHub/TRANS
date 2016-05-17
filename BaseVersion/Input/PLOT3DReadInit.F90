subroutine PLOT3D_Read_And_Init
      use global
      implicit none
      character::well
      character*50 openfilename, str2, bloName,dimensionFile,GridFileName
      integer :: icg,iblock,gblock, kBlock,icon,iboc,NX,NY,NZ,i,j,k,NDim(3), base,NBc,NCon,kc
      type(ConnectivityStruct)   :: AConnectivity
      type(ConnectivityStruct),pointer :: AConnect
      integer :: Prange(6)
      integer::itemp,idtemp,ibtemp
      integer,allocatable::IJKBC(:,:)
      integer::ID_Present,NBoc,i_surf,i_idxa0,i_idxa1,i_idxb0,i_idxb1,Ibgn,Iend,Jbgn,Jend,Kbgn,Kend
      integer::BCtype,LeftOrRight
      integer::TargetNo,TargetConNo
      integer::Imin,Jmin,Kmin,Imax,Jmax,Kmax 
      integer::Transform(3),IJKDirect,revMatrix(3),lreva,lrevb,TransformT(3)

!***********************read global  boundarys & Connectivities*************************
        do iblock=1,NAllBlock
            write (dimensionFile, '(A, I3.3, A)') "grd/dimension_blk", iblock, "_m.dat"
            open(unit=13, file=trim(dimensionFile), MODE='read')
            read(13,*)
            read(13,*)
            read(13,*)well,NX,NY,NZ
            read(13,*)
            read(13,*)
            read(13,*)well,ID_Present,NBoc
            AllBlockName(ID_Present)%NCon=NBoc
            AllBlockName(ID_Present)%NDim=(/NX,NY,NZ/)           
            allocate(AllBlockName(ID_Present)%Prange(6,NBoc)) 
            allocate(AllBlockName(ID_Present)%CPrange(6,NBoc)) 
            !allocate(IJKBC(NBoc,BC_Item))
            allocate(IJKBC(10,NBoc))
            allocate(AllBlockName(ID_Present)%IJKBC(10,NBoc))
            read(13,*)
            read(13,*)
            do iboc=1,NBoc
                read(13,*)well,IJKBC(1,iboc),i_surf,i_idxa0,i_idxa1,i_idxb0,i_idxb1,IJKBC(2,iboc),&
                &   IJKBC(3,iboc),IJKBC(4,iboc),IJKBC(5,iboc),IJKBC(6,iboc),IJKBC(7,iboc),IJKBC(8,iboc),&
                &   IJKBC(9,iboc)
                if(i_surf==1.or.i_surf==2)then
                    if ( i_idxa0 == 0 ) i_idxa0 = 1
                    if ( i_idxa1 == 0 ) i_idxa1 = NY-1
                    if ( i_idxb0 == 0 ) i_idxb0 = 1
                    if ( i_idxb1 == 0 ) i_idxb1 = NZ-1
                    if(i_surf==1)then
                        Ibgn=1
                        Iend=1
                    else
                        Ibgn=NX!-1
                        Iend=NX!-1
                    endif
                    Jbgn=i_idxa0
                    Jend=i_idxa1!+1
                    Kbgn=i_idxb0
                    Kend=i_idxb1!+1
                elseif(i_surf==3.or.i_surf==4)then
                    if ( i_idxa0 == 0 ) i_idxa0 = 1
                    if ( i_idxa1 == 0 ) i_idxa1 = NX-1
                    if ( i_idxb0 == 0 ) i_idxb0 = 1
                    if ( i_idxb1 == 0 ) i_idxb1 = NZ-1
                    if(i_surf==3)then
                        Jbgn=1
                        Jend=1
                    else
                        Jbgn=NY!-1
                        Jend=NY!-1
                    endif
                    Ibgn=i_idxa0
                    Iend=i_idxa1!+1
                    Kbgn=i_idxb0
                    Kend=i_idxb1!+1
                elseif(i_surf==5.or.i_surf==6)then
                    if ( i_idxa0 == 0 ) i_idxa0 = 1
                    if ( i_idxa1 == 0 ) i_idxa1 = NX-1
                    if ( i_idxb0 == 0 ) i_idxb0 = 1
                    if ( i_idxb1 == 0 ) i_idxb1 = NY-1
                    if(i_surf==5)then
                        Kbgn=1
                        Kend=1
                    else
                        Kbgn=NZ!-1
                        Kend=NZ!-1
                    endif
                    Ibgn=i_idxa0
                    Iend=i_idxa1!+1
                    Jbgn=i_idxb0
                    Jend=i_idxb1!+1
                else 
                    write(*,*)"Wrong Dimension Input!"
                endif
                AllBlockName(ID_Present)%Prange(1:6,iboc)=(/Ibgn,Jbgn,Kbgn,Iend,Jend,Kend/)
                if(Ibgn==Iend)then
                    Jend=Jend-1
                    Kend=Kend-1
                    if(Ibgn>1)then
                        Iend=Iend-1
                        Ibgn=Ibgn-1
                    endif
                elseif(Jbgn==Jend)then
                    Iend=Iend-1
                    Kend=Kend-1
                    if(Jbgn>1)then
                        Jbgn=Jbgn-1
                        Jend=Jend-1
                    endif
                elseif(Kbgn==Kend)then
                    Iend=Iend-1
                    Jend=Jend-1
                    if(Kbgn>1)then
                        Kend=Kend-1
                        Kbgn=Kbgn-1
                    endif
                else
                    write(*,*)"Wrong CPrange!",MyProc
                endif
            AllBlockName(ID_Present)%CPrange(1:6,iboc)=(/Ibgn,Jbgn,Kbgn,Iend,Jend,Kend/)
            enddo
            AllBlockName(ID_Present)%IJKBC=IJKBC
            close(13)
            deallocate(IJKBC)
        enddo
!-------------------------***************
        do iblock=1,NBlock
            ibtemp=M_BlockIDs(iblock)
            AllBlocks(iblock)%NBlockGlobal=ibtemp  !block number in total blocks
            NCon=AllBlockName(ibtemp)%NCon
            AllBlocks(iblock)%NConnect= NCon
            allocate(AllBlocks(iblock)%AConnectivity(NCon))
            NDim=AllBlockName(ibtemp)%NDim
            AllBlocks(iblock)%NDim=NDim
        enddo
!******************************************************************************
!******************************set local block Cons and BCs********************
        do iblock=1,NBlock
            NCon=AllBlocks(iblock)%NConnect
            ibtemp=M_BlockIDs(iblock)
            do icon=1,NCon
            !    write(*,*)"Con Start",MyProc,iblock,icon,Ncon
                    Aconnect=>AllBlocks(iblock)%AConnectivity(icon)
                    BCtype=AllBlockName(ibtemp)%IJKBC(2,icon)
                    AConnect%BCtype=BCType
                    AConnect%SubBCType=AllBlockName(ibtemp)%IJKBC(3,icon)
                    AConnect%pointStart=AllBlockName(ibtemp)%Prange(1:3,icon)    
                    AConnect%pointEnd=AllBlockName(ibtemp)%Prange(4:6,icon)    
                    AConnect%CpointStart=AllBlockName(ibtemp)%CPrange(1:3,icon)    
                    AConnect%CpointEnd=AllBlockName(ibtemp)%CPrange(4:6,icon)    
                    TargetNo=AllBlockName(ibtemp)%IJKBC(4,icon)
                    AConnect%TargetBlock=TargetNo
                    TargetConNo=AllBlockName(ibtemp)%IJKBC(5,icon)
                    AConnect%TargetConn=TargetConNo
                    AConnect%lcross=AllBlockName(ibtemp)%IJKBC(6,icon)
                    AConnect%lreva=AllBlockName(ibtemp)%IJKBC(7,icon)
                    AConnect%lrevb=AllBlockName(ibtemp)%IJKBC(8,icon)
                    AConnect%PerType=AllBlockName(ibtemp)%IJKBC(9,icon)
                    
                    AConnect%rotateAngle=(/0,0,0/)
                    if(AConnect%PerType.ne.0) then    !for Periodic Conn
                        AConnect%rotateAngle(1)=real(Npassage*AConnect%PerType)*2*Pi/real(Nblade)     !this angle just for rotor37
                    endif
            enddo
        enddo
!*****************************set target connectivity******************************
        do iblock=1,NBlock
            ibtemp=M_BlockIds(iblock)
!            Ncon=AllBlockName(ibtemp)%NCon
            Ncon=AllBlocks(iblock)%NConnect
            do icon=1,Ncon
                AConnect=>AllBlocks(iblock)%AConnectivity(icon)
                BCtype=AConnect%BCtype
                if(BCtype==9)then
!                    ConNo=AllBlocks(iblock)%AConnectivity(icon)%ConNo
                    TargetNo=AConnect%TargetBlock
                    TargetConNo=AConnect%TargetConn
                    AConnect%TargetPStart(1:3)=AllBlockName(TargetNo)%Prange(1:3,TargetConNo)
                    AConnect%TargetPEnd(1:3)=AllBlockName(TargetNo)%Prange(4:6,TargetConNo)
                    AConnect%TPLijk(1:3)= abs(Aconnect%TargetPStart(1:3)-AConnect%TargetPEnd(1:3))+1
                endif
            enddo
        enddo
        
!***************************
        do iblock=1,Nblock
            Ncon=AllBlocks(iblock)%NConnect
            NDim=AllBlocks(iblock)%NDim
            do icon=1,NCon
                AConnect=>AllBlocks(iblock)%AConnectivity(icon)
                Imin= AConnect%PointStart(1)
                Jmin= AConnect%PointStart(2)
                Kmin= AConnect%PointStart(3)
                Imax= AConnect%PointEnd(1)
                Jmax= AConnect%PointEnd(2)
                Kmax= AConnect%PointEnd(3)
                !! IJKDirect
                if( Imin==Imax )then
                  AConnect%IJKDirect  = 1
                elseif( Jmin==Jmax )then
                  AConnect%IJKDirect  = 2
                elseif( Kmin==Kmax )then
                  AConnect%IJKDirect  = 3
                else
                  write(*,*) 'error, it must be the former three types!'
                  stop
                endif
                !! LeftOrRight
                if((Imin==Imax.and.Imin==1).or.(Jmin==Jmax.and.Jmin==1).or.(Kmin==Kmax.and.Kmin==1))then
                  AConnect%LeftOrRight= 1
                else
                  AConnect%LeftOrRight= 2
                endif

                AConnect%IJKAdd(:)= (/0,0,0/)
                AConnect%IJKAdd1(:)= (/0,0,0/)
                AConnect%IJKAdd2(:)= (/0,0,0/)
                AConnect%CIJKAdd(:)= (/0,0,0/)
                AConnect%CIJKAdd1(:)= (/0,0,0/)
                AConnect%CIJKAdd2(:)= (/0,0,0/)
                AConnect%CEIJKAdd(:)= (/0,0,0/)
                AConnect%CEIJKAdd1(:)= (/0,0,0/)
                AConnect%CEIJKAdd2(:)= (/0,0,0/)
                if( Imin==Imax )then
!                  if( Imin==ib)then
                  if( Imin==1)then
                    AConnect%IJKAdd(1)= 1  !!!??????
                    AConnect%IJKAdd1(1)= 2
                    AConnect%IJKAdd2(1)= 3
                    AConnect%CIJKAdd(1)=0
                    AConnect%CIJKAdd1(1)=1
                    AConnect%CIJKAdd2(1)=2
                    AConnect%CEIJKAdd(1)=-1
                    AConnect%CEIJKAdd1(1)=-2
                    AConnect%CEIJKAdd2(1)=-3
!                  elseif( Imin==im)then
                  elseif( Imin==NDim(1))then
                    AConnect%IJKAdd(1)= -1
                    AConnect%IJKAdd1(1)= -2
                    AConnect%IJKAdd2(1)= -3
                    AConnect%CIJKAdd(1)=-1
                    AConnect%CIJKAdd1(1)=-2
                    AConnect%CIJKAdd2(1)=-3
                    AConnect%CEIJKAdd(1)=0
                    AConnect%CEIJKAdd1(1)=1
                    AConnect%CEIJKAdd2(1)=2
                  else
                    write(*,*) "wrong, quit....,Imin,ib,im", Imin,1,NDim(1)
                    stop
                  endif
                elseif( Jmin==Jmax )then
                  if( Jmin==1)then
                    AConnect%IJKAdd(2)= 1  !!!??????
                    AConnect%IJKAdd1(2)= 2
                    AConnect%IJKAdd2(2)= 3
                    AConnect%CIJKAdd(2)=0 
                    AConnect%CIJKAdd1(2)=1
                    AConnect%CIJKAdd2(2)=2
                    AConnect%CEIJKAdd(2)=-1
                    AConnect%CEIJKAdd1(2)=-2
                    AConnect%CEIJKAdd2(2)=-3
                  elseif( Jmin==NDim(2))then
                    AConnect%IJKAdd(2)= -1
                    AConnect%IJKAdd1(2)= -2
                    AConnect%IJKAdd2(2)= -3
                    AConnect%CIJKAdd(2)=-1
                    AConnect%CIJKAdd1(2)=-2
                    AConnect%CIJKAdd2(2)=-3
                    AConnect%CEIJKAdd(2)=0
                    AConnect%CEIJKAdd1(2)=1
                    AConnect%CEIJKAdd2(2)=2
                  else
                    write(*,*) "wrong, quit....,Jmin,jb,jm", Jmin,1,NDim(2)
                    stop
                  endif
                elseif( Kmin==Kmax )then
                  if( Kmin==1)then
                    AConnect%IJKAdd(3)= 1  !!!??????
                    AConnect%IJKAdd1(3)= 2
                    AConnect%IJKAdd2(3)= 3
                    AConnect%CIJKAdd(3)=0 
                    AConnect%CIJKAdd1(3)=1
                    AConnect%CIJKAdd2(3)=2
                    AConnect%CEIJKAdd(3)=-1
                    AConnect%CEIJKAdd1(3)=-2
                    AConnect%CEIJKAdd2(3)=-3
                  elseif( Kmin==NDim(3))then
                    AConnect%IJKAdd(3)= -1
                    AConnect%IJKAdd1(3)= -2
                    AConnect%IJKAdd2(3)= -3
                    AConnect%CIJKAdd(3)=-1
                    AConnect%CIJKAdd1(3)=-2
                    AConnect%CIJKAdd2(3)=-3
                    AConnect%CEIJKAdd(3)=0
                    AConnect%CEIJKAdd1(3)=1
                    AConnect%CEIJKAdd2(3)=2
                  else
                    write(*,*) "wrong, quit....,Kmin,kb,km", Kmin,1,NDim(3),MyProc,AllBlocks(iblock)%NBlockGlobal,M_BlockIds(iblock)
                    stop
                  endif
                endif
                 AConnect%CIJKAdd(1:3)= (AConnect%IJKAdd(1:3)-1)/2  !! Imin=Imax=ib: 0; =im: -1  inner side
                 AConnect%CEIJKAdd(1:3)=(-AConnect%IJKAdd(1:3)-1)/2 !! Imin=Imax=ib:-1; =im:  0  outer side
            enddo
        enddo
!*******************************Set Transform matrix********************************
    do iblock=1,NBlock
        Ncon=AllBlocks(iblock)%NConnect
        do icon=1,Ncon
            AConnect=>AllBlocks(iblock)%AConnectivity(icon)
            IJKDirect=AConnect%IJKDirect
            lreva=AConnect%lreva
            lrevb=AConnect%lrevb
            Transform(1:3)=(/0,0,0/)
            TransformT=(/0,0,0/)
            revMatrix=(/0,0,0/)

            if(IJKDirect==1)then
                Transform(1)=1
                revMatrix(2)=lreva
                revMatrix(3)=lrevb
            elseif(IJKDirect==2)then
                Transform(2)=1
                revMatrix(1)=lreva
                revMatrix(3)=lrevb
            elseif(IJKDirect==3)then
                Transform(3)=1
                revMatrix(1)=lreva
                revMatrix(2)=lrevb
            else
                write(*,*)"Wrong IJKDirect!",MyProc,iblock,icon
            endif
            AConnect%Transform=Transform
            AConnect%revMatrix=revMatrix
            if(AConnect%BCType==9)then
                Prange(1:3)=AConnect%TargetPStart(1:3)
                Prange(4:6)=AConnect%TargetPEnd(1:3)
!            write(*,*)"tarPE",MyProc,iblock,icon,"range",Prange
                if(Prange(1)==Prange(4))then
                    TransformT(1)=1
                elseif(Prange(2)==Prange(5))then
                    TransformT(2)=1
                elseif(Prange(3)==Prange(6))then
                    TransformT(3)=1
                else
                    write(*,*)"Target Transform error!",MyProc,iblock,icon
                endif
            AConnect%TransformT=TransformT
            endif
        enddo
    enddo        
end subroutine
