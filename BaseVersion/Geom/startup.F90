subroutine IJKStartEnd(I,j,k,Is,Js,Ks,TransformL,TransformT,lcross)
    implicit none
    integer::I,j,k,Is,Js,Ks,lcross
    integer::TransformL(3),TransformT(3)
    integer::gt1,gt2,gt3,ft1,ft2,ft3,temp1,temp2,temp3,lcross1
        
        ft1=TransformL(1)
        ft2=TransformL(2)
        ft3=TransformL(3)
        gt1=TransformT(1)
        gt2=TransformT(2)
        gt3=TransformT(3)
        lcross1=1-lcross

        temp1=ft1*I+ft2*J+ft3*K
        temp2=ft1*(j*lcross1+k*lcross)+ft2*(i*lcross1+k*lcross)+ft3*(i*lcross1+j*lcross)
        temp3=ft1*(j*lcross+k*lcross1)+ft2*(i*lcross+k*lcross1)+ft3*(i*lcross+j*lcross1)

        Is=gt1*temp1+(gt2+gt3)*temp2
        Js=gt2*temp1+gt1*temp2+gt3*temp3
!    write(*,*)"gt",gt1,gt2,gt3,temp1,temp2
        Ks=gt3*temp1+(gt1+gt2)*temp3

    return
end subroutine

subroutine extension
      use Global
      implicit none
!      include 'mpif.h'
      integer  i,j,k, iblock, iL,jL,kL,is,js,ks, NI,NJ,NK,NI1,NJ1,NK1,&
          &     istart, jstart, kstart,  iend, jend, kend,  &
          &     istartT,jstartT,kstartT, iendT,jendT,kendT, &
          &     TargetProc,NTargetB,NTargetConn, &
          &     iadd,jadd,kadd, Transform(3,3), NT(3), icon, &
          &     ibuf,jbuf,kbuf, iLen,jLen,kLen, &
          &     tag,NumGrid,Offset(3), &
          &     ibc,MinMax, PerType,BCType
      real :: dcita,sint,cost,ttx,tty
      real,allocatable :: bufs(:,:,:,:),bufR(:,:,:,:)
      type(ConnectivityStruct),pointer :: AConnect,tempCon
      integer ::  sendq
      integer::revMatrix(3),TransformL(3),TransformT(3),lcross,IIadd,JJadd,KKadd
      character(len=50)::fn
      
      !! init MarkBC
!      do iblock=1, NBlock
!           call pointerSet(iblock)
!        
!           allocate(thisBlock%MarkBCI(2,jq,kq))  !! imin,imax face BC mark
!           allocate(thisBlock%MarkBCJ(2,iq,kq))
!           allocate(thisBlock%MarkBCK(2,iq,jq))
!           thisBlock%MarkBCI(:,:,:)= 1  !! init it to Non-BC
!           thisBlock%MarkBCJ(:,:,:)= 1
!           thisBlock%MarkBCK(:,:,:)= 1
!           do icon=1, ThisBlock%NConnect
!              tempCon=> ThisBlock%AConnectivity(icon)
!              if( tempCon%BCType/=1 )cycle  ! BCType=1 for Wall BC
!                  istart= tempCon%PointStart(1)
!                  jstart= tempCon%PointStart(2)
!                  kstart= tempCon%PointStart(3)
!                  iend  = tempCon%PointEnd(1)
!                  jend  = tempCon%PointEnd(2)
!                  kend  = tempCon%PointEnd(3)
!                 if(tempCon%IJKDirect==1)then
!                    minmax=1
!                   if(istart==im) minmax=2
!                   thisBlock%MarkBCI(minmax,jstart:jend,kstart:kend)=0
!                elseif(tempCon%IJKDirect==2)then
!                    minmax=1
!                if(jstart==jm) minmax=2
!                thisBlock%MarkBCJ(minmax,istart:iend,kstart:kend)=0
!                elseif(tempCon%IJKDirect==3)then
!                   minmax=1
!                  if(kstart==km) minmax=2
!                  thisBlock%MarkBCK(minmax,istart:iend,jstart:jend)=0
!               endif
!         enddo
!      enddo
      
      
      !! set ghost cell grid
!      do iblock=1, NBlock
!      call pointerSet(iblock)
    do iblock=1,Max_Block
        call Global_SetPointers(iblock)
        NI=ThisBlock%NI
        NJ=ThisBlock%NJ
        NK=ThisBlock%NK

      do k=1,NK  !! init the ghost cell grid
      do j=1,NJ
        i=1
        xx(i-1,j,k)= 2*xx(i,j,k)-xx(i+1,j,k)
        yy(i-1,j,k)= 2*yy(i,j,k)-yy(i+1,j,k)
        zz(i-1,j,k)= 2*zz(i,j,k)-zz(i+1,j,k)
!        xx(i-1,j,k)= xx(i,j,k)
!        yy(i-1,j,k)= yy(i,j,k)
!        zz(i-1,j,k)= zz(i,j,k)
        i=NI
        xx(i+1,j,k)= 2*xx(i,j,k)-xx(i-1,j,k)
        yy(i+1,j,k)= 2*yy(i,j,k)-yy(i-1,j,k)
        zz(i+1,j,k)= 2*zz(i,j,k)-zz(i-1,j,k)
!        xx(i+1,j,k)= xx(i,j,k)
!        yy(i+1,j,k)= yy(i,j,k)
!        zz(i+1,j,k)= zz(i,j,k)
      enddo
      enddo
      do k=1,NK
      do i=0,NI+1
        j=1
        xx(i,j-1,k)= 2*xx(i,j,k)-xx(i,j+1,k)
        yy(i,j-1,k)= 2*yy(i,j,k)-yy(i,j+1,k)
        zz(i,j-1,k)= 2*zz(i,j,k)-zz(i,j+1,k)
!        xx(i,j-1,k)= xx(i,j,k)
!        yy(i,j-1,k)= yy(i,j,k)
!        zz(i,j-1,k)= zz(i,j,k)
        j=NJ
        xx(i,j+1,k)= 2*xx(i,j,k)-xx(i,j-1,k)
        yy(i,j+1,k)= 2*yy(i,j,k)-yy(i,j-1,k)
        zz(i,j+1,k)= 2*zz(i,j,k)-zz(i,j-1,k)
!        xx(i,j+1,k)= xx(i,j,k)
!        yy(i,j+1,k)= yy(i,j,k)
!        zz(i,j+1,k)= zz(i,j,k)
      enddo
      enddo
      do j=0,NJ+1
      do i=0,NI+1
        k=1
        xx(i,j,k-1)= 2*xx(i,j,k)-xx(i,j,k+1)
        yy(i,j,k-1)= 2*yy(i,j,k)-yy(i,j,k+1)
        zz(i,j,k-1)= 2*zz(i,j,k)-zz(i,j,k+1)
!        xx(i,j,k-1)=xx(i,j,k)
!        yy(i,j,k-1)=yy(i,j,k)
!        zz(i,j,k-1)=zz(i,j,k)
        k=NK
        xx(i,j,k+1)= 2*xx(i,j,k)-xx(i,j,k-1)
        yy(i,j,k+1)= 2*yy(i,j,k)-yy(i,j,k-1)
        zz(i,j,k+1)= 2*zz(i,j,k)-zz(i,j,k-1)
!        xx(i,j,k+1)=xx(i,j,k)
!        yy(i,j,k+1)=yy(i,j,k)
!        zz(i,j,k+1)=zz(i,j,k)
      enddo
      enddo

      enddo
      !! -----  send grid at ghostcell
!      do iblock=1, NBlock
!      call pointerSet(iblock)
    do iblock=1,Max_Block
        call Global_SetPointers(iblock)

      do icon=1, AllBlocks(iblock)%NConnect
         AConnect=> AllBlocks(iblock)%AConnectivity(icon)
         BCType=AConnect%BCType

        if(BCType==9)then   !for Connectivity 
         NTargetB   = Aconnect%TargetBlock
         NTargetConn= Aconnect%TargetConn
        istart= AConnect%PointStart(1)
        jstart= AConnect%PointStart(2)
        kstart= AConnect%PointStart(3)
        iend  = AConnect%PointEnd(1)
        jend  = AConnect%PointEnd(2)
        kend  = AConnect%PointEnd(3)
!        if(istart/=iend) iend=iend+1  !!?? side range to point range
!        if(jstart/=jend) jend=jend+1
!        if(kstart/=kend) kend=kend+1
!        istart=istart-1; iend=iend+1
!        jstart=jstart-1; jend=jend+1
!        kstart=kstart-1; kend=kend+1
        if(istart/=iend) then
            iend=iend+1  !!?? side range to point range
!            istart=istart-1; iend=iend+1
        endif
        if(jstart/=jend) then
            jend=jend+1
!            jstart=jstart-1; jend=jend+1
        endif
        if(kstart/=kend) then
            kend=kend+1
!            kstart=kstart-1; kend=kend+1
        endif
        iadd  = AConnect%IJKadd(1)
        jadd  = AConnect%IJKadd(2)
        kadd  = AConnect%IJKadd(3)
        
!         write(*,*)"ijk",istart,iend,jstart,jend,kstart,kend
        allocate(bufs(3,iend-istart+1,jend-jstart+1,kend-kstart+1))
        do k=kstart, kend
        do j=jstart, jend
        do i=istart, iend
          is=i+iadd
          js=j+jadd
          ks=k+kadd
          ibuf=i-istart+1
          jbuf=j-jstart+1
          kbuf=k-kstart+1
          bufs(1,ibuf,jbuf,kbuf)= xx(is,js,ks)
          bufs(2,ibuf,jbuf,kbuf)= yy(is,js,ks)
          bufs(3,ibuf,jbuf,kbuf)= zz(is,js,ks)
        enddo
        enddo
        enddo

        TargetProc= PBlock(NTargetB)-1  !! -1 mean actual proc index
        tag= NTargetB*100+ NTargetConn  !! so NTargetConn<100
        numGrid= 3*(iend-istart+1)*(jend-jstart+1)*(kend-kstart+1)
 !        write(*,*) 'Extsend',myProc,'tarProc,TarB,TarCon:',targetProc, NTargetB,NTargetConn,tag,numgrid
        call MPI_IBSend(bufs,numGrid,MPI_REAL,TargetProc,tag, &
                      MPI_COMM_WORLD,sendq, ierr)
        
        call MPI_Wait(sendq,status,ierr)

        deallocate(bufs)
      endif
     enddo
     enddo

     !! receive ghost cell grids, 
!     do iblock=1, NBlock
!      call pointerSet(iblock)
    do iblock=1,Max_Block
        call Global_SetPointers(iblock)
        NI=ThisBlock%NI
        NJ=ThisBlock%NJ
        NK=ThisBlock%NK

     do icon=1, AllBlocks(iblock)%NConnect
        AConnect=> AllBlocks(iblock)%AConnectivity(icon)
        BCType=AConnect%BCType

     if(BCType==9)then      !for Connectivity
        NTargetB= Aconnect%TargetBlock

        istart= AConnect%PointStart(1)
        jstart= AConnect%PointStart(2)
        kstart= AConnect%PointStart(3)
        iend  = AConnect%PointEnd(1)
        jend  = AConnect%PointEnd(2)
        kend  = AConnect%PointEnd(3)
        if(istart/=iend) then
            iend=iend+1  !!?? side range to point range
!            istart=istart-1; iend=iend+1
        endif
        if(jstart/=jend) then
            jend=jend+1
!            jstart=jstart-1; jend=jend+1
        endif
        if(kstart/=kend) then
            kend=kend+1
!            kstart=kstart-1; kend=kend+1
        endif
        iadd  = AConnect%IJKadd(1)
        jadd  = AConnect%IJKadd(2)
        kadd  = AConnect%IJKadd(3)
        OffSet(1:3)= AConnect%Offset(1:3)
        
        PerType= AConnect%PerType
!        dcita  = AConnect%rotateAngle(3)
        dcita  = AConnect%rotateAngle(1)
        sint   = sin(dcita)
        cost   = cos(dcita)
        
!         write(*,*) MyProc,'s,e',istart,jstart,kstart,iend,jend,kend
!         write(*,*) Myproc,'offset', (offset(i),i=1,3)
!         pause
        
        TransformT=AConnect%TransformT
        iLen= AConnect%TPLijk(1)+1*(1-TransformT(1))
        jLen= AConnect%TPLijk(2)+1*(1-TransformT(2))
        kLen= AConnect%TPLijk(3)+1*(1-TransformT(3))
        allocate(bufR(3,iLen,jLen,kLen))
        targetProc= PBlock(NTargetB)-1  !! -1 mean actual proc index
!        tag= (iblock+NblockStart-1)*100+ icon
        tag=AllBlocks(iblock)%NBlockGlobal*100+icon
        numGrid= 3*iLen*jLen*kLen
        ! write(*,*) 'rec',myProc,'tarProc,TarB,TarCon:',targetProc, iblock,icon,tag,numgrid
        call MPI_Recv(bufr,numGrid,MPI_REAL,targetProc,tag, &
                 MPI_COMM_WORLD,status,ierr)

!******************added by ydd******************
        revMatrix=AConnect%revMatrix
        TransformL=AConnect%Transform
        TransformT=AConnect%TransformT
        lcross=AConnect%lcross
        IIadd=Iend-Istart+2
        JJadd=Jend-Jstart+2
        KKadd=Kend-Kstart+2
        do k=kstart, kend
        do j=jstart, jend
        do i=istart, iend
          iL= i-istart+1
          jL= j-jstart+1
          kL= k-kstart+1
            iL=revMatrix(1)*IIadd+(1-2*revMatrix(1))*iL
            jL=revMatrix(2)*JJadd+(1-2*revMatrix(2))*jL
            kL=revMatrix(3)*KKadd+(1-2*revMatrix(3))*kL
            call IJKStartEnd(iL,jL,kL,Is,Js,Ks,TransformL,TransformT,lcross)
            xx(i-iadd,j-jadd,k-kadd)=bufr(1,is,js,ks)
            ttx=bufr(2,is,js,ks)
            tty=bufr(3,is,js,ks)
            yy(i-iadd,j-jadd,k-kadd)=ttx
            zz(i-iadd,j-jadd,k-kadd)=tty
            if(PerType.ne.0)then
                yy(i-iadd,j-jadd,k-kadd)=ttx*cost+tty*sint
                zz(i-iadd,j-jadd,k-kadd)=-ttx*sint+tty*cost
            endif
        enddo
        enddo
        enddo
!********************************************
!        do k=kstart, kend
!        do j=jstart, jend
!        do i=istart, iend
!          iL= i-istart
!          jL= j-jstart
!          kL= k-kstart
!          Transform(:,:)= 0
!          NT(:)= AConnect%Transform(:)
!          Transform(1,abs(NT(1)))= sign(1, NT(1))
!          Transform(2,abs(NT(2)))= sign(1, NT(2))
!          Transform(3,abs(NT(3)))= sign(1, NT(3))

!          is= Transform(1,1)*iL+ Transform(2,1)*jL+ Transform(3,1)*kL+ Offset(1)+1
!          js= Transform(1,2)*iL+ Transform(2,2)*jL+ Transform(3,2)*kL+ Offset(2)+1
!          ks= Transform(1,3)*iL+ Transform(2,3)*jL+ Transform(3,3)*kL+ Offset(3)+1
!		  write(*,*) 'is,js,ks',MyProc,is,js,ks,bufr(1,is,js,ks),bufr(2,is,js,ks),bufr(3,is,js,ks)
      
!          z(i-iadd,j-jadd,k-kadd)= bufr(3,is,js,ks)
!          ttx= bufr(1,is,js,ks)
!          tty= bufr(2,is,js,ks)
!          if( PerType==3 )then
!            x(i-iadd,j-jadd,k-kadd)= ttx*cost+ tty*sint
!            y(i-iadd,j-jadd,k-kadd)=-ttx*sint+ tty*cost
!          else
!            x(i-iadd,j-jadd,k-kadd)= ttx
!            y(i-iadd,j-jadd,k-kadd)= tty
!          endif
!        enddo
!        enddo
!        enddo

        deallocate(bufR)
      endif
      enddo
       !for corner
!       NI1=NI+1;   NJ1=NJ+1;   NK1=NK+1
!       xx(0,0,1:NK)=xx(1,0,1:NK)+xx(0,1,1:NK)-xx(1,1,1:NK)
!       yy(0,0,1:NK)=yy(1,0,1:NK)+yy(0,1,1:NK)-yy(1,1,1:NK)
!       zz(0,0,1:NK)=zz(1,0,1:NK)+zz(0,1,1:NK)-zz(1,1,1:NK)
!       xx(0,NJ1,1:NK)=xx(0,NJ,1:NK)+xx(1,NJ1,1:NK)-xx(1,NJ,1:NK)
!       yy(0,NJ1,1:NK)=yy(0,NJ,1:NK)+yy(1,NJ1,1:NK)-yy(1,NJ,1:NK)
!       zz(0,NJ1,1:NK)=zz(0,NJ,1:NK)+zz(1,NJ1,1:NK)-zz(1,NJ,1:NK)
!       xx(NI1,0,1:NK)=xx(NI,0,1:Nk)+xx(NI1,1,1:NK)-xx(NI,1,1:NK)
!       yy(NI1,0,1:NK)=yy(NI,0,1:Nk)+yy(NI1,1,1:NK)-yy(NI,1,1:NK)
!       zz(NI1,0,1:NK)=zz(NI,0,1:Nk)+zz(NI1,1,1:NK)-zz(NI,1,1:NK)
!       xx(NI1,NJ1,1:NK)=xx(NI,NJ1,1:NK)+xx(NI1,NJ,1:NK)-xx(NI,NJ,1:NK)
!       yy(NI1,NJ1,1:NK)=yy(NI,NJ1,1:NK)+yy(NI1,NJ,1:NK)-yy(NI,NJ,1:NK)
!       zz(NI1,NJ1,1:NK)=zz(NI,NJ1,1:NK)+zz(NI1,NJ,1:NK)-zz(NI,NJ,1:NK)
        

        !! output the block grid
        write(fn,'(I2.2)')ThisBlock%NBlockGlobal
        fn= trim("resu/extGrid"//trim(fn)//".dat")
        open(unit=120,file=fn)
        write(120,*) 'variables="x","y","z"'
        write(120,*) 'zone i=',NI+2, ' j=',NJ+2,' k=',NK+2,' f=point '
        do k=0,NK+1
        do j=0,NJ+1
        do i=0,NI+1
          write(120,*)xx(i,j,k),yy(i,j,k),zz(i,j,k)
        enddo
        enddo
        enddo
        close(120)

        write(fn,'(I2.2)')ThisBlock%NBlockGlobal
        fn= trim("resu/Grid"//trim(fn)//".dat")
        open(unit=120,file=fn)
        write(120,*) 'variables="x","y","z"'
        write(120,*) 'zone i=',NI, ' j=',NJ,' k=',NK,' f=point '
        do k=1,NK
        do j=1,NJ
        do i=1,NI
          write(120,*)xx(i,j,k),yy(i,j,k),zz(i,j,k)
        enddo
        enddo
        enddo
        close(120)
      enddo

      return
end subroutine   

subroutine GeomTransform
    use global
    implicit none
    integer::iblock,icon,istart,iend,jstart,jend,kstart,kend,i,j,k
    integer::iadd,jadd,kadd,iadd1,jadd1,kadd1,iadd2,jadd2,kadd2,Eiadd,&
        &   Ejadd,Ekadd,Eiadd1,Ejadd1,Ekadd1,Eiadd2,Ejadd2,Ekadd2,Siadd,&
        &   Sjadd,Skadd,is,js,ks,is1,js1,ks1,is2,js2,ks2,il,jl,kl,il1,jl1,&
        &   kl1,il2,jl2,kl2,ibuf,jbuf,kbuf,ilen,jlen,klen
    integer::IJKDir,LeftOrRight,NTargetB,NTargetConn, TargetProc,numGrid
    real::out_n,dcita,sint,cost,SDy,SDz
    integer::Transform(3,3),revMatrix(3),TransformL(3),TransformT(3),lcross,IIadd,JJadd,KKadd
    integer::sendq,tag,BCType
    real,allocatable::bufs(:,:,:,:),bufr(:,:,:,:)
    type(ConnectivityStruct),pointer :: AConnect,tempCon
 
! for physical bc
    do iblock=1,Max_Block
        call Global_SetPointers(iblock)
        do icon=1, AllBlocks(iblock)%NConnect
            AConnect=> AllBlocks(iblock)%AConnectivity(icon)
            BCType=AConnect%BCType
            if(BCType.ne.9)then
                istart= AConnect%PointStart(1)
                jstart= AConnect%PointStart(2)
                kstart= AConnect%PointStart(3)
                iend  = AConnect%PointEnd(1)
                jend  = AConnect%PointEnd(2)
                kend  = AConnect%PointEnd(3)
                iadd=AConnect%CIJKAdd(1)
                jadd=AConnect%CIJKAdd(2)
                kadd=AConnect%CIJKAdd(3)
                iadd1=AConnect%CIJKAdd1(1)
                jadd1=AConnect%CIJKAdd1(2)
                kadd1=AConnect%CIJKAdd1(3)
                iadd2=AConnect%CIJKAdd2(1)
                jadd2=AConnect%CIJKAdd2(2)
                kadd2=AConnect%CIJKAdd2(3)
                Eiadd=AConnect%CEIJKAdd(1)
                Ejadd=AConnect%CEIJKAdd(2)
                Ekadd=AConnect%CEIJKAdd(3)
                Eiadd1=AConnect%CEIJKAdd1(1)
                Ejadd1=AConnect%CEIJKAdd1(2)
                Ekadd1=AConnect%CEIJKAdd1(3)
                Eiadd2=AConnect%CEIJKAdd2(1)
                Ejadd2=AConnect%CEIJKAdd2(2)
                Ekadd2=AConnect%CEIJKAdd2(3)
                Siadd=AConnect%IJKAdd(1)
                Sjadd=AConnect%IJKAdd(2)
                Skadd=AConnect%IJKAdd(3)
                IJKDir=AConnect%IJKDirect

                do k=kstart,kend
                do j=jstart,jend
                do i=istart,iend
                    is=i+iadd
                    js=j+jadd
                    ks=k+kadd
                    is1=i+iadd1
                    js1=j+jadd1
                    ks1=k+kadd1
                    is2=i+iadd2
                    js2=j+jadd2
                    ks2=k+kadd2
                    il=i+Eiadd
                    jl=j+Ejadd
                    kl=k+Ekadd
                    il1=i+Eiadd1
                    jl1=j+Ejadd1
                    kl1=k+Ekadd1
                    il2=i+Eiadd2
                    jl2=j+Ejadd2
                    kl2=k+Ekadd2

                    Vol(il,jl,kl)=Vol(is,js,ks)
                    Dst(il,jl,kl)=Dst(is,js,ks)
                    Alagm(IJKDir,il,jl,kl)=Alagm(IJKDir,is,js,ks)
                    Vol(il1,jl1,kl1)=Vol(is,js,ks)
                    Vol(il2,jl2,kl2)=Vol(is,js,ks)
                    Dst(il1,jl1,kl1)=Dst(is,js,ks)
                    Dst(il2,jl2,kl2)=Dst(is,js,ks)
                    Alagm(IJKDir,il1,jl1,kl1)=Alagm(IJKDir,is,js,ks)
                    Alagm(IJKDir,il2,jl2,kl2)=Alagm(IJKDir,is,js,ks)
                                    
                    rad(il,jl,kl)=2.0*rad(is,js,ks)-rad(is1,js1,ks1)
                    rad(il1,jl1,kl1)=2.0*rad(il,jl,kl)-rad(is,js,ks)
                    rad(il2,jl2,kl2)=2.0*rad(il1,jl1,kl1)-rad(il,jl,kl)
                    thtc(il,jl,kl)=2.0*thtc(is,js,ks)-thtc(is1,js1,ks1)
                    thtc(il1,jl1,kl1)=2.0*thtc(il,jl,kl)-thtc(is,js,ks)
                    thtc(il2,jl2,kl2)=2.0*thtc(il1,jl1,kl1)-thtc(il,jl,kl)
                    if(BCType==5)then
                        rad(il,jl,kl)=rad(is,js,ks)
                        rad(il1,jl1,kl1)=rad(is,js,ks)
                        rad(il2,jl2,kl2)=rad(is,js,ks)
                        thtc(il,jl,kl)=thtc(is,js,ks)
                        thtc(il1,jl1,kl1)=thtc(is,js,ks)
                        thtc(il2,jl2,kl2)=thtc(is,js,ks)
                    endif
                    
                    is=i+Siadd
                    js=j+Sjadd
                    ks=k+Skadd
                    il=i-Siadd
                    jl=j-Sjadd
                    kl=k-Skadd
                    SD(IJKDir,1:3,il,jl,kl)=2.0*SD(IJKDir,1:3,i,j,k)-SD(IJKDir,1:3,is,js,ks)
                enddo
                enddo
                enddo
            endif
        enddo
    enddo

endsubroutine 

subroutine GeomTransform1
    use global
    implicit none
    integer::iblock,icon,istart,iend,jstart,jend,kstart,kend,i,j,k
    integer::iadd,jadd,kadd,iadd1,jadd1,kadd1,iadd2,jadd2,kadd2,Eiadd,&
        &   Ejadd,Ekadd,Eiadd1,Ejadd1,Ekadd1,Eiadd2,Ejadd2,Ekadd2,Siadd,&
        &   Sjadd,Skadd,is,js,ks,is1,js1,ks1,is2,js2,ks2,il,jl,kl,il1,jl1,&
        &   kl1,il2,jl2,kl2,ibuf,jbuf,kbuf,ilen,jlen,klen
    integer::IJKDir,LeftOrRight,NTargetB,NTargetConn, TargetProc,numGrid
    real::out_n,dcita,sint,cost,SDy,SDz
    integer::Transform(3,3),revMatrix(3),TransformL(3),TransformT(3),lcross,IIadd,JJadd,KKadd
    integer::sendq,tag,BCType,PerType
    real,allocatable::bufs(:,:,:,:),bufr(:,:,:,:)
    type(ConnectivityStruct),pointer :: AConnect,tempCon
    !send connection
!write(*,*)MyProc,"phy bc start"
    do iblock=1,Max_Block
        call Global_SetPointers(iblock)
        do icon=1, AllBlocks(iblock)%NConnect
            AConnect=> AllBlocks(iblock)%AConnectivity(icon)
            BCType=AConnect%BCType
            if(BCType==9)then    
                istart= AConnect%PointStart(1)
                jstart= AConnect%PointStart(2)
                kstart= AConnect%PointStart(3)
                iend  = AConnect%PointEnd(1)
                jend  = AConnect%PointEnd(2)
                kend  = AConnect%PointEnd(3)
                iadd=AConnect%CIJKAdd(1)
                jadd=AConnect%CIJKAdd(2)
                kadd=AConnect%CIJKAdd(3)
                iadd1=AConnect%CIJKAdd1(1)
                jadd1=AConnect%CIJKAdd1(2)
                kadd1=AConnect%CIJKAdd1(3)
                iadd2=AConnect%CIJKAdd2(1)
                jadd2=AConnect%CIJKAdd2(2)
                kadd2=AConnect%CIJKAdd2(3)
                Eiadd=AConnect%CEIJKAdd(1)
                Ejadd=AConnect%CEIJKAdd(2)
                Ekadd=AConnect%CEIJKAdd(3)
                Eiadd1=AConnect%CEIJKAdd1(1)
                Ejadd1=AConnect%CEIJKAdd1(2)
                Ekadd1=AConnect%CEIJKAdd1(3)
                Eiadd2=AConnect%CEIJKAdd2(1)
                Ejadd2=AConnect%CEIJKAdd2(2)
                Ekadd2=AConnect%CEIJKAdd2(3)
                Siadd=AConnect%IJKAdd(1)
                Sjadd=AConnect%IJKAdd(2)
                Skadd=AConnect%IJKAdd(3)
                IJKDir=AConnect%IJKDirect
                LeftOrRight=AConnect%LeftOrRight

                ilen=iend-istart+1
                jlen=jend-jstart+1
                klen=kend-kstart+1
                allocate(bufs(18,ilen,jlen,klen))

                do k=kstart,kend
                do j=jstart,jend
                do i=istart,iend
                    is=i+iadd
                    js=j+jadd
                    ks=k+kadd
                    is1=i+iadd1
                    js1=j+jadd1
                    ks1=k+kadd1
                    is2=i+iadd2
                    js2=j+jadd2
                    ks2=k+kadd2

                    ibuf=i-istart+1
                    jbuf=j-jstart+1
                    kbuf=k-kstart+1
                    il1=i+Eiadd1
                    jl1=j+Ejadd1
                    kl1=k+Ekadd1
                    il2=i+Eiadd2
                    jl2=j+Ejadd2
                    kl2=k+Ekadd2
                    
                    bufs(1,ibuf,jbuf,kbuf)=Vol(is,js,ks)
                    bufs(2,ibuf,jbuf,kbuf)=Vol(is1,js1,ks1)
                    bufs(3,ibuf,jbuf,kbuf)=Vol(is2,js2,ks2)
                    bufs(4,ibuf,jbuf,kbuf)=Dst(is,js,ks)
                    bufs(5,ibuf,jbuf,kbuf)=Dst(is1,js1,ks1)
                    bufs(6,ibuf,jbuf,kbuf)=Dst(is2,js2,ks2)
                    bufs(7,ibuf,jbuf,kbuf)=Alagm(IJKDir,is,js,ks)
                    bufs(8,ibuf,jbuf,kbuf)=Alagm(IJKDir,is1,js1,ks1)
                    bufs(9,ibuf,jbuf,kbuf)=Alagm(IJKDir,is2,js2,ks2)
                    bufs(10,ibuf,jbuf,kbuf)=rad(is,js,ks)
                    bufs(11,ibuf,jbuf,kbuf)=rad(is1,js1,ks1)
                    bufs(12,ibuf,jbuf,kbuf)=rad(is2,js2,ks2)
                    bufs(13,ibuf,jbuf,kbuf)=thtc(is,js,ks)
                    bufs(14,ibuf,jbuf,kbuf)=thtc(is1,js1,ks1)
                    bufs(15,ibuf,jbuf,kbuf)=thtc(is2,js2,ks2)
                    
                    is=i+Siadd; js=j+Sjadd; ks=k+Skadd
                    out_n=real(2*LeftOrRight-3)
                    bufs(16,ibuf,jbuf,kbuf)=SD(IJKDir,1,is,js,ks)*out_n
                    bufs(17,ibuf,jbuf,kbuf)=SD(IJKDir,2,is,js,ks)*out_n
                    bufs(18,ibuf,jbuf,kbuf)=SD(IJKDir,3,is,js,ks)*out_n
                enddo
                enddo
                enddo
                NTargetB= AConnect%TargetBlock
                NTargetConn= AConnect%TargetConn
                TargetProc= PBlock(NTargetB)-1  !! -1 mean actual proc index
                tag= NTargetB*100+ NTargetConn  !! so NTargetConn<100
                numGrid=18*iLen* jLen* kLen
!                write(*,*) 'iSend',MyProc,icon,targetProc+1,tag,numGrid
                call MPI_IBSend(bufs,numGrid,MPI_REAL,TargetProc,tag,MPI_COMM_WORLD, sendq, ierr)
                call MPI_Wait(sendq,status, ierr)
                deallocate(bufs)
            endif
        enddo
    enddo
    
    !receive connection
    do iblock=1,Max_Block
        call Global_SetPointers(iblock)
        do icon=1, AllBlocks(iblock)%NConnect
            AConnect=> AllBlocks(iblock)%AConnectivity(icon)
            BCType=AConnect%BCType
            PerType=AConnect%PerType
            if(BCType==9)then    
                istart= AConnect%PointStart(1)
                jstart= AConnect%PointStart(2)
                kstart= AConnect%PointStart(3)
                iend  = AConnect%PointEnd(1)
                jend  = AConnect%PointEnd(2)
                kend  = AConnect%PointEnd(3)
                iadd=AConnect%CIJKAdd(1)
                jadd=AConnect%CIJKAdd(2)
                kadd=AConnect%CIJKAdd(3)
                iadd1=AConnect%CIJKAdd1(1)
                jadd1=AConnect%CIJKAdd1(2)
                kadd1=AConnect%CIJKAdd1(3)
                iadd2=AConnect%CIJKAdd2(1)
                jadd2=AConnect%CIJKAdd2(2)
                kadd2=AConnect%CIJKAdd2(3)
                Eiadd=AConnect%CEIJKAdd(1)
                Ejadd=AConnect%CEIJKAdd(2)
                Ekadd=AConnect%CEIJKAdd(3)
                Eiadd1=AConnect%CEIJKAdd1(1)
                Ejadd1=AConnect%CEIJKAdd1(2)
                Ekadd1=AConnect%CEIJKAdd1(3)
                Eiadd2=AConnect%CEIJKAdd2(1)
                Ejadd2=AConnect%CEIJKAdd2(2)
                Ekadd2=AConnect%CEIJKAdd2(3)
                Siadd=AConnect%IJKAdd(1)
                Sjadd=AConnect%IJKAdd(2)
                Skadd=AConnect%IJKAdd(3)
                IJKDir=AConnect%IJKDirect
                LeftOrRight=AConnect%LeftOrRight
                dcita=AConnect%rotateAngle(1)
                sint=sin(dcita)
                cost=cos(dcita)
                NTargetB= AConnect%TargetBlock
                TransformT=AConnect%TransformT

                iLen= AConnect%TPLijk(1)!+TransformT(1)-1
                jLen= AConnect%TPLijk(2)!+TransformT(2)-1
                kLen= AConnect%TPLijk(3)!+TransformT(3)-1
                allocate(bufR(18,iLen,jLen,kLen))
                targetProc= PBlock(NTargetB)-1  !! -1 mean actual proc index
                tag=ThisBlock%NBlockGlobal*100+icon
                numGrid=18*iLen*jLen*kLen
                call MPI_Recv(bufr,numGrid,MPI_REAL,targetProc,tag,MPI_COMM_WORLD,status,ierr)

                revMatrix=AConnect%revMatrix
                TransformL=AConnect%Transform
                TransformT=AConnect%TransformT
                lcross=AConnect%lcross
                IIadd=Iend-Istart+2
                JJadd=Jend-Jstart+2
                KKadd=Kend-Kstart+2

                do k=kstart,kend
                do j=jstart,jend
                do i=istart,iend

                    iL=i-istart+1
                    jL=j-jstart+1
                    kL=k-kstart+1

                    iL=revMatrix(1)*IIadd+(1-2*revMatrix(1))*iL
                    jL=revMatrix(2)*JJadd+(1-2*revMatrix(2))*jL
                    kL=revMatrix(3)*KKadd+(1-2*revMatrix(3))*kL
                    call IJKStartEnd(iL,jL,kL,Is,Js,Ks,TransformL,TransformT,lcross)
                   
                    il=i+Eiadd
                    jl=j+Ejadd
                    kl=k+Ekadd
                    il1=i+Eiadd1
                    jl1=j+Ejadd1
                    kl1=k+Ekadd1
                    il2=i+Eiadd2
                    jl2=j+Ejadd2
                    kl2=k+Ekadd2

                    Vol(il,jl,kl)=bufr(1,is,js,ks) 
                    Vol(il1,jl1,kl1)=bufr(2,is,js,ks) 
                    Vol(il2,jl2,kl2)=bufr(3,is,js,ks) 
                    Dst(il,jl,kl)=bufr(4,is,js,ks) 
                    Dst(il1,jl1,kl1)=bufr(5,is,js,ks) 
                    Dst(il2,jl2,kl2)=bufr(6,is,js,ks) 
                    Alagm(IJKDir,il,jl,kl)=bufr(7,is,js,ks) 
                    Alagm(IJKDir,il1,jl1,kl1)=bufr(8,is,js,ks) 
                    Alagm(IJKDir,il2,jl2,kl2)=bufr(9,is,js,ks) 
                    rad(il,jl,kl)=bufr(10,is,js,ks) 
                    rad(il1,jl1,kl1)=bufr(11,is,js,ks) 
                    rad(il2,jl2,kl2)=bufr(12,is,js,ks) 
                    thtc(il,jl,kl)=bufr(13,is,js,ks)-dcita
                    thtc(il1,jl1,kl1)=bufr(14,is,js,ks)-dcita
                    thtc(il2,jl2,kl2)=bufr(15,is,js,ks)-dcita
                    if(AConnect%PerType.ne.0)then
                        is1=i+iadd; js1=j+jadd; ks1=k+kadd
                        is2=i+iadd1; js2=j+jadd2; ks2=k+kadd1
                        thtc(il,jl,kl)=2.0*thtc(is1,js1,ks1)-thtc(is2,js2,ks2)
                        thtc(il1,jl1,kl1)=2.0*thtc(il,jl,kl)-thtc(is1,js1,ks1)
                        thtc(il2,jl2,kl2)=2.0*thtc(il1,jl1,kl1)-thtc(il,jl,kl)
                        rad(il,jl,kl)=2.0*rad(is1,js1,ks1)-rad(is2,js2,ks2)
                        rad(il1,jl1,kl1)=2.0*rad(il,jl,kl)-rad(is1,js1,ks1)
                        rad(il2,jl2,kl2)=2.0*rad(il1,jl1,kl1)-rad(il,jl,kl)
                    endif
                    
                    il=i-Siadd; jl=j-Sjadd; kl=k-Skadd
                    out_n=real(2*LeftOrRight-3)
                    SD(IJKDir,1,il,jl,kl)=-out_n*bufr(16,is,js,ks)
                    SDy=-out_n*bufr(17,is,js,ks)
                    SDz=-out_n*bufr(18,is,js,ks)
                    SD(IJKDir,2,il,jl,kl)=cost*SDy+sint*SDz
                    SD(IJKDir,3,il,jl,kl)=-sint*SDy+cost*SDz
                    if(PerType.ne.0)then
                        is=i+Siadd;  js=j+Sjadd;  ks=k+Skadd
                        SD(IJKDir,1:3,il,jl,kl)=2.0*SD(IJKDir,1:3,i,j,k)-SD(IJKDir,1:3,is,js,ks)
                    endif

                enddo
                enddo
                enddo
                deallocate(bufr)
            endif
        enddo
    enddo
endsubroutine
        

