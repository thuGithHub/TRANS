SUBROUTINE PARALLEL_BC_sendrecv_blk(iBlock,icon, ImpiTAG, Nmpiproc,sweep )
!subroutine P2PBC(sweep)
    use Global
    implicit none
    INTEGER iBlock,icon, ImpiTAG, Nmpiproc,sweep
    INTEGER i,j,k,l,LL,Ibgn, Iend, Jbgn,Jend, Kbgn, Kend,is,js,ks,il,jl,kl,iadd,jadd,kadd,iadd1,jadd1, &
            &   kadd1,iadd2,jadd2,kadd2,Eiadd,Ejadd,Ekadd,Eiadd1,Ejadd1,Ekadd1,Eiadd2,Ejadd2,Ekadd2
    integer::lcross,IIadd,JJadd,KKadd,ilen,jlen,klen,PerType
    integer::revMatrix(3),TransformL(3),TransformT(3)
    real::spec(3)
    INTEGER:: NSizeI, NSizeJ, NSizeK, NSizeIJK, NIJKmax,NVnum,II
    type(ConnectivityStruct),pointer :: AConnect
    real::dstage,sint,cost,sin2t,cos2t,sct,ugx,ugy,ugz,vgx,vgy,vgz,wgx,wgy,wgz,txy,txz
    integer::STag,Rtag,DestBlock,SourceBlock,DestProc,SourceProc

!    do iblock=1,MAX_Block
      CALL GLOBAL_SetPointers(iBlock)
!    do icon=1,ThisBlock%NConnect
        AConnect=> ThisBlock%AConnectivity(icon)
        DestBlock=AConnect%TargetBlock
        SourceBlock=ThisBlock%NBlockGlobal
        DestProc=PBlock(DestBlock)-1
        SourceProc=PBlock(SourceBlock)-1
!    if(AConnect%BCType==9.and.(SourceProc.ne.DestProc))then
        Ibgn=AConnect%PointStart(1)
        Jbgn=AConnect%PointStart(2)
        Kbgn=AConnect%PointStart(3)
        Iend=AConnect%PointEnd(1)
        Jend=AConnect%PointEnd(2)
        Kend=AConnect%PointEnd(3)
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
        PerType=AConnect%PerType 

        Rtag=100*SourceBlock+icon
        STag=100*DestBlock+AConnect%TargetConn
        

        if (sweep == -1) then
            NVnum = 36
        elseif(sweep==0)then    !for gradient
            NVnum=18        !u,v,w,T,k,Omega
        endif
        NSizeI = Iend-Ibgn+1
        NSizeJ = Jend-Jbgn+1
        NSizeK = Kend-Kbgn+1

        NSizeIJK=NSizeI*NSizeJ*NSizeK
        NIJKmax=NVnum*NSizeIJK
        II = 0
!        dstage=real(PerType)*2*Pi/real(nblade)
        dstage=AConnect%rotateAngle(1)
        sint=sin(dstage)
        cost=cos(dstage)        
        cos2t=cost*cost
        sin2t=sint*sint
        sct=sint*cost
                
        IF (sweep == -1) THEN
            do K=Kbgn,Kend
            do J=Jbgn,Jend
            do I=Ibgn,Iend
                is=i+iadd;  js=j+jadd;  ks=k+kadd
                do L=1,7
                    A(II+L)=V(L,is,js,ks)
                enddo
                A(II+8)=PP(is,js,ks)
                A(II+9)=T(is,js,ks)
                A(II+10)=Rmiu(is,js,ks)
                A(II+11)=Rds(1,is,js,ks)
                A(II+12)=Rds(2,is,js,ks)
                A(II+13)=Rds(3,is,js,ks)
                A(II+14)=shock(is,js,ks)

                is=i+iadd1;  js=j+jadd1;  ks=k+kadd1
                do L=1,7
                    A(II+L+14)=V(L,is,js,ks)
                enddo
                A(II+22)=PP(is,js,ks)
                A(II+23)=T(is,js,ks)
                A(II+24)=Rmiu(is,js,ks)
                A(II+25)=shock(is,js,ks)

                is=i+iadd2;  js=j+jadd2;  ks=k+kadd2
                do L=1,7
                    A(II+L+25)=V(L,is,js,ks)
                enddo
                A(II+33)=PP(is,js,ks)
                A(II+34)=T(is,js,ks)
                A(II+35)=Rmiu(is,js,ks)
                A(II+36)=shock(is,js,ks)
                
                II = II + NVnum
            enddo
            enddo
            enddo
        else if(sweep==0)then    
            do k=kbgn,kend
            do J=jbgn,jend
            do I=ibgn,iend
                is=i+iadd;  js=j+jadd;  ks=k+kadd
                do L=1,18
                    A(II+L)=DqDxyz(L,is,js,ks)
                enddo

                II=II+NVnum
            enddo
            enddo
            enddo
        endif
!	        send & recv
        CALL MPI_SENDRECV(A,NIJKmax,MPI_REAL,Nmpiproc,ImpiTAG,B,NIJKmax,MPI_REAL,Nmpiproc,ImpiTAG,MPI_COMM_WORLD,status,ierr)
!        CALL MPI_SENDRECV(A,NIJKmax,MPI_REAL,DestProc,STag,B,NIJKmax,MPI_REAL,DestProc,Rtag,MPI_COMM_WORLD,status,ierr)
!
        
        revMatrix=AConnect%revMatrix
        TransformL=AConnect%Transform
        TransformT=AConnect%TransformT
        lcross=AConnect%lcross
        IIadd=Iend-Ibgn+2
        JJadd=Jend-Jbgn+2
        KKadd=Kend-Kbgn+2
        ilen=AConnect%TargetPEnd(1)-AConnect%TargetPStart(1)+1
        jlen=AConnect%TargetPEnd(2)-AConnect%TargetPStart(2)+1
        klen=AConnect%TargetPEnd(3)-AConnect%TargetPStart(3)+1
        if(sweep==-1)then
            do K=Kbgn,Kend
            do J=Jbgn,Jend
            do I=Ibgn,Iend
            !    CALL BC_SETMPIrnum(nSizeI,nSizeJ,nSizeK,Ibgn,Jbgn,Kbgn,I,J,K,IJK,IF_lcross,IF_adv_1,IF_adv_2,n_r)  
                il=i-ibgn+1; jl=j-jbgn+1;   kl=k-kbgn+1
                iL=revMatrix(1)*IIadd+(1-2*revMatrix(1))*iL
                jL=revMatrix(2)*JJadd+(1-2*revMatrix(2))*jL
                kL=revMatrix(3)*KKadd+(1-2*revMatrix(3))*kL
                call IJKStartEnd(iL,jL,kL,Is,Js,Ks,TransformL,TransformT,lcross)

                II = ((Ks-1)*Ilen*Jlen+(Js-1)*Ilen+Is-1)*NVnum
!        if(ThisBlock%NBlockGlobal==2.and.PerType.ne.0)then
!            write(*,*)"PerType1",il,jl,kl,DestBlock,AConnect%TargetConn
!        endif
        

                il=i+Eiadd; jl=j+Ejadd; kl=k+Ekadd
                V(1:2,il,jl,kl) = B(II+1:II+2)
                txy=B(II+3)
                txz=B(II+4)
                V(3,il,jl,kl)=cost*txy+sint*txz
                V(4,il,jl,kl)=-sint*txy+cost*txz
                V(5:7,il,jl,kl)=B(II+5:II+7)
                PP(il,jl,kl)=B(II+8)
                T(il,jl,kl)=B(II+9)
                Rmiu(il,jl,kl)=B(II+10)
!        if(ThisBlock%NBlockGlobal==2.and.PerType.ne.0)then
!            write(*,*)"PerType",PerType,il,jl,kl,is,js,ks
!        endif
                spec(1)=B(II+11)
                spec(2)=B(II+12)
                spec(3)=B(II+13)
                call IJKStartEnd(1,2,3,Is,Js,Ks,TransformL,TransformT,lcross)
                Rds(is,il,jl,kl)=spec(1)
                Rds(js,il,jl,kl)=spec(2)
                Rds(ks,il,jl,kl)=spec(3)
                shock(il,jl,kl)=B(II+14)

                il=i+Eiadd1; jl=j+Ejadd1; kl=k+Ekadd1
                V(1:2,il,jl,kl) = B(II+15:II+16)
                txy=B(II+17)
                txz=B(II+18)
                V(3,il,jl,kl)=cost*txy+sint*txz
                V(4,il,jl,kl)=-sint*txy+cost*txz
                V(5:7,il,jl,kl)=B(II+19:II+21)
                PP(il,jl,kl)=B(II+22)
                T(il,jl,kl)=B(II+23)
                Rmiu(il,jl,kl)=B(II+24)
                shock(il,jl,kl)=B(II+25)
                
                il=i+Eiadd2; jl=j+Ejadd2; kl=k+Ekadd2
                V(1:2,il,jl,kl) = B(II+26:II+27)
                txy=B(II+28)
                txz=B(II+29)
                V(3,il,jl,kl)=cost*txy+sint*txz
                V(4,il,jl,kl)=-sint*txy+cost*txz
                V(5:7,il,jl,kl)=B(II+30:II+32)
                PP(il,jl,kl)=B(II+33)
                T(il,jl,kl)=B(II+34)
                Rmiu(il,jl,kl)=B(II+35)
                shock(il,jl,kl)=B(II+36)        

            enddo
            enddo
            enddo
        elseif(sweep == 0)then
            do k=kbgn,kend
            do j=jbgn,jend
            do i=ibgn,iend
                il=i-ibgn+1; jl=j-jbgn+1;   kl=k-kbgn+1
                iL=revMatrix(1)*IIadd+(1-2*revMatrix(1))*iL
                jL=revMatrix(2)*JJadd+(1-2*revMatrix(2))*jL
                kL=revMatrix(3)*KKadd+(1-2*revMatrix(3))*kL
                call IJKStartEnd(iL,jL,kL,Is,Js,Ks,TransformL,TransformT,lcross)

                II = ((Ks-1)*Ilen*Jlen+(Js-1)*Ilen+Is-1)*NVnum

                il=i+Eiadd; jl=j+Ejadd; kl=k+Ekadd
                    ugx=B(II+1)
                    ugy=B(II+2)
                    ugz=B(II+3)
                    vgx=B(II+4)
                    vgy=B(II+5)
                    vgz=B(II+6)
                    wgx=B(II+7)
                    wgy=B(II+8)
                    wgz=B(II+9)
            
                    dQdxyz(1,il,jl,kl)=ugx
                    dQdxyz(2,il,jl,kl)=cost*ugy+sint*ugz
                    dQdxyz(3,il,jl,kl)=-sint*ugy+cost*ugz
                    dQdxyz(4,il,jl,kl)=cost*vgx+sint*wgx
                    dQdxyz(5,il,jl,kl)=cos2t*vgy+sct*(vgz+wgy)+sin2t*wgz
                    dQdxyz(6,il,jl,kl)=cos2t*vgz+sct*(-vgy+wgz)-sin2t*wgy
                    dQdxyz(7,il,jl,kl)=-sint*vgx+cost*wgx
                    dQdxyz(8,il,jl,kl)=-sin2t*vgz+sct*(-vgy+wgz)+cos2t*wgy
                    dQdxyz(9,il,jl,kl)=sin2t*vgy+sct*(vgz+wgy)+cos2t*wgz
            
                    do LL=4,6
                        L=(LL-1)*3+1
                        dQdxyz(L,il,jl,kl) =B(II+L) 
                        ugy=B(II+L+1)
                        ugz=B(II+L+2)
                        dQdxyz(L+1,il,jl,kl)=cost*ugy+sint*ugz
                        dQdxyz(L+2,il,jl,kl)=-sint*ugy+cost*ugz
                    enddo
            enddo
            enddo
            enddo
        end if
!    endif
!            CALL MPI_BARRIER( MPI_COMM_WORLD,ierr )
!    enddo
!    enddo
END SUBROUTINE  ! PARALLEL_BC_sendrecv_blk
