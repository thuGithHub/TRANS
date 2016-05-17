subroutine BC_P2P(ibc,iblock_r,sweep)
    USE Global
    IMPLICIT NONE
    integer::ibc,iblock_r,sweep
    INTEGER i,j,k,l,LL,Ibgn, Iend, Jbgn,Jend, Kbgn, Kend,is,js,ks,il,jl,kl,iadd,jadd,kadd,iadd1,jadd1, &
            &   kadd1,iadd2,jadd2,kadd2,Eiadd,Ejadd,Ekadd,Eiadd1,Ejadd1,Ekadd1,Eiadd2,Ejadd2,Ekadd2
    integer::lcross,IIadd,JJadd,KKadd,PerType,Targetibc,ibuf,jbuf,kbuf,Tibgn,Tjbgn,Tkbgn
    integer::revMatrix(3),TransformL(3),TransformT(3)
    type(ConnectivityStruct),pointer :: AConnect,BConnect
    TYPE(BlockStruct), POINTER:: Block_r
    REAL, POINTER:: V_r(:,:,:,:)
    REAL, POINTER:: PP_r(:,:,:), T_r(:,:,:), Rmiu_r(:,:,:)
    REAL, POINTER:: Rds_r(:,:,:,:)
    REAL, POINTER:: dQdxyz_r(:,:,:,:)
    REAL, POINTER:: Shock_r(:,:,:)
    REAL, POINTER:: VL_r(:,:,:,:,:)
    real::dstage,sint,cost,sin2t,cos2t,sct,ugx,ugy,ugz,vgx,vgy,vgz,wgx,wgy,wgz
    real::spec(3),txy,txz

        Block_r => AllBlocks(iBlock_r)
        V_r => Block_r%V
        PP_r => Block_r%PP
        T_r => Block_r%T
        Rmiu_r => Block_r%Rmiu
        Rds_r => Block_r%Rds
        dQdxyz_r => Block_r%dQdxyz
        shock_r =>  Block_r%shock
        VL_r => Block_r%VL
        AConnect=>ThisBlock%AConnectivity(ibc)
        Targetibc=AConnect%TargetConn
        BConnect=>Block_r%AConnectivity(Targetibc)
        Ibgn=AConnect%PointStart(1)
        Jbgn=AConnect%PointStart(2)
        Kbgn=AConnect%PointStart(3)
        Iend=AConnect%PointEnd(1)
        Jend=AConnect%PointEnd(2)
        Kend=AConnect%PointEnd(3)
        Tibgn=BConnect%PointStart(1)
        Tjbgn=BConnect%PointStart(2)
        Tkbgn=BConnect%PointStart(3)
        iadd=BConnect%CIJKAdd(1)
        jadd=BConnect%CIJKAdd(2)
        kadd=BConnect%CIJKAdd(3)
        iadd1=BConnect%CIJKAdd1(1)
        jadd1=BConnect%CIJKAdd1(2)
        kadd1=BConnect%CIJKAdd1(3)
        iadd2=BConnect%CIJKAdd2(1)
        jadd2=BConnect%CIJKAdd2(2)
        kadd2=BConnect%CIJKAdd2(3)
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
        IIadd=Iend-Ibgn+2
        JJadd=Jend-Jbgn+2
        KKadd=Kend-Kbgn+2
        revMatrix=AConnect%revMatrix
        TransformL=AConnect%Transform
        TransformT=AConnect%TransformT
        lcross=AConnect%lcross

!        dstage=real(PerType)*2.0*Pi/real(nblade)
        dstage=AConnect%rotateAngle(1)
        sint=sin(dstage)
        cost=cos(dstage)
        cos2t=cost*cost
        sin2t=sint*sint
        sct=sint*cost
    
        if(sweep==-1)then        
            DO K=Kbgn,Kend
            DO J=Jbgn,Jend
            DO I=Ibgn,Iend
                iL=i-ibgn+1;    jL=j-jbgn+1;    kl=k-kbgn+1
                iL=revMatrix(1)*IIadd+(1-2*revMatrix(1))*iL
                jL=revMatrix(2)*JJadd+(1-2*revMatrix(2))*jL
                kL=revMatrix(3)*KKadd+(1-2*revMatrix(3))*kL
                call IJKStartEnd(iL,jL,kL,ibuf,jbuf,kbuf,TransformL,TransformT,lcross)
                ibuf=ibuf+Tibgn-1;   jbuf=jbuf+Tjbgn-1;   kbuf=kbuf+Tkbgn-1
        
                is=ibuf+iadd;   js=jbuf+jadd;   ks=kbuf+kadd
                il=i+Eiadd;    jl=j+Ejadd;     kl=k+Ekadd

                V(1:2,il,jl,kl)=V_r(1:2,is,js,ks)
                txy=V_r(3,is,js,ks)
                txz=V_r(4,is,js,ks)                
                V(3,il,jl,kl)=cost*txy+sint*txz
                V(4,il,jl,kl)=-sint*txy+cost*txz
                V(5:7,il,jl,kl)=V_r(5:7,is,js,ks)
                PP(il,jl,kl)=PP_r(is,js,ks)
                T(il,jl,kl)=T_r(is,js,ks)
                Rmiu(il,jl,kl)=Rmiu_r(is,js,ks)
                shock(il,jl,kl)=shock_r(is,js,ks)
                spec(1:3)=Rds_r(1:3,is,js,ks)
                call IJKStartEnd(1,2,3,is,js,ks,TransformL,TransformT,lcross)
                Rds(is,il,jl,kl)=spec(1)
                Rds(js,il,jl,kl)=spec(2)
                Rds(ks,il,jl,kl)=spec(3)
    
                is=ibuf+iadd1;   js=jbuf+jadd1;   ks=kbuf+kadd1
                il=i+Eiadd1;    jl=j+Ejadd1; kl=k+Ekadd1
!        if(ThisBlock%NBlockGlobal==1.and.PerType.ne.0)then
!            write(*,*)"PerType",PerType,il,jl,kl,is,js,ks
!        endif

                V(1:2,il,jl,kl)=V_r(1:2,is,js,ks)
                txy=V_r(3,is,js,ks)
                txz=V_r(4,is,js,ks)
                V(3,il,jl,kl)=cost*txy+sint*txz
                V(4,il,jl,kl)=-sint*txy+cost*txz
                V(5:7,il,jl,kl)=V_r(5:7,is,js,ks)
                PP(il,jl,kl)=PP_r(is,js,ks)
                T(il,jl,kl)=T_r(is,js,ks)
                Rmiu(il,jl,kl)=Rmiu_r(is,js,ks)
                shock(il,jl,kl)=shock_r(is,js,ks)

                is=ibuf+iadd2;   js=jbuf+jadd2;   ks=kbuf+kadd2
                il=i+Eiadd2;    jl=j+Ejadd2; kl=k+Ekadd2

                V(1:2,il,jl,kl)=V_r(1:2,is,js,ks)
                txy=V_r(3,is,js,ks)
                txz=V_r(4,is,js,ks)
                V(3,il,jl,kl)=cost*txy+sint*txz
                V(4,il,jl,kl)=-sint*txy+cost*txz
                V(5:7,il,jl,kl)=V_r(5:7,is,js,ks)
                PP(il,jl,kl)=PP_r(is,js,ks)
                T(il,jl,kl)=T_r(is,js,ks)
                Rmiu(il,jl,kl)=Rmiu_r(is,js,ks)
                shock(il,jl,kl)=shock_r(is,js,ks)        

            enddo
            enddo
            enddo
        elseif(sweep==0)then
            do K=Kbgn,Kend
            do J=Jbgn,Jend
            do I=Ibgn,Iend
                iL=i-ibgn+1;    jL=j-jbgn+1;    kl=k-kbgn+1
                iL=revMatrix(1)*IIadd+(1-2*revMatrix(1))*iL
                jL=revMatrix(2)*JJadd+(1-2*revMatrix(2))*jL
                kL=revMatrix(3)*KKadd+(1-2*revMatrix(3))*kL
                call IJKStartEnd(iL,jL,kL,ibuf,jbuf,kbuf,TransformL,TransformT,lcross)
                ibuf=ibuf+Tibgn-1;   jbuf=jbuf+Tjbgn-1;   kbuf=kbuf+Tkbgn-1
        
                is=ibuf+iadd;   js=jbuf+jadd;   ks=kbuf+kadd
                il=i+Eiadd;    jl=j+Ejadd;     kl=k+Ekadd
                
                ugx=dQdxyz_r(1,is,js,ks)
                ugy=dQdxyz_r(2,is,js,ks)
                ugz=dQdxyz_r(3,is,js,ks)
                vgx=dQdxyz_r(4,is,js,ks)
                vgy=dQdxyz_r(5,is,js,ks)
                vgz=dQdxyz_r(6,is,js,ks)
                wgx=dQdxyz_r(7,is,js,ks)
                wgy=dQdxyz_r(8,is,js,ks)
                wgz=dQdxyz_r(9,is,js,ks)
            
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
                    dQdxyz(L,il,jl,kl) =dQdxyz_r(L,is,js,ks)
                    ugy=dQdxyz_r(L+1,is,js,ks)
                    ugz=dQdxyz_r(L+2,is,js,ks)
                    dQdxyz(L+1,il,jl,kl)=cost*ugy+sint*ugz
                    dQdxyz(L+2,il,jl,kl)=-sint*ugy+cost*ugz
                enddo
            enddo
            enddo
            enddo
        endif
END SUBROUTINE BC_P2P


