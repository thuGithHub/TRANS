!SUBROUTINE BC_Wall(Ibgn,Iend,Jbgn,Jend,Kbgn,Kend,Kindsub_BC ,sweep,IJK,minormax)
subroutine BC_Wall(ibc,sweep)
    USE Global
    IMPLICIT NONE
    
    INTEGER:: Ibgn,Iend,Jbgn,Jend,Kbgn,Kend,i,j,k,l
    INTEGER:: sweep,Kindsub_BC
    integer::iBC,is,js,ks,il,jl,kl,iadd,jadd,kadd,iadd1,jadd1,kadd1,iadd2,jadd2,kadd2,&
        &   Eiadd,Ejadd,Ekadd,Eiadd1,Ejadd1,Ekadd1,Eiadd2,Ejadd2,Ekadd2
    type(ConnectivityStruct),pointer :: Aconnect
    REAL:: Dn,Den,TT,Rmiuw,Omegaw
    real::vibc,vjbc,vkbc,ytemp,ztemp

        AConnect=>ThisBlock%AConnectivity(ibc)
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
        Kindsub_BC=AConnect%SubBCType

        if(sweep==-1)then
          DO K=Kbgn,Kend
          DO J=Jbgn,Jend
          DO I=Ibgn,Iend
            is=i+iadd;  js=j+jadd;  ks=k+kadd
            il=i+Eiadd; jl=j+Ejadd; kl=k+Ekadd

            ztemp=0.5*(Zc(is,js,ks)+Zc(il,jl,kl))
            ytemp=0.5*(Yc(is,js,ks)+Yc(il,jl,kl))
            vibc=0.0
            vjbc=real(Kindsub_BC)*omega(1)*ztemp
            vkbc=-real(Kindsub_BC)*omega(1)*ytemp
!            if(kindsub_BC==1)then
!                vibc=0.0!gridV(IJK,1,I_VR0,J_VR0,K_VR0)
!                vjbc=omega(1)*Zc(i_R0,J_R0,K_R0)!gridV(IJK,2,I_VR0,J_VR0,K_VR0)
!                vkbc=-omega(1)*Yc(I_R0,J_R0,K_R0)!gridV(IJK,3,I_VR0,J_VR0,K_VR0)
!            else
!                vibc=0.0
!                vjbc=0.0
!                vkbc=0.0
!            endif
            PP(il,jl,kl)=PP(is,js,ks)
            T(il,jl,kl)=T(is,js,ks)
            V(1,il,jl,kl)=PP(il,jl,kl)/(T(il,jl,kl)+Tiny)*RXM2 
            den=V(1,il,jl,kl)
            V(2,il,jl,kl)=vibc*den
            V(3,il,jl,kl)=vjbc*den
            V(4,il,jl,kl)=vkbc*den
            Rds(1:3,il,jl,kl)=Rds(1:3,is,js,ks)            
 
            is=i+iadd1;  js=j+jadd1;  ks=k+kadd1
            il=i+Eiadd1; jl=j+Ejadd1; kl=k+Ekadd1

            PP(il,jl,kl)=PP(is,js,ks)
            T(il,jl,kl)=T(is,js,ks)
            V(1,il,jl,kl)=PP(il,jl,kl)/(T(il,jl,kl)+Tiny)*RXM2 
            den=V(1,il,jl,kl)
            V(2,il,jl,kl)=vibc*den
            V(3,il,jl,kl)=vjbc*den
            V(4,il,jl,kl)=vkbc*den
   
            is=i+iadd2;  js=j+jadd2;  ks=k+kadd2
            il=i+Eiadd2; jl=j+Ejadd2; kl=k+Ekadd2
            PP(il,jl,kl)=PP(is,js,ks)
            T(il,jl,kl)=T(is,js,ks)
            V(1,il,jl,kl)=PP(il,jl,kl)/(T(il,jl,kl)+Tiny)*RXM2 
            den=V(1,il,jl,kl)
            V(2,il,jl,kl)=vibc*den
            V(3,il,jl,kl)=vjbc*den
            V(4,il,jl,kl)=vkbc*den

            if(IF_turb)then
                is=i+iadd;  js=j+jadd;  ks=k+kadd
                il=i+Eiadd; jl=j+Ejadd; kl=k+Ekadd
                Dn=Dst(il,jl,kl)
                Den=V(1,il,jl,kl)
                TT=T(il,jl,kl)
                Rmiuw = (1.0+Csthlnd)/(TT*Csth1+Csthlnd)*(Csth1*TT)**1.5
                Omegaw= 60.0*Rmiuw/(Den*0.075*Dn*Dn*ref*ref)
                Rmiuw=Rmiu(is,js,ks)
                
                V(6,il,jl,kl)=0.0
                V(7,il,jl,kl)=Omegaw
                Rmiu(il,jl,kl)=Rmiuw

                is=i+iadd1;  js=j+jadd1;  ks=k+kadd1
                il=i+Eiadd1; jl=j+Ejadd1; kl=k+Ekadd1
                V(6,il,jl,kl)=0.0
                V(7,il,jl,kl)=Omegaw
                Rmiu(il,jl,kl)=Rmiuw

                is=i+iadd2;  js=j+jadd2;  ks=k+kadd2
                il=i+Eiadd2; jl=j+Ejadd2; kl=k+Ekadd2
                V(6,il,jl,kl)=0.0
                V(7,il,jl,kl)=Omegaw
                Rmiu(il,jl,kl)=Rmiuw
            endif
          enddo
          enddo
          enddo
        elseif(sweep==0)then
            do K=Kbgn,Kend
            do J=Jbgn,Jend
            do I=Ibgn,Iend
                is=i+iadd;  js=j+jadd;  ks=k+kadd
                il=i+Eiadd; jl=j+Ejadd; kl=k+Ekadd
                do L=1,18
                    DqDxyz(L,il,jl,kl)=DqDxyz(L,is,js,ks)
                enddo
            enddo
            enddo
            enddo
        endif

END SUBROUTINE BC_Wall
