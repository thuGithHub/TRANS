SUBROUTINE BC_Symm(ibc,sweep)
    USE Global
    IMPLICIT NONE
    INTEGER:: sweep,ibc
    integer::Ibgn,Iend,Jbgn,Jend,Kbgn,Kend,i,j,k,IJKDir
    integer::is,js,ks,il,jl,kl,iadd,jadd,kadd,Eiadd,Ejadd,Ekadd,iadd1,jadd1,kadd1,iadd2,&
        &   jadd2,kadd2,Eiadd1,Ejadd1,Ekadd1,Eiadd2,Ejadd2,Ekadd2,Siadd,Sjadd,Skadd
    type(ConnectivityStruct),pointer :: Aconnect 
    REAL:: Sout1,Sout2,Sout3    
    REAL:: Den,Vxi,Vyi,Vzi,Vc1
    
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
        Siadd=AConnect%IJKAdd(1)
        Sjadd=AConnect%IJKAdd(2)
        Skadd=AConnect%IJKAdd(3)
        IJKDir=AConnect%IJKDirect

      if(sweep==-1)then
        do K=Kbgn,Kend
        do J=Jbgn,Jend
        do I=Ibgn,Iend
            is=i+Siadd; js=j+Sjadd; ks=k+Skadd
            Sout1=SD(IJKDir,1,is,js,ks)/Grad(IJKDir,is,js,ks)
            Sout2=SD(IJKDir,2,is,js,ks)/Grad(IJKDir,is,js,ks)
            Sout3=SD(IJKDir,3,is,js,ks)/Grad(IJKDir,is,js,ks)
       
            is=i+iadd;  js=j+jadd;  ks=k+kadd
            den=V(1,is,js,ks)
            Vxi=V(2,is,js,ks)/den 
            Vyi=V(3,is,js,ks)/den 
            Vzi=V(4,is,js,ks)/den 
            Vc1  = Vxi*Sout1+Vyi*Sout2+Vzi*Sout3
!	     Vc1  = Vxi*Sout1+Vyi*Sout2+Vzi*Sout3  !contravariant velocity
            il=i+Eiadd; jl=j+Ejadd; kl=k+Ekadd
        
            V(1,il,jl,kl ) = V(1,is,js,ks) 
            V(2,il,jl,kl ) = V(2,is,js,ks) -2* Sout1* Vc1* V(1,is,js,ks)
            V(3,il,jl,kl ) = V(3,is,js,ks) -2* Sout2* Vc1* V(1,is,js,ks)
            V(4,il,jl,kl ) = V(4,is,js,ks) -2* Sout3* Vc1* V(1,is,js,ks)
            V(5,il,jl,kl ) = V(5,is,js,ks)
            PP(il,jl,kl ) = PP(is,js,ks)
            T(il,jl,kl ) = T(is,js,ks)

            is=i+iadd1;     js=j+jadd1;     ks=k+kadd1
            il=i+Eiadd1;    jl=j+Ejadd1;    kl=k+Ekadd1
            V(2,il,jl,kl ) = V(2,is,js,ks)-2* Sout1* Vc1* V(1,is,js,ks)
            V(3,il,jl,kl ) = V(3,is,js,ks)-2* Sout2* Vc1* V(1,is,js,ks)
            V(4,il,jl,kl ) = V(4,is,js,ks)-2* Sout3* Vc1* V(1,is,js,ks)
            V(5,il,jl,kl ) = V(5,is,js,ks)
            PP(il,jl,kl ) = PP(is,js,ks)
            T(il,jl,kl ) = T(is,js,ks)
!	       Den=V(1, I_R2 , J_R2 , K_R2)
!           Vxi=V(2, I_R2 , J_R2 , K_R2) / Den
!           Vyi=V(3, I_R2 , J_R2 , K_R2) / Den
!           Vzi=V(4, I_R2 , J_R2 , K_R2) / Den
!           Vc1  = Vxi*Sout1+Vyi*Sout2+Vzi*Sout3  !contravariant velocity
        
            is=i+iadd2;     js=j+jadd2;     ks=k+kadd2
            il=i+Eiadd2;    jl=j+Ejadd2;    kl=k+Ekadd2
            V(1,il,jl,kl ) = V(1,is,js,ks)
            V(2,il,jl,kl ) = V(2,is,js,ks)-2* Sout1* Vc1* V(1,is,js,ks)
            V(3,il,jl,kl ) = V(3,is,js,ks)-2* Sout2* Vc1* V(1,is,js,ks)
            V(4,il,jl,kl ) = V(4,is,js,ks)-2* Sout3* Vc1* V(1,is,js,ks)
            V(5,il,jl,kl ) = V(5,is,js,ks)
            PP(il,jl,kl ) = PP(is,js,ks)
            T(il,jl,kl ) = T(is,js,ks)

            if(IF_turb)then
                is=i+iadd;     js=j+jadd;     ks=k+kadd
                il=i+Eiadd;    jl=j+Ejadd;    kl=k+Ekadd
                V(6:7,il,jl,kl)=V(6:7,is,js,ks)
                Rmiu(il,jl,kl)=Rmiu(is,js,ks)

                is=i+iadd1;     js=j+jadd1;     ks=k+kadd1
                il=i+Eiadd1;    jl=j+Ejadd1;    kl=k+Ekadd1
                V(6:7,il,jl,kl)=V(6:7,is,js,ks)
                Rmiu(il,jl,kl)=Rmiu(is,js,ks)

                is=i+iadd2;     js=j+jadd2;     ks=k+kadd2
                il=i+Eiadd2;    jl=j+Ejadd2;    kl=k+Ekadd2
                V(6:7,il,jl,kl)=V(6:7,is,js,ks)
                Rmiu(il,jl,kl)=Rmiu(is,js,ks)
            endif
        enddo
        enddo
        enddo
    elseif(sweep==0)then
        do K=Kbgn,Kend
        do J=Jbgn,Jend
        do I=Ibgn,Iend
            is=i+iadd;     js=j+jadd;     ks=k+kadd
            il=i+Eiadd;    jl=j+Ejadd;    kl=k+Ekadd
            DqDxyz(1:18,il,jl,kl)=DqDxyz(1:18,is,js,ks)
        enddo
        enddo
        enddo
    endif
!!!!!!
!	IF ((sweep == 0).or.(sweep == 2)) THEN      !by ydd
!		VL(IJK, 1, I_VR0,J_VR0,K_VR0) = Vxi - Vc1*Sout1 !vx
!		VL(IJK, 2, I_VR0,J_VR0,K_VR0) = Vyi - Vc1*Sout2
!		VL(IJK, 3, I_VR0,J_VR0,K_VR0) = Vzi - Vc1*Sout3

!		VL(IJK, 4, I_VR0,J_VR0,K_VR0) = PP(I_R0,J_R0,K_R0)  !p
!		VL(IJK, 5, I_VR0,J_VR0,K_VR0) = T(I_R0,J_R0,K_R0)  !T

END SUBROUTINE BC_Symm



