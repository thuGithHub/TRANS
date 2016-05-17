SUBROUTINE BC_far(ibc,sweep)
    USE Global
    IMPLICIT NONE
    integer::ibc,sweep,LeftOrRight
    REAL::Sout1,Sout2,Sout3,Dnx,Dny,Dnz,out_n
    integer::Ibgn,Iend,Jbgn,Jend,Kbgn,Kend,i,j,k
    integer::is,js,ks,il,jl,kl,iadd,jadd,kadd,Eiadd,Ejadd,Ekadd, &
        &   Eiadd1,Ejadd1,Ekadd1,Eiadd2,Ejadd2,Ekadd2,Siadd,Sjadd,Skadd
    type(ConnectivityStruct),pointer :: Aconnect
    REAL:: Den,Vxi,Vyi,Vzi,P2,PP0,Den0,Vx0,Vy0,Vz0
    REAL:: dir,IJKDir
    real::TT0,TT,RmiuFar

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
        out_n=real(2*LeftOrRight-4)
 
        if(sweep==-1)then
            DO K=Kbgn,Kend
            DO J=Jbgn,Jend
            DO I=Ibgn,Iend
                
                is=i+Siadd; js=j+Sjadd; ks=k+Skadd
                Sout1=SD(IJKDir,1,is,js,ks)/Grad(IJKDir,is,js,ks)
                Sout2=SD(IJKDir,2,is,js,ks)/Grad(IJKDir,is,js,ks)
                Sout3=SD(IJKDir,3,is,js,ks)/Grad(IJKDir,is,js,ks)
                Dnx=Sout1*out_n
                Dny=Sout2*out_n
                Dnz=Sout3*out_n

                is=i+iadd;  js=j+jadd;  ks=k+kadd
                Den=V(1,is,js,ks)
                Vxi=V(2,is,js,ks)/Den
                Vyi=V(3,is,js,ks)/Den
                Vzi=V(4,is,js,ks)/Den
                P2=PP(is,js,ks)
                
                CALL BC_FARBOND(XM,P2,Den,Vxi,Vyi,Vzi,DNx,DNy,DNz,VxF,VyF,VzF,PP0,Den0,VX0,VY0,VZ0,dir)
!                TT0 = RXM2*PP0/Den0
                il=i+Eiadd; jl=j+Ejadd; kl=k+Ekadd
                V(1,il,jl,kl)=2*Den0-V(1,is,js,ks)
                den=V(1,il,jl,kl)
                V(2,il,jl,kl)=2.0*Vx0-Vxi
                V(3,il,jl,kl)=2.0*Vy0-Vyi
                V(4,il,jl,kl)=2.0*Vz0-Vzi
                PP(il,jl,kl)=2.0*PP0-P2
                T(il,jl,kl)=RXM2*PP(il,jl,kl)/den
                
                is=i+Eiadd1;   js=j+Ejadd1;   ks=k+Ekadd1    
!                V(1:5,il1,jl1,kl1,)=2.0*V(1:5,il,jl,kl)-V(1:5,is,js,ks)
!                PP(il1,jl1,kl1)=2.0*PP(il,jl,kl)-PP(is,js,ks)
!                T(il1,jl1,kl1)=2.0*T(il,jl,kl)-T(is,js,ks)
                V(1:5,is,js,ks)=V(1:5,il,jl,kl)
                PP(is,js,ks)=PP(il,jl,kl)
                T(is,js,ks)=T(il,jl,kl)
                
                is=i+Eiadd2;    js=j+Ejadd2;    ks=k+Ekadd2
!                V(1:5,is,js,ks)=2.0*V(1:5,il1,jl1,kl1)-V(1:5,il,jl,kl)
!                PP(is,js,ks)=2.0*PP(il1,jl1,kl1)-PP(il,jl,kl)
!                T(is,js,ks)=2.0*T(il1,jl1,kl1)-T(il,jl,kl)
                V(1:5,is,js,ks)=V(1:5,il,jl,kl)
                PP(is,js,ks)=PP(il,jl,kl)
                T(is,js,ks)=T(il,jl,kl)

                if(IF_turb)then
                    is=i+iadd;  js=j+jadd;  ks=k+kadd
                    il=i+Eiadd; jl=j+Ejadd; kl=k+Ekadd
                    if(dir>0.0)then
                        V(6:7,il,jl,kl)=V(6:7,is,js,ks)      
                        Rmiu(il,jl,kl)=Rmiu(is,js,ks)
                    else
                        TT=T(il,jl,kl)                   !by ydd 20151209
                        Den=V(1,il,jl,kl)
                        RmiuFar=(1.0+Csthlnd)/(TT*Csth1+Csthlnd)*(Csth1*TT)**1.5
                        Rmiu(il,jl,kl)=2.0*VisRatio*RmiuFar-Rmiu(is,js,ks)
                        V(6,il,jl,kl)=2.0*1.5*FSTI*FSTI*Den-V(6,is,js,ks)
                        V(7,il,jl,kl)=V(6,il,jl,kl)/(Rmiu(il,jl,kl)+tiny)*Den
                    endif
                        
                        is=i+Eiadd1;   js=j+Ejadd1;   ks=k+Ekadd1    
                        V(6:7,is,js,ks)=V(6:7,il,jl,kl)
                        Rmiu(is,js,ks)=Rmiu(il,jl,kl)

                        is=i+Eiadd2;   js=j+Ejadd2;   ks=k+Ekadd2   
                        V(6:7,is,js,ks)=V(6:7,il,jl,kl)
                        Rmiu(is,js,ks)=Rmiu(il,jl,kl)
                        !V(6,I_l0,J_l0,K_l0)=Rkf
                        !V(7,I_l0,J_l0,K_l0)=Rof
                        !Rmiu(I_l0,J_l0,K_l0)=Rmiuf
                        !V(6,I_l1,J_l1,K_l1)=Rkf
                        !V(7,I_l1,J_l1,K_l1)=Rof
                        !Rmiu(I_l1,J_l1,K_l1)=Rmiuf
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
                DqDxyz(1:18,il,jl,kl)=DqDxyz(1:18,is,js,ks)
            enddo
            enddo
            enddo
        endif

END SUBROUTINE BC_far



SUBROUTINE BC_FARBOND(XM,P2,RR,Vxi,Vyi,Vzi,DNx,DNy,DNz, VxF,VyF,VzF,PP0,RR0,VX0,VY0,VZ0,dir)
    IMPLICIT NONE
    REAL:: XM
    REAL:: P2,RR,Vxi,Vyi,Vzi,DNx,DNy,DNz
    REAL:: PP0,RR0,VX0,VY0,VZ0
    REAL:: VXF,VYF,VZF
    REAL:: dir
    real:: pGam
    REAL:: AF,AE,QNE,QNF,RF,RE,QN,AS,SON,QTX,QTY,QTZ

        pGam=1.4
        AF=1./XM
        !AE=SQRT(1.4*P2/RR)
        AE=SQRT(pGam*P2/RR)        !speed of sond
        QNE=VXi*DNX+VYi*DNY+VZi*DNZ
        QNF=VXF*DNX+VYF*DNY+VZF*DNZ
        RF=QNF-5.*AF              
        RE=QNE+5.*AE        
        QN=0.5*(RE+RF)
        AS=0.1*(RE-RF)
        IF(QN>-AS.AND.QN<0.)THEN
          !SON=1./1.4/(XM*XM)
          SON=1./pGam/(XM*XM)
          QTX=VXF-QNF*DNX
          QTY=VYF-QNF*DNY
          QTZ=VZF-QNF*DNZ
          VX0=QTX+QN*DNX
          VY0=QTY+QN*DNY
          VZ0=QTZ+QN*DNZ
       !   RR0=(AS*AS/SON/1.4)**2.5
          RR0=(AS*AS/SON/pGam)**2.5
        !  PP0= AS*AS*RR0/1.4
           PP0= AS*AS*RR0/pGam
          
          dir=-1.
        ENDIF

        IF(QN>=0..AND.QN<AS)THEN
         ! SON=P2/RR**1.4
          SON=P2/RR**pGam
          QTX=VXi-QNE*DNX
          QTY=VYi-QNE*DNY
          QTZ=VZi-QNE*DNZ
          VX0=QTX+QN*DNX
          VY0=QTY+QN*DNY
          VZ0=QTZ+QN*DNZ
         ! RR0=(AS*AS/SON/1.4)**2.5
         ! PP0= AS*AS*RR0/1.4
          RR0=(AS*AS/SON/pGam)**2.5
          PP0= AS*AS*RR0/pGam
          
          dir=1.
        ENDIF
        IF(QN<=-AS)THEN
            
          VX0=VXF
          VY0=VYF
          VZ0=VZF
          RR0=1.0
         ! PP0=1./1.4/XM/Xm
          PP0=1./pGam/XM/Xm
          
          dir=-1.
        ENDIF
        IF(QN>=AS)THEN
          VX0=VXi
          VY0=VYi
          VZ0=VZi
          RR0=RR
          PP0=P2
          dir=1.
        ENDIF
END SUBROUTINE BC_FARBOND




