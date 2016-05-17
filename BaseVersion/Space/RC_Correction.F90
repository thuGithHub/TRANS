SUBROUTINE RCCorrection
    USE Global
    IMPLICIT NONE
    
    INTEGER:: I,J,K,L,M,N
    INTEGER:: NI1, NJ1, NK1
    REAL:: volinv
    REAL:: DuDx,DuDy,DuDz,DvDx,DvDy,DvDz,DwDx,DwDy,DwDz,W12,W13,W23,S11,S22,S33,S12,S13,S23
    real::St2(3,3),Wij(3,3)
    integer::II,JJ,KK,MM,NN,LL1,LL2
    real::DsijDt,Sijf1,Sijf2,Sijf3,Sijf4,Sijf5,Sijf6,Sij,vol1
    real::cr1,cr2,cr3,f_rot,DsDt,DsDt2,Ssum,Dsum,Wsum,r_bar,r_star,epsIMN,epsJMN

    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1

!rotate correction, by ydd 20151223
    if(Kind_RotMod==1)then
        cr1=1.0
        cr2=2.0
        cr3=1.0

    do k=1,Nk1
    do j=1,NJ1
    do i=1,NI1

        DuDx= dQdxyz( 1,I,J,K)
        DuDy= dQdxyz( 2,I,J,K)
        DuDz= dQdxyz( 3,I,J,K)
        DvDx= dQdxyz( 4,I,J,K)
        DvDy= dQdxyz( 5,I,J,K)
        DvDz= dQdxyz( 6,I,J,K)
        DwDx= dQdxyz( 7,I,J,K)
        DwDy= dQdxyz( 8,I,J,K)
        DwDz= dQdxyz( 9,I,J,K)
            

        St2(1,1)=dudx
        St2(1,2)=0.5*(dudy+dvdx)
        St2(1,3)=0.5*(dudz+dwdx)
        St2(2,1)=St2(1,2)
        St2(2,2)=dvdy
        St2(2,3)=0.5*(dvdz+dwdy)
        St2(3,1)=St2(1,3)
        St2(3,2)=St2(2,3)

        Wij(1,2)=0.5*(dudy-dvdx)
        Wij(1,3)=0.5*(dudz-dwdx)
        Wij(2,3)=0.5*(dvdz-dwdy)-omega(1)
        Wij(2,1)=-Wij(1,2)
        Wij(3,1)=-Wij(1,3)
        Wij(3,2)=-Wij(2,3)
        SSum=0.0
        WSum=0.0
        do JJ=1,3
        do II=1,3
            SSum=SSum+St2(II,JJ)*St2(II,JJ)
            Wsum=Wsum+Wij(II,JJ)*Wij(II,JJ)
        enddo
        enddo
        SSum=SSum*2.0
        Wsum=Wsum*2.0
        r_bar=0.0
!if(ThisBlock%ID_Present_Blk==0) write(*,*)"begin DsijDt ok",i,j,k
        do JJ=1,3
        do II=1,3
            LL1=(II-1)*3+JJ
            LL2=(JJ-1)*3+II

            Sij=DqDxyz(LL1,i,j,k)+DqDxyz(LL2,i,j,k)

            vol1=vol(i-1,j,k)/(vol(i-1,j,k)+vol(i,j,k)) ! i  
            Sijf1=vol1*(DqDxyz(LL1,i-1,j,k)+DqDxyz(LL2,i-1,j,k))+(1.0-vol1)*Sij    !(DqDxyz(LL1,i,j,k)+DqDxyz(LL2,i,j,k))

            vol1=vol(i+1,j,k)/(vol(i+1,j,k)+vol(i,j,k)) !i+1
            Sijf2=vol1*(DqDxyz(LL1,i+1,j,k)+DqDxyz(LL2,i+1,j,k))+(1.0-vol1)*Sij     !(DqDxyz(LL1,i,j,k)+DqDxyz(LL2,i,j,k))

            vol1=vol(i,j-1,k)/(vol(i,j-1,k)+vol(i,j,k)) !j
            Sijf3=vol1*(DqDxyz(LL1,i,j-1,k)+DqDxyz(LL2,i,j-1,k))+(1.0-vol1)*Sij     !(DqDxyz(LL1,i,j,k)+DqDxyz(LL2,i,j,k))

            vol1=vol(i,j+1,k)/(vol(i,j+1,k)+vol(i,j,k)) !j+1
            Sijf4=vol1*(DqDxyz(LL1,i,j+1,k)+DqDxyz(LL2,i,j+1,k))+(1.0-vol1)*Sij

            vol1=Vol(i,j,k-1)/(vol(i,j,k-1)+vol(i,j,k)) !k
            Sijf5=vol1*(DqDxyz(LL1,i,j,k-1)+DqDxyz(LL2,i,j,k-1))+(1.0-vol1)*Sij

            vol1=Vol(i,j,k+1)/(vol(i,j,k+1)+vol(i,j,k)) !k+1
            Sijf6=vol1*(DqDxyz(LL1,i,j,k+1)+DqDxyz(LL2,i,j,k+1))+(1.0-vol1)*Sij

        volinv=1.0/Vol(I,J,K)
        DsijDt=(VL(1,1,i+1,j,k)*SD(1,1,i+1,j,k)+VL(2,1,i+1,j,k)*SD(2,1,i+1,j,k)+VL(3,1,i+1,j,k)*SD(3,1,i+1,j,k))*Sijf2- &
        &      (VL(1,1,i,j,k)*SD(1,1,i,j,k)+VL(2,1,i,j,k)*SD(2,1,i,j,k)+VL(3,1,i,j,k)*SD(3,1,i,j,k))*Sijf1+    &
        &      (VL(1,2,i,j+1,k)*SD(1,2,i,j+1,k)+VL(2,2,i,j+1,k)*SD(2,2,i,j+1,k)+VL(3,2,i,j+1,k)*SD(3,2,i,j+1,k))*Sijf4- &
        &      (VL(1,2,i,j,k)*SD(1,2,i,j,k)+VL(2,2,i,j,k)*SD(2,2,i,j,k)+VL(3,2,i,j,k)*SD(3,2,i,j,k))*Sijf3+    &
        &      (VL(1,3,i,j,k+1)*SD(1,3,i,j,k+1)+VL(2,3,i,j,k+1)*SD(2,3,i,j,k+1)+VL(3,3,i,j,k+1)*SD(3,3,i,j,k+1))*Sijf6- &
        &      (VL(1,3,i,j,k)*SD(1,3,i,j,k)+VL(2,3,i,j,k)*SD(2,3,i,j,k)+VL(3,3,i,j,k)*SD(3,3,i,j,k))*Sijf5

            DsDt=DsijDt*0.5*volinv
!if(ThisBlock%ID_Present_Blk==0) write(*,*)"DsijDt ok",i,j,k

            DsDt2=0.0
            do MM=1,1
            do NN=1,3
                epsIMN=0.5*(II-MM)*(MM-NN)*(NN-II)
                epsJMN=0.5*(JJ-MM)*(MM-NN)*(NN-JJ)
                DsDt2=DsDt2+(epsIMN*St2(JJ,NN)+epsJMN*St2(II,NN))*omega(M)
            enddo
            enddo
            DsDt=DsDt+DsDt2
!            DsDt=DsDt2
            do kk=1,3
                r_bar=r_bar+Wij(II,KK)*St2(JJ,KK)*DsDt
            enddo
        enddo
        enddo
        Dsum=max(SSum,0.09*(V(7,i,j,k)/V(1,i,j,k))**2.0)
        r_star=sqrt(SSum/WSum)
        r_bar=2.0*r_bar/(sqrt(WSum)*Dsum**1.5)
        f_rot=(1+cr1)*2.0*r_star/(1.0+r_star)*(1-cr3*atan(cr2*r_bar))-cr1
        f_r1(i,j,k)=max(Min(f_rot,1.25),0.0)
    enddo
    enddo
    enddo

    else if(Kind_RotMod==2)then
    do k=1,NK1
    do j=1,NJ1
    do i=1,NI1
        DuDx= dQdxyz( 1,I,J,K)
        DuDy= dQdxyz( 2,I,J,K)
        DuDz= dQdxyz( 3,I,J,K)
        DvDx= dQdxyz( 4,I,J,K)
        DvDy= dQdxyz( 5,I,J,K)
        DvDz= dQdxyz( 6,I,J,K)
        DwDx= dQdxyz( 7,I,J,K)
        DwDy= dQdxyz( 8,I,J,K)
        DwDz= dQdxyz( 9,I,J,K)
        St2(1,1)=dudx
        St2(1,2)=0.5*(dudy+dvdx)
        St2(1,3)=0.5*(dudz+dwdx)
        St2(2,1)=St2(1,2)
        St2(2,2)=dvdy
        St2(2,3)=0.5*(dvdz+dwdy)
        St2(3,1)=St2(1,3)
        St2(3,2)=St2(2,3)

        Wij(1,2)=0.5*(dudy-dvdx)
        Wij(1,3)=0.5*(dudz-dwdx)
        Wij(2,3)=0.5*(dvdz-dwdy)
        Wij(2,1)=-Wij(1,2)
        Wij(3,1)=-Wij(1,3)
        Wij(3,2)=-Wij(2,3)
        SSum=0.0
        WSum=0.0
        do JJ=1,3
        do II=1,3
            SSum=SSum+St2(II,JJ)*St2(II,JJ)
            Wsum=Wsum+Wij(II,JJ)*Wij(II,JJ)
        enddo
        enddo
        SSum=SSum*2.0
        Wsum=Wsum*2.0
        cr1=sqrt(Wsum/SSum)
        cr2=cr1*(cr1-1.0)
        cr3=1.4
        f_r2(i,j,k)=1.0/(1.0+cr2*cr3)
    enddo
    enddo
    enddo
    endif


END SUBROUTINE RCCorrection 



