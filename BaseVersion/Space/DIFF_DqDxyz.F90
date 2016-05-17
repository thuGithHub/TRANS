SUBROUTINE DIFF_DqDxyz
    USE Global
    IMPLICIT NONE
    
    INTEGER:: I,J,K,L,M,N
    INTEGER:: NI, NJ, NK
    INTEGER:: NI1, NJ1, NK1,NI2, NJ2, NK2
    REAL:: volinv, disfe, disfw, disfn, disfs, disft, disfb
   ! for adaptive scheme, added by dzw05,20130410
    REAL:: DuDx,DuDy,DuDz,DvDx,DvDy,DvDz,DwDx,DwDy,DwDz,W12,W13,W23,S11,S22,S33,S12,S13,S23
    REAL:: DkDx,DkDy,DkDz,DoDx,DoDy,DoDz
    REAL:: Vort,Skl,DivU,Den,RmiuT,RmiuL,TT,UUU,VVV,WWW,vK,vOmega,Vc,DLTxyz,Tau,Dist,dTKO,CDKO
    REAL:: arg1,func1,algrd,B0,G0,sok,alturb,A0 !F_limiter,Switchnum

    NI=ThisBlock%NI
    NJ=ThisBlock%NJ
    NK=ThisBlock%NK
    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1
    NI2=ThisBlock%NI2
    NJ2=ThisBlock%NJ2
    NK2=ThisBlock%NK2

      do k=1,NK1 !      do k=1,NK1  by xu
      do j=1,NJ1 !      do j=1,NJ1  by xu
      do i=1,NI1 !      do i=1,NI1  by xu

	   volinv=1.0_8/Vol(I,J,K)
	   !gradient(u)
	     disfe=VL(1,1,i+1,j,k) !disfe=uil(i+1,j,k)
         disfw=VL(1,1,i,j,k)   !disfw=uil(i,j,k)
         disfn=VL(2,1,i,j+1,k) !disfn=ujl(i,j+1,k)
         disfs=VL(2,1,i,j,k)   !disfs=ujl(i,j,k)
         disft=VL(3,1,i,j,k+1) !disft=ukl(i,j,k+1)
         disfb=VL(3,1,i,j,k)   !disfb=ukl(i,j,k)
	   DQdxyz(1,I,J,K)=volinv*(disfe*SD(1,1,i+1,j,k)-disfw*SD(1,1,i,j,k)+disfn*SD(2,1,i,j+1,k)-disfs*SD(2,1,i,j,k)+disft*SD(3,1,i,j,k+1)-disfb*SD(3,1,i,j,k))
	   DQdxyz(2,I,J,K)=volinv*(disfe*SD(1,2,i+1,J,K)-disfw*SD(1,2,I,J,K)+disfn*SD(2,2,i,j+1,k)-disfs*SD(2,2,i,j,k)+disft*SD(3,2,i,j,k+1)-disfb*SD(3,2,i,j,k))
	   DQdxyz(3,I,J,K)=volinv*(disfe*SD(1,3,I+1,J,K)-disfw*SD(1,3,I,J,K)+disfn*SD(2,3,i,j+1,k)-disfs*SD(2,3,i,j,k)+disft*SD(3,3,i,j,k+1)-disfb*SD(3,3,i,j,k))

	     disfe=VL(1,2,i+1,j,k) !disfe=uil(i+1,j,k)
         disfw=VL(1,2,i,j,k)   !disfw=uil(i,j,k)
         disfn=VL(2,2,i,j+1,k) !disfn=ujl(i,j+1,k)
         disfs=VL(2,2,i,j,k)   !disfs=ujl(i,j,k)
         disft=VL(3,2,i,j,k+1) !disft=ukl(i,j,k+1)
         disfb=VL(3,2,i,j,k)   !disfb=ukl(i,j,k)
	   DQdxyz(4,I,J,K)=volinv*(disfe*SD(1,1,i+1,j,k)-disfw*SD(1,1,i,j,k)+disfn*SD(2,1,i,j+1,k)-disfs*SD(2,1,i,j,k)+disft*SD(3,1,i,j,k+1)-disfb*SD(3,1,i,j,k))
	   DQdxyz(5,I,J,K)=volinv*(disfe*SD(1,2,I+1,J,K)-disfw*SD(1,2,I,J,K)+disfn*SD(2,2,i,j+1,k)-disfs*SD(2,2,i,j,k)+disft*SD(3,2,i,j,k+1)-disfb*SD(3,2,i,j,k))
	   DQdxyz(6,I,J,K)=volinv*(disfe*SD(1,3,I+1,J,K)-disfw*SD(1,3,I,J,K)+disfn*SD(2,3,i,j+1,k)-disfs*SD(2,3,i,j,k)+disft*SD(3,3,i,j,k+1)-disfb*SD(3,3,i,j,k))
	   !gradient(W)
	     disfe=VL(1,3,i+1,j,k) !disfe=uil(i+1,j,k)
         disfw=VL(1,3,i,j,k)   !disfw=uil(i,j,k)
         disfn=VL(2,3,i,j+1,k) !disfn=ujl(i,j+1,k)
         disfs=VL(2,3,i,j,k)   !disfs=ujl(i,j,k)
         disft=VL(3,3,i,j,k+1) !disft=ukl(i,j,k+1)
         disfb=VL(3,3,i,j,k)   !disfb=ukl(i,j,k)
	   DQdxyz(7,I,J,K)=volinv*(disfe*SD(1,1,i+1,j,k)-disfw*SD(1,1,i,j,k)+disfn*SD(2,1,i,j+1,k)-disfs*SD(2,1,i,j,k)+disft*SD(3,1,i,j,k+1)-disfb*SD(3,1,i,j,k))
	   DQdxyz(8,I,J,K)=volinv*(disfe*SD(1,2,I+1,J,K)-disfw*SD(1,2,I,J,K)+disfn*SD(2,2,i,j+1,k)-disfs*SD(2,2,i,j,k)+disft*SD(3,2,i,j,k+1)-disfb*SD(3,2,i,j,k))
	   DQdxyz(9,I,J,K)=volinv*(disfe*SD(1,3,I+1,J,K)-disfw*SD(1,3,I,J,K)+disfn*SD(2,3,i,j+1,k)-disfs*SD(2,3,i,j,k)+disft*SD(3,3,i,j,k+1)-disfb*SD(3,3,i,j,k))
	   !gradient(T)
	     disfe=VL(1,5,i+1,j,k) !disfe=uil(i+1,j,k)
         disfw=VL(1,5,i,j,k)   !disfw=uil(i,j,k)
         disfn=VL(2,5,i,j+1,k) !disfn=ujl(i,j+1,k)
         disfs=VL(2,5,i,j,k)   !disfs=ujl(i,j,k)
         disft=VL(3,5,i,j,k+1) !disft=ukl(i,j,k+1)
         disfb=VL(3,5,i,j,k)   !disfb=ukl(i,j,k)
	   DQdxyz(10,I,J,K)=volinv*(disfe*SD(1,1,i+1,j,k)-disfw*SD(1,1,i,j,k)+disfn*SD(2,1,i,j+1,k)-disfs*SD(2,1,i,j,k)+disft*SD(3,1,i,j,k+1)-disfb*SD(3,1,i,j,k))
	   DQdxyz(11,I,J,K)=volinv*(disfe*SD(1,2,I+1,J,K)-disfw*SD(1,2,I,J,K)+disfn*SD(2,2,i,j+1,k)-disfs*SD(2,2,i,j,k)+disft*SD(3,2,i,j,k+1)-disfb*SD(3,2,i,j,k))
	   DQdxyz(12,I,J,K)=volinv*(disfe*SD(1,3,I+1,J,K)-disfw*SD(1,3,I,J,K)+disfn*SD(2,3,i,j+1,k)-disfs*SD(2,3,i,j,k)+disft*SD(3,3,i,j,k+1)-disfb*SD(3,3,i,j,k))

!if(ThisBlock%ID_Present_Blk==2) write(*,*)i,j,k,DQDxyz(1:3,i,j,k)
	end do
	end do
	end do

! gradient in ghost cells -1 & NI/NJ/NK
    do k=1,NK1
    do j=1,NJ1
    do i=0,NI,NI

	   volinv=1.0_8/Vol(I,J,K)
	   !gradient(u)
	     disfe=VL(1,1,i+1,j,k) !disfe=uil(i+1,j,k)
         disfw=VL(1,1,i,j,k)   !disfw=uil(i,j,k)
         disfn=VL(2,1,i,j+1,k) !disfn=ujl(i,j+1,k)
         disfs=VL(2,1,i,j,k)   !disfs=ujl(i,j,k)
         disft=VL(3,1,i,j,k+1) !disft=ukl(i,j,k+1)
         disfb=VL(3,1,i,j,k)   !disfb=ukl(i,j,k)
	   DQdxyz(1,I,J,K)=volinv*(disfe*SD(1,1,i+1,j,k)-disfw*SD(1,1,i,j,k)+disfn*SD(2,1,i,j+1,k)-disfs*SD(2,1,i,j,k)+disft*SD(3,1,i,j,k+1)-disfb*SD(3,1,i,j,k))
	   DQdxyz(2,I,J,K)=volinv*(disfe*SD(1,2,i+1,J,K)-disfw*SD(1,2,I,J,K)+disfn*SD(2,2,i,j+1,k)-disfs*SD(2,2,i,j,k)+disft*SD(3,2,i,j,k+1)-disfb*SD(3,2,i,j,k))
	   DQdxyz(3,I,J,K)=volinv*(disfe*SD(1,3,I+1,J,K)-disfw*SD(1,3,I,J,K)+disfn*SD(2,3,i,j+1,k)-disfs*SD(2,3,i,j,k)+disft*SD(3,3,i,j,k+1)-disfb*SD(3,3,i,j,k))
        !gradient(V)
	     disfe=VL(1,2,i+1,j,k) !disfe=uil(i+1,j,k)
         disfw=VL(1,2,i,j,k)   !disfw=uil(i,j,k)
         disfn=VL(2,2,i,j+1,k) !disfn=ujl(i,j+1,k)
         disfs=VL(2,2,i,j,k)   !disfs=ujl(i,j,k)
         disft=VL(3,2,i,j,k+1) !disft=ukl(i,j,k+1)
         disfb=VL(3,2,i,j,k)   !disfb=ukl(i,j,k)
	   DQdxyz(4,I,J,K)=volinv*(disfe*SD(1,1,i+1,j,k)-disfw*SD(1,1,i,j,k)+disfn*SD(2,1,i,j+1,k)-disfs*SD(2,1,i,j,k)+disft*SD(3,1,i,j,k+1)-disfb*SD(3,1,i,j,k))
	   DQdxyz(5,I,J,K)=volinv*(disfe*SD(1,2,I+1,J,K)-disfw*SD(1,2,I,J,K)+disfn*SD(2,2,i,j+1,k)-disfs*SD(2,2,i,j,k)+disft*SD(3,2,i,j,k+1)-disfb*SD(3,2,i,j,k))
	   DQdxyz(6,I,J,K)=volinv*(disfe*SD(1,3,I+1,J,K)-disfw*SD(1,3,I,J,K)+disfn*SD(2,3,i,j+1,k)-disfs*SD(2,3,i,j,k)+disft*SD(3,3,i,j,k+1)-disfb*SD(3,3,i,j,k))
	   !gradient(W)
	     disfe=VL(1,3,i+1,j,k) !disfe=uil(i+1,j,k)
         disfw=VL(1,3,i,j,k)   !disfw=uil(i,j,k)
         disfn=VL(2,3,i,j+1,k) !disfn=ujl(i,j+1,k)
         disfs=VL(2,3,i,j,k)   !disfs=ujl(i,j,k)
         disft=VL(3,3,i,j,k+1) !disft=ukl(i,j,k+1)
         disfb=VL(3,3,i,j,k)   !disfb=ukl(i,j,k)
	   DQdxyz(7,I,J,K)=volinv*(disfe*SD(1,1,i+1,j,k)-disfw*SD(1,1,i,j,k)+disfn*SD(2,1,i,j+1,k)-disfs*SD(2,1,i,j,k)+disft*SD(3,1,i,j,k+1)-disfb*SD(3,1,i,j,k))
	   DQdxyz(8,I,J,K)=volinv*(disfe*SD(1,2,I+1,J,K)-disfw*SD(1,2,I,J,K)+disfn*SD(2,2,i,j+1,k)-disfs*SD(2,2,i,j,k)+disft*SD(3,2,i,j,k+1)-disfb*SD(3,2,i,j,k))
	   DQdxyz(9,I,J,K)=volinv*(disfe*SD(1,3,I+1,J,K)-disfw*SD(1,3,I,J,K)+disfn*SD(2,3,i,j+1,k)-disfs*SD(2,3,i,j,k)+disft*SD(3,3,i,j,k+1)-disfb*SD(3,3,i,j,k))
	   !gradient(T)
	     disfe=VL(1,5,i+1,j,k) !disfe=uil(i+1,j,k)
         disfw=VL(1,5,i,j,k)   !disfw=uil(i,j,k)
         disfn=VL(2,5,i,j+1,k) !disfn=ujl(i,j+1,k)
         disfs=VL(2,5,i,j,k)   !disfs=ujl(i,j,k)
         disft=VL(3,5,i,j,k+1) !disft=ukl(i,j,k+1)
         disfb=VL(3,5,i,j,k)   !disfb=ukl(i,j,k)
	   DQdxyz(10,I,J,K)=volinv*(disfe*SD(1,1,i+1,j,k)-disfw*SD(1,1,i,j,k)+disfn*SD(2,1,i,j+1,k)-disfs*SD(2,1,i,j,k)+disft*SD(3,1,i,j,k+1)-disfb*SD(3,1,i,j,k))
	   DQdxyz(11,I,J,K)=volinv*(disfe*SD(1,2,I+1,J,K)-disfw*SD(1,2,I,J,K)+disfn*SD(2,2,i,j+1,k)-disfs*SD(2,2,i,j,k)+disft*SD(3,2,i,j,k+1)-disfb*SD(3,2,i,j,k))
	   DQdxyz(12,I,J,K)=volinv*(disfe*SD(1,3,I+1,J,K)-disfw*SD(1,3,I,J,K)+disfn*SD(2,3,i,j+1,k)-disfs*SD(2,3,i,j,k)+disft*SD(3,3,i,j,k+1)-disfb*SD(3,3,i,j,k))

!if(ThisBlock%ID_Present_Blk==2) write(*,*)i,j,k,DQDxyz(1:3,i,j,k)
    enddo
    enddo
    enddo


    do k=1,NK1
    do j=0,NJ,NJ
    do i=1,NI1

	   volinv=1.0_8/Vol(I,J,K)
	   !gradient(u)
	     disfe=VL(1,1,i+1,j,k) !disfe=uil(i+1,j,k)
         disfw=VL(1,1,i,j,k)   !disfw=uil(i,j,k)
         disfn=VL(2,1,i,j+1,k) !disfn=ujl(i,j+1,k)
         disfs=VL(2,1,i,j,k)   !disfs=ujl(i,j,k)
         disft=VL(3,1,i,j,k+1) !disft=ukl(i,j,k+1)
         disfb=VL(3,1,i,j,k)   !disfb=ukl(i,j,k)
	   DQdxyz(1,I,J,K)=volinv*(disfe*SD(1,1,i+1,j,k)-disfw*SD(1,1,i,j,k)+disfn*SD(2,1,i,j+1,k)-disfs*SD(2,1,i,j,k)+disft*SD(3,1,i,j,k+1)-disfb*SD(3,1,i,j,k))
	   DQdxyz(2,I,J,K)=volinv*(disfe*SD(1,2,i+1,J,K)-disfw*SD(1,2,I,J,K)+disfn*SD(2,2,i,j+1,k)-disfs*SD(2,2,i,j,k)+disft*SD(3,2,i,j,k+1)-disfb*SD(3,2,i,j,k))
	   DQdxyz(3,I,J,K)=volinv*(disfe*SD(1,3,I+1,J,K)-disfw*SD(1,3,I,J,K)+disfn*SD(2,3,i,j+1,k)-disfs*SD(2,3,i,j,k)+disft*SD(3,3,i,j,k+1)-disfb*SD(3,3,i,j,k))
        !gradient(V)
	     disfe=VL(1,2,i+1,j,k) !disfe=uil(i+1,j,k)
         disfw=VL(1,2,i,j,k)   !disfw=uil(i,j,k)
         disfn=VL(2,2,i,j+1,k) !disfn=ujl(i,j+1,k)
         disfs=VL(2,2,i,j,k)   !disfs=ujl(i,j,k)
         disft=VL(3,2,i,j,k+1) !disft=ukl(i,j,k+1)
         disfb=VL(3,2,i,j,k)   !disfb=ukl(i,j,k)
	   DQdxyz(4,I,J,K)=volinv*(disfe*SD(1,1,i+1,j,k)-disfw*SD(1,1,i,j,k)+disfn*SD(2,1,i,j+1,k)-disfs*SD(2,1,i,j,k)+disft*SD(3,1,i,j,k+1)-disfb*SD(3,1,i,j,k))
	   DQdxyz(5,I,J,K)=volinv*(disfe*SD(1,2,I+1,J,K)-disfw*SD(1,2,I,J,K)+disfn*SD(2,2,i,j+1,k)-disfs*SD(2,2,i,j,k)+disft*SD(3,2,i,j,k+1)-disfb*SD(3,2,i,j,k))
	   DQdxyz(6,I,J,K)=volinv*(disfe*SD(1,3,I+1,J,K)-disfw*SD(1,3,I,J,K)+disfn*SD(2,3,i,j+1,k)-disfs*SD(2,3,i,j,k)+disft*SD(3,3,i,j,k+1)-disfb*SD(3,3,i,j,k))
	   !gradient(W)
	     disfe=VL(1,3,i+1,j,k) !disfe=uil(i+1,j,k)
         disfw=VL(1,3,i,j,k)   !disfw=uil(i,j,k)
         disfn=VL(2,3,i,j+1,k) !disfn=ujl(i,j+1,k)
         disfs=VL(2,3,i,j,k)   !disfs=ujl(i,j,k)
         disft=VL(3,3,i,j,k+1) !disft=ukl(i,j,k+1)
         disfb=VL(3,3,i,j,k)   !disfb=ukl(i,j,k)
	   DQdxyz(7,I,J,K)=volinv*(disfe*SD(1,1,i+1,j,k)-disfw*SD(1,1,i,j,k)+disfn*SD(2,1,i,j+1,k)-disfs*SD(2,1,i,j,k)+disft*SD(3,1,i,j,k+1)-disfb*SD(3,1,i,j,k))
	   DQdxyz(8,I,J,K)=volinv*(disfe*SD(1,2,I+1,J,K)-disfw*SD(1,2,I,J,K)+disfn*SD(2,2,i,j+1,k)-disfs*SD(2,2,i,j,k)+disft*SD(3,2,i,j,k+1)-disfb*SD(3,2,i,j,k))
	   DQdxyz(9,I,J,K)=volinv*(disfe*SD(1,3,I+1,J,K)-disfw*SD(1,3,I,J,K)+disfn*SD(2,3,i,j+1,k)-disfs*SD(2,3,i,j,k)+disft*SD(3,3,i,j,k+1)-disfb*SD(3,3,i,j,k))
	   !gradient(T)
	     disfe=VL(1,5,i+1,j,k) !disfe=uil(i+1,j,k)
         disfw=VL(1,5,i,j,k)   !disfw=uil(i,j,k)
         disfn=VL(2,5,i,j+1,k) !disfn=ujl(i,j+1,k)
         disfs=VL(2,5,i,j,k)   !disfs=ujl(i,j,k)
         disft=VL(3,5,i,j,k+1) !disft=ukl(i,j,k+1)
         disfb=VL(3,5,i,j,k)   !disfb=ukl(i,j,k)
	   DQdxyz(10,I,J,K)=volinv*(disfe*SD(1,1,i+1,j,k)-disfw*SD(1,1,i,j,k)+disfn*SD(2,1,i,j+1,k)-disfs*SD(2,1,i,j,k)+disft*SD(3,1,i,j,k+1)-disfb*SD(3,1,i,j,k))
	   DQdxyz(11,I,J,K)=volinv*(disfe*SD(1,2,I+1,J,K)-disfw*SD(1,2,I,J,K)+disfn*SD(2,2,i,j+1,k)-disfs*SD(2,2,i,j,k)+disft*SD(3,2,i,j,k+1)-disfb*SD(3,2,i,j,k))
	   DQdxyz(12,I,J,K)=volinv*(disfe*SD(1,3,I+1,J,K)-disfw*SD(1,3,I,J,K)+disfn*SD(2,3,i,j+1,k)-disfs*SD(2,3,i,j,k)+disft*SD(3,3,i,j,k+1)-disfb*SD(3,3,i,j,k))

!if(ThisBlock%ID_Present_Blk==2) write(*,*)i,j,k,DQDxyz(1:3,i,j,k)
    enddo
    enddo
    enddo


    do k=0,NK,NK
    do j=1,NJ1
    do i=1,NI1

	   volinv=1.0_8/Vol(I,J,K)
	   !gradient(u)
	     disfe=VL(1,1,i+1,j,k) !disfe=uil(i+1,j,k)
         disfw=VL(1,1,i,j,k)   !disfw=uil(i,j,k)
         disfn=VL(2,1,i,j+1,k) !disfn=ujl(i,j+1,k)
         disfs=VL(2,1,i,j,k)   !disfs=ujl(i,j,k)
         disft=VL(3,1,i,j,k+1) !disft=ukl(i,j,k+1)
         disfb=VL(3,1,i,j,k)   !disfb=ukl(i,j,k)
	   DQdxyz(1,I,J,K)=volinv*(disfe*SD(1,1,i+1,j,k)-disfw*SD(1,1,i,j,k)+disfn*SD(2,1,i,j+1,k)-disfs*SD(2,1,i,j,k)+disft*SD(3,1,i,j,k+1)-disfb*SD(3,1,i,j,k))
	   DQdxyz(2,I,J,K)=volinv*(disfe*SD(1,2,i+1,J,K)-disfw*SD(1,2,I,J,K)+disfn*SD(2,2,i,j+1,k)-disfs*SD(2,2,i,j,k)+disft*SD(3,2,i,j,k+1)-disfb*SD(3,2,i,j,k))
	   DQdxyz(3,I,J,K)=volinv*(disfe*SD(1,3,I+1,J,K)-disfw*SD(1,3,I,J,K)+disfn*SD(2,3,i,j+1,k)-disfs*SD(2,3,i,j,k)+disft*SD(3,3,i,j,k+1)-disfb*SD(3,3,i,j,k))
        !gradient(V)
	     disfe=VL(1,2,i+1,j,k) !disfe=uil(i+1,j,k)
         disfw=VL(1,2,i,j,k)   !disfw=uil(i,j,k)
         disfn=VL(2,2,i,j+1,k) !disfn=ujl(i,j+1,k)
         disfs=VL(2,2,i,j,k)   !disfs=ujl(i,j,k)
         disft=VL(3,2,i,j,k+1) !disft=ukl(i,j,k+1)
         disfb=VL(3,2,i,j,k)   !disfb=ukl(i,j,k)
	   DQdxyz(4,I,J,K)=volinv*(disfe*SD(1,1,i+1,j,k)-disfw*SD(1,1,i,j,k)+disfn*SD(2,1,i,j+1,k)-disfs*SD(2,1,i,j,k)+disft*SD(3,1,i,j,k+1)-disfb*SD(3,1,i,j,k))
	   DQdxyz(5,I,J,K)=volinv*(disfe*SD(1,2,I+1,J,K)-disfw*SD(1,2,I,J,K)+disfn*SD(2,2,i,j+1,k)-disfs*SD(2,2,i,j,k)+disft*SD(3,2,i,j,k+1)-disfb*SD(3,2,i,j,k))
	   DQdxyz(6,I,J,K)=volinv*(disfe*SD(1,3,I+1,J,K)-disfw*SD(1,3,I,J,K)+disfn*SD(2,3,i,j+1,k)-disfs*SD(2,3,i,j,k)+disft*SD(3,3,i,j,k+1)-disfb*SD(3,3,i,j,k))
	   !gradient(W)
	     disfe=VL(1,3,i+1,j,k) !disfe=uil(i+1,j,k)
         disfw=VL(1,3,i,j,k)   !disfw=uil(i,j,k)
         disfn=VL(2,3,i,j+1,k) !disfn=ujl(i,j+1,k)
         disfs=VL(2,3,i,j,k)   !disfs=ujl(i,j,k)
         disft=VL(3,3,i,j,k+1) !disft=ukl(i,j,k+1)
         disfb=VL(3,3,i,j,k)   !disfb=ukl(i,j,k)
	   DQdxyz(7,I,J,K)=volinv*(disfe*SD(1,1,i+1,j,k)-disfw*SD(1,1,i,j,k)+disfn*SD(2,1,i,j+1,k)-disfs*SD(2,1,i,j,k)+disft*SD(3,1,i,j,k+1)-disfb*SD(3,1,i,j,k))
	   DQdxyz(8,I,J,K)=volinv*(disfe*SD(1,2,I+1,J,K)-disfw*SD(1,2,I,J,K)+disfn*SD(2,2,i,j+1,k)-disfs*SD(2,2,i,j,k)+disft*SD(3,2,i,j,k+1)-disfb*SD(3,2,i,j,k))
	   DQdxyz(9,I,J,K)=volinv*(disfe*SD(1,3,I+1,J,K)-disfw*SD(1,3,I,J,K)+disfn*SD(2,3,i,j+1,k)-disfs*SD(2,3,i,j,k)+disft*SD(3,3,i,j,k+1)-disfb*SD(3,3,i,j,k))
	   !gradient(T)
	     disfe=VL(1,5,i+1,j,k) !disfe=uil(i+1,j,k)
         disfw=VL(1,5,i,j,k)   !disfw=uil(i,j,k)
         disfn=VL(2,5,i,j+1,k) !disfn=ujl(i,j+1,k)
         disfs=VL(2,5,i,j,k)   !disfs=ujl(i,j,k)
         disft=VL(3,5,i,j,k+1) !disft=ukl(i,j,k+1)
         disfb=VL(3,5,i,j,k)   !disfb=ukl(i,j,k)
	   DQdxyz(10,I,J,K)=volinv*(disfe*SD(1,1,i+1,j,k)-disfw*SD(1,1,i,j,k)+disfn*SD(2,1,i,j+1,k)-disfs*SD(2,1,i,j,k)+disft*SD(3,1,i,j,k+1)-disfb*SD(3,1,i,j,k))
	   DQdxyz(11,I,J,K)=volinv*(disfe*SD(1,2,I+1,J,K)-disfw*SD(1,2,I,J,K)+disfn*SD(2,2,i,j+1,k)-disfs*SD(2,2,i,j,k)+disft*SD(3,2,i,j,k+1)-disfb*SD(3,2,i,j,k))
	   DQdxyz(12,I,J,K)=volinv*(disfe*SD(1,3,I+1,J,K)-disfw*SD(1,3,I,J,K)+disfn*SD(2,3,i,j+1,k)-disfs*SD(2,3,i,j,k)+disft*SD(3,3,i,j,k+1)-disfb*SD(3,3,i,j,k))

!if(ThisBlock%ID_Present_Blk==2) write(*,*)i,j,k,DQDxyz(1:3,i,j,k)
    enddo
    enddo
    enddo

	if (IF_turb) then

        do k=1,NK1 !      do k=1,NK1  by xu
        do j=1,NJ1 !      do j=1,NJ1  by xu
        do i=1,NI1 !      do i=1,NI1  by xu

	     volinv=1.0_8/Vol(I,J,K)
           ! gradient(TKO)
	       disfe=VL(1,6,i+1,j,k) !disfe=teil(i+1,j,k)
           disfw=VL(1,6,i,j,k)   !disfw=teil(i,j,k)
           disfn=VL(2,6,i,j+1,k) !disfn=tejl(i,j+1,k)
           disfs=VL(2,6,i,j,k)   !disfs=tejl(i,j,k)
           disft=VL(3,6,i,j,k+1) !disft=tekl(i,j,k+1)
           disfb=VL(3,6,i,j,k)   !disfb=tekl(i,j,k)

	   DQdxyz(13,I,J,K)=volinv*(disfe*SD(1,1,i+1,j,k)-disfw*SD(1,1,i,j,k)+disfn*SD(2,1,i,j+1,k)-disfs*SD(2,1,i,j,k)+disft*SD(3,1,i,j,k+1)-disfb*SD(3,1,i,j,k))
	   DQdxyz(14,I,J,K)=volinv*(disfe*SD(1,2,i+1,J,K)-disfw*SD(1,2,I,J,K)+disfn*SD(2,2,i,j+1,k)-disfs*SD(2,2,i,j,k)+disft*SD(3,2,i,j,k+1)-disfb*SD(3,2,i,j,k))
	   DQdxyz(15,I,J,K)=volinv*(disfe*SD(1,3,I+1,J,K)-disfw*SD(1,3,I,J,K)+disfn*SD(2,3,i,j+1,k)-disfs*SD(2,3,i,j,k)+disft*SD(3,3,i,j,k+1)-disfb*SD(3,3,i,j,k))
           ! gradient(TOO)
	       disfe=VL(1,7,i+1,j,k) !disfe=teil(i+1,j,k)
           disfw=VL(1,7,i,j,k)   !disfw=teil(i,j,k)
           disfn=VL(2,7,i,j+1,k) !disfn=tejl(i,j+1,k)
           disfs=VL(2,7,i,j,k)   !disfs=tejl(i,j,k)
           disft=VL(3,7,i,j,k+1) !disft=tekl(i,j,k+1)
           disfb=VL(3,7,i,j,k)   !disfb=tekl(i,j,k)
!dtedx=volinv*(disfe*aix(i+1,j,k)-disfw*aix(i,j,k)+disfn*ajx(i,j+1,k)-disfs*ajx(i,j,k)+disft*akx(i,j,k+1)-disfb*akx(i,j,k))
	   DQdxyz(16,I,J,K)=volinv*(disfe*SD(1,1,i+1,j,k)-disfw*SD(1,1,i,j,k)+disfn*SD(2,1,i,j+1,k)-disfs*SD(2,1,i,j,k)+disft*SD(3,1,i,j,k+1)-disfb*SD(3,1,i,j,k))
	   DQdxyz(17,I,J,K)=volinv*(disfe*SD(1,2,i+1,J,K)-disfw*SD(1,2,I,J,K)+disfn*SD(2,2,i,j+1,k)-disfs*SD(2,2,i,j,k)+disft*SD(3,2,i,j,k+1)-disfb*SD(3,2,i,j,k))
	   DQdxyz(18,I,J,K)=volinv*(disfe*SD(1,3,I+1,J,K)-disfw*SD(1,3,I,J,K)+disfn*SD(2,3,i,j+1,k)-disfs*SD(2,3,i,j,k)+disft*SD(3,3,i,j,k+1)-disfb*SD(3,3,i,j,k))

	  end do
	  end do
	  end do
    endif
! gradient in ghost cells, -1,NI/NJ/NK, by ydd
!        do k=1,NK1 !      
!        do j=1,NJ1 !      
!        do i=0,NI,NI
!
!	     volinv=1.0_8/Vol(I,J,K)
!           ! gradient(TKO)
!	       disfe=VL(1,6,i+1,j,k) !disfe=teil(i+1,j,k)
!           disfw=VL(1,6,i,j,k)   !disfw=teil(i,j,k)
!           disfn=VL(2,6,i,j+1,k) !disfn=tejl(i,j+1,k)
!           disfs=VL(2,6,i,j,k)   !disfs=tejl(i,j,k)
!           disft=VL(3,6,i,j,k+1) !disft=tekl(i,j,k+1)
!           disfb=VL(3,6,i,j,k)   !disfb=tekl(i,j,k)
!
!	   DQdxyz(13,I,J,K)=volinv*(disfe*SD(1,1,i+1,j,k)-disfw*SD(1,1,i,j,k)+disfn*SD(2,1,i,j+1,k)-disfs*SD(2,1,i,j,k)+disft*SD(3,1,i,j,k+1)-disfb*SD(3,1,i,j,k))
!	   DQdxyz(14,I,J,K)=volinv*(disfe*SD(1,2,i+1,J,K)-disfw*SD(1,2,I,J,K)+disfn*SD(2,2,i,j+1,k)-disfs*SD(2,2,i,j,k)+disft*SD(3,2,i,j,k+1)-disfb*SD(3,2,i,j,k))
!	   DQdxyz(15,I,J,K)=volinv*(disfe*SD(1,3,I+1,J,K)-disfw*SD(1,3,I,J,K)+disfn*SD(2,3,i,j+1,k)-disfs*SD(2,3,i,j,k)+disft*SD(3,3,i,j,k+1)-disfb*SD(3,3,i,j,k))
!           ! gradient(TOO)
!	       disfe=VL(1,7,i+1,j,k) !disfe=teil(i+1,j,k)
!           disfw=VL(1,7,i,j,k)   !disfw=teil(i,j,k)
!           disfn=VL(2,7,i,j+1,k) !disfn=tejl(i,j+1,k)
!           disfs=VL(2,7,i,j,k)   !disfs=tejl(i,j,k)
!           disft=VL(3,7,i,j,k+1) !disft=tekl(i,j,k+1)
!           disfb=VL(3,7,i,j,k)   !disfb=tekl(i,j,k)
!!dtedx=volinv*(disfe*aix(i+1,j,k)-disfw*aix(i,j,k)+disfn*ajx(i,j+1,k)-disfs*ajx(i,j,k)+disft*akx(i,j,k+1)-disfb*akx(i,j,k))
!	   DQdxyz(16,I,J,K)=volinv*(disfe*SD(1,1,i+1,j,k)-disfw*SD(1,1,i,j,k)+disfn*SD(2,1,i,j+1,k)-disfs*SD(2,1,i,j,k)+disft*SD(3,1,i,j,k+1)-disfb*SD(3,1,i,j,k))
!	   DQdxyz(17,I,J,K)=volinv*(disfe*SD(1,2,i+1,J,K)-disfw*SD(1,2,I,J,K)+disfn*SD(2,2,i,j+1,k)-disfs*SD(2,2,i,j,k)+disft*SD(3,2,i,j,k+1)-disfb*SD(3,2,i,j,k))
!	   DQdxyz(18,I,J,K)=volinv*(disfe*SD(1,3,I+1,J,K)-disfw*SD(1,3,I,J,K)+disfn*SD(2,3,i,j+1,k)-disfs*SD(2,3,i,j,k)+disft*SD(3,3,i,j,k+1)-disfb*SD(3,3,i,j,k))

!if(ThisBlock%ID_Present_Blk==2) write(*,*)i,j,k,"Tur",DQDxyz(13,i,j,k),DQDxyz(16,i,j,k)
!	  end do
!	  end do
!	  end do

       DO K=1,NK1
	   DO J=1,NJ1
	   DO I=1,NI1 
         
	    DuDx= dQdxyz( 1,I,J,K)
	    DuDy= dQdxyz( 2,I,J,K)
	    DuDz= dQdxyz( 3,I,J,K)
		DvDx= dQdxyz( 4,I,J,K)
	    DvDy= dQdxyz( 5,I,J,K)
	    DvDz= dQdxyz( 6,I,J,K)
	    DwDx= dQdxyz( 7,I,J,K)
	    DwDy= dQdxyz( 8,I,J,K)
	    DwDz= dQdxyz( 9,I,J,K)

         W12 =(Dudy-Dvdx)/2.
	     W13 =(Dudz-Dwdx)/2.
	     W23 =(Dvdz-Dwdy)/2.
	     S11 =(Dudx+Dudx)/2.
	     S22 =(Dvdy+Dvdy)/2.
	     S33 =(Dwdz+Dwdz)/2.
	     S12 =(Dudy+Dvdx)/2.
	     S13 =(Dudz+Dwdx)/2.
	     S23 =(Dwdy+Dvdz)/2.

           Vort=SQRT(2.*(2.*W12*W12+2.*W13*W13+2.*W23*W23))
           Skl =SQRT(2.*(1.*S11*S11+1.*S22*S22+1.*S33*S33 &
                     & +2.*S12*S12+2.*S13*S13+2.*S23*S23))
	     
		 DivU=sqrt(DuDx*DuDx+DvDy*DvDy+DwDz*DwDz)


		 Den   = V(1,I,J,K)
	     RmiuT =Rmiu(I,J,K)
         TT    =   T(I,J,K)
         RMIUL = (1.+Csthlnd)/(TT*Csth1+Csthlnd)*(Csth1*TT)**1.5


	     UUU   =V(2,I,J,K)/Den
	     VVV   =V(3,I,J,K)/Den
	     WWW   =V(4,I,J,K)/Den
	     vK    =V(6,I,J,K)/Den
	     vOmega=V(7,I,J,K)/Den
	     
         Vc = (UUU*UUU+VVV*VVV+WWW*WWW)**0.5_8
		 DLTxyz= aLeng(I,J,K)
         tau = 1.0_8


         dkdx=DQdxyz(13,I,J,K)
		 dkdy=DQdxyz(14,I,J,K)
		 dkdz=DQdxyz(15,I,J,K)
		 dodx=DQdxyz(16,I,J,K)
		 dody=DQdxyz(17,I,J,K)
		 dodz=DQdxyz(18,I,J,K)
	     Dist=Dst(I,J,K)




	    dTKO=dkdx*dodx+dkdy*dody+dkdz*dodz
	    CDKO=max(2.0_8*DEN*sigo2/vomega*(dTKO),1.0d-10)

		arg1=max(sqrt(vK)/Beta_star/vOmega/Dist/Ref,&
                 & 500.0_8*RMIUL/DEN/Dist/Dist/vOmega/Ref/Ref)

           arg1=min(arg1,4.0_8*DEN*sigo2*vK/CDKO/Dist/Dist)

           FunC1=tanh(arg1**4)

	       aLgrd=((1.-FunC1)*Cdes_epsi+FunC1*Cdes_omeg)*DLTxyz
	     
           B0=CH3*Vort*max(Vort,Skl)/max((Vort**2.+Skl**2.)/2.,alscale*0.1)!near shock,S & Omega is large
      
	      g0=tanh(B0*B0*B0*B0)
	     
		 sok=max(((Vort**2.+Skl**2.)/2.)**0.5,0.1/tau)  !0.1/tau
         alturb=((RmiuT+RMIUL)/Den/(Ref*sok*cmu**1.5))**0.5
	    
         A0=CH2*max(((aLgrd/alturb/g0)-0.1),0.0) 
!	     Fdes_rans(I,J,K)=aLgrd/alturb
      Fscheme(I,J,K)= max(1.*tanh(A0**CH1),F_limiter1)
	  shock(i,j,k) =  DivU*DivU/max(DivU*DivU+Vort*Vort,1._8)
    enddo
	enddo
	enddo

    
   
   do L=1,3
    shock(1-L,:,:) = shock(1,:,:)   
    shock(:,1-L,:) = shock(:,1,:)   
    shock(:,:,1-L) = shock(:,:,1)

    shock(NI1+L,:,:) = shock(NI1,:,:)   
    shock(:,Nj1+L,:) = shock(:,NJ1,:)   
    shock(:,:,NK1+L) = shock(:,:,NK1)
  enddo


    Fscheme( 0,:,:) = Fscheme(1,:,:)   
    Fscheme(:, 0,:) = Fscheme(:,1,:)   
    Fscheme(:,:, 0) = Fscheme(:,:,1)

    Fscheme( NI,:,:) = Fscheme(NI1,:,:)   
    Fscheme(:, NJ,:) = Fscheme(:,NJ1,:)   
    Fscheme(:,:, NK) = Fscheme(:,:,NK1)


END SUBROUTINE DIFF_DqDxyz



