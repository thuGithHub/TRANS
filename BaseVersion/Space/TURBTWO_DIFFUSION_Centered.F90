SUBROUTINE TURBTWO_DIFFUSION_Centered
    USE Global
    IMPLICIT NONE
    
    INTEGER:: I,J,K,L,Lxyz,LL
    INTEGER:: I1,J1,K1
    INTEGER:: NI, NJ, NK
    INTEGER:: NI1, NJ1, NK1

    INTEGER:: Iswall0, Iswallm
    REAL:: NMD=-2./3.

    REAL:: alsinv, all, alr
    REAL:: disfpp,disfnn,disfp,disfn

    REAL:: dTkdx1,dTkdy1,dTkdz1
    REAL:: dTkdx,dTkdy,dTkdz
    REAL:: dTkdxL,dTkdyL,dTkdzL

    REAL:: dTOdx1,dTOdy1,dTOdz1
    REAL:: dTOdx,dTOdy,dTOdz
    REAL:: dTOdxL,dTOdyL,dTOdzL

    REAL:: TTL,RmiuL,RmiuT,RmiuT1,RmiuTL

    REAL:: DIFFTki,DIFFToi
    REAL:: DIFFTkj,DIFFToj
    REAL:: DIFFTkk,DIFFTok

    REAL:: SigkAB, SigoAB
    REAL:: DEN0,TKK0,TOO0
    REAL:: DEN1,TKK1,TOO1
    REAL:: Den,TKK,TOO
    REAL:: dTK0,CDK0
    REAL:: distfc,arg1,fc1c
    REAL:: vmute,vmuoq
    
    REAL:: TT,TT1,dTKO,CDKO
    REAL:: Dx,Dy,Dz

   
    NI=ThisBlock%NI
    NJ=ThisBlock%NJ
    NK=ThisBlock%NK
    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1






!              diffusion for SST turbulence mode
!     I direction 

	DO J=1,NJ1
      DO K=1,NK1
      DO I=1,NI
		Iswall0=0
		Iswallm=0
		if ((I.eq. 1).and.(MarkwallI0(J,K).eq.1)) Iswall0=1
		if ((I.eq.NI).and.(MarkwallIm(J,K).eq.1)) Iswallm=1

          alsinv=1./(vol(i,j,k)+vol(i-1,j,k))
          all=vol(i-1,j,k)*alsinv
          alr=vol(i,j,k)*alsinv

	    if (Iswall0 .or. Iswallm) alsinv=0.0  !wall diff

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    DTkDx1= dQdxyz( 13,I-1,J,K)
	    DTkDy1= dQdxyz( 14,I-1,J,K)
	    DTkDz1= dQdxyz( 15,I-1,J,K)

	    DTkDx= dQdxyz( 13,I,J,K)
	    DTkDy= dQdxyz( 14,I,J,K)
	    DTkDz= dQdxyz( 15,I,J,K)

          dTkdxL=all*DTkDx1+alr*DTkDx   !         VL(3,9,-2:MI,-2:MJ,-2:MK)  ! U V W P T k omega muT gm 
          dTkdyL=all*DTkDy1+alr*DTkDy
          dTkdzL=all*DTkDz1+alr*DTkDz

!	    Den    = V(1,I,J,K) 	         !UU     = V(2,I1,J1,K1)/Den  !VV     = V(3,I1,J1,K1)/Den

          disfpp=VL(1,6,i+1,j,k)            !
          disfnn=VL(1,6,i-1,j,k)            !
          disfp=V(6,I,J,K)/ V(1,I,J,K)      !
          disfn=V(6,I-1,J,K)/ V(1,I-1,J,K)  !


!          dtedx=alsinv*(-aix(i+1,j,k)*(disfpp-disfp)-aix(i-1,j,k)*(disfn-disfnn)+aix(i,j,k)*(disfp-disfn))+dtedx
           dTkdxL=alsinv*(-SD(1,1,i+1,j,k)*(disfpp-disfp)   &
     &                    -SD(1,1,i-1,j,k)*(disfn-disfnn)   &
     &                    +SD(1,1,i,j,k)*(disfp-disfn)   )  &
     &          +dTkdxL

!          dtedy=alsinv*(-aiy(i+1,j,k)*(disfpp-disfp)-aiy(i-1,j,k)*(disfn-disfnn)+aiy(i,j,k)*(disfp-disfn))+dtedy
           dTkdyL=alsinv*(-SD(1,2,i+1,j,k)*(disfpp-disfp)   &
     &                    -SD(1,2,i-1,j,k)*(disfn-disfnn)   &
     &                    +SD(1,2,i,j,k)*(disfp-disfn)   )  &
     &          +dTkdyL

!          dtedz=alsinv*(-aiz(i+1,j,k)*(disfpp-disfp)-aiz(i-1,j,k)*(disfn-disfnn)+aiz(i,j,k)*(disfp-disfn))+dtedz
           dTkdzL=alsinv*(-SD(1,3,i+1,j,k)*(disfpp-disfp)   &
     &                    -SD(1,3,i-1,j,k)*(disfn-disfnn)   &
     &                    +SD(1,3,i,j,k)*(disfp-disfn)   )  &
     &          +dTkdzL


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    DTODx1= dQdxyz( 16,I-1,J,K)
	    DTODy1= dQdxyz( 17,I-1,J,K)
	    DTODz1= dQdxyz( 18,I-1,J,K)

	    DTODx= dQdxyz( 16,I,J,K)
	    DTODy= dQdxyz( 17,I,J,K)
	    DTODz= dQdxyz( 18,I,J,K)

          dTOdxL=all*DTODx1+alr*DTODx   !         VL(ML,-2:MI,-2:MJ,-2:MK)  ! U V W P T k omega muT gm 
          dTOdyL=all*DTODy1+alr*DTODy
          dTOdzL=all*DTODz1+alr*DTODz

          disfpp=VL(1,7,i+1,j,k)           !disfpp=vil(i+1,j,k)
          disfnn=VL(1,7,i-1,j,k)           !disfnn=vil(i-1,j,k)
          disfp=V(7,I,J,K)/ V(1,I,J,K)     !vy(i,j,k)
          disfn=V(7,I-1,J,K)/ V(1,I-1,J,K) !vy(i-1,j,k)


!         doqdx=alsinv*(-aix(i+1,j,k)*(disfpp-disfp)-aix(i-1,j,k)*(disfn-disfnn)+aix(i,j,k)*(disfp-disfn))+doqdx
           dTOdxL=alsinv*(-SD(1,1,i+1,j,k)*(disfpp-disfp)   &
     &                    -SD(1,1,i-1,j,k)*(disfn-disfnn)   &
     &                    +SD(1,1,i,j,k)*(disfp-disfn)   )  &
     &           +dTOdxL

!          doqdy=alsinv*(-aiy(i+1,j,k)*(disfpp-disfp)-aiy(i-1,j,k)*(disfn-disfnn)+aiy(i,j,k)*(disfp-disfn))+doqdy
           dTOdyL=alsinv*(-SD(1,2,i+1,j,k)*(disfpp-disfp)   &  
     &                    -SD(1,2,i-1,j,k)*(disfn-disfnn)   &
     &                    +SD(1,2,i,j,k)*(disfp-disfn)   )  &
     &           +dTOdyL

!          doqdz=alsinv*(-aiz(i+1,j,k)*(disfpp-disfp)-aiz(i-1,j,k)*(disfn-disfnn)+aiz(i,j,k)*(disfp-disfn))+doqdz

           dTOdzL=alsinv*(-SD(1,3,i+1,j,k)*(disfpp-disfp)   &
     &                    -SD(1,3,i-1,j,k)*(disfn-disfnn)   &
     &                    +SD(1,3,i,j,k)*(disfp-disfn)   )  &
     &          +dTOdzL


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          alsinv=1.0/(Alagm(1,i,j,k)+Alagm(1,i-1,j,k))
          all=Alagm(1,i,j,k)*alsinv
          alr=Alagm(1,i-1,j,k)*alsinv

	    TT=T(I,J,K)
	    TT1=T(I-1,J,K)
	    TTL=all*TT1+alr*TT
!		TTL = VL(1,5,I,J,K)  !interface
          RmiuL=(1.0+Csthlnd)/(TTL+Csthlnd)*TTL**1.5   !vmuc=all*vmu(i-1,j,k)+alr*vmu(i,j,k)

         RmiuT=Rmiu(I,J,K)
         RmiuT1=Rmiu(I-1,J,K)
	    RmiuTL=all*RmiuT1+alr*RmiuT  !vmut=vmuc-vmulc
	    if (Iswall0 .or. Iswallm) RmiuTL=0.0  !wall diff
!		RmiuTL = VL(1,8,I,J,K)

	  IF (Kind_model == 2) then !WD+06
		SigkAB = sigk1
		SigoAB = sigo1
	  ELSEIF (Kind_model == 1) then
		DEN0 = V(1,I,J,K)
		TKK0 = V(6,I,J,K)/DEN0
		TOO0 = V(7,I,J,K)/DEN0
		DEN1 = V(1,I-1,J,K)
		TKK1 = V(6,I-1,J,K)/DEN1
		TOO1 = V(7,I-1,J,K)/DEN1
		
		Den = all*Den1+alr*Den0
		TKK = all*TKK1+alr*TKK0
		TOO = all*TOO1+alr*TOO0

!		Den=VL(1,4,I,J,K)/PPF/VL(1,5,I,J,K) !Den=p/rcpcv/T
!		TKK=VL(1,6,I,J,K)
!		TOO=VL(1,7,I,J,K)
      
 !       cdko=max(2.*rofc/(oqfc*prgo)*(dtedx*doqdx+dtedy*doqdy+dtedz*doqdz),1.0d-20)
	   dTKO=DTKDXL*DTODXL+DTKDyL*DTODyL+DTKDzL*DTODZL
	   CDKO=max(2.0*DEN*sigo2/TOO*(dTKO),1.0d-10)

         distfc=alr*Dst(i,j,k)+all*Dst(i-1,j,k)  !! ?????
!        arg1=max(sqrt(tefc)/(0.09*oqfc*distfc),500.*vmulc/(rofc*distfc*distfc*oqfc))
         arg1=max(sqrt(TKK)/Beta_star/TOO/distfc/Ref,           &
     &           500.0*RMIUL/DEN/distfc/distfc/TOO/Ref/Ref)
!         arg1=min(arg1,4.*rofc*tefc/(cdko*distfc*distfc*prgo))
          arg1=min(arg1,4.0*DEN*sigo2*TKK/CDKO/distfc/distfc)

          fc1c=tanh(arg1**4)

          SigkAB=fc1c*sigk1+(1.-fc1c)*sigk2
          SigoAB=fc1c*sigo1+(1.-fc1c)*sigo2
	  ENDIF

          vmute=RMIUL+RMIUTL*SigkAB
          vmuoq=RMIUL+RMIUTL*SigoAB

          Dx=SD(1,1,i,j,k) !aix(i,j,k)
          Dy=SD(1,2,i,j,k) !aiy(i,j,k)
          Dz=SD(1,3,i,j,k) !aiz(i,j,k)

          DIFFTKi=(vmute*(dTKdxL*DX+dTkdyL*DY+dTkdzL*DZ)) / Ref
      
          DIFFTOi=(vmuoq*(dTOdxL*Dx+dTOdyL*Dy+dTOdzL*Dz)) / Ref
 
	    D(6,I,J,K)=D(6,I,J,K)-DIFFTKi
          D(6,I-1,J,K)=D(6,I-1,J,K)+DIFFTKi
	
	    D(7,I,J,K)=D(7,I,J,K)-DIFFTOi
          D(7,I-1,J,K)=D(7,I-1,J,K)+DIFFTOi


	end do
	end do
	end do




!     J direction 

      do k=1,NK1
      do j=1,NJ
      do i=1,NI1    
		Iswall0=0
		Iswallm=0
		if ((J.eq. 1).and.(MarkwallJ0(I,K).eq.1)) Iswall0=1
		if ((J.eq.NJ).and.(MarkwallJm(I,K).eq.1)) Iswallm=1

	    alsinv=1.0/(vol(i,j,k)+vol(i,j-1,k))
          all=vol(i,j-1,k)*alsinv
          alr=vol(i,j,k)*alsinv

	    if (Iswall0 .or. Iswallm) alsinv=0.0  !wall diff

	    DTkDx1= dQdxyz(13,I,J-1,K)  !uir(i,j-1,k)
	    DTkDy1= dQdxyz(14,I,J-1,K)  !ujr(i,j-1,k)
	    DTkDz1= dQdxyz(15,I,J-1,K)  !ukr(i,j-1,k)

	    DTkDx= dQdxyz( 13,I,J,K)     !uir(i,j,k)
	    DTkDy= dQdxyz( 14,I,J,K)     !ujr(i,j,k)
	    DTkDz= dQdxyz( 15,I,J,K)     !ukr(i,j,k)

          dTkdxL=all*DTkdx1+alr*DTkDx  !all*uir(i,j-1,k)+alr*uir(i,j,k)
          dTkdyL=all*DTkdy1+alr*DTkDy  !all*ujr(i,j-1,k)+alr*ujr(i,j,k)
          dTkdzL=all*DTkdz1+alr*DTkDz  !all*ukr(i,j-1,k)+alr*ukr(i,j,k)

          disfpp=VL(2,6,i,j+1,k)       !ujl(i,j+1,k)
          disfnn=VL(2,6,i,j-1,k)       !ujl(i,j-1,k)
          disfp=V(6,I,J,K)/ V(1,I,J,K) !vx(i,j,k)
          disfn=V(6,I,J-1,K)/ V(1,I,J-1,K) !vx(i,j-1,k)

!          dudx=alsinv*(-ajx(i,j+1,k)*(disfpp-disfp)-ajx(i,j-1,k)*(disfn-disfnn)+ajx(i,j,k)*(disfp-disfn))+dudx

           dTkdxL=alsinv*(-SD(2,1,i,j+1,k)*(disfpp-disfp)   &
     &                    -SD(2,1,i,j-1,k)*(disfn-disfnn)   &
     &                    +SD(2,1,i,j,k)*(disfp-disfn)   )  &
     &          +dTkdxL

!          dudy=alsinv*(-ajy(i,j+1,k)*(disfpp-disfp)-ajy(i,j-1,k)*(disfn-disfnn)+ajy(i,j,k)*(disfp-disfn))+dudy
           dTkdyL=alsinv*(-SD(2,2,i,j+1,k)*(disfpp-disfp)   &
     &                    -SD(2,2,i,j-1,k)*(disfn-disfnn)   &
     &                    +SD(2,2,i,j,k)*(disfp-disfn)   )  &
     &          +dTkdyL
!          dudz=alsinv*(-ajz(i,j+1,k)*(disfpp-disfp)-ajz(i,j-1,k)*(disfn-disfnn)&+ajz(i,j,k)*(disfp-disfn))+dudz
           dTkdzL=alsinv*(-SD(2,3,i,j+1,k)*(disfpp-disfp)   &
     &                    -SD(2,3,i,j-1,k)*(disfn-disfnn)   &
     &                    +SD(2,3,i,j,k)*(disfp-disfn)   )  &
     &          +dTkdzL


	    DTODx1= dQdxyz( 16,I,J-1,K)
	    DTODy1= dQdxyz( 17,I,J-1,K)
	    DTODz1= dQdxyz( 18,I,J-1,K)

	    DTODx= dQdxyz( 16,I,J,K)
	    DTODy= dQdxyz( 17,I,J,K)
	    DTODz= dQdxyz( 18,I,J,K)

          dTOdxL=all*DTODx1+alr*DTODx   !         VL(ML,-2:MI,-2:MJ,-2:MK)  ! U V W PP T k omega muT gm 
          dTOdyL=all*DTODy1+alr*DTODy
          dTOdzL=all*DTODz1+alr*DTODz

          disfpp=VL(2,7,i,j+1,k)           !disfpp=vil(i+1,j,k)
          disfnn=VL(2,7,i,j-1,k)           !disfnn=vil(i-1,j,k)
          disfp=V(7,I,J,K)/ V(1,I,J,K)     !vy(i,j,k)
          disfn=V(7,I,J-1,K)/ V(1,I,J-1,K) !vy(i-1,j,k)

!          dvdx=alsinv*(-ajx(i,j+1,k)*(disfpp-disfp)-ajx(i,j-1,k)*(disfn-disfnn)+ajx(i,j,k)*(disfp-disfn))+dvdx
           dTOdxL=alsinv*(-SD(2,1,i,j+1,k)*(disfpp-disfp)   &
     &                    -SD(2,1,i,j-1,k)*(disfn-disfnn)   &
     &                    +SD(2,1,i,j,k)*(disfp-disfn)   )  &
     &           +dTOdxL

!          dvdy=alsinv*(-ajy(i,j+1,k)*(disfpp-disfp)-ajy(i,j-1,k)*(disfn-disfnn)+ajy(i,j,k)*(disfp-disfn))+dvdy
           dTOdyL=alsinv*(-SD(2,2,i,j+1,k)*(disfpp-disfp)   &
     &                    -SD(2,2,i,j-1,k)*(disfn-disfnn)   &
     &                    +SD(2,2,i,j,k)*(disfp-disfn)   )  &
     &           +dTOdyL

!          dvdz=alsinv*(-ajz(i,j+1,k)*(disfpp-disfp)-ajz(i,j-1,k)*(disfn-disfnn)+ajz(i,j,k)*(disfp-disfn))+dvdz
           dTOdzL=alsinv*(-SD(2,3,i,j+1,k)*(disfpp-disfp)   &
     &                    -SD(2,3,i,j-1,k)*(disfn-disfnn)   &
     &                    +SD(2,3,i,j,k)*(disfp-disfn)   )  &
     &           +dTOdzL
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          alsinv=1.0/(Alagm(2,i,j,k)+Alagm(2,i,j-1,k))
          all=Alagm(2,i,j,k)*alsinv
          alr=Alagm(2,i,j-1,k)*alsinv


	    TT=T(I,J,K)
	    TT1=T(I,J-1,K)
	    TTL=all*TT1+alr*TT
!		TTL = VL(2,5,I,J,K)  !interface
          RmiuL=(1.0+Csthlnd)/(TTL+Csthlnd)*TTL**1.5          !vmuc=all*vmu(i-1,j,k)+alr*vmu(i,j,k)

          RmiuT=Rmiu(I,J,K)
          RmiuT1=Rmiu(I,J-1,K)
	    RmiuTL=all*RmiuT1+alr*RmiuT                       !vmut=vmuc-vmulc
	    if (Iswall0 .or. Iswallm) RmiuTL=0.0  !wall diff
!		RmiuTL = VL(2,8,I,J,K)
        
	  IF (Kind_model == 2) then !k-omega
		SigkAB = sigk1
		SigoAB = sigo1
	  ELSEIF (Kind_model == 1 ) then
		DEN0 = V(1,I,J,K)
		TKK0 = V(6,I,J,K)/DEN0
		TOO0 = V(7,I,J,K)/DEN0
		DEN1 = V(1,I,J-1,K)
		TKK1 = V(6,I,J-1,K)/DEN1
		TOO1 = V(7,I,J-1,K)/DEN1
		
		Den = all*Den1+alr*Den0
		TKK = all*TKK1+alr*TKK0
		TOO = all*TOO1+alr*TOO0
!		Den=VL(2,4,I,J,K)/PPF/VL(2,5,I,J,K) !Den=p/rcpcv/T
!		TKK=VL(2,6,I,J,K)
!		TOO=VL(2,7,I,J,K)
      
 !       cdko=max(2.*rofc/(oqfc*prgo)*(dtedx*doqdx+dtedy*doqdy+dtedz*doqdz),1.0d-20)
	   dTKO=DTKDXL*DTODXL+DTKDyL*DTODyL+DTKDzL*DTODZL
	   CDKO=max(2.0*DEN*sigo2/TOO*(dTKO),1.0d-10)

         distfc=alr*Dst(i,j,k)+all*Dst(i,j-1,k)  !! ?????
!        arg1=max(sqrt(tefc)/(0.09*oqfc*distfc),500.*vmulc/(rofc*distfc*distfc*oqfc))
         arg1=max(sqrt(TKK)/Beta_star/TOO/distfc/Ref,           &
     &           500.0*RMIUL/DEN/distfc/distfc/TOO/Ref/Ref)
!         arg1=min(arg1,4.*rofc*tefc/(cdko*distfc*distfc*prgo))
          arg1=min(arg1,4.0*DEN*sigo2*TKK/CDKO/distfc/distfc)

          fc1c=tanh(arg1**4)

          SigkAB=fc1c*sigk1+(1.-fc1c)*sigk2
          SigoAB=fc1c*sigo1+(1.-fc1c)*sigo2
	  ENDIF

         vmute=RMIUL+RMIUTL*SigkAB
         vmuoq=RMIUL+RMIUTL*SigoAB

        
          Dx=SD(2,1,i,j,k)
		Dy=SD(2,2,i,j,k)
          Dz=SD(2,3,i,j,k)

          DIFFTKj=(vmute*(dTKdxL*DX+dTkdyL*DY+dTkdzL*DZ)) / Ref
      
          DIFFTOj=(vmuoq*(dTOdxL*Dx+dTOdyL*Dy+dTOdzL*Dz)) / Ref


	    D(6,I,J,K)=D(6,I,J,K)-DIFFTkj
          D(6,I,J-1,K)=D(6,I,J-1,K)+DIFFTkj
	
	    D(7,I,J,K)=D(7,I,J,K)-DIFFTOj
          D(7,I,J-1,K)=D(7,I,J-1,K)+DIFFTOj

	end do 
	end do
	end do

!     K direction 

	do k=1,NK
      do j=1,NJ1
      do i=1,NI1
		Iswall0=0
		Iswallm=0
		if ((K.eq. 1).and.(MarkwallK0(I,J).eq.1)) Iswall0=1
		if ((K.eq.NK).and.(MarkwallKm(I,J).eq.1)) Iswallm=1

          alsinv=1.0/(vol(i,j,k)+vol(i,j,k-1))
          all=vol(i,j,k-1)*alsinv
          alr=vol(i,j,k)*alsinv

	    if (Iswall0 .or. Iswallm) alsinv=0.0  !wall diff

	    DTkDx1=dQdxyz( 13,I,J,K-1)         !uir(i,j,k-1)
	    DTkDy1=dQdxyz( 14,I,J,K-1)         !ujr(i,j,k-1)
	    DTkDz1=dQdxyz( 15,I,J,K-1)         !ukr(i,j,k-1)

	    DTkDx=dQdxyz( 13,I,J,K)            !uir(i,j,k)
	    DTkDy=dQdxyz( 14,I,J,K)            !ujr(i,j,k)
	    DTkDz=dQdxyz( 15,I,J,K)            !ukr(i,j,k)


	    dTkdxL=all*DTkDx1+alr*DTkDx         !all*uir(i,j,k-1)+alr*uir(i,j,k)
          dTkdyL=all*DTkDy1+alr*DTkDy         !all*ujr(i,j,k-1)+alr*ujr(i,j,k)
          dTkdzL=all*DTkDz1+alr*DTkDz         !all*ukr(i,j,k-1)+alr*ukr(i,j,k)

          disfpp=VL(3,6,i,j,k+1)           !ukl(i,j,k+1)
          disfnn=VL(3,6,i,j,k-1)           !ukl(i,j,k-1)
          disfp=V(6,I,J,K)/ V(1,I,J,K)     !vx(i,j,k)
          disfn=V(6,I,J,K-1)/ V(1,I,J,K-1) !vx(i,j,k-1)

!         dudx=alsinv*(-akx(i,j,k+1)*(disfpp-disfp)-akx(i,j,k-1)*(disfn-disfnn)+akx(i,j,k)*(disfp-disfn))+dudx
           dTkdxL=alsinv*(-SD(3,1,i,j,k+1)*(disfpp-disfp)   &
     &                    -SD(3,1,i,j,k-1)*(disfn-disfnn)   &
     &                    +SD(3,1,i,j,k)*(disfp-disfn)   )  &
     &          +dTkdxL

!         dudy=alsinv*(-aky(i,j,k+1)*(disfpp-disfp)&-aky(i,j,k-1)*(disfn-disfnn)+aky(i,j,k)*(disfp-disfn))+dudy
           dTkdyL=alsinv*(-SD(3,2,i,j,k+1)*(disfpp-disfp)   &
     &                   -SD(3,2,i,j,k-1)*(disfn-disfnn)    &
     &                   +SD(3,2,i,j,k)*(disfp-disfn)   )   &
     &          +dTkdyL

!         dudz=alsinv*(-akz(i,j,k+1)*(disfpp-disfp)&-akz(i,j,k-1)*(disfn-disfnn)+akz(i,j,k)*(disfp-disfn))+dudz
           dTkdzL=alsinv*(-SD(3,3,i,j,k+1)*(disfpp-disfp)   &
     &                    -SD(3,3,i,j,k-1)*(disfn-disfnn)   &
     &                    +SD(3,3,i,j,k)*(disfp-disfn)   )  &
     &          +dTkdzL !dTodzl
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	    DTODx1=dQdxyz( 16,I,J,K-1)
	    DTODy1=dQdxyz( 17,I,J,K-1)
	    DTODz1=dQdxyz( 18,I,J,K-1)

	    DTODx= dQdxyz( 16,I,J,K)
	    DTODy= dQdxyz( 17,I,J,K)
	    DTODz= dQdxyz( 18,I,J,K)

          dTOdxL=all*DTODx1+alr*DTODx   !         VL(ML,-2:MI,-2:MJ,-2:MK)  ! U V W PP T k omega muT gm 
          dTOdyL=all*DTODy1+alr*DTODy
          dTOdzL=all*DTODz1+alr*DTODz
          
          disfpp=VL(3,7,i,j,k+1)           !vkl(i,j,k+1)
          disfnn=VL(3,7,i,j,k-1)           !vkl(i,j,k-1)
          disfp=V(7,I,J,K)/ V(1,I,J,K)     !vy(i,j,k)
          disfn=V(7,I,J,K-1)/ V(1,I,J,K-1) !vy(i,j,k-1)

 !        dvdx=alsinv*(-akx(i,j,k+1)*(disfpp-disfp)-akx(i,j,k-1)*(disfn-disfnn)+akx(i,j,k)*(disfp-disfn))+dvdx
           dTOdxL=alsinv*(-SD(3,1,i,j,k+1)*(disfpp-disfp)   &
     &                   -SD(3,1,i,j,k-1)*(disfn-disfnn)    &
     &                   +SD(3,1,i,j,k)*(disfp-disfn)   )   &
     &           +dTOdxL
 !        dvdy=alsinv*(-aky(i,j,k+1)*(disfpp-disfp)-aky(i,j,k-1)*(disfn-disfnn)+aky(i,j,k)*(disfp-disfn))+dvdy
            dTOdyL=alsinv*(-SD(3,2,i,j,k+1)*(disfpp-disfp)  &
     &                    -SD(3,2,i,j,k-1)*(disfn-disfnn)   &
     &                    +SD(3,2,i,j,k)*(disfp-disfn)   )  &
     &            +dTOdyL
 !        dvdz=alsinv*(-akz(i,j,k+1)*(disfpp-disfp)&-akz(i,j,k-1)*(disfn-disfnn)&+akz(i,j,k)*(disfp-disfn))+dvdz
           dTOdzL=alsinv*(-SD(3,3,i,j,k+1)*(disfpp-disfp)   &
     &                   -SD(3,3,i,j,k-1)*(disfn-disfnn)    &
     &                   +SD(3,3,i,j,k)*(disfp-disfn)   )   &
     &            +dTOdzL
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
 	   
		
	    alsinv=1.0/(alagm(3,i,j,k)+alagm(3,i,j,k-1))
          all=alagm(3,i,j,k)*alsinv
          alr=alagm(3,i,j,k-1)*alsinv

	    TT=T(I,J,K)
	    TT1=T(I,J,K-1)
	    TTL=all*TT1+alr*TT
!		TTL = VL(3,5,I,J,K)  !interface
          RmiuL=(1.0+Csthlnd)/(TTL+Csthlnd)*TTL**1.5   !vmuc=all*vmu(i-1,j,k)+alr*vmu(i,j,k)

          RmiuT=Rmiu(I,J,K)
          RmiuT1=Rmiu(I,J,K-1)
	    RmiuTL=all*RmiuT1+alr*RmiuT  !vmut=vmuc-vmulc
	    !RMIUTL=0.
	    if (Iswall0 .or. Iswallm) RmiuTL=0.0  !wall diff
!		RmiuTL = VL(3,8,I,J,K)

	  IF (Kind_model == 2) then !k-omega
		SigkAB = sigk1
		SigoAB = sigo1
	  ELSEIF (Kind_model == 1 ) then
		DEN0 = V(1,I,J,K)
		TKK0 = V(6,I,J,K)/DEN0
		TOO0 = V(7,I,J,K)/DEN0
		DEN1 = V(1,I,J,K-1)
		TKK1 = V(6,I,J,K-1)/DEN1
		TOO1 = V(7,I,J,K-1)/DEN1
		
		Den = all*Den1+alr*Den0
		TKK = all*TKK1+alr*TKK0
		TOO = all*TOO1+alr*TOO0
!		Den=VL(3,4,I,J,K)/PPF/VL(3,5,I,J,K) !Den=p/rcpcv/T
!		TKK=VL(3,6,I,J,K)
!		TOO=VL(3,7,I,J,K)
      
 !       cdko=max(2.*rofc/(oqfc*prgo)*(dtedx*doqdx+dtedy*doqdy+dtedz*doqdz),1.0d-20)
	   dTKO=DTKDXL*DTODXL+DTKDyL*DTODyL+DTKDzL*DTODZL
	   CDKO=max(2.0*DEN*sigo2/TOO*(dTKO),1.0d-10)

         distfc=alr*Dst(i,j,k)+all*Dst(i,j,k-1)  !! ?????
!        arg1=max(sqrt(tefc)/(0.09*oqfc*distfc),500.*vmulc/(rofc*distfc*distfc*oqfc))
         arg1=max(sqrt(TKK)/Beta_star/TOO/distfc/Ref,           &
     &           500.0*RMIUL/DEN/distfc/distfc/TOO/Ref/Ref)
!         arg1=min(arg1,4.*rofc*tefc/(cdko*distfc*distfc*prgo))
          arg1=min(arg1,4.0*DEN*sigo2*TKK/CDKO/distfc/distfc)

          fc1c=tanh(arg1**4)

          SigkAB=fc1c*sigk1+(1.-fc1c)*sigk2
          SigoAB=fc1c*sigo1+(1.-fc1c)*sigo2
	  ENDIF

         vmute=RMIUL+RMIUTL*SigkAB
         vmuoq=RMIUL+RMIUTL*SigoAB

          Dx=SD(3,1,i,j,k)
          Dy=SD(3,2,i,j,k)
          Dz=SD(3,3,i,j,k)

          DIFFTKk=(vmute*(dTKdxL*DX+dTkdyL*DY+dTkdzL*DZ)) / Ref
      
          DIFFTOk=(vmuoq*(dTOdxL*Dx+dTOdyL*Dy+dTOdzL*Dz)) / Ref



	    D(6,I,J,K)=D(6,I,J,K)-DIFFTKk
          D(6,i,j,k-1)=D(6,i,j,k-1)+DIFFTKk
	
	    D(7,I,J,K)=D(7,I,J,K)-DIFFTOk
          D(7,i,j,k-1)=D(7,i,j,k-1)+DIFFTOk

	end do
	end do
	end do

END SUBROUTINE TURBTWO_DIFFUSION_Centered