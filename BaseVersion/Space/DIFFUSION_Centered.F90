SUBROUTINE DIFFUSION_Centered
    USE Global
    IMPLICIT NONE
    
    INTEGER:: I,J,K,L,Lxyz,LL
    INTEGER:: I1,J1,K1
    INTEGER:: NI, NJ, NK
    INTEGER:: NI1, NJ1, NK1
    INTEGER:: II,JJ


    INTEGER:: Iswall0, Iswallm
    REAL:: NMD=-2./3.

    REAL:: alsinv, all, alr
    REAL:: disfpp,disfnn,disfp,disfn

    REAL:: dudx1,dudy1,dudz1
    REAL:: dudx,dudy,dudz
    REAL:: dudxL,dudyL,dudzL

    REAL:: dvdx1,dvdy1,dvdz1
    REAL:: dvdx,dvdy,dvdz
    REAL:: dvdxL,dvdyL,dvdzL

    REAL:: dwdx1,dwdy1,dwdz1
    REAL:: dwdx,dwdy,dwdz
    REAL:: dwdxL,dwdyL,dwdzL
    
    REAL:: dTdx1,dTdy1,dTdz1
    REAL:: dTdx,dTdy,dTdz
    REAL:: dTdxL,dTdyL,dTdzL
    
    REAL:: YL(6)



    REAL:: TTL,RmiuL,RmiuT,RmiuT1,RmiuTL
    REAL:: TxxLam,TyyLam,TzzLam,TxyLam,TyzLam,TzxLam
    REAL:: TxxTur,TyyTur,TzzTur,TxyTur,TyzTur,TzxTur
    REAL:: Txx,Tyy,Tzz,Txy,Tyz,Tzx
    REAL:: prl,prte,cp,akmu
    REAL:: ucc,vcc,wcc
    REAL:: Phix,Phiy,Phiz
    REAL:: Dx,Dy,Dz
    REAL:: DIFFroi,DIFFXi,DIFFYi,DIFFZi,DIFFEi
    REAL:: DIFFroj,DIFFXj,DIFFYj,DIFFZj,DIFFEj
    REAL:: DIFFrok,DIFFXk,DIFFYk,DIFFZk,DIFFEk

        real::ss3,ss4,Row


    NI=ThisBlock%NI
    NJ=ThisBlock%NJ
    NK=ThisBlock%NK
    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1



!     I direction 

	  DO J=1,NJ1
      DO K=1,NK1
      DO I=1,NI
		Iswall0=0
		Iswallm=0
		if ((I.eq. 1).and.(MarkwallI0(J,K).eq.1)) Iswall0=1
		if ((I.eq.NI).and.(MarkwallIm(J,K).eq.1)) Iswallm=1

          alsinv=1.0/(vol(i,j,k)+vol(i-1,j,k))
          all=vol(i-1,j,k)*alsinv
          alr=vol(i,j,k)*alsinv

	    if (Iswall0 .or. Iswallm) alsinv=0.0  !wall diff

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    DuDx1= dQdxyz( 1,I-1,J,K)
	    DuDy1= dQdxyz( 2,I-1,J,K)
	    DuDz1= dQdxyz( 3,I-1,J,K)

	    DuDx= dQdxyz( 1,I,J,K)
	    DuDy= dQdxyz( 2,I,J,K)
	    DuDz= dQdxyz( 3,I,J,K)
        
          dudxL=all*DuDx1+alr*DuDx   !         VL(3,ML,-2:MI,-2:MJ,-2:MK)  ! U V W P T k omega muT gm 
          dudyL=all*DuDy1+alr*DuDy
          dudzL=all*DuDz1+alr*DuDz
          
!          open(unit=111,file='rho record.txt')
 !         write(111,*) I, J, K, V(1,I,J,K)
  !        continue
          
         
!	    Den    = V(1,I,J,K) 	         !UU     = V(2,I1,J1,K1)/Den  !VV     = V(3,I1,J1,K1)/Den

          disfpp=VL(1,1,i+1,j,k)            !uil(i+1,j,k)
          disfnn=VL(1,1,i-1,j,k)            !uil(i-1,j,k)
          disfp=V(2,I,J,K)/ V(1,I,J,K)      !vx(i,j,k)
          disfn=V(2,I-1,J,K)/ V(1,I-1,J,K)  !vx(i-1,j,k)
          
          
         
!if(thisBlock%ID_Present_blk==4.and.i==1.and.k==1)then
!    write(*,*)"before DuDx",i,j,k,DudxL
!endif
          
!          dudx=alsinv*(-aix(i+1,j,k)*(disfpp-disfp)-aix(i-1,j,k)*(disfn-disfnn)&+aix(i,j,k)*(disfp-disfn))+dudx
           dudxL=alsinv*(-SD(1,1,i+1,j,k)*(disfpp-disfp)    &
     &                   -SD(1,1,i-1,j,k)*(disfn-disfnn)    &
     &                   +SD(1,1,i,j,k)*(disfp-disfn)   )   &
     &          +dudxL
           

!if(thisBlock%ID_Present_blk==4.and.i==1.and.k==1)then
!    write(*,*)"after Dudx",i,j,k,DudxL,alsinv,SD(1,1,i-1:i+1,j,k),disfpp,disfp,disfn,disfnn
!endif
!          dudy=alsinv*(-aiy(i+1,j,k)*(disfpp-disfp)-aiy(i-1,j,k)*(disfn-disfnn)+aiy(i,j,k)*(disfp-disfn))+dudy
           dudyL=alsinv*(-SD(1,2,i+1,j,k)*(disfpp-disfp)    &
     &                   -SD(1,2,i-1,j,k)*(disfn-disfnn)    &
     &                   +SD(1,2,i,j,k)*(disfp-disfn)   )   &
     &          +dudyL

!          dudz=alsinv*(-aiz(i+1,j,k)*(disfpp-disfp)-aiz(i-1,j,k)*(disfn-disfnn)+aiz(i,j,k)*(disfp-disfn))+dudz
           dudzL=alsinv*(-SD(1,3,i+1,j,k)*(disfpp-disfp)    &
     &                   -SD(1,3,i-1,j,k)*(disfn-disfnn)    &
     &                   +SD(1,3,i,j,k)*(disfp-disfn)   )   &
     &          +dudzL
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    DvDx1= dQdxyz( 4,I-1,J,K)
	    DvDy1= dQdxyz( 5,I-1,J,K)
	    DvDz1= dQdxyz( 6,I-1,J,K)

	    DvDx= dQdxyz( 4,I,J,K)
	    DvDy= dQdxyz( 5,I,J,K)
	    DvDz= dQdxyz( 6,I,J,K)

          dvdxL=all*DvDx1+alr*DvDx   !         VL(ML,-2:MI,-2:MJ,-2:MK)  ! U V W PP T k omega muT gm 
          dvdyL=all*DvDy1+alr*DvDy
          dvdzL=all*DvDz1+alr*DvDz

          disfpp=VL(1,2,i+1,j,k)           !disfpp=vil(i+1,j,k)
          disfnn=VL(1,2,i-1,j,k)           !disfnn=vil(i-1,j,k)
          disfp=V(3,I,J,K)/ V(1,I,J,K)     !vy(i,j,k)
          disfn=V(3,I-1,J,K)/ V(1,I-1,J,K) !vy(i-1,j,k)

!          dvdx=alsinv*(-aix(i+1,j,k)*(disfpp-disfp)-aix(i-1,j,k)*(disfn-disfnn)+aix(i,j,k)*(disfp-disfn))+dvdx
           dvdxL=alsinv*(-SD(1,1,i+1,j,k)*(disfpp-disfp)    &
     &                   -SD(1,1,i-1,j,k)*(disfn-disfnn)    &
     &                   +SD(1,1,i,j,k)*(disfp-disfn)   )   &
     &          +dvdxL
!          dvdy=alsinv*(-aiy(i+1,j,k)*(disfpp-disfp)-aiy(i-1,j,k)*(disfn-disfnn)+aiy(i,j,k)*(disfp-disfn))+dvdy
           dvdyL=alsinv*(-SD(1,2,i+1,j,k)*(disfpp-disfp)    &
     &                   -SD(1,2,i-1,j,k)*(disfn-disfnn)    &
     &                   +SD(1,2,i,j,k)*(disfp-disfn)   )   &
     &          +dvdyL
!          dvdz=alsinv*(-aiz(i+1,j,k)*(disfpp-disfp)-aiz(i-1,j,k)*(disfn-disfnn)+aiz(i,j,k)*(disfp-disfn))+dvdz
           dvdzL=alsinv*(-SD(1,3,i+1,j,k)*(disfpp-disfp)    &
     &                   -SD(1,3,i-1,j,k)*(disfn-disfnn)    &
     &                   +SD(1,3,i,j,k)*(disfp-disfn)   )   &
     &          +dvdzL
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    DwDx1= dQdxyz( 7,I-1,J,K)
	    DwDy1= dQdxyz( 8,I-1,J,K)
	    DwDz1= dQdxyz( 9,I-1,J,K)

	    DwDx= dQdxyz( 7,I,J,K)
	    DwDy= dQdxyz( 8,I,J,K)
	    DwDz= dQdxyz( 9,I,J,K)

          dwdxL=all*DwDx1+alr*DwDx   !         VL(ML,-2:MI,-2:MJ,-2:MK)  ! U V W PP T k omega muT gm 
          dwdyL=all*DwDy1+alr*DwDy
          dwdzL=all*DwDz1+alr*DwDz

          disfpp=VL(1,3,i+1,j,k)            !disfpp=wil(i+1,j,k)
          disfnn=VL(1,3,i-1,j,k)            !disfnn=wil(i-1,j,k)
          disfp=V(4,I,J,K)/V(1,I,J,K)       !disfp=vz(i,j,k)
          disfn=V(4,I-1,J,K)/V(1,I-1,J,K)   !disfn=vz(i-1,j,k)

!         dwdx=alsinv*(-aix(i+1,j,k)*(disfpp-disfp)-aix(i-1,j,k)*(disfn-disfnn)+aix(i,j,k)*(disfp-disfn))+dwdx
          dwdxL=alsinv*(-SD(1,1,i+1,j,k)*(disfpp-disfp)     &
     &                  -SD(1,1,i-1,j,k)*(disfn-disfnn)     &
     &                  +SD(1,1,i,j,k)*(disfp-disfn)   )    &
     &          +dwdxL
!          dwdy=alsinv*(-aiy(i+1,j,k)*(disfpp-disfp)-aiy(i-1,j,k)*(disfn-disfnn)+aiy(i,j,k)*(disfp-disfn))+dwdy
           dwdyL=alsinv*(-SD(1,2,i+1,j,k)*(disfpp-disfp)    &
     &                   -SD(1,2,i-1,j,k)*(disfn-disfnn)    &
     &                   +SD(1,2,i,j,k)*(disfp-disfn)   )   &
     &          +dwdyL
!          dwdz=alsinv*(-aiz(i+1,j,k)*(disfpp-disfp)-aiz(i-1,j,k)*(disfn-disfnn)+aiz(i,j,k)*(disfp-disfn))+dwdz
           dwdzL=alsinv*(-SD(1,3,i+1,j,k)*(disfpp-disfp)    &
     &                   -SD(1,3,i-1,j,k)*(disfn-disfnn)    &
     &                   +SD(1,3,i,j,k)*(disfp-disfn)   )   &
     &          +dwdzL
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    DTDx1= dQdxyz( 10,I-1,J,K)
	    DTDy1= dQdxyz( 11,I-1,J,K)
	    DTDz1= dQdxyz( 12,I-1,J,K)

	    DTDx= dQdxyz( 10,I,J,K)
	    DTDy= dQdxyz( 11,I,J,K)
	    DTDz= dQdxyz( 12,I,J,K)

          dTdxL=all*DTDx1+alr*DTDx   !         VL(ML,-2:MI,-2:MJ,-2:MK)  ! U V W PP T k omega muT gm 
          dTdyL=all*DTDy1+alr*DTDy
          dTdzL=all*DTDz1+alr*DTDz

          disfpp=VL(1,5,i+1,j,k)            !disfpp=til(i+1,j,k)
          disfnn=VL(1,5,i-1,j,k)            !disfnn=til(i-1,j,k)
          disfp=T(I,J,K)                    !disfp=tp(i,j,k)
          disfn=T(I-1,J,K)                  !disfn=tp(i-1,j,k)  
       
!          dtdx=alsinv*(-aix(i+1,j,k)*(disfpp-disfp)&-aix(i-1,j,k)*(disfn-disfnn)+aix(i,j,k)*(disfp-disfn))+dtdx
          dTdxL=alsinv*(-SD(1,1,i+1,j,k)*(disfpp-disfp)     &
     &                  -SD(1,1,i-1,j,k)*(disfn-disfnn)     &
     &                  +SD(1,1,i,j,k)*(disfp-disfn)   )    &
     &          +dTdxL
!          dtdy=alsinv*(-aiy(i+1,j,k)*(disfpp-disfp)&-aiy(i-1,j,k)*(disfn-disfnn)+aiy(i,j,k)*(disfp-disfn))+dtdy
           dTdyL=alsinv*(-SD(1,2,i+1,j,k)*(disfpp-disfp)    &
     &                   -SD(1,2,i-1,j,k)*(disfn-disfnn)    &
     &                   +SD(1,2,i,j,k)*(disfp-disfn)   )   &
     &          +dTdyL
!          dtdz=alsinv*(-aiz(i+1,j,k)*(disfpp-disfp)&-aiz(i-1,j,k)*(disfn-disfnn)+aiz(i,j,k)*(disfp-disfn))+dtdz
           dTdzL=alsinv*(-SD(1,3,i+1,j,k)*(disfpp-disfp)    &
     &                   -SD(1,3,i-1,j,k)*(disfn-disfnn)    &
     &                   +SD(1,3,i,j,k)*(disfp-disfn)   )   &
     &          +dTdzL
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          alsinv=1.0/(Alagm(1,i,j,k)+Alagm(1,i-1,j,k))
          all=Alagm(1,i,j,k)*alsinv
          alr=Alagm(1,i-1,j,k)*alsinv

!	    TT=T(I,J,K)
!	    TT1=T(I-1,J,K)
!	    TTL=all*TT1+alr*TT

		TTL = VL(1,5,I,J,K)  !interface
          RmiuL=(1.0+Csthlnd)/(TTL*Csth1+Csthlnd)*(Csth1*TTL)**1.5   !vmuc=all*vmu(i-1,j,k)+alr*vmu(i,j,k)
          RmiuT=Rmiu(I,J,K)
          RmiuT1=Rmiu(I-1,J,K)
	    RmiuTL=all*RmiuT1+alr*RmiuT  !vmut=vmuc-vmulc
        
	    if (Iswall0 .or. Iswallm) RmiuTL=0.0  !wall diff
!		RmiuTL = VL(1,8,I,J,K)
!         write(*,*) RMIUTL

	    TxxLam=RmiuL*(NMD*(DuDxL+DvDyL+DwDzL)+2.0*DuDxL)   !txx=2.0*vmuc*dudx+vlac*div
	    TyyLam=RmiuL*(NMD*(DuDxL+DvDyL+DwDzL)+2.0*DvDyL)   !tyy=2.0*vmuc*dvdy+vlac*div
	    TzzLam=RmiuL*(NMD*(DuDxL+DvDyL+DwDzL)+2.0*DwDzL)   !tzz=2.0*vmuc*dwdz+vlac*div
	    TxyLam=RmiuL*(     DuDyL+DvDxL)                   !txy=vmuc*(dudy+dvdx)
	    TyzLam=RmiuL*(     DvDzL+DwDyL)                   !tyz=vmuc*(dvdz+dwdy)
	    TzxLam=RmiuL*(     DuDzL+DwDxL)                   !tzx=vmuc*(dwdx+dudz)

	    TxxTur=RmiuTL*(NMD*(DuDxL+DvDyL+DwDzL)+2.0*DuDxL)
	    TyyTur=RmiuTL*(NMD*(DuDxL+DvDyL+DwDzL)+2.0*DvDyL)
	    TzzTur=RmiuTL*(NMD*(DuDxL+DvDyL+DwDzL)+2.0*DwDzL)
	    TxyTur=RmiuTL*(     DuDyL+DvDxL)
	    TyzTur=RmiuTL*(     DvDzL+DwDyL)
		TzxTur=RmiuTL*(     DuDzL+DwDxL)

!if(thisBlock%ID_Present_blk==4)then
!    write(*,*)"before D",KI,i,j,k,DudxL,DuDyL,DuDzL,DvDxL,DvDyL,DvDzL,DwDxL,DwDyL,DwDzL
!endif
		
		 
	    Txx=TxxLam+TxxTur
	    Tyy=TyyLam+TyyTur
	    Tzz=TzzLam+TzzTur
	    Txy=TxyLam+TxyTur
	    Tyz=TyzLam+TyzTur
	    Tzx=TzxLam+TzxTur

!!!!!!!!!??????? Prantl?
		prl=0.72
		prte=0.9
	    cp= 1.4*PPF/(1.4-1.0)   !gam(I,J,K)*PPF/(gam(I,J,K)-1.0)
	    akmu=cp*(RmiuL/prl+RmiuTL/prte)


          ucc=VL(1,1,i,j,k)  !ucc=uil(i,j,k)
          vcc=VL(1,2,i,j,k)  !vcc=vil(i,j,k)
          wcc=VL(1,3,i,j,k)  !wcc=wil(i,j,k)

          Phix=akmu*dTdxL+ucc*Txx+vcc*Txy+wcc*Tzx
          Phiy=akmu*dTdyL+ucc*Txy+vcc*Tyy+wcc*Tyz
          Phiz=akmu*dTdzL+ucc*Tzx+vcc*Tyz+wcc*Tzz
        
          Dx=SD(1,1,i,j,k)
          Dy=SD(1,2,i,j,k)
          Dz=SD(1,3,i,j,k)

	    DIFFroi =   0.0                                              !droi= 0.0
	    DIFFXi =(  Txx * Dx +  Txy * Dy +  Tzx * Dz ) / Ref        !dxi = aix(i  ,j,k)*txx+aiy(i  ,j,k)*txy+aiz(i  ,j,k)*tzx
	    DIFFYi =(  Txy * Dx +  Tyy * Dy +  Tyz * Dz ) / Ref        !dyi = aix(i  ,j,k)*txy+aiy(i  ,j,k)*tyy+aiz(i  ,j,k)*tyz
	    DIFFZi =(  Tzx * Dx +  Tyz * Dy +  Tzz * Dz ) / Ref        !dzi = aix(i  ,j,k)*tzx+aiy(i  ,j,k)*tyz+aiz(i  ,j,k)*tzz
	    DIFFEi =( Phix * Dx + Phiy * Dy + Phiz * Dz ) / Ref        !drei= aix(i  ,j,k)*btx+aiy(i  ,j,k)*bty+aiz(i  ,j,k)*btz
	    D(1,I,J,K)=D(1,I,J,K)-DIFFroi
          D(1,I-1,J,K)=D(1,I-1,J,K)+DIFFroi
	
	    D(2,I,J,K)=D(2,I,J,K)-DIFFXi
          D(2,I-1,J,K)=D(2,I-1,J,K)+DIFFXi

	    D(3,I,J,K)=D(3,I,J,K)-DIFFYi
          D(3,I-1,J,K)=D(3,I-1,J,K)+DIFFYi

	    D(4,I,J,K)=D(4,I,J,K)-DIFFZi
          D(4,I-1,J,K)=D(4,I-1,J,K)+DIFFZi

	    D(5,I,J,K)=D(5,I,J,K)-DIFFEi
          D(5,I-1,J,K)=D(5,I-1,J,K)+DIFFEi

!if(thisBlock%ID_Present_blk==4)then
!    write(*,*)"after D",KI,i,j,k,D(1:5,i,j,k)
!endif
	end do
	end do
	end do

     continue

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

	    if (Iswall0 .or. Iswallm) alsinv=0.  !wall diff

	    DuDx1= dQdxyz( 1,I,J-1,K)  !uir(i,j-1,k)
	    DuDy1= dQdxyz( 2,I,J-1,K)  !ujr(i,j-1,k)
	    DuDz1= dQdxyz( 3,I,J-1,K)  !ukr(i,j-1,k)

	    DuDx= dQdxyz( 1,I,J,K)     !uir(i,j,k)
	    DuDy= dQdxyz( 2,I,J,K)     !ujr(i,j,k)
	    DuDz= dQdxyz( 3,I,J,K)     !ukr(i,j,k)

          dudxL=all*Dudx1+alr*DuDx  !all*uir(i,j-1,k)+alr*uir(i,j,k)
          dudyL=all*Dudy1+alr*DuDy  !all*ujr(i,j-1,k)+alr*ujr(i,j,k)
          dudzL=all*Dudz1+alr*DuDz  !all*ukr(i,j-1,k)+alr*ukr(i,j,k)

          disfpp=VL(2,1,i,j+1,k)       !ujl(i,j+1,k)
          disfnn=VL(2,1,i,j-1,k)       !ujl(i,j-1,k)
          disfp=V(2,I,J,K)/ V(1,I,J,K) !vx(i,j,k)
          disfn=V(2,I,J-1,K)/ V(1,I,J-1,K) !vx(i,j-1,k)

!          dudx=alsinv*(-ajx(i,j+1,k)*(disfpp-disfp)-ajx(i,j-1,k)*(disfn-disfnn)+ajx(i,j,k)*(disfp-disfn))+dudx

           dudxL=alsinv*(-SD(2,1,i,j+1,k)*(disfpp-disfp)    &
     &                   -SD(2,1,i,j-1,k)*(disfn-disfnn)    &
     &                   +SD(2,1,i,j,k)*(disfp-disfn)   )   &
     &          +dudxL

!          dudy=alsinv*(-ajy(i,j+1,k)*(disfpp-disfp)-ajy(i,j-1,k)*(disfn-disfnn)+ajy(i,j,k)*(disfp-disfn))+dudy
           dudyL=alsinv*(-SD(2,2,i,j+1,k)*(disfpp-disfp)    &
     &                   -SD(2,2,i,j-1,k)*(disfn-disfnn)    &
     &                   +SD(2,2,i,j,k)*(disfp-disfn)   )   &
     &          +dudyL
!          dudz=alsinv*(-ajz(i,j+1,k)*(disfpp-disfp)-ajz(i,j-1,k)*(disfn-disfnn)&+ajz(i,j,k)*(disfp-disfn))+dudz
           dudzL=alsinv*(-SD(2,3,i,j+1,k)*(disfpp-disfp)    &
     &                   -SD(2,3,i,j-1,k)*(disfn-disfnn)    &
     &                   +SD(2,3,i,j,k)*(disfp-disfn)   )   &
     &          +dudzL


	    DvDx1= dQdxyz( 4,I,J-1,K)
	    DvDy1= dQdxyz( 5,I,J-1,K)
	    DvDz1= dQdxyz( 6,I,J-1,K)

	    DvDx= dQdxyz( 4,I,J,K)
	    DvDy= dQdxyz( 5,I,J,K)
	    DvDz= dQdxyz( 6,I,J,K)

          dvdxL=all*DvDx1+alr*DvDx   !         VL(ML,-2:MI,-2:MJ,-2:MK)  ! U V W PP T k omega muT gm 
          dvdyL=all*DvDy1+alr*DvDy
          dvdzL=all*DvDz1+alr*DvDz

          disfpp=VL(2,2,i,j+1,k)           !disfpp=vil(i+1,j,k)
          disfnn=VL(2,2,i,j-1,k)           !disfnn=vil(i-1,j,k)
          disfp=V(3,I,J,K)/ V(1,I,J,K)     !vy(i,j,k)
          disfn=V(3,I,J-1,K)/ V(1,I,J-1,K) !vy(i-1,j,k)

!          dvdx=alsinv*(-ajx(i,j+1,k)*(disfpp-disfp)-ajx(i,j-1,k)*(disfn-disfnn)+ajx(i,j,k)*(disfp-disfn))+dvdx
           dvdxL=alsinv*(-SD(2,1,i,j+1,k)*(disfpp-disfp)    &
     &                   -SD(2,1,i,j-1,k)*(disfn-disfnn)    &
     &                   +SD(2,1,i,j,k)*(disfp-disfn)   )   &
     &          +dvdxL

!          dvdy=alsinv*(-ajy(i,j+1,k)*(disfpp-disfp)-ajy(i,j-1,k)*(disfn-disfnn)+ajy(i,j,k)*(disfp-disfn))+dvdy
           dvdyL=alsinv*(-SD(2,2,i,j+1,k)*(disfpp-disfp)    &
     &                   -SD(2,2,i,j-1,k)*(disfn-disfnn)    &
     &                   +SD(2,2,i,j,k)*(disfp-disfn)   )   &
     &          +dvdyL

!          dvdz=alsinv*(-ajz(i,j+1,k)*(disfpp-disfp)-ajz(i,j-1,k)*(disfn-disfnn)+ajz(i,j,k)*(disfp-disfn))+dvdz
           dvdzL=alsinv*(-SD(2,3,i,j+1,k)*(disfpp-disfp)    &
     &                   -SD(2,3,i,j-1,k)*(disfn-disfnn)    &
     &                   +SD(2,3,i,j,k)*(disfp-disfn)   )   &
     &          +dvdzL
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    DwDx1= dQdxyz( 7,I,J-1,K)
	    DwDy1= dQdxyz( 8,I,J-1,K)
	    DwDz1= dQdxyz( 9,I,J-1,K)

	    DwDx= dQdxyz( 7,I,J,K)
	    DwDy= dQdxyz( 8,I,J,K)
	    DwDz= dQdxyz( 9,I,J,K)

          dwdxL=all*DwDx1+alr*DwDx   !         VL(ML,-2:MI,-2:MJ,-2:MK)  ! U V W PP T k omega muT gm 
          dwdyL=all*DwDy1+alr*DwDy
          dwdzL=all*DwDz1+alr*DwDz

          disfpp=VL(2,3,i,j+1,k)            !disfpp=wil(i+1,j,k)
          disfnn=VL(2,3,i,j-1,k)            !disfnn=wil(i-1,j,k)
          disfp=V(4,I,J,K)/V(1,I,J,K)       !disfp=vz(i,j,k)
          disfn=V(4,I,J-1,K)/V(1,I,J-1,K)   !disfn=vz(i-1,j,k)


!         dwdx=alsinv*(-ajx(i,j+1,k)*(disfpp-disfp)-ajx(i,j-1,k)*(disfn-disfnn)+ajx(i,j,k)*(disfp-disfn))+dwdx

          dwdxL=alsinv*(-SD(2,1,i,j+1,k)*(disfpp-disfp)     &
     &                  -SD(2,1,i,j-1,k)*(disfn-disfnn)     &
     &                  +SD(2,1,i,j,k)*(disfp-disfn)   )    &
     &          +dwdxL


!          dwdy=alsinv*(-ajy(i,j+1,k)*(disfpp-disfp)-ajy(i,j-1,k)*(disfn-disfnn)+ajy(i,j,k)*(disfp-disfn))+dwdy

           dwdyL=alsinv*(-SD(2,2,i,j+1,k)*(disfpp-disfp)    &
     &                   -SD(2,2,i,j-1,k)*(disfn-disfnn)    &
     &                   +SD(2,2,i,j,k)*(disfp-disfn)   )   &
     &          +dwdyL


!          dwdz=alsinv*(-ajz(i,j+1,k)*(disfpp-disfp)-ajz(i,j-1,k)*(disfn-disfnn)&+ajz(i,j,k)*(disfp-disfn))+dwdz
           dwdzL=alsinv*(-SD(2,3,i,j+1,k)*(disfpp-disfp)    &
     &                   -SD(2,3,i,j-1,k)*(disfn-disfnn)    &
     &                   +SD(2,3,i,j,k)*(disfp-disfn)   )   &
     &          +dwdzL

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    DTDx1= dQdxyz( 10,I,J-1,K)
	    DTDy1= dQdxyz( 11,I,J-1,K)
	    DTDz1= dQdxyz( 12,I,J-1,K)

	    DTDx= dQdxyz( 10,I,J,K)
	    DTDy= dQdxyz( 11,I,J,K)
	    DTDz= dQdxyz( 12,I,J,K)

          dTdxL=all*DTDx1+alr*DTDx   !         VL(ML,-2:MI,-2:MJ,-2:MK)  ! U V W PP T k omega muT gm 
          dTdyL=all*DTDy1+alr*DTDy
          dTdzL=all*DTDz1+alr*DTDz

          disfpp=VL(2,5,i,j+1,k)            ! disfpp=tjl(i,j+1,k)
          disfnn=VL(2,5,i,j-1,k)            ! disfnn=tjl(i,j-1,k)
          disfp=T(I,J,K)                    ! disfp=tp(i,j,k)
          disfn=T(I,J-1,K)                  ! disfn=tp(i,j-1,k)  
       

!         dtdx=alsinv*(-ajx(i,j+1,k)*(disfpp-disfp)-ajx(i,j-1,k)*(disfn-disfnn)+ajx(i,j,k)*(disfp-disfn))+dtdx
          dTdxL=alsinv*(-SD(2,1,i,j+1,k)*(disfpp-disfp)     &
     &                  -SD(2,1,i,j-1,k)*(disfn-disfnn)     &
     &                  +SD(2,1,i,j,k)*(disfp-disfn)   )    &
     &          +dTdxL


!          dtdy=alsinv*(-ajy(i,j+1,k)*(disfpp-disfp)-ajy(i,j-1,k)*(disfn-disfnn)&+ajy(i,j,k)*(disfp-disfn))+dtdy

           dTdyL=alsinv*(-SD(2,2,i,j+1,k)*(disfpp-disfp)    &
     &                   -SD(2,2,i,j-1,k)*(disfn-disfnn)    &
     &                   +SD(2,2,i,j,k)*(disfp-disfn)   )   &
     &          +dTdyL

!          dtdz=alsinv*(-ajz(i,j+1,k)*(disfpp-disfp)&-ajz(i,j-1,k)*(disfn-disfnn)&+ajz(i,j,k)*(disfp-disfn))+dtdz
           dTdzL=alsinv*(-SD(2,3,i,j+1,k)*(disfpp-disfp)    &
     &                   -SD(2,3,i,j-1,k)*(disfn-disfnn)    &
     &                   +SD(2,3,i,j,k)*(disfp-disfn)   )   &
     &          +dTdzL
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          alsinv=1.0/(Alagm(2,i,j,k)+Alagm(2,i,j-1,k))
          all=Alagm(2,i,j,k)*alsinv
          alr=Alagm(2,i,j-1,k)*alsinv

!	    TT=T(I,J,K)
!	    TT1=T(I,J-1,K)
!	    TTL=all*TT1+alr*TT
		TTL = VL(2,5,I,J,K)  !interface
          RmiuL=(1.0+Csthlnd)/(TTL*Csth1+Csthlnd)*(Csth1*TTL)**1.5   !vmuc=all*vmu(i-1,j,k)+alr*vmu(i,j,k)

          RmiuT=Rmiu(I,J,K)
          RmiuT1=Rmiu(I,J-1,K)
          RmiuTL=all*RmiuT1+alr*RmiuT                       !vmut=vmuc-vmulc
          
	    if (Iswall0 .or. Iswallm) RmiuTL=0.0  !wall diff
!		RmiuTL = VL(2,8,I,J,K)

	    TxxLam=RmiuL*(NMD*(DuDxL+DvDyL+DwDzL)+2.0*DuDxL)   !txx=2.0*vmuc*dudx+vlac*div
	    TyyLam=RmiuL*(NMD*(DuDxL+DvDyL+DwDzL)+2.0*DvDyL)   !tyy=2.0*vmuc*dvdy+vlac*div
	    TzzLam=RmiuL*(NMD*(DuDxL+DvDyL+DwDzL)+2.0*DwDzL)   !tzz=2.0*vmuc*dwdz+vlac*div
	    TxyLam=RmiuL*(     DuDyL+DvDxL)                   !txy=vmuc*(dudy+dvdx)
	    TyzLam=RmiuL*(     DvDzL+DwDyL)                   !tyz=vmuc*(dvdz+dwdy)
	    TzxLam=RmiuL*(     DuDzL+DwDxL)                   !tzx=vmuc*(dwdx+dudz)

	    TxxTur=RmiuTL*(NMD*(DuDxL+DvDyL+DwDzL)+2.0*DuDxL)
	    TyyTur=RmiuTL*(NMD*(DuDxL+DvDyL+DwDzL)+2.0*DvDyL)
	    TzzTur=RmiuTL*(NMD*(DuDxL+DvDyL+DwDzL)+2.0*DwDzL)
	    TxyTur=RmiuTL*(     DuDyL+DvDxL)
	    TyzTur=RmiuTL*(     DvDzL+DwDyL)
		TzxTur=RmiuTL*(     DuDzL+DwDxL)


	    Txx=TxxLam+TxxTur
	    Tyy=TyyLam+TyyTur
	    Tzz=TzzLam+TzzTur
	    Txy=TxyLam+TxyTur
	    Tyz=TyzLam+TyzTur
	    Tzx=TzxLam+TzxTur
        
	    cp= 1.4*PPF/(1.4-1.0)   !gam(I,J,K)*PPF/(gam(I,J,K)-1.0)
	    akmu=cp*(RmiuL/prl+RmiuTL/prte)


          ucc=VL(2,1,i,j,k)  !ucc=ujl(i,j,k)
          vcc=VL(2,2,i,j,k)  !vcc=vjl(i,j,k)
          wcc=VL(2,3,i,j,k)  !wcc=wjl(i,j,k)

          Phix=akmu*dTdxL+ucc*Txx+vcc*Txy+wcc*Tzx
          Phiy=akmu*dTdyL+ucc*Txy+vcc*Tyy+wcc*Tyz
          Phiz=akmu*dTdzL+ucc*Tzx+vcc*Tyz+wcc*Tzz
        
          Dx=SD(2,1,i,j,k)
          Dy=SD(2,2,i,j,k)
          Dz=SD(2,3,i,j,k)

	    DIFFroj =   0.                                             !droj= 0.
	    DIFFXj =(  Txx * Dx +  Txy * Dy +  Tzx * Dz ) / Ref        !dxi = ajx(i  ,j,k)*txx+ajy(i  ,j,k)*txy+ajz(i  ,j,k)*tzx
	    DIFFYj =(  Txy * Dx +  Tyy * Dy +  Tyz * Dz ) / Ref        !dyj = ajx(i  ,j,k)*txy+ajy(i  ,j,k)*tyy+ajz(i  ,j,k)*tyz
	    DIFFZj =(  Tzx * Dx +  Tyz * Dy +  Tzz * Dz ) / Ref        !dzj = ajx(i  ,j,k)*tzx+ajy(i  ,j,k)*tyz+ajz(i  ,j,k)*tzz
	    DIFFEj =( Phix * Dx + Phiy * Dy + Phiz * Dz ) / Ref        !drej= ajx(i  ,j,k)*btx+ajy(i  ,j,k)*bty+ajz(i  ,j,k)*btz

      
	    D(1,I,J,K)=D(1,I,J,K)-DIFFroj
          D(1,I,J-1,K)=D(1,I,J-1,K)+DIFFroj
	
	    D(2,I,J,K)=D(2,I,J,K)-DIFFXj
          D(2,I,J-1,K)=D(2,I,J-1,K)+DIFFXj

	    D(3,I,J,K)=D(3,I,J,K)-DIFFYj
          D(3,I,J-1,K)=D(3,I,J-1,K)+DIFFYj

	    D(4,I,J,K)=D(4,I,J,K)-DIFFZj
          D(4,I,J-1,K)=D(4,I,J-1,K)+DIFFZj

	    D(5,I,J,K)=D(5,I,J,K)-DIFFEj
          D(5,I,J-1,K)=D(5,I,J-1,K)+DIFFEj
          
       
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

	    DuDx1=dQdxyz( 1,I,J,K-1)         !uir(i,j,k-1)
	    DuDy1=dQdxyz( 2,I,J,K-1)         !ujr(i,j,k-1)
	    DuDz1=dQdxyz( 3,I,J,K-1)         !ukr(i,j,k-1)

	    DuDx=dQdxyz( 1,I,J,K)            !uir(i,j,k)
	    DuDy=dQdxyz( 2,I,J,K)            !ujr(i,j,k)
	    DuDz=dQdxyz( 3,I,J,K)            !ukr(i,j,k)


	    dudxL=all*DuDx1+alr*DuDx         !all*uir(i,j,k-1)+alr*uir(i,j,k)
          dudyL=all*DuDy1+alr*DuDy         !all*ujr(i,j,k-1)+alr*ujr(i,j,k)
          dudzL=all*DuDz1+alr*DuDz         !all*ukr(i,j,k-1)+alr*ukr(i,j,k)

          disfpp=VL(3,1,i,j,k+1)           !ukl(i,j,k+1)
          disfnn=VL(3,1,i,j,k-1)           !ukl(i,j,k-1)
          disfp=V(2,I,J,K)/ V(1,I,J,K)     !vx(i,j,k)
          disfn=V(2,I,J,K-1)/ V(1,I,J,K-1) !vx(i,j,k-1)

        
 !if(ThisBlock%ID_Present_Blk==4)write(*,*)"KDuDxL",KI,i,j,k,SD(3,1,i,j,k-1:k+1),DuDxL,disfpp,disfp,disfnn,disfn
!         dudx=alsinv*(-akx(i,j,k+1)*(disfpp-disfp)-akx(i,j,k-1)*(disfn-disfnn)+akx(i,j,k)*(disfp-disfn))+dudx
           dudxL=alsinv*(-SD(3,1,i,j,k+1)*(disfpp-disfp)    &
     &                   -SD(3,1,i,j,k-1)*(disfn-disfnn)    &
     &                   +SD(3,1,i,j,k)*(disfp-disfn)   )   &
     &          +dudxL

!         dudy=alsinv*(-aky(i,j,k+1)*(disfpp-disfp)&-aky(i,j,k-1)*(disfn-disfnn)+aky(i,j,k)*(disfp-disfn))+dudy
           dudyL=alsinv*(-SD(3,2,i,j,k+1)*(disfpp-disfp)    &
     &                   -SD(3,2,i,j,k-1)*(disfn-disfnn)    &
     &                   +SD(3,2,i,j,k)*(disfp-disfn)   )   &
     &          +dudyL

!         dudz=alsinv*(-akz(i,j,k+1)*(disfpp-disfp)&-akz(i,j,k-1)*(disfn-disfnn)+akz(i,j,k)*(disfp-disfn))+dudz
           dudzL=alsinv*(-SD(3,3,i,j,k+1)*(disfpp-disfp)    &
     &                   -SD(3,3,i,j,k-1)*(disfn-disfnn)    &
     &                   +SD(3,3,i,j,k)*(disfp-disfn)   )   &
     &          +dudzL
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	    DvDx1= dQdxyz( 4,I,J,K-1)
	    DvDy1= dQdxyz( 5,I,J,K-1)
	    DvDz1= dQdxyz( 6,I,J,K-1)

	    DvDx= dQdxyz( 4,I,J,K)
	    DvDy= dQdxyz( 5,I,J,K)
	    DvDz= dQdxyz( 6,I,J,K)

          dvdxL=all*DvDx1+alr*DvDx   !         VL(ML,-2:MI,-2:MJ,-2:MK)  ! U V W PP T k omega muT gm 
          dvdyL=all*DvDy1+alr*DvDy
          dvdzL=all*DvDz1+alr*DvDz
          
          disfpp=VL(3,2,i,j,k+1)           !vkl(i,j,k+1)
          disfnn=VL(3,2,i,j,k-1)           !vkl(i,j,k-1)
          disfp=V(3,I,J,K)/ V(1,I,J,K)     !vy(i,j,k)
          disfn=V(3,I,J,K-1)/ V(1,I,J,K-1) !vy(i,j,k-1)

 !        dvdx=alsinv*(-akx(i,j,k+1)*(disfpp-disfp)-akx(i,j,k-1)*(disfn-disfnn)+akx(i,j,k)*(disfp-disfn))+dvdx
           dvdxL=alsinv*(-SD(3,1,i,j,k+1)*(disfpp-disfp)    &
     &                   -SD(3,1,i,j,k-1)*(disfn-disfnn)    &
     &                   +SD(3,1,i,j,k)*(disfp-disfn)   )   &
     &          +dvdxL
 !        dvdy=alsinv*(-aky(i,j,k+1)*(disfpp-disfp)-aky(i,j,k-1)*(disfn-disfnn)+aky(i,j,k)*(disfp-disfn))+dvdy
            dvdyL=alsinv*(-SD(3,2,i,j,k+1)*(disfpp-disfp)   &
     &                    -SD(3,2,i,j,k-1)*(disfn-disfnn)   &
     &                    +SD(3,2,i,j,k)*(disfp-disfn)   )  &
     &          +dvdyL
 !        dvdz=alsinv*(-akz(i,j,k+1)*(disfpp-disfp)&-akz(i,j,k-1)*(disfn-disfnn)&+akz(i,j,k)*(disfp-disfn))+dvdz
           dvdzL=alsinv*(-SD(3,3,i,j,k+1)*(disfpp-disfp)    &
     &                   -SD(3,3,i,j,k-1)*(disfn-disfnn)    &
     &                   +SD(3,3,i,j,k)*(disfp-disfn)   )   &
     &          +dvdzL
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
 	    DwDx1= dQdxyz( 7,I,J,K-1)
	    DwDy1= dQdxyz( 8,I,J,K-1)
	    DwDz1= dQdxyz( 9,I,J,K-1)

	    DwDx= dQdxyz( 7,I,J,K)
	    DwDy= dQdxyz( 8,I,J,K)
	    DwDz= dQdxyz( 9,I,J,K)

          dwdxL=all*DwDx1+alr*DwDx   !         VL(ML,-2:MI,-2:MJ,-2:MK)  ! U V W PP T k omega muT gm 
          dwdyL=all*DwDy1+alr*DwDy
          dwdzL=all*DwDz1+alr*DwDz

          disfpp=VL(3,3,i,j,k+1)            !disfpp=wkl(i,j,k+1)
          disfnn=VL(3,3,i,j,k-1)            !disfnn=wkl(i,j,k-1)
          disfp=V(4,I,J,K)/V(1,I,J,K)       !disfp=vz(i,j,k)
          disfn=V(4,I,J,K-1)/V(1,I,J,K-1)   !disfn=vz(i,j,k-1)
 
!         dwdx=alsinv*(-akx(i,j,k+1)*(disfpp-disfp)-akx(i,j,k-1)*(disfn-disfnn)+akx(i,j,k)*(disfp-disfn))+dwdx
          dwdxL=alsinv*(-SD(3,1,i,j,k+1)*(disfpp-disfp)     &
     &                  -SD(3,1,i,j,k-1)*(disfn-disfnn)     &
     &                  +SD(3,1,i,j,k)*(disfp-disfn)   )    &
     &          +dwdxL
!          dwdy=alsinv*(-aky(i,j,k+1)*(disfpp-disfp)-aky(i,j,k-1)*(disfn-disfnn)+aky(i,j,k)*(disfp-disfn))+dwdy
           dwdyL=alsinv*(-SD(3,2,i,j,k+1)*(disfpp-disfp)    &
     &                   -SD(3,2,i,j,k-1)*(disfn-disfnn)    &
     &                   +SD(3,2,i,j,k)*(disfp-disfn)   )   &
     &          +dwdyL
!          dwdz=alsinv*(-akz(i,j,k+1)*(disfpp-disfp)-akz(i,j,k-1)*(disfn-disfnn)+akz(i,j,k)*(disfp-disfn))+dwdz
           dwdzL=alsinv*(-SD(3,3,i,j,k+1)*(disfpp-disfp)    &
     &                   -SD(3,3,i,j,k-1)*(disfn-disfnn)    &
     &                   +SD(3,3,i,j,k)*(disfp-disfn)   )   &
     &          +dwdzL
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		DTDx1= dQdxyz( 10,I,J,K-1)
	    DTDy1= dQdxyz( 11,I,J,K-1)
	    DTDz1= dQdxyz( 12,I,J,K-1)

	    DTDx= dQdxyz( 10,I,J,K)
	    DTDy= dQdxyz( 11,I,J,K)
	    DTDz= dQdxyz( 12,I,J,K)

          dTdxL=all*DTDx1+alr*DTDx   !         VL(3,ML,-2:MI,-2:MJ,-2:MK)  ! U V W P T k omega muT gm 
          dTdyL=all*DTDy1+alr*DTDy
          dTdzL=all*DTDz1+alr*DTDz

          disfpp=VL(3,5,i,j,k+1)            !disfpp=til(i,j,k+1)
          disfnn=VL(3,5,i,j,k-1)            !disfnn=til(i,j,k-1)
          disfp=T(I,J,K)                    !disfp=tp(i,j,k)
          disfn=T(I,J,K-1)                  !disfn=tp(i,j,k-1)  

!         dtdx=alsinv*(-akx(i,j,k+1)*(disfpp-disfp)-akx(i,j,k-1)*(disfn-disfnn)+akx(i,j,k)*(disfp-disfn))+dtdx
          dTdxL=alsinv*(-SD(3,1,i,j,k+1)*(disfpp-disfp)     &
     &                  -SD(3,1,i,j,k-1)*(disfn-disfnn)     &
     &                  +SD(3,1,i,j,k)*(disfp-disfn)   )    &
     &          +dTdxL
!         dtdy=alsinv*(-aky(i,j,k+1)*(disfpp-disfp)-aky(i,j,k-1)*(disfn-disfnn)+aky(i,j,k)*(disfp-disfn))+dtdy
           dTdyL=alsinv*(-SD(3,2,i,j,k+1)*(disfpp-disfp)    &
     &                   -SD(3,2,i,j,k-1)*(disfn-disfnn)    &
     &                   +SD(3,2,i,j,k)*(disfp-disfn)   )   &
     &          +dTdyL
!         dtdz=alsinv*(-akz(i,j,k+1)*(disfpp-disfp)-akz(i,j,k-1)*(disfn-disfnn)+akz(i,j,k)*(disfp-disfn))+dtdz
           dTdzL=alsinv*(-SD(3,3,i,j,k+1)*(disfpp-disfp)    &
     &                   -SD(3,3,i,j,k-1)*(disfn-disfnn)    &
     &                   +SD(3,3,i,j,k)*(disfp-disfn)   )   &
     &          +dTdzL
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    alsinv=1.0/(alagm(3,i,j,k)+alagm(3,i,j,k-1))
          all=alagm(3,i,j,k)*alsinv
          alr=alagm(3,i,j,k-1)*alsinv

!	    TT=T(I,J,K)
!	    TT1=T(I,J,K-1)
!	    TTL=all*TT1+alr*TT
		TTL = VL(3,5,I,J,K)  !interface
          RmiuL=(1.0+Csthlnd)/(TTL*Csth1+Csthlnd)*(Csth1*TTL)**1.5   !vmuc=all*vmu(i-1,j,k)+alr*vmu(i,j,k)

          RmiuT=Rmiu(I,J,K)
          RmiuT1=Rmiu(I,J,K-1)
	    RmiuTL=all*RmiuT1+alr*RmiuT  !vmut=vmuc-vmulc
!	    RMIUTL=0.
	    if (Iswall0 .or. Iswallm) RmiuTL=0.0  !wall diff

	    TxxLam=RmiuL*(NMD*(DuDxL+DvDyL+DwDzL)+2.0*DuDxL)   !txx=2.0*vmuc*dudx+vlac*div
	    TyyLam=RmiuL*(NMD*(DuDxL+DvDyL+DwDzL)+2.0*DvDyL)   !tyy=2.0*vmuc*dvdy+vlac*div
	    TzzLam=RmiuL*(NMD*(DuDxL+DvDyL+DwDzL)+2.0*DwDzL)   !tzz=2.0*vmuc*dwdz+vlac*div
	    TxyLam=RmiuL*(     DuDyL+DvDxL)                   !txy=vmuc*(dudy+dvdx)
	    TyzLam=RmiuL*(     DvDzL+DwDyL)                   !tyz=vmuc*(dvdz+dwdy)
	    TzxLam=RmiuL*(     DuDzL+DwDxL)                   !tzx=vmuc*(dwdx+dudz)

	    TxxTur=RmiuTL*(NMD*(DuDxL+DvDyL+DwDzL)+2.0*DuDxL)
	    TyyTur=RmiuTL*(NMD*(DuDxL+DvDyL+DwDzL)+2.0*DvDyL)
	    TzzTur=RmiuTL*(NMD*(DuDxL+DvDyL+DwDzL)+2.0*DwDzL)
	    TxyTur=RmiuTL*(     DuDyL+DvDxL)
	    TyzTur=RmiuTL*(     DvDzL+DwDyL)
		TzxTur=RmiuTL*(     DuDzL+DwDxL)


	    Txx=TxxLam+TxxTur
	    Tyy=TyyLam+TyyTur
	    Tzz=TzzLam+TzzTur
	    Txy=TxyLam+TxyTur
	    Tyz=TyzLam+TyzTur
	    Tzx=TzxLam+TzxTur
        
! if(ThisBlock%ID_Present_Blk==3.and.KI>8)write(*,*)"KDUFF=",KI,i,j,k,DuDxL,DuDyL,DuDzL,Txx,Txy
	    cp= 1.4*PPF/(1.4-1.0)   !gam(I,J,K)*PPF/(gam(I,J,K)-1.0)
	    akmu=cp*(RmiuL/prl+RmiuTL/prte)

          ucc=VL(3,1,i,j,k)  !ucc=ukl(i,j,k)
          vcc=VL(3,2,i,j,k)  !vcc=vkl(i,j,k)
          wcc=VL(3,3,i,j,k)  !wcc=wkl(i,j,k)

          Phix=akmu*dTdxL+ucc*Txx+vcc*Txy+wcc*Tzx
          Phiy=akmu*dTdyL+ucc*Txy+vcc*Tyy+wcc*Tyz
          Phiz=akmu*dTdzL+ucc*Tzx+vcc*Tyz+wcc*Tzz
        
          Dx=SD(3,1,i,j,k)
          Dy=SD(3,2,i,j,k)
          Dz=SD(3,3,i,j,k)

	    DIFFrok =   0.0                                             !drok= 0.
	    DIFFXk =(  Txx * Dx +  Txy * Dy +  Tzx * Dz ) / Ref        !dxk = akx(i  ,j,k)*txx+aky(i  ,j,k)*txy+akz(i  ,j,k)*tzx
	    DIFFYk =(  Txy * Dx +  Tyy * Dy +  Tyz * Dz ) / Ref        !dyk = akx(i  ,j,k)*txy+aky(i  ,j,k)*tyy+akz(i  ,j,k)*tyz
	    DIFFZk =(  Tzx * Dx +  Tyz * Dy +  Tzz * Dz ) / Ref        !dzk = akx(i  ,j,k)*tzx+aky(i  ,j,k)*tyz+akz(i  ,j,k)*tzz
	    DIFFEk =( Phix * Dx + Phiy * Dy + Phiz * Dz ) / Ref        !drek= akx(i  ,j,k)*btx+aky(i  ,j,k)*bty+akz(i  ,j,k)*btz


! if(ThisBlock%ID_Present_Blk==4)write(*,*)"DIFFUSIOND12345=",KI,i,j,k,D(1:5,i,j,k),Txx,Txy,Tzx,Dx,Dy,Dz
	    D(1,I,J,K)=D(1,I,J,K)-DIFFrok
          D(1,i,j,k-1)=D(1,i,j,k-1)+DIFFrok
	
	    D(2,I,J,K)=D(2,I,J,K)-DIFFXk
          D(2,i,j,k-1)=D(2,i,j,k-1)+DIFFXk

	    D(3,I,J,K)=D(3,I,J,K)-DIFFYk
          D(3,i,j,k-1)=D(3,i,j,k-1)+DIFFYk

	    D(4,I,J,K)=D(4,I,J,K)-DIFFZk
          D(4,i,j,k-1)=D(4,i,j,k-1)+DIFFZk

	    D(5,I,J,K)=D(5,I,J,K)-DIFFEk
          D(5,i,j,k-1)=D(5,i,j,k-1)+DIFFEk
          
        
          
!if(ThisBlock%ID_Present_Blk==4) write(*,*)"DIFF",KI,i,j,k,D(1:2,i,j,k),DIFFrok,DIFFXk,DIFFYk,DIFFZk,DIFFEK
        
	end do
	end do
	end do
	
!---------------!ydd add rotating source terms 20150909-------------
    do k=1,NK1
    do j=1,NJ1
    do i=1,NI1
        RoW=omega(1)
        ss3=2.0*Row*V(4,i,j,k)+Row*Yc(i,j,k)*omega(1)*V(1,i,j,k)
        ss4=-2.0*Row*V(3,i,j,k)+Row*Zc(i,j,k)*omega(1)*V(1,i,j,k)
        D(3,i,j,k)=D(3,i,j,k)+ss3*Vol(i,j,k)
        D(4,i,j,k)=D(4,i,j,k)+ss4*Vol(i,j,k)
    enddo
    enddo
    enddo
!----------------------------------------------------------------------------          



END SUBROUTINE DIFFUSION_Centered
