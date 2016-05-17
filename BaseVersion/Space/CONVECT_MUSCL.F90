SUBROUTINE CONVECT_MUSCL

    USE Global
    IMPLICIT NONE
    
    INTEGER:: I,J,K,L
    INTEGER:: NI, NJ, NK
    INTEGER:: NI1, NJ1, NK1
    
    REAL:: UU,VV,WW,PP0,TT
    REAL:: DR,DL
    real::rrr

    INTERFACE
        REAL FUNCTION alimiter(var1,var2)  !muscl3

        END FUNCTION
    END INTERFACE
    
    NI=ThisBlock%NI
    NJ=ThisBlock%NJ
    NK=ThisBlock%NK
    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1



!        VL(3,ML,-2:MI,-2:MJ,-2:MK)  ! U V W PP T k omega muT gm         DO i=0,NI+1
       DO i=0,NI+1
       DO j=0,NJ+1
	   DO k=0,NK+1
	   DO L=1,5
		  VR(1,L,i,j,k)=VL(1,L,i,j,k)
		  VR(2,L,i,j,k)=VL(2,L,i,j,k)
		  VR(3,L,i,j,k)=VL(3,L,i,j,k)
	   end do
	   end do
	   end do
       end do


!        i direction
         Do k=1,NK1
         Do j=1,NJ1
         Do i=0,NI

!        limite U

         UU=V(2,i,j,k)/V(1,i,j,k)


         DR=VL(1,1,i+1,j,k)-UU                 ! dr=uil(i+1,j,k)-  vx(i,j,k)
         DL=UU-VR(1,1,i,j,k)                   ! dl=vx(i,j,k)    -uir(i,j,k)

         VL(1,1,i+1,j,k)=UU+alimiter(DR,DL)    ! uil(i+1,j,k)=vx(i,j,k)+alimiter(dr,dl)

	   VR(1,1,i,j,k)=UU-alimiter(DL,DR)      ! uir(i,j,k)  =vx(i,j,k)-alimiter(dl,dr)


!        limite V

         VV=V(3,i,j,k)/V(1,i,j,k)


         DR=VL(1,2,i+1,j,k)-VV                 ! dr=vil(i+1,j,k)-  vy(i,j,k)
         DL=VV-VR(1,2,i,j,k)                   ! dl=vy(i,j,k)    -vir(i,j,k) 

         VL(1,2,i+1,j,k)=VV+alimiter(DR,DL)    ! vil(i+1,j,k)=vy(i,j,k)+alimiter(dr,dl)

	   VR(1,2,i,j,k)=VV-alimiter(DL,DR)      ! vir(i,j,k)  =vy(i,j,k)-alimiter(dl,dr)


!        limite W

         WW=V(4,i,j,k)/V(1,i,j,k)


         DR=VL(1,3,i+1,j,k)-WW                 ! dr=(wil(i+1,j,k)-  vz(i,j,k))
         DL=WW-VR(1,3,i,j,k)                   ! dl=(vz(i,j,k)    -wir(i,j,k))

         VL(1,3,i+1,j,k)=WW+alimiter(DR,DL)    ! wil(i+1,j,k)=vz(i,j,k)+alimiter(dr,dl)

	   VR(1,3,i,j,k)=WW-alimiter(DL,DR)      ! wir(i,j,k)  =vz(i,j,k)-alimiter(dl,dr)



!         limite P


         PP0=PP(i,j,k)

         DR=VL(1,4,i+1,j,k)-PP0                 ! dr=(pil(i+1,j,k)-  p(i,j,k))
         DL=PP0-VR(1,4,i,j,k)                   ! dl=(p(i,j,k)    -pir(i,j,k))

         VL(1,4,i+1,j,k)=PP0+alimiter(DR,DL)    ! pil(i+1,j,k)=p(i,j,k)+alimiter(dr,dl)

	   VR(1,4,i,j,k)=PP0-alimiter(DL,DR)      ! pir(i,j,k)  =p(i,j,k)-alimiter(dl,dr)


         !!!!!!!!!!!!!!!!!!!!!!!!!!
         TT=T(i,j,k)


         DR=VL(1,5,i+1,j,k)-TT                 ! dr=(til(i+1,j,k)-  tp(i,j,k))
         DL=TT-VR(1,5,i,j,k)                   ! dl=(tp(i,j,k)    -tir(i,j,k))

         VL(1,5,i+1,j,k)=TT+alimiter(DR,DL)    ! til(i+1,j,k)=tp(i,j,k)+alimiter(dr,dl)

	   VR(1,5,i,j,k)=TT-alimiter(DL,DR)      ! tir(i,j,k)  =tp(i,j,k)-alimiter(dl,dr)
    !limit rcc, by ydd
        rrr=rad(i,j,k)
        DR=rccl(1,i+1,j,k)-rrr
        DL=rrr-rccr(1,i,j,k)
        rccl(1,i+1,j,k)=rrr+alimiter(DR,DL)
        rccr(1,i,j,k)=rrr-alimiter(DL,DR)
	   end do
	   end do
	   end do



!        j direction

	   Do k=1,NK1
         Do j=0,NJ
         Do i=1,NI1

!        limite U

         UU=V(2,i,j,k)/V(1,i,j,k)


         DR=VL(2,1,i,j+1,k)-UU                 !  dr=ujl(i,j+1,k)-  vx(i,j,k)
         DL=UU-VR(2,1,i,j,k)                   !  dl=vx(i,j,k)    -ujr(i,j,k))

         VL(2,1,i,j+1,k)=UU+alimiter(DR,DL)    ! ujl(i,j+1,k)=vx(i,j,k)+alimiter(dr,dl)

	   VR(2,1,i,j,k)=UU-alimiter(DL,DR)      ! ujr(i,j,k)  =vx(i,j,k)-alimiter(dl,dr)


!        limite V

         VV=V(3,i,j,k)/V(1,i,j,k)


         DR=VL(2,2,i,j+1,k)-VV                 !  dr=(vjl(i,j+1,k)-  vy(i,j,k))
         DL=VV-VR(2,2,i,j,k)                   !  dl=(vy(i,j,k)    -vjr(i,j,k))

         VL(2,2,i,j+1,k)=VV+alimiter(DR,DL)    ! vjl(i,j+1,k)=vy(i,j,k)+alimiter(dr,dl)

	   VR(2,2,i,j,k)=VV-alimiter(DL,DR)      ! vjr(i,j,k)  =vy(i,j,k)-alimiter(dl,dr)


!        limite W

         WW=V(4,i,j,k)/V(1,i,j,k)


         DR=VL(2,3,i,j+1,k)-WW                 ! dr=(wjl(i,j+1,k)-  vz(i,j,k))
         DL=WW-VR(2,3,i,j,k)                   ! dl=(vz(i,j,k)    -wjr(i,j,k))

         VL(2,3,i,j+1,k)=WW+alimiter(DR,DL)    ! wjl(i,j+1,k)=vz(i,j,k)+alimiter(dr,dl)

	   VR(2,3,i,j,k)=WW-alimiter(DL,DR)      ! wjr(i,j,k)  =vz(i,j,k)-alimiter(dl,dr)

    
!         limite P


         PP0=PP(i,j,k)
         

         DR=VL(2,4,i,j+1,k)-PP0                 ! dr=(pjl(i,j+1,k)-  p(i,j,k))
         DL=PP0-VR(2,4,i,j,k)                   ! dl=(p(i,j,k)    -pjr(i,j,k))

         VL(2,4,i,j+1,k)=PP0+alimiter(DR,DL)    ! pjl(i,j+1,k)=p(i,j,k)+alimiter(dr,dl)

	   VR(2,4,i,j,k)=PP0-alimiter(DL,DR)      ! pjr(i,j,k)  =p(i,j,k)-alimiter(dl,dr)


         !!!!!!!!!!!!!!!!!!!!!!!!!!
         TT=T(i,j,k)


         DR=VL(2,5,i,j+1,k)-TT                 ! dr=(tjl(i,j+1,k)-  tp(i,j,k))
         DL=TT-VR(2,5,i,j,k)                   ! dl=(tp(i,j,k)    -tjr(i,j,k))

         VL(2,5,i,j+1,k)=TT+alimiter(DR,DL)    ! tjr(i,j,k)  =tp(i,j,k)-alimiter(dl,dr)

	   VR(2,5,i,j,k)=TT-alimiter(DL,DR)      ! tjl(i,j+1,k)=tp(i,j,k)+alimiter(dr,dl)
!limit rcc
        rrr=rad(i,j,k)
        DR=rccl(2,i,j+1,k)-rrr
        DL=rrr-rccr(2,i,j,k)
        rccl(2,i,j+1,k)=rrr+alimiter(DR,DL)
        rccr(2,i,j,k)=rrr-alimiter(DL,DR)
	   end do
	   end do
	   end do


!        k direction
         Do k=0,NK
         Do j=1,NJ1
         Do i=1,NI1

!        limite U

         UU=V(2,i,j,k)/V(1,i,j,k)


         DR=VL(3,1,i,j,k+1)-UU                 ! dr=ukl(i,j,k+1)-  vx(i,j,k)
         DL=UU-VR(3,1,i,j,k)                   ! dl=vx(i,j,k)    -ukr(i,j,k)

         VL(3,1,i,j,k+1)=UU+alimiter(DR,DL)    ! ukr(i,j,k)  =vx(i,j,k)-alimiter(dl,dr)

	   VR(3,1,i,j,k)=UU-alimiter(DL,DR)      ! ukl(i,j,k+1)=vx(i,j,k)+alimiter(dr,dl)


!        limite V

         VV=V(3,i,j,k)/V(1,i,j,k)


         DR=VL(3,2,i,j,k+1)-VV                 ! dr=(vkl(i,j,k+1)-  vy(i,j,k))
         DL=VV-VR(3,2,i,j,k)                   ! dl=(vy(i,j,k)    -vkr(i,j,k))

         VL(3,2,i,j,k+1)=VV+alimiter(DR,DL)    ! vkr(i,j,k)  =vy(i,j,k)-alimiter(dl,dr)

	   VR(3,2,i,j,k)=VV-alimiter(DL,DR)      ! vkl(i,j,k+1)=vy(i,j,k)+alimiter(dr,dl)


!        limite W

         WW=V(4,i,j,k)/V(1,i,j,k)


         DR=VL(3,3,i,j,k+1)-WW                 ! dr=uil(i+1,j,k)-  vx(i,j,k)
         DL=WW-VR(3,3,i,j,k)                   ! dl=vx(i,j,k)    -uir(i,j,k)

         VL(3,3,i,j,k+1)=WW+alimiter(DR,DL)    ! uil(i+1,j,k)=vx(i,j,k)+alimiter(dr,dl)

	   VR(3,3,i,j,k)=WW-alimiter(DL,DR)      ! uir(i,j,k)  =vx(i,j,k)-alimiter(dl,dr)



!         limite P


         PP0=PP(i,j,k)

         DR=VL(3,4,i,j,k+1)-PP0                 ! dr=uil(i+1,j,k)-  vx(i,j,k)
         DL=PP0-VR(3,4,i,j,k)                   ! dl=vx(i,j,k)    -uir(i,j,k)

         VL(3,4,i,j,k+1)=PP0+alimiter(DR,DL)    ! uil(i+1,j,k)=vx(i,j,k)+alimiter(dr,dl)

	   VR(3,4,i,j,k)=PP0-alimiter(DL,DR)      ! uir(i,j,k)  =vx(i,j,k)-alimiter(dl,dr)


         !!!!!!!!!!!!!!!!!!!!!!!!!!
         TT=T(i,j,k)


         DR=VL(3,5,i,j,k+1)-TT                 ! dr=uil(i+1,j,k)-  vx(i,j,k)
         DL=TT-VR(3,5,i,j,k)                   ! dl=vx(i,j,k)    -uir(i,j,k)

         VL(3,5,i,j,k+1)=TT+alimiter(DR,DL)    ! uil(i+1,j,k)=vx(i,j,k)+alimiter(dr,dl)

	   VR(3,5,i,j,k)=TT-alimiter(DL,DR)      ! uir(i,j,k)  =vx(i,j,k)-alimiter(dl,dr)
    !limit rcc, by ydd
        rrr=rad(i,j,k)
        DR=rccl(3,i,j,k+1)-rrr
        DL=rrr-rccr(3,i,j,k)
        rccl(3,i,j,k+1)=rrr+alimiter(DR,DL)
        rccr(3,i,j,k)=rrr-alimiter(DL,DR)
!if(ThisBlock%ID_Present_Blk==2) write(*,*)"VLR",i,j,k,VL(3,1:3,i,j,k),VR(3,1:3,i,j,k)

	   end do
	   end do
	   end do


END SUBROUTINE CONVECT_MUSCL
