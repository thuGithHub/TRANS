SUBROUTINE GEOM_CalcDSTs
    USE Global
    IMPLICIT NONE

    INTEGER iBlock
    TYPE(BlockStruct), POINTER:: Block

    DO iBlock = 1, Max_Block

        Block => AllBlocks(iBlock)
        CALL GLOBAL_SetPointers(iBlock)
        ! 赋值文件名变量Block%FLgeo为grd/grid_blk***.geo
        ! 赋值文件名变量Block%FLdst为grd/grid_blk***.dst
        WRITE(Block%FLgeo, '(A,A, I3.3, A)') trim(FilePathPrefix), "grd/grid_blk", Block%ID_Present_Blk, ".geo"
                WRITE(Block%FLdst, '(A,A, I3.3, A)') trim(FilePathPrefix), "grd/grid_blk", Block%ID_Present_Blk, ".dst"
                !WRITE(Block%FLgeo, '(A,A, I3.3, A)') trim(FilePathPrefix), "grd/grid_blk", iBlock, ".geo"
                !WRITE(Block%FLdst, '(A,A, I3.3, A)') trim(FilePathPrefix), "grd/grid_blk", iBlock, ".dst"

                CALL GEOM_CalcDST(iBlock)      !网格几何量计算
                CALL GEOM_OutputGEO(iBlock)    !网格几何量输出

            ENDDO
                  
        END SUBROUTINE GEOM_CalcDSTs



        SUBROUTINE GEOM_CalcDST(iBlock)
            USE Global
            IMPLICIT NONE
            
            INTEGER:: iBlock
            TYPE(BlockStruct), POINTER:: Block

            
            INTEGER:: I,J,K,I1,J1,K1,II,JJ,KK,L
            INTEGER:: N,NI,NJ,NK,NI1,NJ1,NK1,NI2,NJ2,NK2
            REAL:: Xcnt, Ycnt, Zcnt

            REAL:: Dst_minN !(Nblk_sur)       
            REAL:: dst_minimum                
            
            REAL:: Xsurf, Ysurf, Zsurf         
            REAL:: Dx, Dy, Dz
            REAL:: DDD
            
            REAL:: dl,xv,yv,zv
            REAL:: xsv,ysv,zsv
            REAL:: R1X,R1Y,R1Z,R2X,R2Y,R2Z
            REAL:: RS
            REAL:: X1,Y1,Z1,X2,Y2,Z2,T1,T2

            REAL:: Xcs, Ycs, Zcs

            REAL:: DxW,DyW,DzW
            REAL:: DxE,DyE,DzE
            REAL:: DxS,DyS,DzS
            REAL:: DxN,DyN,DzN
            REAL:: DxL,DyL,DzL
            REAL:: DxU,DyU,DzU
            REAL:: DDw,DDe,DDs,DDn,DDl,DDu
            REAL:: DLTxyz, DSTxyz

            REAL:: C13,XR,YR,ZR,OMGv
            
            REAL:: Xw,Yw,Zw
            REAL:: Xe,Ye,Ze
            REAL:: Xs,Ys,Zs
            REAL:: Xn,Yn,Zn
            REAL:: Xb,Yb,Zb
            REAL:: Xt,Yt,Zt
            
            real::txc,tyc,tzc,xVec,yVec,zVec
            real::W(3)
            real::normalS,velocityS


            Block => AllBlocks(iBlock)

            NI=Block%NI
            NJ=Block%NJ
            NK=Block%NK
            NI1=Block%NI1
            NJ1=Block%NJ1
            NK1=Block%NK1
            NI2=Block%NI2
            NJ2=Block%NJ2
            NK2=Block%NK2

            !added by ydd
            W(1:3)=omega(1:3)
            
            !!!!!!
            !	  XYZc
            !!!!!!
                !计算网格中心点坐标
                DO I=1,NI1
                   I1=I+1
                   DO J=1,NJ1
                      J1=J+1
                      DO K=1,NK1
                         K1=K+1
                         Zc(I,J,K)=(ZZ(I1,J1,K1)+ZZ(I1,J1,K)+ZZ(I1,J,K1) &
             &                     +ZZ(I1,J,K)  +ZZ(I,J1,K1)+ZZ(I,J1,K) &
             &                     +ZZ(I,J,K1)  +ZZ(I,J,K))/8.
                         Yc(I,J,K)=(YY(I1,J1,K1)+YY(I1,J1,K)+YY(I1,J,K1) &
             &                     +YY(I1,J,K)  +YY(I,J1,K1)+YY(I,J1,K) &
             &                     +YY(I,J,K1)  +YY(I,J,K))/8.
                         Xc(I,J,K)=(XX(I1,J1,K1)+XX(I1,J1,K)+XX(I1,J,K1) &
             &                     +XX(I1,J,K)  +XX(I,J1,K1)+XX(I,J1,K) &
             &                     +XX(I,J,K1)  +XX(I,J,K))/8.

                    thtc(I,J,K)=atan2(Zc(I,j,k),Yc(i,j,k))
                    rad(i,j,k)=sqrt(Yc(i,j,k)**2.0+Zc(i,j,k)**2.0)
                        ENDDO
                     ENDDO
                  ENDDO




            !!!!!!
            !	  minDST
            !!!!!!
                  IF(IF_dst == 0)THEN

                     DO KK=1,NK1
                        WRITE(*,*) "K lyaer", KK
                        !! MYID write(*,*)"K layer",KK, "myid is ", Block%myid
                     DO JJ=1,NJ1
                     DO II=1,NI1

                        Xcnt=Xc(II,JJ,KK)
                        Ycnt=Yc(II,JJ,KK)
                        Zcnt=Zc(II,JJ,KK)

                                dst_minimum = 1.E20
                    ! 网格中心点到物面距离
                            do N=1,Num_surface
                           Dst_minN = 1.e+20 !(N)=1.E+20
                           DO J=1, AllSurfs(N)%MSur-1  !Index2(N)-1
                               DO I=1, AllSurfs(N)%NSur-1  !Index1(N)-1
                                ! 物面网格中心点坐标
                            Xsurf=( AllSurfs(N)%X(I,J) + AllSurfs(N)%X(I+1,J) + AllSurfs(N)%X(I,J+1) + AllSurfs(N)%X(I+1,J+1))/4.    !(Xsur(N,I,J)+Xsur(N,I+1,J)+Xsur(N,I,J+1)+Xsur(N,I+1,J+1))/4.
                                Ysurf=( AllSurfs(N)%Y(I,J) + AllSurfs(N)%Y(I+1,J) + AllSurfs(N)%Y(I,J+1) + AllSurfs(N)%Y(I+1,J+1))/4. !(Ysur(N,I,J)+Ysur(N,I+1,J)+Ysur(N,I,J+1)+Ysur(N,I+1,J+1))/4.
                                Zsurf=( AllSurfs(N)%Z(I,J) + AllSurfs(N)%Z(I+1,J) + AllSurfs(N)%Z(I,J+1) + AllSurfs(N)%Z(I+1,J+1))/4. !(Zsur(N,I,J)+Zsur(N,I+1,J)+Zsur(N,I,J+1)+Zsur(N,I+1,J+1))/4.

                              Dx=Xcnt-Xsurf
                              Dy=Ycnt-Ysurf
                              Dz=Zcnt-Zsurf

                              DDD=Dx*Dx+Dy*Dy+Dz*Dz
                              if(DDD < Dst_minN )Dst_minN = DDD
                           ENDDO
                           ENDDO

                                  dst_minimum = amax1(dst_minimum, tiny)
                           if(Dst_minN < dst_minimum ) dst_minimum = Dst_minN
                        enddo
                    ! 第II,JJ,KK网格到物面距离
                        Dst(II,JJ,KK)= sqrt( dst_minimum)
                     ENDDO
                     ENDDO
                     ENDDO
               
                !SAVE DST
                OPEN(UNIT=14, FILE=Block%FLdst, MODE='WRITE', FORM='UNFORMATTED') !, SHARE='DENYWR')
!                OPEN(UNIT=14, FILE=Block%FLdst)
                     DO K=1,NK1
                     DO J=1,NJ1
                 DO I=1,NI1
                      WRITE(14) Dst(I,J,K)
                     ENDDO
                 ENDDO
                     ENDDO
                     CLOSE(14)
            !!!!!!
            !	  Read DST
            ! 已经计算过物面距离读入数据 
                  ELSE
                OPEN(UNIT=14, FILE=Block%FLdst, MODE='READ', FORM='UNFORMATTED') !, SHARE='DENYWR')
!                OPEN(UNIT=14, FILE=Block%FLdst)
                     DO K=1,NK1
                     DO J=1,NJ1
                 DO I=1,NI1
                      READ(14)Dst(I,J,K)
                     ENDDO
                     ENDDO
                     ENDDO
                 CLOSE(14)
                  ENDIF


            !!!!!! SD & Grad
            !????NEEDED?	  Delta=delta_BL
            !
            !!!!!!
            
                dl=1.
                  xv=xc(2,1,1)-xc(1,1,1)
                  yv=yc(2,1,1)-yc(1,1,1)
                  zv=zc(2,1,1)-zc(1,1,1)
                ! I面面积矢量、面积计算
                DO K=1,NK1
                   K1=K+1
                   DO J=1,NJ1
                      J1=J+1
                      DO I=1,NI
                         ! R1、R2分别为I面两条对角线矢量
                         R1X=XX(I,J1,K1)-XX(I,J,K)
                         R1Y=YY(I,J1,K1)-YY(I,J,K)
                         R1Z=ZZ(I,J1,K1)-ZZ(I,J,K)
                         R2X=XX(I,J,K1) -XX(I,J1,K)
                         R2Y=YY(I,J,K1) -YY(I,J1,K)
                         R2Z=ZZ(I,J,K1) -ZZ(I,J1,K)

                           if(i == 1.and.j == 1.and.k == 1) then
                           xsv=(R1Y*R2Z-R1Z*R2Y)*0.5
                           ysv=(R1Z*R2X-R1X*R2Z)*0.5
                           zsv=(R1X*R2Y-R1Y*R2X)*0.5
                           if(xv*xsv+yv*ysv+zv*zsv < 0) dl=-1.
                                !	 dl=-1.
                             if(dl == -1.) print*,'i reversed'
                           endif
                         ! SD为网格面积矢量
                         SD(1,1,I,J,K)=(R1Y*R2Z-R1Z*R2Y)*0.5*dl
                         SD(1,2,I,J,K)=(R1Z*R2X-R1X*R2Z)*0.5*dl
                         SD(1,3,I,J,K)=(R1X*R2Y-R1Y*R2X)*0.5*dl
            !                 SD(1,4,I,J,K)=0.
                        ! Grad 为网格面积
                           Grad(1,I,J,K)=SQRT(SD(1,1,I,J,K)*SD(1,1,I,J,K)+ &
             &                              SD(1,2,I,J,K)*SD(1,2,I,J,K)+ &
             &                              SD(1,3,I,J,K)*SD(1,3,I,J,K)+tiny)
                        ENDDO
                        SD(1,1,0,J,K)=SD(1,1,1,J,K)
                        SD(1,2,0,J,K)=SD(1,2,1,J,K)
                        SD(1,3,0,J,K)=SD(1,3,1,J,K)
                        Grad(1,0,J,K)=Grad(1,1,J,K)
                     ENDDO
                  ENDDO
                ! J面面积矢量、面积计算
                DO K=1,NK1
                   K1=K+1
                   DO I=1,NI1
                      I1=I+1
                      DO J=1,NJ
                         R1X=XX(I,J,K1)-XX(I1,J,K)
                         R1Y=YY(I,J,K1)-YY(I1,J,K)
                         R1Z=ZZ(I,J,K1)-ZZ(I1,J,K)
                         R2X=XX(I1,J,K1)-XX(I,J,K)
                         R2Y=YY(I1,J,K1)-YY(I,J,K)
                         R2Z=ZZ(I1,J,K1)-ZZ(I,J,K)

                           if(i == 1.and.j == 1.and.k == 1) then
                           xsv=(R1Y*R2Z-R1Z*R2Y)*0.5
                           ysv=(R1Z*R2X-R1X*R2Z)*0.5
                           zsv=(R1X*R2Y-R1Y*R2X)*0.5
                             if(dl == -1.) print*,'j reversed'
                           endif

                         SD(2,1,I,J,K)=(R1Y*R2Z-R1Z*R2Y)*0.5*dl
                         SD(2,2,I,J,K)=(R1Z*R2X-R1X*R2Z)*0.5*dl
                         SD(2,3,I,J,K)=(R1X*R2Y-R1Y*R2X)*0.5*dl
            !                 SD(2,4,I,J,K)=0.

                           Grad(2,I,J,K)=SQRT(SD(2,1,I,J,K)*SD(2,1,I,J,K)+ &
             &                              SD(2,2,I,J,K)*SD(2,2,I,J,K)+ &
             &                              SD(2,3,I,J,K)*SD(2,3,I,J,K)+tiny)
                        ENDDO
                        SD(2,1,I,0,K)=SD(2,1,I,1,K)
                        SD(2,2,I,0,K)=SD(2,2,I,1,K)
                        SD(2,3,I,0,K)=SD(2,3,I,1,K)
                        Grad(2,I,0,K)=Grad(2,I,1,K)
                     ENDDO
                  ENDDO
                ! K面面积矢量、面积计算
                DO J=1,NJ1
                   J1=J+1
                   DO I=1,NI1
                      I1=I+1
                      DO K=1,NK
                         R1X=XX(I1,J,K)-XX(I,J1,K)
                         R1Y=YY(I1,J,K)-YY(I,J1,K)
                         R1Z=ZZ(I1,J,K)-ZZ(I,J1,K)
                         R2X=XX(I1,J1,K)-XX(I,J,K)
                         R2Y=YY(I1,J1,K)-YY(I,J,K)
                         R2Z=ZZ(I1,J1,K)-ZZ(I,J,K)

                           if(i == 1.and.j == 1.and.k == 1) then
                           xsv=(R1Y*R2Z-R1Z*R2Y)*0.5
                           ysv=(R1Z*R2X-R1X*R2Z)*0.5
                           zsv=(R1X*R2Y-R1Y*R2X)*0.5
                             if(dl == -1.) print*,'k reversed'
                           endif

                         SD(3,1,I,J,K)=(R1Y*R2Z-R1Z*R2Y)*0.5*dl
                         SD(3,2,I,J,K)=(R1Z*R2X-R1X*R2Z)*0.5*dl
                         SD(3,3,I,J,K)=(R1X*R2Y-R1Y*R2X)*0.5*dl
!if(ThisBlock%ID_Present_blk==6.and.K==17.and.J==64.and.(i>38.and.i<48))then
!write(*,*)ThisBlock%ID_Present_blk,SD(3,:,i,j,k),R1x,R1y,R1z,R2x,R2y,R2z
!endif
            !                 SD(3,4,I,J,K)=0.

                           Grad(3,I,J,K)=SQRT(SD(3,1,I,J,K)*SD(3,1,I,J,K)+ &
             &                              SD(3,2,I,J,K)*SD(3,2,I,J,K)+ &
             &                              SD(3,3,I,J,K)*SD(3,3,I,J,K)+tiny)
                        ENDDO
                        SD(3,1,I,J,0)=SD(3,1,I,J,1)
                        SD(3,2,I,J,0)=SD(3,2,I,J,1)
                        SD(3,3,I,J,0)=SD(3,3,I,J,1)
                        Grad(3,I,J,0)=Grad(3,I,J,1)
                     ENDDO
                  ENDDO







        !!!!!!!!!!!              by xu interpolation coefficient
            DO K=1,NK1
                DO J=1,NJ1
            DO I=1,NI1
              ! I面 几何中心坐标
              xw=0.25_8*(xx(i,j,k)+xx(i,j+1,k)+xx(i,j,k+1)+xx(i,j+1,k+1))
              yw=0.25_8*(yy(i,j,k)+yy(i,j+1,k)+yy(i,j,k+1)+yy(i,j+1,k+1))
              zw=0.25_8*(zz(i,j,k)+zz(i,j+1,k)+zz(i,j,k+1)+zz(i,j+1,k+1))
              ! I+1面几何中心坐标
              xe=0.25_8*(xx(i+1,j,k)+xx(i+1,j+1,k)+xx(i+1,j,k+1)+xx(i+1,j+1,k+1))
              ye=0.25_8*(yy(i+1,j,k)+yy(i+1,j+1,k)+yy(i+1,j,k+1)+yy(i+1,j+1,k+1))
              ze=0.25_8*(zz(i+1,j,k)+zz(i+1,j+1,k)+zz(i+1,j,k+1)+zz(i+1,j+1,k+1))
              ! I+1面和I面几何中心距离
              Alagm(1,i,j,k)=sqrt((xw-xe)**2.0_8+(yw-ye)**2.0_8+(zw-ze)**2.0_8)
              ! J面 几何中心坐标
              xs=0.25_8*(xx(i,j,k)+xx(i+1,j,k)+xx(i,j,k+1)+xx(i+1,j,k+1))
              ys=0.25_8*(yy(i,j,k)+yy(i+1,j,k)+yy(i,j,k+1)+yy(i+1,j,k+1))
              zs=0.25_8*(zz(i,j,k)+zz(i+1,j,k)+zz(i,j,k+1)+zz(i+1,j,k+1))
              ! J+1面几何中心坐标
              xn=0.25_8*(xx(i,j+1,k)+xx(i+1,j+1,k)+xx(i,j+1,k+1)+xx(i+1,j+1,k+1)) 
              yn=0.25_8*(yy(i,j+1,k)+yy(i+1,j+1,k)+yy(i,j+1,k+1)+yy(i+1,j+1,k+1)) 
              zn=0.25_8*(zz(i,j+1,k)+zz(i+1,j+1,k)+zz(i,j+1,k+1)+zz(i+1,j+1,k+1)) 
              ! J+1面和I面几何中心距离
              Alagm(2,i,j,k)=sqrt((xs-xn)**2.0_8+(ys-yn)**2.0_8+(zs-zn)**2.0_8)
              ! K面 几何中心坐标
              xb=0.25_8*(xx(i,j,k)+xx(i+1,j,k)+xx(i,j+1,k)+xx(i+1,j+1,k)) 
              yb=0.25_8*(yy(i,j,k)+yy(i+1,j,k)+yy(i,j+1,k)+yy(i+1,j+1,k)) 
              zb=0.25_8*(zz(i,j,k)+zz(i+1,j,k)+zz(i,j+1,k)+zz(i+1,j+1,k)) 
              ! K+1面几何中心坐标
              xt=0.25_8*(xx(i,j,k+1)+xx(i+1,j,k+1)+xx(i,j+1,k+1)+xx(i+1,j+1,k+1))
              yt=0.25_8*(yy(i,j,k+1)+yy(i+1,j,k+1)+yy(i,j+1,k+1)+yy(i+1,j+1,k+1))
              zt=0.25_8*(zz(i,j,k+1)+zz(i+1,j,k+1)+zz(i,j+1,k+1)+zz(i+1,j+1,k+1))
               ! K+1面和K面几何中心距离
              Alagm(3,i,j,k)=sqrt((xb-xt)**2.0_8+(yb-yt)**2.0_8+(zb-zt)**2.0_8)
              
            ENDDO
                ENDDO
                ENDDO







            !!!!!!
            !	  BC Grids
            !!!!!!
            ! 网格块边界面几何中心计算
                  DO K=1,NK1    !! I direction
                  DO J=1,NJ1
                   Xcs=(XX(1,J,K)+XX(1,J+1,K)+XX(1,J,K+1)+XX(1,J+1,K+1))/4.
                   Ycs=(YY(1,J,K)+YY(1,J+1,K)+YY(1,J,K+1)+YY(1,J+1,K+1))/4.
                   Zcs=(ZZ(1,J,K)+ZZ(1,J+1,K)+ZZ(1,J,K+1)+ZZ(1,J+1,K+1))/4.
                     Xc(0,J,K)=Xcs
                     Yc(0,J,K)=Ycs
                     Zc(0,J,K)=Zcs

                   Xcs=(XX(NI,J,K)+XX(NI,J+1,K)+XX(NI,J,K+1)+XX(NI,J+1,K+1))/4.
                   Ycs=(YY(NI,J,K)+YY(NI,J+1,K)+YY(NI,J,K+1)+YY(NI,J+1,K+1))/4.
                   Zcs=(ZZ(NI,J,K)+ZZ(NI,J+1,K)+ZZ(NI,J,K+1)+ZZ(NI,J+1,K+1))/4.
                     Xc(NI,J,K)=Xcs
                     Yc(NI,J,K)=Ycs
                     Zc(NI,J,K)=Zcs
                  ENDDO
                  ENDDO

                  DO K=1,NK1    !! J direction
                  DO I=1,NI1
                   Xcs=(XX(I,1,K)+XX(I+1,1,K)+XX(I,1,K+1)+XX(I+1,1,K+1))/4.
                   Ycs=(YY(I,1,K)+YY(I+1,1,K)+YY(I,1,K+1)+YY(I+1,1,K+1))/4.
                   Zcs=(ZZ(I,1,K)+ZZ(I+1,1,K)+ZZ(I,1,K+1)+ZZ(I+1,1,K+1))/4.
                     Xc(i,0,K)=Xcs
                     Yc(i,0,K)=Ycs
                     Zc(i,0,K)=Zcs

                   Xcs=(XX(I,NJ,K)+XX(I+1,NJ,K)+XX(I,NJ,K+1)+XX(I+1,NJ,K+1))/4.
                   Ycs=(YY(I,NJ,K)+YY(I+1,NJ,K)+YY(I,NJ,K+1)+YY(I+1,NJ,K+1))/4.
                   Zcs=(ZZ(I,NJ,K)+ZZ(I+1,NJ,K)+ZZ(I,NJ,K+1)+ZZ(I+1,NJ,K+1))/4.
                     Xc(I,NJ,K)=Xcs
                     Yc(I,NJ,K)=Ycs
                     Zc(I,NJ,K)=Zcs
                  ENDDO
                  ENDDO

                  DO J=1,NJ1    !! K direction
                  DO I=1,NI1
                   Xcs=(XX(I,J,1)+XX(I+1,J,1)+XX(I,J+1,1)+XX(I+1,J+1,1))/4.
                   Ycs=(YY(I,J,1)+YY(I+1,J,1)+YY(I,J+1,1)+YY(I+1,J+1,1))/4.
                   Zcs=(ZZ(I,J,1)+ZZ(I+1,J,1)+ZZ(I,J+1,1)+ZZ(I+1,J+1,1))/4.
                     Xc(i,j,0)=Xcs
                     Yc(i,j,0)=Ycs
                     Zc(i,j,0)=Zcs

                   Xcs=(XX(I,J,NK)+XX(I+1,J,NK)+XX(I,J+1,NK)+XX(I+1,J+1,NK))/4.
                   Ycs=(YY(I,J,NK)+YY(I+1,J,NK)+YY(I,J+1,NK)+YY(I+1,J+1,NK))/4.
                   Zcs=(ZZ(I,J,NK)+ZZ(I+1,J,NK)+ZZ(I,J+1,NK)+ZZ(I+1,J+1,NK))/4.
                     Xc(I,J,NK)=Xcs
                     Yc(I,J,NK)=Ycs
                     Zc(I,J,NK)=Zcs
                  ENDDO
                  ENDDO







                  DO K=1,NK1
                     DO J=1,NJ1
                        DO I=1,NI1
                           DxW=Xc(I-1,J,K)-Xc(I,J,K)
                           DyW=Yc(I-1,J,K)-Yc(I,J,K)
                           DzW=Zc(I-1,J,K)-Zc(I,J,K)

                           DxE=Xc(I+1,J,K)-Xc(I,J,K)
                           DyE=Yc(I+1,J,K)-Yc(I,J,K)
                           DzE=Zc(I+1,J,K)-Zc(I,J,K)

                           DxS=Xc(I,J-1,K)-Xc(I,J,K)
                           DyS=Yc(I,J-1,K)-Yc(I,J,K)
                           DzS=Zc(I,J-1,K)-Zc(I,J,K)

                           DxN=Xc(I,J+1,K)-Xc(I,J,K)
                           DyN=Yc(I,J+1,K)-Yc(I,J,K)
                           DzN=Zc(I,J+1,K)-Zc(I,J,K)

                           DxL=Xc(I,J,K-1)-Xc(I,J,K)
                           DyL=Yc(I,J,K-1)-Yc(I,J,K)
                           DzL=Zc(I,J,K-1)-Zc(I,J,K)

                           DxU=Xc(I,J,K+1)-Xc(I,J,K)
                           DyU=Yc(I,J,K+1)-Yc(I,J,K)
                           DzU=Zc(I,J,K+1)-Zc(I,J,K)
        !              ! 网格中心与相邻6个网格中心的距离
                           DDw=SQRT(DxW*DxW+DyW*DyW+DzW*DzW)
                           DDe=SQRT(DxE*DxE+DyE*DyE+DzE*DzE)
                           DDs=SQRT(DxS*DxS+DyS*DyS+DzS*DzS)
                           DDn=SQRT(DxN*DxN+DyN*DyN+DzN*DzN)
                           DDl=SQRT(DxL*DxL+DyL*DyL+DzL*DzL)
                           DDu=SQRT(DxU*DxU+DyU*DyU+DzU*DzU)
        !              ! 网格中心与相邻6个网格中心的最大距离和最小距离
                           DLTxyz=amax1(DDw,DDe,DDs,DDn,DDl,DDu)
                           DSTxyz=amin1(DDw,DDe,DDs,DDn,DDl,DDu)

                           if(NI1 == 2) DLTxyz=amax1(DDs,DDn,DDl,DDu)
                           if(NJ1 == 2) DLTxyz=amax1(DDw,DDe,DDl,DDu)
                           if(NK1 == 2) DLTxyz=amax1(DDw,DDe,DDs,DDn)

                               if(NI1 == 2) DSTxyz=amin1(DDs,DDn,DDl,DDu)
                           if(NJ1 == 2) DSTxyz=amin1(DDw,DDe,DDl,DDu)
                           if(NK1 == 2) DSTxyz=amin1(DDw,DDe,DDs,DDn)


                           aLeng(I,J,K)=DLTxyz            !! for Delta: the scale of grids
                         aSeng(I,J,K)=DSTxyz            
                        ENDDO
                     ENDDO
                  ENDDO



            !!!!!!
            !	  Vol
            !!!!!!
                ! 计算网格体积
                C13=1./3.
                DO K=1,NK1
                   K1=K+1
                   DO J=1,NJ1
                      J1=J+1
                      DO I=1,NI1
                         I1=I+1
                         XR=XX(I1,J1,K1)-XX(I,J1,K1)
                         YR=YY(I1,J1,K1)-YY(I,J1,K1)
                         ZR=ZZ(I1,J1,K1)-ZZ(I,J1,K1)
                         OMGv=XR*SD(1,1,I,J,K)+YR*SD(1,2,I,J,K)+ZR*SD(1,3,I,J,K)

                         XR=XX(I1,J1,K1)-XX(I1,J,K1)
                         YR=YY(I1,J1,K1)-YY(I1,J,K1)
                         ZR=ZZ(I1,J1,K1)-ZZ(I1,J,K1)
                         OMGv=OMGv+XR*SD(2,1,I,J,K)+YR*SD(2,2,I,J,K)+ZR*SD(2,3,I,J,K)

                         XR=XX(I1,J1,K1)-XX(I1,J1,K)
                         YR=YY(I1,J1,K1)-YY(I1,J1,K)
                         ZR=ZZ(I1,J1,K1)-ZZ(I1,J1,K)
                         OMGv=OMGv+XR*SD(3,1,I,J,K)+YR*SD(3,2,I,J,K)+ZR*SD(3,3,I,J,K)

                         Vol(I,J,K)=C13* OMGv

                           if(Vol(I,J,K) < Tiny)then
                          write(*,*) iBlock
                                      write(*,*)I,J,K
                                          write(*,*)"The volume is negative !!!! Please test "
                              pause
                           endif

                        ENDDO
                     ENDDO
                  ENDDO





        !  by hjb, center points on the corner interfaces

        !XC
                do j=1,NJ1
                        Xc(0,j,0) = 1.5 * Xc(0,j,1) - Xc(0,j,2) * 0.5
                        Xc(0,j,NK) = 1.5 * Xc(0,j,NK1) - Xc(0,j,NK1-1) * 0.5
                        Xc(NI,j,0) = 1.5 * Xc(NI,j,1) - Xc(NI,j,2) * 0.5
                        Xc(NI,j,NK) = 1.5 * Xc(NI,j,NK1) - Xc(NI,j,NK1-1) * 0.5

                        Yc(0,j,0) = 1.5 * Yc(0,j,1) - Yc(0,j,2) * 0.5
                        Yc(0,j,NK) = 1.5 * Yc(0,j,NK1) - Yc(0,j,NK1-1) * 0.5
                        Yc(NI,j,0) = 1.5 * Yc(NI,j,1) - Yc(NI,j,2) * 0.5
                        Yc(NI,j,NK) = 1.5 * Yc(NI,j,NK1) - Yc(NI,j,NK1-1) * 0.5

                        Zc(0,j,0) = 1.5 * Zc(0,j,1) - Zc(0,j,2) * 0.5
                        Zc(0,j,NK) = 1.5 * Zc(0,j,NK1) - Zc(0,j,NK1-1) * 0.5
                        Zc(NI,j,0) = 1.5 * Zc(NI,j,1) - Zc(NI,j,2) * 0.5
                        Zc(NI,j,NK) = 1.5 * Zc(NI,j,NK1) - Zc(NI,j,NK1-1) * 0.5

                enddo

                do k=0,NK
                        Xc(0,0,k) = 1.5 * Xc(0,1,k) - Xc(0,2,k) * 0.5
                        Xc(0,NJ,k) = 1.5 * Xc(0,NJ1,k) - Xc(0,NJ1-1,k) * 0.5
                        Xc(NI,0,k) = 1.5 * Xc(NI,1,k) - Xc(NI,2,k) * 0.5
                        Xc(NI,NJ,k) = 1.5 * Xc(NI,NJ1,k) - Xc(NI,NJ1-1,k) * 0.5

                        Yc(0,0,k) = 1.5 * Yc(0,1,k) - Yc(0,2,k) * 0.5
                        Yc(0,NJ,k) = 1.5 * Yc(0,NJ1,k) - Yc(0,NJ1-1,k) * 0.5
                        Yc(NI,0,k) = 1.5 * Yc(NI,1,k) - Yc(NI,2,k) * 0.5
                        Yc(NI,NJ,k) = 1.5 * Yc(NI,NJ1,k) - Yc(NI,NJ1-1,k) * 0.5

                        Zc(0,0,k) = 1.5 * Zc(0,1,k) - Zc(0,2,k) * 0.5
                        Zc(0,NJ,k) = 1.5 * Zc(0,NJ1,k) - Zc(0,NJ1-1,k) * 0.5
                        Zc(NI,0,k) = 1.5 * Zc(NI,1,k) - Zc(NI,2,k) * 0.5
                        Zc(NI,NJ,k) = 1.5 * Zc(NI,NJ1,k) - Zc(NI,NJ1-1,k) * 0.5
                enddo

                do i=0,NI
                        Xc(i,0,0) = 1.5 * Xc(i,1,0) - Xc(i,2,0) * 0.5
                        Xc(i,NJ,0) = 1.5 * Xc(i,NJ1,0) - Xc(i,NJ1-1,0) * 0.5
                        Xc(i,0,NK) = 1.5 * Xc(i,1,NK) - Xc(i,2,NK) * 0.5
                        Xc(i,NJ,NK) = 1.5 * Xc(i,NJ1,NK) - Xc(i,NJ1-1,NK) * 0.5

                        Yc(i,0,0) = 1.5 * Yc(i,1,0) - Yc(i,2,0) * 0.5
                        Yc(i,NJ,0) = 1.5 * Yc(i,NJ1,0) - Yc(i,NJ1-1,0) * 0.5
                        Yc(i,0,NK) = 1.5 * Yc(i,1,NK) - Yc(i,2,NK) * 0.5
                        Yc(i,NJ,NK) = 1.5 * Yc(i,NJ1,NK) - Yc(i,NJ1-1,NK) * 0.5


                        Zc(i,0,0) = 1.5 * Zc(i,1,0) - Zc(i,2,0) * 0.5
                        Zc(i,NJ,0) = 1.5 * Zc(i,NJ1,0) - Zc(i,NJ1-1,0) * 0.5
                        Zc(i,0,NK) = 1.5 * Zc(i,1,NK) - Zc(i,2,NK) * 0.5
                        Zc(i,NJ,NK) = 1.5 * Zc(i,NJ1,NK) - Zc(i,NJ1-1,NK) * 0.5
                enddo



              DO k=1,NK1
                DO j=1,NJ1
                 Vol( 0,j,k)=Vol(  1,j,k)
                 Vol(NI,j,k)=Vol(NI1,j,k)
                 Dst( 0,j,k)=Dst(  1,j,k)
                 Dst(NI,j,k)=Dst(NI1,j,k)
                 Alagm(1,0,j,k)=Alagm(1,1,j,k)
                 Alagm(1,NI,j,k)=Alagm(1,NI1,j,k)
              END DO
                END DO

              DO k=1,NK1
                DO i=1,NI1
                 Vol( i, 0,k)=Vol(i,  1,k)
                 Vol( i,NJ,k)=Vol(i,NJ1,k)
                 Dst( i, 0,k)=Dst(i,  1,k)
                 Dst( i,NJ,k)=Dst(i,NJ1,k)
                 Alagm(2,i,0,k)=Alagm(2,i,1,k)
                 Alagm(2,i,NJ,k)=Alagm(2,i,NJ1,k)
               END DO
                END DO

              DO j=1,NJ1
                DO i=1,NI1
                 Vol( i,J, 0)=Vol(i,j,  1)
                 Vol( i,J,NK)=Vol(i,j,NK1)
                 Dst( i,J, 0)=Dst(i,j,  1)
                 Dst( i,J,NK)=Dst(i,j,NK1)
                 Alagm(3,i,j,0)=Alagm(3,i,j,1)
                 Alagm(3,i,j,NK)=Alagm(3,i,j,NK1)
              END DO
                END DO


                DO L=1,3
              DO k=1,NK
                DO j=1,NJ
                   SD(1,L,   0,J,K)=SD(1,L,  2,J,K)
                   SD(1,L,NI+1,J,K)=SD(1,L,NI1,J,K)
              END DO
                END DO

              DO k=1,NK
                DO i=1,NI
                   SD(2,L,I,   0,K)=SD(2,L,I,  2,K)
                   SD(2,L,I,NJ+1,K)=SD(2,L,I,NJ1,K)
              END DO
                END DO

              DO j=1,NJ
                DO i=1,NI
                   SD(3,L,I,J,0   ) =SD(3,L,I,J,  2)
                   SD(3,L,I,J,NK+1) =SD(3,L,I,J,NK1)
              END DO
                END DO
                ENDDO

            do k=1,NK1
            do j=1,NJ1
            do i=1,NI
                txc=0.25*(XX(i,j,k)+XX(i,j,k+1)+XX(i,j+1,k)+XX(i,j+1,k+1))
                tyc=0.25*(YY(i,j,k)+YY(i,j,k+1)+YY(i,j+1,k)+YY(i,j+1,k+1))
                tzc=0.25*(ZZ(i,j,k)+ZZ(i,j,k+1)+ZZ(i,j+1,k)+ZZ(i,j+1,k+1))
                radSurf(1,i,j,k)=sqrt(tyc**2.0+tzc**2.0)
                xVec=txc-Rot_Ori(1)
                yVec=tyc-Rot_Ori(2)
                zVec=tzc-Rot_Ori(3)
                velocityS=-(-(w(2)*zVec-w(3)*yVec))
                normalS=sqrt(SD(1,1,i,j,k)**2.0+SD(1,2,i,j,k)**2.0+SD(1,3,i,j,k)**2.0)
             !   gridV(1,1,i,j,k)=(w(2)*zVec-w(3)*yVec)*SD(1,1,i,j,k)/normalS
              !  gridV(1,2,i,j,k)=(w(3)*xVec-w(1)*zVec)*SD(1,2,i,j,k)/normalS
               ! gridV(1,3,i,j,k)=(w(1)*yVec-w(2)*xVec)*SD(1,3,i,j,k)/normalS
                !vib=(wxr)*n
              !  vibn(i,j,k)=-(w(2)*zVec-w(3)*yVec)*SD(1,1,i,j,k)+(w(1)*zVec-w(3)*xVec)*SD(1,2,i,j,k)- &
              !  &   (w(1)*yVec-w(2)*xVec)*SD(1,3,i,j,k)

                thtf(1,i,j,k)=atan2(tzc,tyc)
            enddo
            enddo
            enddo
            
            do K=1,NK1
            do J=1,NJ
            do I=1,NI1
                txc=0.25*(XX(i,j,k)+XX(i,j,k+1)+XX(i+1,j,k)+XX(i+1,j,k+1))
                tyc=0.25*(YY(i,j,k)+YY(i,j,k+1)+YY(i+1,j,k)+YY(i+1,j,k+1))
                tzc=0.25*(ZZ(i,j,k)+ZZ(i,j,k+1)+ZZ(i+1,j,k)+ZZ(i+1,j,k+1))
                radSurf(2,i,j,k)=sqrt(tyc**2.0+tzc**2.0)
                xVec=txc-Rot_Ori(1)
                yVec=tyc-Rot_Ori(2)
                zVec=tzc-Rot_Ori(3)
                velocityS=-(w(1)*zVec-w(3)*xVec)
                normalS=sqrt(SD(2,1,i,j,k)**2.0+SD(2,2,i,j,k)**2.0+SD(2,3,i,j,k)**2.0)
          !      gridV(2,1,i,j,k)=(w(2)*zVec-w(3)*yVec)*SD(2,1,i,j,k)/normalS
           !     gridV(2,2,i,j,k)=(w(3)*xVec-w(1)*zVec)*SD(2,2,i,j,k)/normalS
            !    gridV(2,3,i,j,k)=(w(1)*yVec-w(2)*xVec)*SD(2,3,i,j,k)/normalS
                !vib=(wxr)*n
                !vibn(i,j,k)=-(w(2)*zVec-w(3)*yVec)*SD(1,1,i,j,k)+(w(1)*zVec-w(3)*xVec)*SD(1,2,i,j,k)- &
                !vib=(wxr)*n
        !        vjbn(i,j,k)=-(w(2)*zVec-w(3)*yVec)*SD(2,1,i,j,k)+(w(1)*zVec-w(3)*xVec)*SD(2,2,i,j,k)- &
         !       &   (w(1)*yVec-w(2)*xVec)*SD(2,3,i,j,k)

                thtf(2,i,j,k)=atan2(tzc,tyc)
            enddo
            enddo
            enddo

            do k=1,NK
            do j=1,NJ1
            do i=1,NI1
                txc=0.25*(XX(i,j,k)+XX(i,j+1,k)+XX(i+1,j,k)+XX(i+1,j+1,k))
                tyc=0.25*(YY(i,j,k)+YY(i,j+1,k)+YY(i+1,j,k)+YY(i+1,j+1,k))
                tzc=0.25*(ZZ(i,j,k)+ZZ(i,j+1,k)+ZZ(i+1,j,k)+ZZ(i+1,j+1,k))
                radSurf(3,i,j,k)=sqrt(tyc**2.0+tzc**2.0)
                xVec=txc-Rot_Ori(1)
                yVec=tyc-Rot_Ori(2)
                zVec=tzc-Rot_Ori(3)
                velocityS=-(-w(1)*yVec+w(2)*xVec)
                normalS=sqrt(SD(3,1,i,j,k)**2.0+SD(3,2,i,j,k)**2.0+SD(3,3,i,j,k)**2.0)
               ! gridV(3,1,i,j,k)=(w(2)*zVec-w(3)*yVec)*SD(3,1,i,j,k)/normalS
               ! gridV(3,2,i,j,k)=(w(3)*xVec-w(1)*zVec)*SD(3,2,i,j,k)/normalS
               ! gridV(3,3,i,j,k)=(w(1)*yVec-w(2)*xVec)*SD(3,3,i,j,k)/normalS
        !vib=(wxr)*n
        !vibn(i,j,k)=-(w(2)*zVec-w(3)*yVec)*SD(1,1,i,j,k)+(w(1)*zVec-w(3)*xVec)*SD(1,2,i,j,k)- &
        !vibn=(wxr)*n
      !  vkbn(i,j,k)=-(w(2)*zVec-w(3)*yVec)*SD(3,1,i,j,k)+(w(1)*zVec-w(3)*xVec)*SD(3,2,i,j,k)- &
       ! &   (w(1)*yVec-w(2)*xVec)*SD(3,3,i,j,k)
        
        thtf(3,i,j,k)=atan2(tzc,tyc)
!write(*,*)i,j,k,tzc,tyc,thtf(3,i,j,k),atan2(tzc,tyc)
    enddo
    enddo
    enddo

END SUBROUTINE GEOM_CalcDST





SUBROUTINE GEOM_OutputGEO(iBlock)
    USE Global
    IMPLICIT NONE
    INTEGER:: iBlock
    TYPE(BlockStruct), POINTER:: Block

    include "../tecio.F90"  !how?
    
    INTEGER*4:: I,J,K
    INTEGER*4:: MaxI,MaxJ,MaxK
    CHARACTER(LEN=1):: NULLCHR
    INTEGER*4:: Debug,II,III,VIsDouble
    REAL*4:: var(100)
    INTEGER*4:: Nul(100)
    INTEGER*4:: Nul1(1)
    POINTER (NullPtr,Nulp)
    INTEGER*4:: Nulp(*)
    REAL*4,ALLOCATABLE:: VARs(:,:,:)
    
    NullPtr=0
    Nul = 0

    Block => AllBlocks(iBlock)

        NULLCHR   = CHAR(0)
        Debug     = 0						 
        VIsDouble = 0						 
        
        MaxI=Block%NI1
        MaxJ=Block%NJ1
        MaxK=Block%NK1
    !!!!!
    !	  Block format
    !!!!!
        I = TecIni100('SIMPLE DATASET'//NULLCHR, &
     &  'X Y Z Dst aL aS Vol thtc' &
     &                  //NULLCHR, &
     &             trim(Block%FLgeo)//NULLCHR, &
     &             '.'  //NULLCHR, &
     &             Debug, &
     &             VIsDouble)

!        I = TecZne100('Simple Zone'//NULLCHR,&
!     &             MaxI, &
!     &             MaxJ, &
!     &             MaxK, &
!     &             'POINT'//NULLCHR, &
!     &             NULLCHR) 

        I = TecZne100('Simple Zone'//NULLCHR,&
     &             0, &
     &             MaxI, &
     &             MaxJ, &
     &             MaxK, &

     &             0, &
     &             0, &
     &             0, &

     &             1, &
     &   0,0,nulp,nulp,0)

!     &             'POINT'//NULLCHR, &
!     &             NULLCHR) 

    ALLOCATE(VARs(MaxI,MaxJ,MaxK))
    
      III=MaxI*MaxJ*MaxK

      do K=1,MaxK
	  do J=1,MaxJ
	  do I=1,MaxI
!	     var(1)=   Xc(I,J,K)
!	     var(2)=   Yc(I,J,K)
!	     var(3)=   Zc(I,J,K)
!	     var(4)=  Dst(I,J,K)
!	     var(5)=aLeng(I,J,K)
 !          var(6)=aSeng(I,J,K)
	!     var(7)=  Vol(I,J,K)
!		 II   = TecDat100(III,var,VIsDouble)
          VARs(I,J,K) = Xc(I,J,K)           !输出网格中心x坐标
	  enddo
	  enddo
      enddo
      II=TecDat100(III,VARs,VIsDouble)
      
      do K=1,MaxK
	  do J=1,MaxJ
	  do I=1,MaxI
          VARs(I,J,K) = Yc(I,J,K)           !输出网格中心y坐标
	  enddo
	  enddo
      enddo
      II=TecDat100(III,VARs,VIsDouble)

      do K=1,MaxK
	  do J=1,MaxJ
	  do I=1,MaxI
          VARs(I,J,K) = Zc(I,J,K)           !输出网格中心k坐标
	  enddo
	  enddo
      enddo
      II=TecDat100(III,VARs,VIsDouble)
      
      do K=1,MaxK
	  do J=1,MaxJ
	  do I=1,MaxI
          VARs(I,J,K) = Dst(I,J,K)          !输出网格中心到物面距离
	  enddo
	  enddo
      enddo
      II=TecDat100(III,VARs,VIsDouble)
      
      do K=1,MaxK
	  do J=1,MaxJ
	  do I=1,MaxI
          VARs(I,J,K) = Aleng(I,J,K)        !输出网格中心到邻接网格中心点最大距离
	  enddo
	  enddo
      enddo
      II=TecDat100(III,VARs,VIsDouble)
      
      do K=1,MaxK
	  do J=1,MaxJ
	  do I=1,MaxI
          VARs(I,J,K) = ASeng(I,J,K)        !输出网格中心到邻接网格中心点最小距离
	  enddo
	  enddo
      enddo
      II=TecDat100(III,VARs,VIsDouble)
      
      do K=1,MaxK
	  do J=1,MaxJ
	  do I=1,MaxI
          VARs(I,J,K) = Vol(I,J,K)          !输出网格体积
	  enddo
	  enddo
      enddo
      II=TecDat100(III,VARs,VIsDouble)
      
      do k=1,MaxK
      do j=1,MaxJ
      do i=1,MaxI
        VARs(I,J,K)=thtc(i,j,k)
      enddo
      enddo
      enddo
       II=TecDat100(III,VARs,VIsDouble)

      I = TecEnd100()
      
      DEALLOCATE(VARs)

END SUBROUTINE GEOM_OutputGEO


subroutine Ghost_FaceNormal !by ydd
use global

implicit none

    integer::iblock,i,j,k,NI,NJ,NK,NI1,NJ1,Nk1,I1,J1,K1
    real::tyc,tzc,R1X,R1Y,R1Z,R2X,R2Y,R2Z
    

    do iblock=1,Max_block
        CALL GLOBAL_SetPointers(iBlock)
        NI=ThisBlock%NI
        NJ=ThisBlock%NJ
        NK=ThisBlock%NK
        NI1=NI-1
        NJ1=NJ-1
        NK1=NK-1
        
        !i direction
        do k=1,NK1
            K1=k+1
        do j=1,NJ1
            J1=J+1
        do i=0,NI+1,NI+1
            
                         R1X=XX(I,J1,K1)-XX(I,J,K)
                         R1Y=YY(I,J1,K1)-YY(I,J,K)
                         R1Z=ZZ(I,J1,K1)-ZZ(I,J,K)
                         R2X=XX(I,J,K1) -XX(I,J1,K)
                         R2Y=YY(I,J,K1) -YY(I,J1,K)
                         R2Z=ZZ(I,J,K1) -ZZ(I,J1,K)
                         SD(1,1,I,J,K)=(R1Y*R2Z-R1Z*R2Y)
                         SD(1,2,I,J,K)=(R1Z*R2X-R1X*R2Z)
                         SD(1,3,I,J,K)=(R1X*R2Y-R1Y*R2X)
                            Grad(1,i,j,k)=sqrt(SD(1,1,i,j,k)**2.0+SD(1,2,i,j,k)**2.0+SD(1,3,i,j,k)**2.0)
                        tyc=0.25*(YY(I,J,K)+YY(I,J,K1)+YY(I,J1,K)+YY(I,J1,K1))
                        tzc=0.25*(ZZ(I,J,K)+ZZ(I,J,K1)+ZZ(I,J1,K)+ZZ(I,J1,K1))
                        thtf(1,i,j,k)=atan2(tzc,tyc)

        enddo
        enddo
        enddo
        
        do k=1,NK1
            k1=k+1
        do j=1,NJ
        do i=0,NI,NI
            I1=I+1
                         R1X=XX(I,J,K1)-XX(I1,J,K)
                         R1Y=YY(I,J,K1)-YY(I1,J,K)
                         R1Z=ZZ(I,J,K1)-ZZ(I1,J,K)
                         R2X=XX(I1,J,K1)-XX(I,J,K)
                         R2Y=YY(I1,J,K1)-YY(I,J,K)
                         R2Z=ZZ(I1,J,K1)-ZZ(I,J,K)
                         SD(2,1,I,J,K)=(R1Y*R2Z-R1Z*R2Y)
                         SD(2,2,I,J,K)=(R1Z*R2X-R1X*R2Z)
                         SD(2,3,I,J,K)=(R1X*R2Y-R1Y*R2X)
                    Grad(2,I,J,K)=SQRT(SD(2,1,I,J,K)*2.0+SD(2,2,I,J,K)*2.0+SD(2,3,I,J,K)*2.0+tiny)
                        tyc=0.25*(YY(I,J,K)+YY(I,J,K1)+YY(I1,J,K)+YY(I1,J,K1))
                        tzc=0.25*(ZZ(I,J,K)+ZZ(I,J,K1)+ZZ(I1,J,K)+ZZ(I1,J,K1))
                        thtf(2,i,j,k)=atan2(tzc,tyc)
        enddo
        enddo
        enddo

        do K=1,NK
        do J=1,NJ1
            J1=J+1
        do i=0,NI,NI
            I1=I+1
                         R1X=XX(I1,J,K)-XX(I,J1,K)
                         R1Y=YY(I1,J,K)-YY(I,J1,K)
                         R1Z=ZZ(I1,J,K)-ZZ(I,J1,K)
                         R2X=XX(I1,J1,K)-XX(I,J,K)
                         R2Y=YY(I1,J1,K)-YY(I,J,K)
                        SD(3,1,I,J,K)=(R1Y*R2Z-R1Z*R2Y)
                        SD(3,2,I,J,K)=(R1Z*R2X-R1X*R2Z)
                        SD(3,3,I,J,K)=(R1X*R2Y-R1Y*R2X)
                    Grad(3,i,j,k)=sqrt(SD(3,1,i,j,k)**2.0+SD(3,2,i,j,k)**2.0+SD(3,3,i,j,k)**2.0)
                        tyc=0.25*(YY(I,J,K)+YY(I,J1,K)+YY(I1,J,K)+YY(I1,J1,K))
                        tzc=0.25*(ZZ(I,J,K)+ZZ(I,J1,K)+ZZ(I1,J,K)+ZZ(I1,J1,K))
                        thtf(3,i,j,k)=atan2(tzc,tyc)
        enddo
        enddo
        enddo


        !j direction
        do k=1,NK1
            K1=k+1
        do j=0,NJ,NJ
            J1=J+1
        do i=1,NI
            
                         R1X=XX(I,J1,K1)-XX(I,J,K)
                         R1Y=YY(I,J1,K1)-YY(I,J,K)
                         R1Z=ZZ(I,J1,K1)-ZZ(I,J,K)
                         R2X=XX(I,J,K1) -XX(I,J1,K)
                         R2Y=YY(I,J,K1) -YY(I,J1,K)
                         R2Z=ZZ(I,J,K1) -ZZ(I,J1,K)
                         SD(1,1,I,J,K)=(R1Y*R2Z-R1Z*R2Y)
                         SD(1,2,I,J,K)=(R1Z*R2X-R1X*R2Z)
                         SD(1,3,I,J,K)=(R1X*R2Y-R1Y*R2X)
                            Grad(1,i,j,k)=sqrt(SD(1,1,i,j,k)**2.0+SD(1,2,i,j,k)**2.0+SD(1,3,i,j,k)**2.0)
                        tyc=0.25*(YY(I,J,K)+YY(I,J,K1)+YY(I,J1,K)+YY(I,J1,K1))
                        tzc=0.25*(ZZ(I,J,K)+ZZ(I,J,K1)+ZZ(I,J1,K)+ZZ(I,J1,K1))
                        thtf(1,i,j,k)=atan2(tzc,tyc)

        enddo
        enddo
        enddo
        
        do k=1,NK1
            k1=k+1
        do j=0,NJ+1,NJ+1
        do i=1,NI1
            I1=I+1
                         R1X=XX(I,J,K1)-XX(I1,J,K)
                         R1Y=YY(I,J,K1)-YY(I1,J,K)
                         R1Z=ZZ(I,J,K1)-ZZ(I1,J,K)
                         R2X=XX(I1,J,K1)-XX(I,J,K)
                         R2Y=YY(I1,J,K1)-YY(I,J,K)
                         R2Z=ZZ(I1,J,K1)-ZZ(I,J,K)
                         SD(2,1,I,J,K)=(R1Y*R2Z-R1Z*R2Y)
                         SD(2,2,I,J,K)=(R1Z*R2X-R1X*R2Z)
                         SD(2,3,I,J,K)=(R1X*R2Y-R1Y*R2X)
                    Grad(2,I,J,K)=SQRT(SD(2,1,I,J,K)*2.0+SD(2,2,I,J,K)*2.0+SD(2,3,I,J,K)*2.0+tiny)
                        tyc=0.25*(YY(I,J,K)+YY(I,J,K1)+YY(I1,J,K)+YY(I1,J,K1))
                        tzc=0.25*(ZZ(I,J,K)+ZZ(I,J,K1)+ZZ(I1,J,K)+ZZ(I1,J,K1))
                        thtf(2,i,j,k)=atan2(tzc,tyc)
        enddo
        enddo
        enddo

        do K=1,NK
        do J=0,NJ,NJ
            J1=J+1
        do i=1,NI1
            I1=I+1
                         R1X=XX(I1,J,K)-XX(I,J1,K)
                         R1Y=YY(I1,J,K)-YY(I,J1,K)
                         R1Z=ZZ(I1,J,K)-ZZ(I,J1,K)
                         R2X=XX(I1,J1,K)-XX(I,J,K)
                         R2Y=YY(I1,J1,K)-YY(I,J,K)
                        SD(3,1,I,J,K)=(R1Y*R2Z-R1Z*R2Y)
                        SD(3,2,I,J,K)=(R1Z*R2X-R1X*R2Z)
                        SD(3,3,I,J,K)=(R1X*R2Y-R1Y*R2X)
                    Grad(3,i,j,k)=sqrt(SD(3,1,i,j,k)**2.0+SD(3,2,i,j,k)**2.0+SD(3,3,i,j,k)**2.0)
                        tyc=0.25*(YY(I,J,K)+YY(I,J1,K)+YY(I1,J,K)+YY(I1,J1,K))
                        tzc=0.25*(ZZ(I,J,K)+ZZ(I,J1,K)+ZZ(I1,J,K)+ZZ(I1,J1,K))
                        thtf(3,i,j,k)=atan2(tzc,tyc)
        enddo
        enddo
        enddo


        !K direction
        do k=0,NK,NK
            K1=k+1
        do j=1,NJ1
            J1=J+1
        do i=1,NI
            
                         R1X=XX(I,J1,K1)-XX(I,J,K)
                         R1Y=YY(I,J1,K1)-YY(I,J,K)
                         R1Z=ZZ(I,J1,K1)-ZZ(I,J,K)
                         R2X=XX(I,J,K1) -XX(I,J1,K)
                         R2Y=YY(I,J,K1) -YY(I,J1,K)
                         R2Z=ZZ(I,J,K1) -ZZ(I,J1,K)
                         SD(1,1,I,J,K)=(R1Y*R2Z-R1Z*R2Y)
                         SD(1,2,I,J,K)=(R1Z*R2X-R1X*R2Z)
                         SD(1,3,I,J,K)=(R1X*R2Y-R1Y*R2X)
                            Grad(1,i,j,k)=sqrt(SD(1,1,i,j,k)**2.0+SD(1,2,i,j,k)**2.0+SD(1,3,i,j,k)**2.0)
                        tyc=0.25*(YY(I,J,K)+YY(I,J,K1)+YY(I,J1,K)+YY(I,J1,K1))
                        tzc=0.25*(ZZ(I,J,K)+ZZ(I,J,K1)+ZZ(I,J1,K)+ZZ(I,J1,K1))
                        thtf(1,i,j,k)=atan2(tzc,tyc)

        enddo
        enddo
        enddo
        
        do k=0,NK,NK
            k1=k+1
        do j=1,NJ
        do i=1,NI1
            I1=I+1
                         R1X=XX(I,J,K1)-XX(I1,J,K)
                         R1Y=YY(I,J,K1)-YY(I1,J,K)
                         R1Z=ZZ(I,J,K1)-ZZ(I1,J,K)
                         R2X=XX(I1,J,K1)-XX(I,J,K)
                         R2Y=YY(I1,J,K1)-YY(I,J,K)
                         R2Z=ZZ(I1,J,K1)-ZZ(I,J,K)
                         SD(2,1,I,J,K)=(R1Y*R2Z-R1Z*R2Y)
                         SD(2,2,I,J,K)=(R1Z*R2X-R1X*R2Z)
                         SD(2,3,I,J,K)=(R1X*R2Y-R1Y*R2X)
                    Grad(2,I,J,K)=SQRT(SD(2,1,I,J,K)*2.0+SD(2,2,I,J,K)*2.0+SD(2,3,I,J,K)*2.0+tiny)
                        tyc=0.25*(YY(I,J,K)+YY(I,J,K1)+YY(I1,J,K)+YY(I1,J,K1))
                        tzc=0.25*(ZZ(I,J,K)+ZZ(I,J,K1)+ZZ(I1,J,K)+ZZ(I1,J,K1))
                        thtf(2,i,j,k)=atan2(tzc,tyc)
        enddo
        enddo
        enddo

        do K=0,NK+1,NK+1
        do J=1,NJ1
            J1=J+1
        do i=1,NI1
            I1=I+1
                         R1X=XX(I1,J,K)-XX(I,J1,K)
                         R1Y=YY(I1,J,K)-YY(I,J1,K)
                         R1Z=ZZ(I1,J,K)-ZZ(I,J1,K)
                         R2X=XX(I1,J1,K)-XX(I,J,K)
                         R2Y=YY(I1,J1,K)-YY(I,J,K)
                        SD(3,1,I,J,K)=(R1Y*R2Z-R1Z*R2Y)
                        SD(3,2,I,J,K)=(R1Z*R2X-R1X*R2Z)
                        SD(3,3,I,J,K)=(R1X*R2Y-R1Y*R2X)
                    Grad(3,i,j,k)=sqrt(SD(3,1,i,j,k)**2.0+SD(3,2,i,j,k)**2.0+SD(3,3,i,j,k)**2.0)
                        tyc=0.25*(YY(I,J,K)+YY(I,J1,K)+YY(I1,J,K)+YY(I1,J1,K))
                        tzc=0.25*(ZZ(I,J,K)+ZZ(I,J1,K)+ZZ(I1,J,K)+ZZ(I1,J1,K))
                        thtf(3,i,j,k)=atan2(tzc,tyc)
        enddo
        enddo
        enddo
    enddo

end
