      subroutine svm(n2,m2,mval2,sigmasq2,Num_ctilde2,ctil2,point2,y2,
     &     miss)

c==========================================================================c
c                                                                          c
c        Main Program for NPA Algorithm for Dense Kernels.                 c
c                                                                          c
c==========================================================================c
c
c  Date of Last Revision: March 15, 1999
c
c  This is a Fortran program that implements NPA, the Nearest Point
c  Algorithm described in the following Technical Report:
c
c  S.S. Keerthi, S.K. Shevade, C. Bhattacharya and K.R.K. Murthy,
c  A Fast Iterative Nearest Point Algorithm for Support Vector Machine
c  Classifier Design, Technical Report TR-ISL-99-03, Intelligent Systems
c  Lab, Dept. of Computer Science and Automation, Indian Institute of
c  Science, Bangalore, India, 1999.
c
c  The following parameters have to be supplied in a file named, "innpa".
c
c  Line 1         :      n, the number of input variables
c  Line 2         :      m, the number of training vectors
c  Line 3         :      mval, the number of validation vectors
c  Line 4         :      sigmasq, the square of sigma in Gaussian kernel
c  Line 5         :      Num_ctilde, the number of ctilde values
c  Line 6         :      The ctilde values, one by one in the same line
c  Next m lines   :      The training pairs, one pair in each line.
c                        One training pair consists of an input vector
c                        and a target value (+1 for class 1, -1 for class 2)
c  Next mval lines:      The validation pairs, one pair in each line.
c                        One validation pair consists of an input vector
c                        and a target value (+1 for class 1, -1 for class 2)
c                        The number of misclassifications and the sum of
c                        squared violations on the validation set will be 
c                        printed in the output file mentioned below.
c
c  * Currently the program has been written for n limited to 14. If you want
c  to use this code for n > 14, simply replace all occurences of 14 in this 
c  code with the number that you need.
c
c  * Similarly, m+mval is limited to 2700. If you want to use this code for
c  m+mval > 2700, simply replace all occurences of 2700 in this code with the
c  number that you need.
c
c  * The output will be returned in a file named, "outnpa". The program has
c  been developed and tested using g77. Please report any problems in using
c  this code, by e-mail to ssk@csa.iisc.ernet.in.
c
c  * This code is being made available free only for non-commercial use. Any
c  work done using this code or its modifications should make a reference
c  to Indian Institute of Science and the Technical Report mentioned above.
c
c==========================================================================c
      implicit none
c
      integer n2, m2, mval2, y2(2700),Num_ctilde2
      double precision sigmasq2, point2(100,2700),ctil2(100)
      integer miss(2700)
c
      integer n, m, mval, y(2700), alind(2700), k, l, i, j, nsupp
      integer Examine_NonSV, SV_Optimality, NonSV_Optimality
      integer Num_Type1_Updates, Num_Type2_Updates, success
      integer nsupp_old, Max_Updates, Loop_Completed
      integer imax, jmax, kmax, Num_ctilde, kc, val_miss
c
      double precision ctilde, sigmasq, point(100,2700)
      double precision eu(2700), ev(2700), beta(2700), ctil(100)
      double precision deluu, delvv, deluv, delzz, Numk, zkzl
      double precision zdk, zdu, zdv, tol, eps, gU, gV, di, dj
      double precision h_U, h_V, bhat, gamma, wdzk, data_error, xik
c
      common /main1/ n, m, ctilde, sigmasq, point, y, Numk
      common /main2/ deluu, delvv, deluv, delzz, eu, ev
      common /main3/ alind, beta
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      n=n2
      m=m2
      mval=mval2
      sigmasq=sigmasq2
      Num_ctilde=Num_ctilde2
      do i=1,Num_ctilde
         ctil(i)=ctil2(i)
      enddo
      do i=1,m+mval
         do j=1,n
            point(j,i)=point2(j,i)
         enddo
         y(i)=y2(i)
      enddo
c  Open write file...
c
c      OPEN(6, FILE = 'outnpa', ERR = 500)
c
c  Apply the NPA algm, one ctilde at a time...
c
      do 400 kc=1,Num_ctilde
c
      ctilde = ctil(kc)
c      write(6,*)'Solution for ctilde:',ctilde
c
      eps = 0.001d0
c  Initialization...
      i = 0
      j = 0
      do 20 k=1,m
         if (i.eq.0 .and. y(k).eq.1) i = k
         if (j.eq.0 .and. y(k).eq.-1) j = k
         alind(k) = 0
         beta(k) = 0.0d0
   20 continue
      if (i.eq.0) then
c         write(6,*)'Class 1 is empty.'
         stop
      end if
      if (j.eq.0) then
c         write(6,*)'Class 2 is empty.'
         stop
      end if
c
      alind(i) = 1
      alind(j) = 1
c
      beta(i) = 1.0d0
      beta(j) = 1.0d0
c
      nsupp = 2
c
      call kernel(i,i,deluu)
      call kernel(j,j,delvv)
      call kernel(i,j,deluv)
      delzz = deluu + delvv - 2.0d0*deluv
      eu(i) = deluu
      eu(j) = deluv
      ev(i) = deluv
      ev(j) = delvv
c
      Examine_NonSV = 1
      SV_Optimality = 1
      NonSV_Optimality = 0
      Num_Type1_Updates = 0
      Num_Type2_Updates = 0
c  Initialize some counters...
      Numk = 0
c
  100 continue
      if (SV_Optimality .ne. 0 .and. NonSV_Optimality .ne. 0) go to 200
c
c  Cleaning up any zero multipliers from the support vector list...
      do 105 k=1,m
         if (alind(k) .eq. 1 .and. beta(k) .le. 0.0d0) then
            alind(k) = 0
            beta(k) = 0.0d0
         end if
  105 continue
c
c      write(6,*)'-----------------------------------------------------'
c      write(6,*)'After a loop...'
c      write(6,*)'Wolfe-Dual Obj fn:', 2.0d0/delzz
c      write(6,*)'Num of Kernel Evaluations:', Numk
c
      if (Examine_NonSV .eq. 1) then
         Num_Type1_Updates = 0
         NonSV_Optimality = 1
         Examine_NonSV = 0
         do 150 k=1,m
            if (alind(k) .eq. 1) go to 150
c  Compute eu and ev for k...
            eu(k) = 0.0d0
            ev(k) = 0.0d0
            do 110 l=1,m
               if (alind(l) .eq. 0) go to 110
               call kernel(k,l,zkzl)
               if (y(l) .eq. 1) then
                  eu(k) = eu(k) + beta(l)*zkzl
               else
                  ev(k) = ev(k) + beta(l)*zkzl
               end if
  110       continue 
            zdk = eu(k) - ev(k)
            zdu = deluu - deluv
            zdv = deluv - delvv
            tol = 0.5d0*eps*delzz
c  Check if ii will give sufficient improvement...
            if ((y(k).eq. 1 .and. zdu-zdk .ge. tol) .or.
     1          (y(k).eq.-1 .and. zdk-zdv .ge. tol)) then
               NonSV_Optimality = 0
               beta(k) = 0.0d0
               alind(k) = 1
               call Take_Step(k,success)
               if (success .eq. 0) then
                  alind(k) = 0
                  beta(k) = 0.0d0
               else
                  Num_Type1_Updates = Num_Type1_Updates + 1
               end if
            end if
  150    continue
      else
         Examine_NonSV = 1
         nsupp_old = nsupp
         nsupp = 0
         do 160 k=1,m
            if (alind(k) .eq. 0) go to 160
            nsupp = nsupp + 1
  160    continue
         if (iabs(nsupp_old-nsupp) .ge. (nsupp/20)) then
            Max_Updates = m
         else
            Max_Updates = 10*m
         end if
         Loop_Completed = 0
         Num_Type2_Updates = 0
  170    continue
         if (Loop_Completed .eq. 1) go to 190
            zdu = deluu - deluv
            zdv = deluv - delvv
            imax = 0
            gU = 0.0d0
            jmax = 0
            gV = 0.0d0
            do 180 k=1,m
               if (alind(k) .eq. 0) go to 180
               if (y(k) .eq. 1) then
                  i = k
                  di = zdu - eu(i) + ev(i) 
                  if (di .gt. gU) then
                     gU = di
                     imax = i
                  end if
               else
                  j = k
                  dj = eu(j) - ev(j) - zdv
                  if (dj .gt. gV) then
                     gV = dj
                     jmax = j
                  end if
               end if
  180       continue
            tol = 0.5d0*eps*delzz
            if (gU .le. tol .and. gV .le. tol) then
               Loop_Completed = 1
               SV_Optimality = 1
            else
               if (gU .ge. gV) then
                  kmax = imax
               else
                  kmax = jmax
               end if
               call Take_Step(kmax,success)
               if (success .eq. 0) then
                  Loop_Completed = 1
                  SV_Optimality = 0
               else
                  Num_Type2_Updates = Num_Type2_Updates + 1
               end if
               if (Num_Type2_Updates .gt. Max_Updates) then
                  Loop_Completed = 1
                  SV_Optimality = 0
               end if
            end if
            go to 170
  190    continue
      end if
      if (Num_Type1_Updates .eq. 0 .and. Num_Type2_Updates .eq. 0) then
c         write(6,*)'Optimality Criteria not satisfied.'
c         write(6,*)'Algorithm is unable to make any more progress.'
         go to 200
      end if
      go to 100
c
  200 continue
c  Compute the number of support vectors...
      nsupp = 0
      do 210 k=1,m
         if (alind(k) .eq. 0) go to 210
         nsupp = nsupp + 1
  210 continue
c
c  Print final output...
c
c      write(6,*)'-----------------------------------------------------'
c      write(6,*)'Final Output:'
c      if (SV_Optimality .eq. 1 .and. NonSV_Optimality .eq. 1) 
c     1   write(6,*)'Optimality Criteria satisfied.'
c      write(6,*)'Wolfe-Dual Obj fn:', 2.0d0/delzz
c      write(6,*)'Total Num of Kernel Evalns for getting Solution:',Numk
c      write(6,*)'Num of Support Vectors:',nsupp
c
c  FOR INFERENCE: 
c
c  Set h_U = ev(imax)-eu(imax), h_V = eu(jmax)-ev(jmax), 
c  compute gamma and bhat using equations (19) and (20) of the Tech Rept
c  and print alind, beta, gamma and bhat to some output file.
c
c  Then, for later inference, use the sign of 
c     [gamma sum_{l\in alind} y_l beta_l kernel(x,x_l)] + bhat 
c  to decide if a given x belongs to Class 1 or Class 2.
c
c  If you are used to standard SVM notations, then define b = bhat,
c  alpha(l) = gamma*beta(l) for all l, and use the sign of
c     sum_{l\in alind} y_l alpha_l kernel(x,x_l) + b 
c  to decide if a given x belongs to Class 1 or Class 2.
c
      h_U = ev(imax) - eu(imax)
      h_V = eu(jmax) - ev(jmax)
      gamma = 2.0d0 / (-h_U - h_V)
      bhat  = (h_V - h_U) / (h_V + h_U)
      val_miss = 0
      data_error = 0.0d0
      if (mval .eq. 0) go to 360
c  Compute eu and ev for each k in the validation set...
      do k=1,mval
         miss(k)=0
      enddo
      do 350 k=m+1, m+mval
         eu(k) = 0.0d0
         ev(k) = 0.0d0
         do 310 l=1,m
            if (alind(l) .eq. 0) go to 310
            call kernel(k,l,zkzl)
            if (y(l) .eq. 1) then
               eu(k) = eu(k) + beta(l)*zkzl
            else
               ev(k) = ev(k) + beta(l)*zkzl
            end if
  310    continue 
         wdzk = gamma * (eu(k) - ev(k))
         if ( (y(k) .eq.  1 .and. wdzk+bhat .lt. 0.0d0) .or.
     1        (y(k) .eq. -1 .and. wdzk+bhat .gt. 0.0d0) ) then
c            write(6,*) (point(l,k),l=1,n), y(k)
            val_miss = val_miss + 1
            miss(k-m)=1
         end if
         xik = y(k)*(wdzk+bhat) - 1.0d0
         if (xik .ge. 0.0d0) xik = 0.0d0
         data_error = data_error + xik**2
  350 continue
  360 continue
c      write(6,*)'Num of MisClassifications on Validation Set:',val_miss
c      write(6,*)'Sum Squared Violations over Validation Set:',data_error
c      write(6,*)'====================================================='
c
  400 continue
      return
c
c  500 write(6,*)'Error in reading.'
c      close(6)
      return
      end
c-----------------------------------------------------------------------
      subroutine kernel(i,j,kxixj)
c
      implicit none
c
      integer n, m, i, j, k, y(2700)
c
      double precision ctilde, sigmasq, point(100,2700), Numk
c
      double precision kxixj, normsq
c
      common /main1/ n, m, ctilde, sigmasq, point, y, Numk
c
      normsq = 0.0d0
      do 10 k=1,n
         normsq = normsq + (point(k,i)-point(k,j))**2
   10 continue
      kxixj = dexp(-normsq/(2.0d0*sigmasq))
      if (i .eq. j) kxixj = kxixj + (1.0d0/ctilde)
      Numk = Numk + 1.0
      return
      end
c-----------------------------------------------------------------------
      subroutine Line_Segment(d11,d22,d12,d,lambda,vert)
c
      implicit none
c
      integer vert
c
      double precision d11, d22, d12, d, lambda
      double precision num, den, dtil, ltil
c
      d = d11
      lambda = 0.0d0
      vert = 1
      if (d22 .lt. d) then
         d = d22
         lambda = 1.0d0
         vert = 2
      end if
      den = d11 + d22 - 2.0d0*d12
      if (den .le. 0.0d0) return
      num = (d11*d22 - d12*d12)
      dtil = num/den
      ltil = (d11 - d12) / den
      if (0.0d0 .lt. ltil .and. ltil .lt. 1.0d0 .and. 
     1    dtil .lt. d) then
         d = dtil
         lambda = ltil
         vert = 0
      end if
      return
      end
c-----------------------------------------------------------------------
      subroutine Triangle(d11,d22,d33,d12,d13,d23,d,la1,la2,la3,flag)
c
      implicit none
c
      integer flag, vert
c
      double precision d11, d22, d33, d12, d13, d23, d, la1, la2, la3
      double precision e11, e12, e22, f1, f2, den, dtil
      double precision lt1, lt2, lt3
c
      call Line_Segment(d11,d22,d12,d,la2,vert)
      la1 = 1.0d0 - la2
      la3 = 0.0d0
      if (vert .eq. 1) then
         flag = 1
      else if (vert .eq. 2) then
         flag = 2
      else
         flag = 4
      end if
c
      call Line_Segment(d22,d33,d23,dtil,lt3,vert)
      if (dtil .lt. d) then
         d = dtil
         la1 = 0.0d0
         la2 = 1.0d0 - lt3
         la3 = lt3
         if (vert .eq. 1) then
            flag = 2
         else if (vert .eq. 2) then
            flag = 3
         else
            flag = 5
         end if
      end if
c
      call Line_Segment(d33,d11,d13,dtil,lt1,vert)
      if (dtil .lt. d) then
         d = dtil
         la1 = lt1
         la2 = 0.0d0
         la3 = 1.0d0 - lt1 
         if (vert .eq. 1) then
            flag = 3
         else if (vert .eq. 2) then
            flag = 1
         else
            flag = 6
         end if
      end if
c
      e11 = d22 + d11 - 2.0d0*d12
      e22 = d33 + d11 - 2.0d0*d13
      e12 = d23 - d12 - d13 + d11
      den = e11*e22 - e12*e12
      if (den .le. 0.0d0) return
      f1 = d11 - d12
      f2 = d11 - d13
      lt2 = (e22*f1 - e12*f2) / den
      lt3 = (-e12*f1 + e11*f2) / den
      lt1 = 1.0d0 - lt2 - lt3
      dtil = d11 - lt2*f1 - lt3*f2
      if (lt1 .gt. 0.0d0 .and. lt2 .gt. 0.0d0 .and. 
     1    lt3 .gt. 0.0d0 .and. dtil .lt. d) then
         d = dtil
         la1 = lt1
         la2 = lt2
         la3 = lt3
         flag = 0
      end if
c
      return
      end
c-----------------------------------------------------------------------
      subroutine Two_Lines(r,s1,s2,d1,d2,dold,d,
     1                     la1,la2,flag,vert1,vert2)
c
      implicit none
c
      integer flag, vert1, vert2
c
      double precision r, s1, s2, d1, d2, dold, d, la1, la2, den
c
      if (d1 .le. 0.0d0 .or. d2 .le. 0.0d0) then
         flag = 0
         return
      end if
      den = d1*d2 - r*r
      if (den .le. 0.0d0) then
         flag = 0
         return
      end if
c
      la1 = (s1*d2 - s2*r) / den
      if (la1 .lt. 0.0d0) then
         la1 = 0.0d0
         vert1 = 1
      else if (la1 .gt. 1.0d0) then
         la1 = 1.0d0
         vert1 = 2
      else
         vert1 = 0
      end if
c
      la2 = (la1*r - s2) / d2
      if (la2 .lt. 0.0d0) then
         la2 = 0.0d0
         vert2 = 1
      else if (la2 .gt. 1.0d0) then
         la2 = 1.0d0
         vert2 = 2
      else
         d = dold + d1*la1*la1 + d2*la2*la2 - 2.0d0*r*la1*la2 + 
     1              2.0d0*s2*la2 - 2.0d0*s1*la1
         flag = 1
         vert2 = 0
         return
      end if
c
      la1 = (la2*r + s1) / d1
      if (la1 .lt. 0.0d0) then
         la1 = 0.0d0
         vert1 = 1
      else if (la1 .gt. 1.0d0) then
         la1 = 1.0d0
         vert1 = 2
      else
         vert1 = 0
      end if
c
      d = dold + d1*la1*la1 + d2*la2*la2 - 2.0d0*r*la1*la2 + 
     1           2.0d0*s2*la2 - 2.0d0*s1*la1
      flag = 1
      return
      end
c-----------------------------------------------------------------------
      subroutine Gilbert_Step(kmax,success)
c
      implicit none
c
      integer n, m, kmax, success, y(2700), alind(2700)
      integer imax, jmax, vert, k
c
      double precision ctilde, sigmasq, deluu, delvv, deluv, delzz, Numk 
      double precision eu(2700), ev(2700), beta(2700), point(100,2700)
      double precision d11, dimax, djmax, d22, d12, d, lambda, dot, oml
c
      common /main1/ n, m, ctilde, sigmasq, point, y, Numk
      common /main2/ deluu, delvv, deluv, delzz, eu, ev
      common /main3/ alind, beta
c
      if (y(kmax) .eq. 1) then
         imax = kmax
         d11 = delzz
         call kernel(imax,imax,dimax)
         d22 = dimax + delvv - 2.0d0*ev(imax)
         d12 = eu(imax) - ev(imax) - deluv + delvv
         call Line_Segment(d11,d22,d12,d,lambda,vert)
         if (d .ge. delzz) then
            success = 0
            return
         end if
         success = 1
         delzz = d
         oml = 1.0d0 - lambda
         deluu = oml*oml*deluu + lambda*lambda*dimax +
     1              2.0d0*lambda*oml*eu(imax)
         deluv = oml*deluv + lambda*ev(imax)
         do 10 k=1,m
            if (alind(k) .eq. 0) go to 10
            call kernel(imax,k,dot)
            eu(k) = oml*eu(k) + lambda*dot
            if (y(k) .eq. 1) beta(k) = oml*beta(k)
   10    continue
         beta(imax) = beta(imax) + lambda
      else
         jmax = kmax
         d11 = delzz
         call kernel(jmax,jmax,djmax)
         d22 = djmax + deluu - 2.0d0*eu(jmax)
         d12 = ev(jmax) - eu(jmax) - deluv + deluu
         call Line_Segment(d11,d22,d12,d,lambda,vert)
         if (d .ge. delzz) then
            success = 0
            return
         end if
         success = 1
         delzz = d
         oml = 1.0d0 - lambda
         delvv = oml*oml*delvv + lambda*lambda*djmax +
     1              2.0d0*lambda*oml*ev(jmax)
         deluv = oml*deluv + lambda*eu(jmax)
         do 20 k=1,m
            if (alind(k) .eq. 0) go to 20
            call kernel(jmax,k,dot)
            ev(k) = oml*ev(k) + lambda*dot
            if (y(k) .eq. -1) beta(k) = oml*beta(k)
   20    continue
         beta(jmax) = beta(jmax) + lambda
      end if
      return
      end
c-----------------------------------------------------------------------
      subroutine Take_Step(kmax,success)
c
      implicit none
c
      integer n, m, kmax, success, y(2700), alind(2700)
      integer imin, jmin, imax, jmax, kmin, k, flag
      integer vert1, vert2, i, j
c
      double precision ctilde, sigmasq, deluu, delvv, deluv, delzz, Numk 
      double precision eu(2700), ev(2700), beta(2700), point(100,2700)
      double precision zdu, zdv, eumin, evmin, eumax, evmax
      double precision mu, psi, psi1, psi2, d11, t1, d22, d33
      double precision d12, d13, d23, d, la1, la2, la3
      double precision la1b, la2b, la3b, zimin, zimax, zjmin, zjmax
      double precision r, s1, s2, d1, d2, dold, dimin, djmin, di, dj
c
      common /main1/ n, m, ctilde, sigmasq, point, y, Numk
      common /main2/ deluu, delvv, deluv, delzz, eu, ev
      common /main3/ alind, beta
c
c  Compute kmin.
c
      zdu = deluu - deluv
      zdv = deluv - delvv
      imin = 0
      dimin = 0.0d0
      jmin = 0
      djmin = 0.0d0
      do 10 k=1,m
         if (alind(k) .eq. 0 .or. beta(k) .le. 0.0d0) go to 10
         if (y(k) .eq. 1) then
            i = k
            di = zdu - eu(i) + ev(i) 
            if (di .le. dimin) then
               dimin = di
               imin = i
            end if
         else
            j = k
            dj = -zdv + eu(j) - ev(j) 
            if (dj .le. djmin) then
               djmin = dj
               jmin = j
            end if
         end if
   10 continue
      if (dimin .lt. djmin) then
         kmin = imin
      else
         kmin = jmin
      end if
c
c  Do a modified Gilbert step if necessary.
c
      if (kmin .eq. 0 .or. beta(kmin) .ge. 1.0d0) then
         call Gilbert_Step(kmax,success)
         return
      end if
c
      mu = beta(kmin) / (1.0d0 - beta(kmin))
      call kernel(kmin,kmax,psi)
      call kernel(kmax,kmax,psi1)
      call kernel(kmin,kmin,psi2)
      eumax = eu(kmax)
      evmax = ev(kmax)
      eumin = eu(kmin)
      evmin = ev(kmin)
c
c  Case 1. Triangle of class 1...
c
      if (y(kmax) .eq. 1 .and. y(kmin) .eq. 1) then
         imax = kmax
         imin = kmin
         d11 = delzz
         t1 = deluu - eumin - deluv + evmin 
         d22 = delzz + mu*mu*(deluu + psi2 - 2.0d0*eumin) + 
     1                 2.0d0*mu*t1
         d33 = psi1 + delvv - 2.0d0*evmax
         d12 = delzz + mu*t1
         d13 = eumax - deluv - evmax + delvv
         d23 = d13 + mu*(eumax - deluv - psi + evmin)
         call Triangle(d11,d22,d33,d12,d13,d23,d,la1,la2,la3,flag)
         if (d .ge. delzz) then
            success = 0
            return
         end if
         success = 1
         delzz = d
         la1b = la1 + la2 + la2*mu
         la2b = -la2*mu
         la3b = la3
         deluu = la1b*la1b*deluu + la2b*la2b*psi2 + la3b*la3b*psi1 +
     1           2.0d0*la1b*la2b*eumin + 2.0d0*la1b*la3b*eumax +
     2           2.0d0*la2b*la3b*psi
         deluv = la1b*deluv + la2b*evmin + la3b*evmax
         do 20 k=1,m
            if (alind(k) .eq. 0) go to 20 
            call kernel(imax,k,zimax)
            call kernel(imin,k,zimin)
            eu(k) = la1b*eu(k) + la2b*zimin + la3b*zimax
            if (y(k) .eq. 1) beta(k) = beta(k)*la1b
   20    continue
         beta(imax) = beta(imax) + la3b
         if (flag .eq. 2 .or. flag .eq. 5 .or. flag .eq. 3) then
            alind(imin) = 0
            beta(imin) = 0.0d0
         else
            beta(imin) = beta(imin) + la2b
         end if
c
c  Case 2. Triangle of class 2...
c
      else if (y(kmax) .eq. -1 .and. y(kmin) .eq. -1) then
         jmax = kmax
         jmin = kmin
         d11 = delzz
         t1 = delvv - evmin - deluv + eumin 
         d22 = delzz + mu*mu*(delvv + psi2 - 2.0d0*evmin) + 
     1                 2.0d0*mu*t1
         d33 = psi1 + deluu - 2.0d0*eumax
         d12 = delzz + mu*t1
         d13 = evmax - deluv - eumax + deluu
         d23 = d13 + mu*(eumin - deluv - psi + evmax)
         call Triangle(d11,d22,d33,d12,d13,d23,d,la1,la2,la3,flag)
         if (d .ge. delzz) then
            success = 0
            return
         end if
         success = 1
         delzz = d
         la1b = la1 + la2 + la2*mu
         la2b = -la2*mu
         la3b = la3
         delvv = la1b*la1b*delvv + la2b*la2b*psi2 + la3b*la3b*psi1 +
     1           2.0d0*la1b*la2b*evmin + 2.0d0*la1b*la3b*evmax +
     2           2.0d0*la2b*la3b*psi
         deluv = la1b*deluv + la2b*eumin + la3b*eumax
         do 30 k=1,m
            if (alind(k) .eq. 0) go to 30 
            call kernel(jmax,k,zjmax)
            call kernel(jmin,k,zjmin)
            ev(k) = la1b*ev(k) + la2b*zjmin + la3b*zjmax
            if (y(k) .eq. -1) beta(k) = beta(k)*la1b
   30    continue
         beta(jmax) = beta(jmax) + la3b
         if (flag .eq. 2 .or. flag .eq. 5 .or. flag .eq. 3) then
            alind(jmin) = 0
            beta(jmin) = 0.0d0
         else
            beta(jmin) = beta(jmin) + la2b
         end if
c
c  Case 3. Gilbert line of Class1 vs Reflected line of Class2 ...
c
      else if (y(kmax) .eq. 1 .and. y(kmin) .eq. -1) then
         imax = kmax
         jmin = kmin
         r = mu*(evmax - psi - deluv + eumin)
         s1 = (evmax - eumax - deluv + deluu)
         s2 = mu*(delvv - deluv - evmin + eumin)
         d1 = psi1 + deluu - 2.0d0*eumax
         d2 = mu*mu*(delvv + psi2 - 2.0d0*evmin)
         dold = delzz
         call Two_Lines(r,s1,s2,d1,d2,dold,d,
     1                  la1,la2,flag,vert1,vert2)
         if (flag .eq. 0 .or. d .ge. delzz) then
            success = 0
            return
         end if
         success = 1
         delzz = d
         deluu = (1.0d0 - la1)*(1.0d0 -la1)*deluu + 
     1           la1*la1*psi1 + 2.0d0*(1.0d0 - la1)*la1*eumax
         delvv = (1.0d0 + la2*mu)*(1.0d0 + la2*mu)*delvv +
     1           la2*la2*mu*mu*psi2 - 
     2           2.0d0*(1.0d0 + la2*mu)*la2*mu*evmin
         deluv = (1.0d0 - la1)*(1.0d0 + la2*mu)*deluv -
     1           (1.0d0 - la1)*la2*mu*eumin + 
     2           (1.0d0 + la2*mu)*la1*evmax -
     3           la1*la2*mu*psi
         do 40 k=1,m
            if (alind(k) .eq. 0) go to 40
            call kernel(imax,k,zimax)
            call kernel(jmin,k,zjmin)
            eu(k) = (1.0d0 - la1)*eu(k) + la1*zimax
            ev(k) = (1.0d0 + la2*mu)*ev(k) - la2*mu*zjmin
            if (y(k) .eq. 1) beta(k) = beta(k)*(1.0d0 - la1)
            if (y(k) .eq. -1) beta(k) = beta(k)*(1.0d0 + la2*mu)
   40    continue
         beta(imax) = beta(imax) + la1 
         if (vert2 .eq. 2) then
            alind(jmin) = 0
            beta(jmin) = 0.0d0
         else
            beta(jmin) = beta(jmin) - la2*mu
         end if
      else
c
c  Case 4. Gilbert line of Class2 vs Reflected line of Class1 ...
c
         jmax = kmax
         imin = kmin
         r = mu*(eumax - psi - deluv + evmin)
         s1 = (eumax - evmax - deluv + delvv)
         s2 = mu*(deluu - deluv - eumin + evmin)
         d1 = psi1 + delvv - 2.0d0*evmax
         d2 = mu*mu*(deluu + psi2 - 2.0d0*eumin)
         dold = delzz
         call Two_Lines(r,s1,s2,d1,d2,dold,d,
     1                  la1,la2,flag,vert1,vert2)
         if (flag .eq. 0 .or. d .ge. delzz) then
            success = 0
            return
         end if
         success = 1
         delzz = d
         delvv = (1.0d0 - la1)*(1.0d0 -la1)*delvv + 
     1           la1*la1*psi1 + 2.0d0*(1.0d0 - la1)*la1*evmax
         deluu = (1.0d0 + la2*mu)*(1.0d0 + la2*mu)*deluu +
     1           la2*la2*mu*mu*psi2 - 
     2           2.0d0*(1.0d0 + la2*mu)*la2*mu*eumin
         deluv = (1.0d0 - la1)*(1.0d0 + la2*mu)*deluv -
     1           (1.0d0 - la1)*la2*mu*evmin + 
     2           (1.0d0 + la2*mu)*la1*eumax -
     3           la1*la2*mu*psi
         do 50 k=1,m
            if (alind(k) .eq. 0) go to 50
            call kernel(jmax,k,zjmax)
            call kernel(imin,k,zimin)
            ev(k) = (1.0d0 - la1)*ev(k) + la1*zjmax
            eu(k) = (1.0d0 + la2*mu)*eu(k) - la2*mu*zimin
            if (y(k) .eq. -1) beta(k) = beta(k)*(1.0d0 - la1)
            if (y(k) .eq. 1) beta(k) = beta(k)*(1.0d0 + la2*mu)
   50    continue
         beta(jmax) = beta(jmax) + la1 
         if (vert2 .eq. 2) then
            alind(imin) = 0
            beta(imin) = 0.0d0
         else
            beta(imin) = beta(imin) - la2*mu
         end if

      end if
      return
      end
c









