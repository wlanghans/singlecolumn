;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; get the benchmark profile
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
f=addfile(benchdir+"/benchmark_SAM_BOMEXhighlarge_qt.nc","r")
zmean=f->height1
zflx=f->height2
nzmean = dimsizes(zmean)
nzflx = dimsizes(zflx)

qt =  f->QT
thl = f->THL
thvd=f->THVD
qtd=f->QTD

wthl=f->THLflx
wqt=f->QTflx

if (casename.eq."BOMEX") then

qc  = f->QC

else

thgrad=new(2,float)
do i =0,1
tmp=thl(1+i,:)
zlow=1200.+i*200.
zhigh=zlow+400.
zlow=700.
tmp=where(zmean.gt.zhigh.or.zmean.lt.zlow,-999.,tmp)
tmp@_FillValue=-999.
thgrad(i)=regline(zmean,tmp)
delete(tmp)
end do

end if



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;get the 1D model results 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
f=addfile(pathin+"/"+filein,"r")

tavg=-1.
if (casename.eq."CPBL2".or.casename.eq."CPBL4") then
tavg=600.
t_ind=(/3.,5./)*3600.
end if
if (casename.eq."BOMEX") then
tavg=1800.
t_ind=(/2.5/)*3600.
end if

if (tavg.eq.-1.) then
  print("ERROR: Case not specified in NCL postproc!")
  exit
end if

z=(f->z)
zi=(f->z2)
time=f->time
nzm=dimsizes(z)

thl1D = new((/dimsizes(t_ind),nzm/),float)
thv1D = new((/dimsizes(t_ind),nzm/),float)
qt1D = new((/dimsizes(t_ind),nzm/),float)
thvd1D = new((/dimsizes(t_ind),nzm/),float)
thvd1D@_FillValue=-999.
thvd1D=0.0
qtd1D = new((/dimsizes(t_ind),nzm/),float)
qtd1D@_FillValue=-999.
qtd1D=0.0

wthl1D = new((/dimsizes(t_ind),nzm+1/),float)
wqt1D = new((/dimsizes(t_ind),nzm+1/),float)

if (casename.eq."BOMEX") then
qc1D = new(nzm,float)
;rhoi1D= new(nzm+1,float)
else
thgrad1D=new((/dimsizes(t_ind)/),float)
end if

do t=0,dimsizes(t_ind)-1
    tmin_vec = ind(abs(time-t_ind(t)+tavg).eq.min(abs(time-t_ind(t)+tavg)))
    tmax_vec = ind(abs(time-t_ind(t)-tavg).eq.min(abs(time-t_ind(t)-tavg)))
    tmin=tmin_vec(0)
    tmax=tmax_vec(0)
    nt  =tmax-tmin+1
    avgind = ispan(tmin,tmax,1)
    thl1D(t,:)=(dim_avg_n(f->thl(tmin:tmax,:),0))
    thv1D(t,:)=(dim_avg_n(f->thr(tmin:tmax,:),0))
    qt1D(t,:)=(dim_avg_n(f->qt(tmin:tmax,:),0))

    upthv=f->upthv(tmin:tmax,:,:)
    upqt=f->upqt(tmin:tmax,:,:)
    fraca=f->upa(tmin:tmax,:,:)
    fractot=dim_sum_n(fraca,1)

    do k=0,nzm-1
    navg=0
    do indi=0,nt-1
      if (fractot(indi,k).ge.0.01) then
         fracsmall=0.
         fracsmallold=0.
         np2=39
         lfrac=True
         do while (np2.ge.0.and.lfrac)
          fracsmall=fracsmallold+fraca(indi,np2,k)
          if (fracsmall.gt.0.01) then
            wgt = (0.01-fracsmallold)/0.01
          else
            wgt = fraca(indi,np2,k)/0.01
          end if
          thvd1D(t,k) = thvd1D(t,k) + upthv(indi,np2,k) * wgt
          qtd1D(t,k) = qtd1D(t,k) + upqt(indi,np2,k) * wgt
          fracsmallold=fracsmall
          if (fracsmallold.gt.0.01) then
            lfrac=False
          end if
          np2=np2-1
         end do
         navg=navg+1
      end if
    end do
    if (navg.gt.0) then
       thvd1D(t,k) = thvd1D(t,k) / int2flt(navg)
       qtd1D(t,k) = qtd1D(t,k) / int2flt(navg)
    else
       thvd1D(t,k) = -999.
       qtd1D(t,k) = -999.
    end if
    end do
  

    do k=0,nzm-2
        if ((.not.ismissing(thvd1D(t,k))).and.(.not.ismissing(thvd1D(t,k+1)))) then
           tmp = 0.5*(thvd1D(t,k)+thvd1D(t,k+1))
           thvd1D(t,k)= tmp - thv1D(t,k)
           tmp = 0.5*(qtd1D(t,k)+qtd1D(t,k+1))
           qtd1D(t,k)= tmp - qt1D(t,k)
           delete(tmp)
        else
           thvd1D(t,k)= -999.
           qtd1D(t,k)= -999.
        end if
    end do
    thvd1D(t,nzm-1)= -999.
    qtd1D(t,nzm-1)= -999.

  
    wthl1D(t,:)=(dim_avg_n(f->tflx(avgind,:),0) )
    wqt1D(t,:)=(dim_avg_n(f->qtflx(avgind,:),0)) 

    if (casename.eq."BOMEX") then
    qc1D(:)=(dim_avg_n(f->qcl(avgind,:),0))
    else
      tmp=thl1D(t,:)
      zlow=1200.+t*200.
      zhigh=zlow+400.
      zlow=700.
      tmp=where(z.gt.zhigh.or.z.lt.zlow,-999.,tmp)
      tmp@_FillValue=-999.
      thgrad1D(t)=regline(doubletofloat(z),tmp)
      delete(tmp)
    end if
    delete(avgind)
    delete(upthv)
    delete(upqt)
    delete(fraca)
    delete(fractot)
    delete(tmin_vec)
    delete(tmax_vec)
end do


rmse_thl = new(dimsizes(t_ind),float)
rmse_qt  = new(dimsizes(t_ind),float)
rmse_thl = 0.
rmse_qt  = 0.
rmse_thvd  = new(dimsizes(t_ind),float)
rmse_thvd=0.
rmse_qtd  = new(dimsizes(t_ind),float)
rmse_qtd=0.

rmse_thgrad  = new(dimsizes(t_ind),float)
rmse_thgrad=0.

rmse_wthl = new(dimsizes(t_ind),float)
rmse_wqt  = new(dimsizes(t_ind),float)
rmse_wthl = 0.
rmse_wqt  = 0.


do n=0,dimsizes(t_ind) - 1 

if (casename.ne."BOMEX") then
rmse_thgrad(n)=abs(thgrad1D(n) - thgrad(n))
else
rmse_qc   = 0.
end if

ncount = 0
ncount_thvd = 0
do k=0,nzmean-1

  if (casename.eq."BOMEX") then
    if (zmean(k).lt.2900.) then
      ;if (zmean(k).gt.10.) then
        ncount=ncount+1
        rmse_thl(n) = rmse_thl(n) + (thl(0,k)-thl1D(n,k))^2. 
        rmse_qt(n)  = rmse_qt(n)  + (qt(0,k)-qt1D(n,k))^2. 
        if (.not.ismissing(thvd1D(n,k))) then
          rmse_thvd(n)  = rmse_thvd(n)  + (thvd(0,k)-thvd1D(n,k))^2. 
          rmse_qtd(n)  = rmse_qtd(n)  + (qtd(0,k)-qtd1D(n,k))^2. 
          ncount_thvd=ncount_thvd+1
        end if
      ;end if
    end if
  end if
  if (casename.eq."CPBL2") then
    if (zmean(k).gt.10..and.zmean(k).lt.2900.) then
      ncount=ncount+1
      rmse_thl(n) = rmse_thl(n) + (thl(1+n,k)-thl1D(n,k))^2. 
      rmse_qt(n)  = rmse_qt(n)  + (qt(1+n,k)-qt1D(n,k))^2. 
      if (.not.ismissing(thvd1D(n,k))) then
        rmse_thvd(n)  = rmse_thvd(n)  + (thvd(1+n,k)-thvd1D(n,k))^2. 
        rmse_qtd(n)  = rmse_qtd(n)  + (qtd(1+n,k)-qtd1D(n,k))^2. 
        ncount_thvd=ncount_thvd+1
      end if
    end if
  end if
  if (casename.eq."CPBL4") then
    if (zmean(k).gt.10..and.zmean(k).lt.2900.) then
      ncount=ncount+1
      rmse_thl(n) = rmse_thl(n) + (thl(3+n,k)-thl1D(n,k))^2. 
      rmse_qt(n)  = rmse_qt(n)  + (qt(3+n,k)-qt1D(n,k))^2. 
      if (.not.ismissing(thvd1D(n,k))) then
        rmse_thvd(n)  = rmse_thvd(n)  + (thvd(3+n,k)-thvd1D(n,k))^2. 
        rmse_qtd(n)  = rmse_qtd(n)  + (qtd(3+n,k)-qtd1D(n,k))^2. 
        ncount_thvd=ncount_thvd+1
      end if
    end if
  end if

end do

  rmse_thl(n) = sqrt(rmse_thl(n)/int2flt(ncount))
  rmse_qt(n)  = sqrt(rmse_qt(n) /int2flt(ncount))
  if (ncount_thvd.gt.0) then
    rmse_thvd(n)= sqrt(rmse_thvd(n) /int2flt(ncount_thvd))
    rmse_qtd(n)= sqrt(rmse_qtd(n) /int2flt(ncount_thvd))
  else
    rmse_thvd(n)=0. 
    rmse_qtd(n)=0. 
  end if

ncount=0
ncount_thvd=0
do k=0,nzmean-1

  if (zmean(k).lt.2900.) then
  ncount=ncount+1
  if (casename.eq."BOMEX") then
    rmse_wthl(n) = rmse_wthl(n) + (wthl(0,k)-wthl1D(n,k+1))^2. 
    rmse_wqt(n)  = rmse_wqt(n)  + (wqt(0,k)-wqt1D(n,k+1))^2. 
    rmse_qc   = rmse_qc   + (qc(0,k)-qc1D(k))^2.
    ncount_thvd=ncount_thvd+1
  else
    rmse_wthl(n) = rmse_wthl(n) + (wthl(1+n,k)-wthl1D(n,k+1))^2. 
    rmse_wqt(n)  = rmse_wqt(n)  + (wqt(1+n,k)-wqt1D(n,k+1))^2. 
  end if
  end if

end do
  rmse_wthl(n)=sqrt(rmse_wthl(n)/int2flt(ncount))
  rmse_wqt(n)=sqrt(rmse_wqt(n)/int2flt(ncount))
  if (casename.eq."BOMEX") then
  rmse_qc=sqrt(rmse_qc/int2flt(ncount_thvd))
  end if

end do

lusethvd=(/True,False/)
do n=0,dimsizes(lusethvd)-1

if (casename.eq."BOMEX") then
  thl_0 = 0.298500898501127
  qt_0  = 4.709051607813674d-4
  qc_0  = 2.291507888805414d-6
  wthl_0=8.29152631298648
  wqt_0 =37.5370361230699
  thvd0 = 0.110751867145755
  qtd0 = 4.754294124718259d-4 
  if (lusethvd(n)) then
    rmse= 1./7.*(rmse_thl(0)/thl_0 + rmse_qt(0)/qt_0 + \
    rmse_wqt(0)/wqt_0 + rmse_wthl(0)/wthl_0+rmse_qc/qc_0 + rmse_thvd(0)/thvd0+rmse_qtd(0)/qtd0)
  else
    rmse= 1./5.*(rmse_thl(0)/thl_0 + rmse_qt(0)/qt_0 + \
    rmse_wqt(0)/wqt_0 + rmse_wthl(0)/wthl_0+rmse_qc/qc_0)
  end if

else
  thl_0 = 0.142648591929376
  thgrad0=9.579204282494366e-005
  qt_0  = 1.561694507888805e-004
  wthl_0=5.55608093415176
  wqt_0 =18.0438175345439
  thvd0 = 0.119297832160030
  qtd0 = 1.220929376408715d-4
  if (lusethvd(n)) then
     rmse = 1./7. * (avg(rmse_wthl)/wthl_0+avg(rmse_wqt)/wqt_0+avg(rmse_thl)/thl_0 + avg(rmse_qt)/qt_0 + avg(rmse_thvd)/thvd0 + avg(rmse_thgrad)/thgrad0+avg(rmse_qtd)/qtd0)
  else
     rmse = 1./5. * (avg(rmse_wthl)/wthl_0+avg(rmse_wqt)/wqt_0+avg(rmse_thl)/thl_0 + avg(rmse_qt)/qt_0 + avg(rmse_thgrad)/thgrad0) 
  end if
end if


; write case rmse's to run directory
if (casename.eq."CPBL2".or.casename.eq."CPBL4") then
  rmse_qc  =-999.
else
  rmse_thgrad=-999.
end if

if (.not.lusethvd(n)) then
  rmse_thvd2=rmse_thvd
  rmse_thvd2=-999.
  rmse_qtd2=rmse_qtd
  rmse_qtd2=-999.

  system("echo "+ sprintf("%13.8f",rmse) + " "+sprintf("%13.8f",avg(rmse_thl(:))) + " "+ sprintf("%13.8f",avg(rmse_thgrad(:))) + " "+sprintf("%13.8f",avg(rmse_qt(:))) \
  +" "+sprintf("%13.8f",rmse_qc)+" "+sprintf("%13.8f",avg(rmse_wthl(:)))+" "+sprintf("%13.8f",avg(rmse_wqt(:))) +" "+sprintf("%13.8f",avg(rmse_thvd2(:)))+" "+sprintf("%13.8f",avg(rmse_qtd2(:)))+" > " +pathin+"/rmse_"+casename)

else

  system("echo "+ sprintf("%13.8f",rmse) + " "+sprintf("%13.8f",avg(rmse_thl(:))) + " "+sprintf("%13.8f",avg(rmse_thgrad(:))) + " "+ sprintf("%13.8f",avg(rmse_qt(:))) \
  +" "+sprintf("%13.8f",rmse_qc)+" "+sprintf("%13.8f",avg(rmse_wthl(:)))+" "+sprintf("%13.8f",avg(rmse_wqt(:))) +" "+sprintf("%13.8f",avg(rmse_thvd(:)))+" "+sprintf("%13.8f",avg(rmse_qtd(:)))+" > " +pathin+"/rmse_thvd_"+casename)

end if


end do
