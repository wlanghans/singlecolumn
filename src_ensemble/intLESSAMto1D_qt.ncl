;create target grid
nzmean2  = floattointeger((3010.-10.)/20.) + 1
zmean2=fspan(10.,3010.,nzmean2)
nzflx2  = floattointeger((3000.-20.)/20.) + 1
zflx2=fspan(20.,3000.,nzflx2)

; create target variables 
;BOMEX: (qt,thl,wthl,wqt,ql,bu)
;CPBL:  (qt,thl,wthl,wqt,ql,bu)
mean=new((/25,nzmean2/),float)
flx=new((/10,nzflx2/),float)
mean=0.0
flx=0.0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; get the profiles from BOMEX
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
bomexdir="/global/scratch/langhans/OUT_STAT/"
;"/global/home/users/langhans/BOMEX_LES/"  
;fin=addfile(bomexdir+"/BOMEXmean.nc","r")
fin=addfile(bomexdir+"/BOMEX_LES_REFHL.nc","r")

plt_times=2.5

z=(fin->z)
zi=z
zi(0)=0.
zi(1:dimsizes(z)-1)=0.5*(z(0:dimsizes(z)-2) +z(1:dimsizes(z)-1))
time=fin->time
;day to sec
time=time *24.*3600.
nzm=dimsizes(z)

thl = new((/dimsizes(plt_times),nzm/),float)
qt = new((/dimsizes(plt_times),nzm/),float)
ql = new((/dimsizes(plt_times),nzm/),float)
wthl = new((/dimsizes(plt_times),nzm/),float)
wqt = new((/dimsizes(plt_times),nzm/),float)
thvd = new((/dimsizes(plt_times),nzm/),float)
qtd = new((/dimsizes(plt_times),nzm/),float)

plt_times=plt_times*3600.

do t=0,dimsizes(plt_times)-1
    tminh = ind(abs(time-plt_times(t)+1800.).eq.min(abs(time-plt_times(t)+1800.)))
    if (dimsizes(tminh).gt.1) then
      tmin = tminh(0)
    else
      tmin=tminh
    end if
    delete(tminh)

    tmaxh = ind(abs(time-plt_times(t)-1800.).eq.min(abs(time-plt_times(t)-1800.)))
    if (dimsizes(tmaxh).gt.1) then
      tmax = tmaxh(0)
    else
      tmax=tmaxh
    end if
    delete(tmaxh)

    nt  =tmax-tmin+1
    print("Using "+nt+" timesteps for averaging the SAM statistics")
    avgind = ispan(tmin,tmax,1)
    thl(t,:)=dim_avg_n(fin->THETAL(avgind,:),0)
    qt(t,:)=dim_avg_n(fin->QT(avgind,:),0)/1000.
    ql(t,:)=dim_avg_n(fin->QCL(avgind,:),0)/1000.
    wthl(t,:)=dim_avg_n(fin->TLFLUX(avgind,:),0)
    wqt(t,:)=dim_avg_n(fin->QTFLUX(avgind,:),0)
    thvd(t,:)=dim_avg_n(fin->TVw1_A(avgind,:),0)
    qtd(t,:)=dim_avg_n(fin->QTw1_(avgind,:),0)/1000. - qt(t,:)
end do
delete(plt_times)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;interpolate to target grid
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
zmean=z
zflx=zi
nzmean=dimsizes(zmean)
nzflx=dimsizes(zflx)

do k=0,nzmean2-1

  lfound=0
  do l=0,nzmean-2
    if (zmean(l+1).ge.zmean2(k).and.zmean(l).le.zmean2(k)) then
      mean(0,k) = thl(0,l) + (zmean2(k)-zmean(l))*(thl(0,l+1)-thl(0,l))/(zmean(l+1)-zmean(l))
      mean(1,k) = qt(0,l) + (zmean2(k)-zmean(l))*(qt(0,l+1)-qt(0,l))/(zmean(l+1)-zmean(l))
      mean(2,k) = ql(0,l) + (zmean2(k)-zmean(l))*(ql(0,l+1)-ql(0,l))/(zmean(l+1)-zmean(l))
      mean(3,k) = thvd(0,l) + (zmean2(k)-zmean(l))*(thvd(0,l+1)-thvd(0,l))/(zmean(l+1)-zmean(l))
      mean(20,k)= qtd(0,l) + (zmean2(k)-zmean(l))*(qtd(0,l+1)-qtd(0,l))/(zmean(l+1)-zmean(l))
      lfound=1
    end if
  end do
  if (lfound.eq.0) then
      mean(0,k) = -999.
      mean(1,k) = -999.
      mean(2,k) = -999.
      mean(3,k) = -999.
      mean(20,k) = -999.
  end if

end do

do k=0,nzflx2-1

  lfound=0
  do l=0,nzflx-2
    if (zflx(l+1).ge.zflx2(k).and.zflx(l).le.zflx2(k)) then
      flx(0,k) = wthl(0,l) + (zflx2(k)-zflx(l))*(wthl(0,l+1)-wthl(0,l))/(zflx(l+1)-zflx(l))
      flx(1,k) = wqt(0,l) + (zflx2(k)-zflx(l))*(wqt(0,l+1)-wqt(0,l))/(zflx(l+1)-zflx(l))
      lfound=1
    end if
  end do
  if (lfound.eq.0) then
      flx(0,k)=-999.
      flx(1,k)=-999.
  end if

end do


delete(z)
delete(zi)
delete(thl)
delete(qt)
delete(ql)
delete(wthl)
delete(wqt)
delete(thvd)
delete(qtd)
delete(time)
delete(zmean)
delete(zflx)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; get the profiles from CPBL cases
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ncpbl=1
dir="/global/scratch/langhans/OUT_STAT/"
cpblname=(/"CPBL2_LES_REF","CPBL4_LES_REF"/)

do n=0,ncpbl-1

fin=addfile(dir+"/"+cpblname(n)+".nc","r")
plt_times=(/3.,5./)


z=(fin->z)
zi=z
zi(0)=0.
zi(1:dimsizes(z)-1)=0.5*(z(0:dimsizes(z)-2) +z(1:dimsizes(z)-1))
time=fin->time
;day to sec
time=time *24.*3600.
nzm=dimsizes(z)

thl = new((/dimsizes(plt_times),nzm/),float)
qt = new((/dimsizes(plt_times),nzm/),float)
ql = new((/dimsizes(plt_times),nzm/),float)
wthl = new((/dimsizes(plt_times),nzm/),float)
wqt = new((/dimsizes(plt_times),nzm/),float)
thvd = new((/dimsizes(plt_times),nzm/),float)
qtd = new((/dimsizes(plt_times),nzm/),float)

plt_times=plt_times*3600.

do t=0,dimsizes(plt_times)-1

    delete(avgind)
    tminh = ind(abs(time-plt_times(t)+600.).eq.min(abs(time-plt_times(t)+600.)))
    if (dimsizes(tminh).gt.1) then
      tmin = tminh(0)
    else
      tmin=tminh
    end if
    delete(tminh)

    tmaxh = ind(abs(time-plt_times(t)-600.).eq.min(abs(time-plt_times(t)-600.)))
    if (dimsizes(tmaxh).gt.1) then
      tmax = tmaxh(0)
    else
      tmax=tmaxh
    end if
    delete(tmaxh)

    nt  =tmax-tmin+1
    print("Using "+nt+" timesteps for averaging the SAM statistics")
    avgind = ispan(tmin,tmax,1)
    thl(t,:)=dim_avg_n(fin->THETAL(avgind,:),0)
    qt(t,:)=dim_avg_n(fin->QT(avgind,:),0)/1000.
    ql(t,:)=dim_avg_n(fin->QCL(avgind,:),0)/1000.
    wthl(t,:)=dim_avg_n(fin->TLFLUX(avgind,:),0)
    wqt(t,:)=dim_avg_n(fin->QTFLUX(avgind,:),0)
    thvd(t,:)=dim_avg_n(fin->TVw1_A(avgind,:),0)
    qtd(t,:)=dim_avg_n(fin->QTw1_(avgind,:),0)/1000. - qt(t,:)

end do

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;interpolate to target grid
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
zmean=z
zflx=zi
nzmean=dimsizes(zmean)
nzflx=dimsizes(zflx)

do t=0,dimsizes(plt_times)-1

do k=0,nzmean2-1

  lfound=0
  do l=0,nzmean-2
    if (zmean(l+1).ge.zmean2(k).and.zmean(l).le.zmean2(k)) then
      mean(4+4*t+8*n,k) = thl(t,l) + (zmean2(k)-zmean(l))*(thl(t,l+1)-thl(t,l))/(zmean(l+1)-zmean(l))
      mean(5+4*t+8*n,k) = qt(t,l) + (zmean2(k)-zmean(l))*(qt(t,l+1)-qt(t,l))/(zmean(l+1)-zmean(l))
      mean(6+4*t+8*n,k) = ql(t,l) + (zmean2(k)-zmean(l))*(ql(t,l+1)-ql(t,l))/(zmean(l+1)-zmean(l))
      mean(7+4*t+8*n,k) = thvd(t,l) + (zmean2(k)-zmean(l))*(thvd(t,l+1)-thvd(t,l))/(zmean(l+1)-zmean(l))
      mean(21+t+2*n,k) = qtd(t,l) + (zmean2(k)-zmean(l))*(qtd(t,l+1)-qtd(t,l))/(zmean(l+1)-zmean(l))
      lfound=1
    end if
  end do
  if (lfound.eq.0) then
      mean(4+4*t+8*n,k) = -999.
      mean(5+4*t+8*n,k) = -999.
      mean(6+4*t+8*n,k) = -999.
      mean(7+4*t+8*n,k) = -999.
      mean(21+t+2*n,k) = -999.
  end if

end do

do k=0,nzflx2-1

  lfound=0
  do l=0,nzflx-2
    if (zflx(l+1).ge.zflx2(k).and.zflx(l).le.zflx2(k)) then
      flx(2+2*t+4*n,k) = wthl(t,l) + (zflx2(k)-zflx(l))*(wthl(t,l+1)-wthl(t,l))/(zflx(l+1)-zflx(l))
      flx(3+2*t+4*n,k) = wqt(t,l) + (zflx2(k)-zflx(l))*(wqt(t,l+1)-wqt(t,l))/(zflx(l+1)-zflx(l))
      lfound=1
    end if
  end do
  if (lfound.eq.0) then
      flx(2+2*t+4*n,k)=-999.
      flx(3+2*t+4*n,k)=-999.
  end if

end do



end do

delete(plt_times)
delete(avgind)
delete(z)
delete(zi)
delete(thl)
delete(qt)
delete(ql)
delete(wthl)
delete(wqt)

end do


       system("/bin/rm -f benchmark_SAM_BOMEXhighlarge_qt.nc")   ; remove any pre-existing file
       ncdf = addfile("benchmark_SAM_BOMEXhighlarge_qt.nc" ,"c")  ; open output netCDF file

       fAtt               = True    
       fAtt@title         = "Benchmark profiles from LES"
       fAtt@source_file   =  "SAM LES"
       fAtt@Conventions   = "None"
       fAtt@creation_date = systemfunc ("date")
       fileattdef( ncdf, fAtt )            ; copy file attributes

       outmean=(/mean((/0,4,8,12,16/),:)/)

       outmean!0 = "case"
       outmean!1 = "height1"
       outmean&case = (/0,1,2,3,4/)
       outmean&height1  = zmean2
       outmean@long_name = "liquid potential temperature"           ; assign attributes 
       outmean@_FillValue = -999.           ; assign attributes 
       ncdf->THL  = outmean                                         

       outmean2=(/mean((/2,6,10,14,18/),:)/)
       outmean2!0 = "case"
       outmean2!1 = "height1"
       outmean2&case = (/0,1,2,3,4/)
       outmean2&height1  = zmean2
       outmean2@long_name = "liquid water mixing ratio"           ; assign attributes 
       outmean2@_FillValue = -999.           ; assign attributes 
       ncdf->QC    = outmean2                         
      
       delete(outmean)

       outmean=(/mean((/1,5,9,13,17/),:)/)
       outmean!0 = "case"
       outmean!1 = "height1"
       outmean&case = (/0,1,2,3,4/)
       outmean&height1  = zmean2
       outmean@long_name = "total water mixing ratio"           ; assign attributes 
       outmean@_FillValue = -999.           ; assign attributes 
       ncdf->QT    =  outmean                        

       outmean3 = (/flx((/1,3,5,7,9/),:)/)
       outmean3!0 = "case"
       outmean3!1 = "height2"
       outmean3&case = (/0,1,2,3,4/)
       outmean3&height2  = zflx2
       outmean3@long_name = "qt flux"           
       outmean3@_FillValue = -999.           ; assign attributes 
       ncdf->QTflx =  outmean3                       

       delete(outmean3)

       outmean3 = (/flx((/0,2,4,6,8/),:)/)
       outmean3!0 = "case"
       outmean3!1 = "height2"
       outmean3&case = (/0,1,2,3,4/)
       outmean3&height2  = zflx2
       outmean3@long_name = "theta_l flux"           
       outmean3@_FillValue = -999.           ; assign attributes 
       ncdf->THLflx = outmean3                         

       delete(outmean3)

       outmean3 = (/mean((/3,7,11,15,19/),:)/)
       outmean3!0 = "case"
       outmean3!1 = "height1"
       outmean3&case = (/0,1,2,3,4/)
       outmean3&height1  = zmean2
       outmean3@long_name = "average virt pot temp anomaly of w>99%"           
       outmean3@_FillValue = -999.           ; assign attributes 
       ncdf->THVD = outmean3                         

       delete(outmean3)

       outmean3 = (/mean((/20,21,22,23,24/),:)/)
       outmean3!0 = "case"
       outmean3!1 = "height1"
       outmean3&case = (/0,1,2,3,4/)
       outmean3&height1  = zmean2
       outmean3@long_name = "average qt anomaly of w>99%"           
       outmean3@_FillValue = -999.           ; assign attributes 
       ncdf->QTD = outmean3                         

       zmean2!0 = "height1"
       ncdf->height1 = (/zmean2/)
       zflx2!0 = "height2"
       ncdf->height2 = (/zflx2/)
