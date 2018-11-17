load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

npara=6
pmin=(/0.,0.,0.,0.,100.,0./)
pmax=(/4.,4.,2000.,2.,1000.,2./)
;pmin=(/0.,0.,0.,100.,0./)
;pmax=(/4.,4.,2.,1000.,2./)


filebase=(/"x_srbf","x_dycors"/)
colors=(/"red","blue","blue","red"/)

ofs=asciiread("./opti_new/of_srbf.txt",(/300,20/),"float")
ofd=asciiread("./opti_new/of_dy.txt",(/300,20/),"float")

stat_s=stat_dispersion(ofs,False)
stat_d=stat_dispersion(ofd,False)

ym=new((/5,npara/),float)

;ym(3,:) =  (/0.34,  0.82,  0. ,0.37 , 434. , 0.47 /)
;ym(4,:) =  (/0.25,  0.63,  0.,  0.55,  445.,  0.50 /)


masks=new((/300,20/),float)
masks=-999.
masks@_FillValue=-999.
maskd=new((/300,20/),float)
maskd=-999.
maskd@_FillValue=-999.
;get index of single-best configurations
smin=999.
smin2=999.
smin3=999.
dmin=999.
ofs_ind=new(20,integer)
ofd_ind=new(20,integer)
sbins=new(2,integer)
sbins2=new(2,integer)
sbind=new(2,integer)

  filename="./opti_new/"+filebase(0) +"_"+"2"+".txt"
  tmp= asciiread(filename,(/300,20/),"float")
  tmp@_FillValue=-999.
  tmp=pmin(2) + tmp * (pmax(2)-pmin(2))

  filename="./opti_new/"+filebase(0) +"_"+"3"+".txt"
  tmp3= asciiread(filename,(/300,20/),"float")
  tmp3@_FillValue=-999.
  tmp3=pmin(3) + tmp3 * (pmax(3)-pmin(3))

  mask3=masks
  mask3=-999.
  mask3=where(tmp.eq.0.0.and.tmp3.lt.0.3.and.ofs.lt.0.43,1.,-999.)
  print("number from Stoch RBF:" + sum(mask3))
  sbin3=sbins

do i=0,19
  ofs_ind(i) = ind(ofs(:,i).eq. min(ofs(:,i)) )
  ofd_ind(i) = ind(ofd(:,i).eq. min(ofd(:,i)) )
  masks(ofs_ind(i),i) = 1.
  maskd(ofd_ind(i),i) = 1.
  if (ofs(ofs_ind(i),i).lt.smin) then
    sbins(0) = ofs_ind(i)
    sbins(1) = i
    smin=ofs(ofs_ind(i),i)
  end if
  if (ofd(ofd_ind(i),i).lt.dmin) then
    sbind(0) = ofd_ind(i)
    sbind(1) = i
    dmin=ofd(ofd_ind(i),i)
  end if
  if (ofs(ofs_ind(i),i).lt.smin2.and.tmp(ofs_ind(i),i).lt.300.0.and.tmp3(ofs_ind(i),i).lt.0.3) then
    sbins2(0) = ofs_ind(i)
    sbins2(1) = i
    smin2=ofs(ofs_ind(i),i)
  end if
  if (num(.not.ismissing(mask3(:,i))).ge.1) then
     if (min(mask3(:,i)*ofs(:,i)).le.smin3 ) then
         print("Setting sbin3")
         smin3 = min(mask3(:,i)*ofs(:,i))
         sbin3(0) = ind((mask3(:,i)*ofs(:,i)).eq.smin3)
         sbin3(1) = i
     end if
  end if
end do


masks=-999.
maskd=-999.
masks= where(ofs.lt.stat_s(23),1.,-999.)
maskd= where(ofd.lt.stat_d(23),1.,-999.)
print(stat_s(23) + " "+ stat_d(23))


data=new((/12,5/),float)

do i=0,npara-1
  filename="./opti_new/"+filebase(0) +"_"+i+".txt"
  tmp= asciiread(filename,(/300,20/),"float")
  tmp@_FillValue=-999.
  tmp=pmin(i) + tmp * (pmax(i)-pmin(i))
  ym(1,i) = tmp(sbins(0),sbins(1))
  ym(3,i) = tmp(sbins2(0),sbins2(1))
  print("SB Srbf p"+i+"="+tmp(sbin3(0),sbin3(1)) + " "+smin3 )
  tmp2=tmp*masks
  tmp2@_FillValue=-999.
  tmp3=tmp3*masks
  tmp3@_FillValue=-999.
  if (i.eq.2) then
    numlt100=num(tmp2.eq.0..and.tmp3.lt.0.2)
    print("numer of simulations with tl<100 "+ numlt100)
     
  end if
  stat=stat_dispersion(tmp2,False)
  data(2*i,0) =stat(2) 
  data(2*i,1) =stat(6) 
  data(2*i,2) =stat(8) 
  data(2*i,3) =stat(10) 
  data(2*i,4) =stat(14) 
  filename="./opti_new/"+filebase(1) +"_"+i+".txt"
  tmp= asciiread(filename,(/300,20/),"float")
  tmp@_FillValue=-999.
  tmp=pmin(i) + tmp * (pmax(i)-pmin(i))
  ym(2,i) = tmp(sbind(0),sbind(1))
  print("SB Dycors p"+i+"="+tmp(sbind(0),sbind(1)))
  tmp2=tmp*maskd
  tmp2@_FillValue=-999.
  if (i.eq.2) then
    numlt100=num(tmp2.eq.0.)
    print("numer of simulations with tl<100 "+ numlt100)
  end if
  stat=stat_dispersion(tmp2,False)
  data(2*i+1,0) =stat(2) 
  data(2*i+1,1) =stat(6) 
  data(2*i+1,2) =stat(8) 
  data(2*i+1,3) =stat(10) 
  data(2*i+1,4) =stat(14) 
end do


pnames=(/"a!BF8!q!F!v!N!","a!B!qt!N!","!F8!t!FB!l!NF!","w!B!min!N!","!F8!t!F!!BF8!e!NF!","b!BF8!e!NF!"/)


x=(/1.,1.15/)
plot=new(6,graphic)

wks = gsn_open_wks("eps","box_parameter")               ; send graphics to PNG file

  res            = True                         ; plot mods desired
  res@tmXBLabels = (/"StochRBF","DYCORS"/) ; labels for each box
  res@tmXBLabelsOn   = False
  res@tiMainFontHeightF=0.06
  res@tmXBOn=False
;***********************************************
; the function boxplot will accept three different
; resource lists. In this default example, we set
; two of them to False.
;**********************************************

 llres                   = True
 llres@gsLineThicknessF  = 1.                 ; line thickness 
;**********************************************
; resources that control color and width of boxes
;**********************************************  
  opti          = True			
  opti@boxWidth = .09				; Width of box (x units)
  opti@boxColors = (/"red","blue"/) 

  res@trXMinF         = 0.92          ; MÃ­nimo valor do eixo y.      
  res@trXMaxF         = 1.23         ; MÃ¡ximo valor do eixo y.

  mres               = True                     ; marker mods desired
  mres@gsMarkerIndex = 1                        ; polymarker style
  mres@gsMarkerSizeF = 20.                      ; polymarker size
  mres@gsMarkerColor = "darkgreen"                    ; polymarker color
  xm=avg(x)
  ym(0,0) = 0.58
  ym(0,1) = 0.32
  ym(0,2) = 400.
  ym(0,3) = 1.28
  ym(0,4) = 400.
  ym(0,5) = 0.5

  dum = new(6,graphic)
  dums = new(6,graphic)
  dums2 = new(6,graphic)
  dumd = new(6,graphic)
  dumd2 = new(6,graphic)

  do i=0,5
  mres@gsMarkerSizeF = 20.                      ; polymarker size
  mres@gsMarkerIndex = 1                        ; polymarker style
    res@tiMainString =pnames(i) 
    res@trYMinF         = pmin(i)          ; MÃ­nimo valor do eixo y.      
    res@trYMaxF         = pmax(i)         ; MÃ¡ximo valor do eixo y.
    res@tiYAxisString = pnames(i)
    plot(i) = boxplot(wks,x,data(2*i:2*i+1,:),opti,res,llres)	
  mres@gsMarkerColor = "darkgreen"                    ; polymarker color
    dum(i) = gsn_add_polymarker(wks,plot(i),xm,ym(0,i),mres) 
  mres@gsMarkerColor = "red"                    ; polymarker color
    dums(i) = gsn_add_polymarker(wks,plot(i),x(0)-0.02,ym(1,i),mres) 
  mres@gsMarkerColor = "blue"                    ; polymarker color
    dumd(i) = gsn_add_polymarker(wks,plot(i),x(1)-0.02,ym(2,i),mres) 
  mres@gsMarkerColor = "red"                    ; polymarker color
  mres@gsMarkerIndex = 4                   ; polymarker style
  mres@gsMarkerSizeF = 7.                      ; polymarker size
    dums2(i) = gsn_add_polymarker(wks,plot(i),x(0)+0.02,ym(3,i),mres) 
  mres@gsMarkerColor = "blue"                    ; polymarker color
  ;  dumd2(i) = gsn_add_polymarker(wks,plot(i),x(1)+0.02,ym(4,i),mres) 
  
  end do

  resp=True
  ;resp@gsnPanelXWhiteSpacePercent=5.

  gsn_panel(wks,plot,(/1,6/),resp)



print(ym(3,:))
print(ofs(sbins(0),sbins(1)))






