;This script cycles through all the files in a "data" directory and calculates the standard devations of the timeseries and their correlations/regressions

;History:
; 1. It was originally used to create the supplementary table for Brierley & Wainer (2018) - although that was then tidied up in Excel to make it presentable
; 2. It was then modified to output amplitude of EOF-based modes, as well as index-based ones. The output of this script was then plotted as Fig 3 in Rehfeld et al (2020,Earth Sys. Dynam.)
; 3. Now it's being overwritten to provide some tables as inputs for my PMIP4_past2future_analyzer repo

; run with:
; conda activate ncl-nco; ncl -n compile_tidy_numbers.ncl

load "cvdp_data.functions.ncl"

data_dir="cvdp_output_by_experiment/"
files=systemfunc("ls "+data_dir+"*/*.nc")
files = str_sub_str(files, data_dir, "")
;files=(/"INM-CM4-8_historical.cvdp_data.1850-2014.nc"/) ;select a single file for testing...

;define the filename for the outputting
fname="tidy_numbers_per_simulations.csv"
system("rm -f "+fname)


; loop over the various CVDP files
do fi = 0, dimsizes(files)-1
  print("Working on "+files(fi))
  ;set up what to look over
  variables=new(50, string)
  headers=new(50, string)
  headers(0)="model"
  headers(1)="experiment"
  headers(2)="years"
  file_parts = str_split(files(fi),"/")
  expt=file_parts(0)
  filename=file_parts(1)
  name_parts = str_split(filename,".")
  name_strs0 = str_split(name_parts(0),"_")
  name_strs2 = str_split(name_parts(2),"-")
  variables(0)=name_strs0(0) ;model
  if dimsizes(name_strs0).eq.1 then
    variables(1)="observations" ;experiment
  else
    variables(1)=name_strs0(1) ;experiment
  end if
  variables(2)=toint(name_strs2(1))-toint(name_strs2(0))+1 ;years
  delete([/file_parts,expt,filename,name_parts,name_strs0,name_strs2/])

  ts_opt=True
  ts_opt@NO_EXIT_ON_ERROR=True
  if any(variables(1).eq.(/"1pctCO2","abrupt4xCO2","rcp26","rcp85"/)) then
    ts_opt@subset=True
    ts_opt@subset_syear=-30
    ts_opt@subset_eyear=0
    ts_opt@DETREND=True
  end if

  ;Global Mean TAS
  headers(3)="GMTemp"
  fil=addfile(data_dir+files(fi),"r")
  if .not.isfilevar(fil,"tas_spatialmean_ann") then
    variables(3)=variables@_FillValue
  else
    variables(3)=extract_latlon_areastat(data_dir+files(fi),"tas_spatialmean_ann",(/-90.,90./),(/0.,360./),"mean")
  end if
  delete(fil)

  ;IOD mean gradient
  ts_opt_abs=ts_opt
  ts_opt_abs@make_absolute=True
  ts_opt_abs@DETREND=False
  headers(4)="ATL3 mean"
  variables(4)=stat_ts_var(data_dir+files(fi),"atl3_timeseries_mon","mean",ts_opt_abs)

  ;IOD Amplitude
  headers(5)="ATL3 stddev"
  variables(5)=stat_ts_var(data_dir+files(fi),"atl3_timeseries_mon","stddev",ts_opt)

  ;Nino3.4 Amplitude
  headers(6)="Nino34 stddev"
  variables(6)=stat_ts_var(data_dir+files(fi),"nino34","stddev",ts_opt)

  ;atlantic_nino mean
  headers(7)="AtlNino mean"
  variables(7)=stat_ts_var(data_dir+files(fi),"atlantic_nino","mean",ts_opt_abs)

  ;atlantic_nino Mode
  headers(8)="AtlNino stddev"
  variables(8)=stat_ts_var(data_dir+files(fi),"atlantic_nino","stddev",ts_opt)

  ;atlantic_nino Annual Cycle
  headers(9)="AtlNino AnnCycAmp"
  variables(9)=stat_ts_var(data_dir+files(fi),"atlantic_nino","AnnCycAmp",ts_opt)

  ;Atlantic Meridional Mode
  headers(10)="AMM mean"
  variables(10)=stat_ts_var(data_dir+files(fi),"atlantic_meridional_mode","mean",ts_opt_abs)

  ;Atlantic Meridional Mode
  headers(11)="AMM stddev"
  variables(11)=stat_ts_var(data_dir+files(fi),"atlantic_meridional_mode","stddev",ts_opt)

  ;AMM Annual Cycle
  headers(9)="AMM AnnCycAmp"
  variables(9)=stat_ts_var(data_dir+files(fi),"atlantic_meridional_mode","AnnCycAmp",ts_opt)

  ;South Asian Monsoon Area
  headers(12)="SAMSmonsoon area mean"
  variables(12)=stat_ts_var(data_dir+files(fi),"monsoon_area_SAMS","mean",ts_opt)

  ;South Asian Monsoon Rainrate
  headers(13)="SAMSmonsoon rain mean"
  variables(13)=stat_ts_var(data_dir+files(fi),"monsoon_rain_SAMS","mean",ts_opt)

  ;South Asian Monsoon Area
  headers(14)="SAMSmonsoon area stddev"
  variables(14)=stat_ts_var(data_dir+files(fi),"monsoon_area_SAMS","stddev",ts_opt)

  ;South Asian Monsoon Rainrate
  headers(15)="SAMSmonsoon rain stddev"
  variables(15)=stat_ts_var(data_dir+files(fi),"monsoon_rain_SAMS","stddev",ts_opt)

  ;NE Brazil annual rain
  headers(16)="NEB rain ann"
  variables(16)=stat_ts_var(data_dir+files(fi),"ipcc_NEB_pr","mean",ts_opt_abs)

  ;NE Brazil summer rain
  ts_opt_jja_abs=ts_opt_abs
  ts_opt_jja_abs@SEASON="JJA"
  headers(17)="NEB rain jja"
  variables(17)=stat_ts_var(data_dir+files(fi),"ipcc_NEB_pr","mean",ts_opt_jja_abs)

  ;NE Brazil winter rain
  ts_opt_djf_abs=ts_opt_abs
  ts_opt_djf_abs@SEASON="DJF"
  headers(18)="NEB rain djf"
  variables(18)=stat_ts_var(data_dir+files(fi),"ipcc_NEB_pr","mean",ts_opt_djf_abs)

  ;NE Brazil summer rain stddev
  ts_opt_jja=ts_opt
  ts_opt_jja@SEASON="JJA"
  headers(19)="NEB rain jja stddev"
  variables(19)=stat_ts_var(data_dir+files(fi),"ipcc_NEB_pr","stddev",ts_opt_jja)

  ;NE Brazil winter rain stddev
  ts_opt_djf=ts_opt
  ts_opt_djf@SEASON="DJF"
  headers(20)="NEB rain djf stddev"
  variables(20)=stat_ts_var(data_dir+files(fi),"ipcc_NEB_pr","stddev",ts_opt_djf)

  ;compute correlations
  ts_opt_smooth=ts_opt
  ts_opt_smooth@SMOOTH=True
  ts_opt_smooth@nave=5
  nino34=read_ts(data_dir+files(fi),"nino34",ts_opt_smooth)
  atl3=read_ts(data_dir+files(fi),"atlantic_nino",ts_opt_smooth) 
  amo=read_ts(data_dir+files(fi),"amo_timeseries_lowpass_mon",ts_opt_smooth)

  headers(21)="r(atl3.nino34)"
  if (.not.all(ismissing(nino34))).and.(.not.all(ismissing(atl3))).and.(.not.(stddev(atl3).lt.0.001)) then
    variables(21)=escorc(nino34,atl3)   
  else
    variables(21)=variables@_FillValue   
  end if

  headers(29)="r(atl3.amo)"
  if (.not.all(ismissing(amo))).and.(.not.all(ismissing(atl3)).and.(.not.(stddev(amo).lt.0.001))) then
    variables(29)=escorc(amo,atl3)   
  else
    variables(29)=variables@_FillValue   
  end if  
  delete([/nino34,atl3,amo/])

  ;write out to file
  exists=ind(.not.ismissing(headers))
  if fi.eq.0 then
    hlist=[/str_join(headers(exists),",")/]
    write_table(fname, "w", hlist, "%s")
  end if 
  vlist=[/str_join(variables(exists),",")/]
  write_table(fname, "a", vlist, "%s")

  delete([/headers,variables,ts_opt,ts_opt_abs,ts_opt_jja,ts_opt_djf,ts_opt_jja_abs,ts_opt_djf_abs,ts_opt_smooth/])
end do
