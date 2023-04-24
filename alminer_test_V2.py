'''
cd /Users/garethjones/Desktop/ALminer_Enhanced
conda activate py39
python alminer_test_V2.py 
'''

#What should we do?
make_spec_plots=True
make_cont_plots=True
make_spat_plots=True
make_tables=True
make_line_plots=True
download_products=False
infile='INPUTS/IFS_Table.txt'

#------------------------

import alminer
import pandas
from astropy.io import ascii
from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import os
import tarfile
import astropy.io.fits as fits
import astropy.wcs as wcs

fs=13; fs2=25
font = {'family' : 'sans-serif','weight' : 'bold','size' : fs}
matplotlib.rc('font', **font)

almabands=[]
almabands.append([3,84.0,116.0])
almabands.append([4,125.0,163.0])
almabands.append([5,158.0,211.0])
almabands.append([6,211.0,275.0])
almabands.append([7,275.0,373.0])
almabands.append([8,385.0,500.0])
almabands.append([9,602.0,720.0])
almabands.append([10,787.0,950.0])

colorscheme=["#f42457",
"#71de5f",
"#ca5ee5",
"#bac400",
"#0060c7",
"#458e00",
"#d823a4",
"#00803b",
"#c7003f",
"#009575",
"#dd7800",
"#85a3ff",
"#ff8d44",
"#2f5095",
"#f1be64",
"#d9baf8",
"#884e00",
"#ffa5a1",
"#8c9b65",
"#833b54"]

#Finds right subplot (2x4) based on frequency
def getsubplotnum(freq):
	tempnum=[-1,-1]
	for sbi in range(len(almabands)):
		if almabands[sbi][1]<freq and almabands[sbi][2]>freq:
			if almabands[sbi][0]==3: tempnum=[0,0]
			if almabands[sbi][0]==4: tempnum=[0,1]
			if almabands[sbi][0]==5: tempnum=[0,2]
			if almabands[sbi][0]==6: tempnum=[0,3]
			if almabands[sbi][0]==7: tempnum=[1,0]
			if almabands[sbi][0]==8: tempnum=[1,1]
			if almabands[sbi][0]==9: tempnum=[1,2]
			if almabands[sbi][0]==10: tempnum=[1,3]
	return tempnum

#Finds right subplot (2x4) based on number
def getsubplotnum2(num):
	tempnum=[-1,-1]
	if num==0: tempnum=[0,0]
	if num==1: tempnum=[0,1]
	if num==2: tempnum=[0,2]
	if num==3: tempnum=[0,3]
	if num==4: tempnum=[1,0]
	if num==5: tempnum=[1,1]
	if num==6: tempnum=[1,2]
	if num==7: tempnum=[1,3]
	return tempnum

#Finds right subplot (3x4) based on number
def getsubplotnum3(num):
	tempnum=[-1,-1]
	if num==0: tempnum=[0,0]
	if num==1: tempnum=[0,1]
	if num==2: tempnum=[0,2]
	if num==3: tempnum=[0,3]
	if num==4: tempnum=[1,0]
	if num==5: tempnum=[1,1]
	if num==6: tempnum=[1,2]
	if num==7: tempnum=[1,3]
	if num==8: tempnum=[2,0]
	if num==9: tempnum=[2,1]
	if num==10: tempnum=[2,2]
	if num==11: tempnum=[2,3]	
	return tempnum

#Checks if two values are within a delta of each other
def iswithin(value1, value2, delta):
	if value1+delta>value2 and value1-delta<value2:
		return True
	else:
		return False

#Returns indicies of sorted list
def sortpl(projlist):
	temp_pl=[]
	for spli in range(len(projlist)):
		spli2=projlist[spli].split('.')
		temp_value_split=float(spli2[0])+float(spli2[2])*(1E-5)
		temp_pl.append(temp_value_split)
	index_split=np.argsort(temp_pl)
	return index_split

#Gets rid of repeated entries for overlapping SPWs
def spwoverlap(entry1,entry2):
	same_proj=(entry1['project_code']==entry2['project_code'])
	same_obj=(entry1['ALMA_source_name']==entry2['ALMA_source_name'])
	same_ra=(entry1['RAJ2000']==entry2['RAJ2000'])
	same_dec=(entry1['DEJ2000']==entry2['DEJ2000'])
	same_AR=(entry1['ang_res_arcsec']==entry2['ang_res_arcsec'])
	same_ET=(entry1['t_exptime']==entry2['t_exptime'])
	if same_proj and same_obj and same_ra and same_dec and same_AR and same_ET:
		return True
	else:
		return False

#Plot circle
#Circle ICRS 53.157917 -27.773611 0.007827
def plotCircle(RA_C,DEC_C,RAD_C,FILL_C,linestyle='-'):
	axes.add_patch(plt.Circle((RA_C,DEC_C), RAD_C, edgecolor=colorscheme[numcolor], alpha=0.8, fill=FILL_C, facecolor=colorscheme[numcolor],linestyle=linestyle))
				
#-------

f=open(infile,'r')
ff=f.readlines()

#Get atmospheric transmission
atm_file='INPUTS/Atmosphere.txt'
f_atm=open(atm_file,'r')
ff_atm=f_atm.readlines()
atm_data_3=[[],[]]
atm_data_4=[[],[]]
atm_data_5=[[],[]]
atm_data_6=[[],[]]
atm_data_7=[[],[]]
atm_data_8=[[],[]]
atm_data_9=[[],[]]
atm_data_10=[[],[]]
for i in range(5,len(ff_atm)):
	atm_temp=ff_atm[i].split('  ')
	if float(atm_temp[0])>=float(almabands[0][1]) and float(atm_temp[0])<=float(almabands[0][2]):
		atm_data_3[0].append(float(atm_temp[0]))
		atm_data_3[1].append(float(atm_temp[1].replace(' \n','')))
	if float(atm_temp[0])>=float(almabands[1][1]) and float(atm_temp[0])<=float(almabands[1][2]):
		atm_data_4[0].append(float(atm_temp[0]))
		atm_data_4[1].append(float(atm_temp[1].replace(' \n','')))
	if float(atm_temp[0])>=float(almabands[2][1]) and float(atm_temp[0])<=float(almabands[2][2]):
		atm_data_5[0].append(float(atm_temp[0]))
		atm_data_5[1].append(float(atm_temp[1].replace(' \n','')))
	if float(atm_temp[0])>=float(almabands[3][1]) and float(atm_temp[0])<=float(almabands[3][2]):
		atm_data_6[0].append(float(atm_temp[0]))
		atm_data_6[1].append(float(atm_temp[1].replace(' \n','')))
	if float(atm_temp[0])>=float(almabands[4][1]) and float(atm_temp[0])<=float(almabands[4][2]):
		atm_data_7[0].append(float(atm_temp[0]))
		atm_data_7[1].append(float(atm_temp[1].replace(' \n','')))
	if float(atm_temp[0])>=float(almabands[5][1]) and float(atm_temp[0])<=float(almabands[5][2]):
		atm_data_8[0].append(float(atm_temp[0]))
		atm_data_8[1].append(float(atm_temp[1].replace(' \n','')))
	if float(atm_temp[0])>=float(almabands[6][1]) and float(atm_temp[0])<=float(almabands[6][2]):
		atm_data_9[0].append(float(atm_temp[0]))
		atm_data_9[1].append(float(atm_temp[1].replace(' \n','')))
	if float(atm_temp[0])>=float(almabands[7][1]) and float(atm_temp[0])<=float(almabands[7][2]):
		atm_data_10[0].append(float(atm_temp[0]))
		atm_data_10[1].append(float(atm_temp[1].replace(' \n','')))

#Reset line table
if make_tables:
	line_table_txt=open('OUTPUTS/LineTable.txt','w');line_table_txt.close()
	#
	cont_table_txt=open('OUTPUTS/ContTable.txt','w');cont_table_txt.close()
	#
	mous_table_txt=open('OUTPUTS/MOUSTable.txt','w');mous_table_txt.close()

for i in range(len(ff)):

	#Get RA and DEC of target
	temp=ff[i].split(' ')
	name=temp[0]
	RA=(180./12.)*(float(temp[1])+(float(temp[2])/60.)+(float(temp[3])/3600.))
	decfact=1.
	if float(temp[4])<0.:
		decfact=-1.
		temp[4]=abs(float(temp[4]))
	DEC=decfact*(float(temp[4])+(float(temp[5])/60.)+(float(temp[6])/3600.))
	zred=float(temp[7])

	#Get lines
	lines=[]
	linefile='INPUTS/splatalogue_1.tsv'
	g=open(linefile)
	gg=g.readlines()
	for j in range(1,len(gg)):
		templine=gg[j].split('\t')
		templine_name=templine[0]+' '+templine[3]
		templine_freq=float(templine[2].split(', ')[0])
		lines.append([templine_name,templine_freq])
	g.close()


	#Make sure the object directory exists
	try:
		os.mkdir('OUTPUTS/'+name)
	except FileExistsError:
		pass

	#Query ALMA archive for RA, DEC
	print('->',name)
	myquery = alminer.conesearch(ra=RA, dec=DEC, public=True, point=True, tap_service='NRAO',print_query=True)

	#Get basic details of each returned observation
	EXPLORATION=alminer.explore(myquery, allcols=True, allrows=True)

	master_MOUS_list=[]
	try:
		for e_i in EXPLORATION.index:
			#Only get data for projects that aren't mosaics
			if EXPLORATION.loc[e_i]['is_mosaic']=='F':
				RA0=RA*(np.pi/180.)
				DEC0=DEC*(np.pi/180.)
				RA1=EXPLORATION.loc[e_i]['RAJ2000']*(np.pi/180.)
				DEC1=EXPLORATION.loc[e_i]['DEJ2000']*(np.pi/180.)
				temp_DEL=np.sin(DEC0)*np.sin(DEC1) + np.cos(DEC0)*np.cos(DEC1)*np.cos(RA0 - RA1)
				DEL=np.arccos(temp_DEL)*(180./np.pi)*3600.
				#DEL=np.sqrt((RA0-RA1)**2+(DEC0-DEC1)**2)*(180./np.pi)*3600.
				temp_MOUS_list={}
				temp_MOUS_list['project_code']=EXPLORATION.loc[e_i]['project_code']
				temp_MOUS_list['ALMA_source_name']=EXPLORATION.loc[e_i]['ALMA_source_name']
				temp_MOUS_list['RAJ2000']=EXPLORATION.loc[e_i]['RAJ2000']
				temp_MOUS_list['DEJ2000']=EXPLORATION.loc[e_i]['DEJ2000']
				temp_MOUS_list['ang_res_arcsec']=EXPLORATION.loc[e_i]['ang_res_arcsec']
				temp_MOUS_list['min_freq_GHz']=EXPLORATION.loc[e_i]['min_freq_GHz']
				temp_MOUS_list['max_freq_GHz']=EXPLORATION.loc[e_i]['max_freq_GHz']
				temp_MOUS_list['freq_res_kHz']=EXPLORATION.loc[e_i]['freq_res_kHz']
				temp_MOUS_list['FoV_arcsec']=EXPLORATION.loc[e_i]['FoV_arcsec'] 
				temp_MOUS_list['data_rights']=EXPLORATION.loc[e_i]['data_rights'] 
				temp_MOUS_list['t_exptime']=EXPLORATION.loc[e_i]['t_exptime'] 
				temp_MOUS_list['central_freq_GHz']=EXPLORATION.loc[e_i]['central_freq_GHz'] 
				temp_MOUS_list['cont_sens_bandwidth']=EXPLORATION.loc[e_i]['cont_sens_bandwidth'] 
				temp_MOUS_list['member_ous_uid']=EXPLORATION.loc[e_i]['member_ous_uid']
				master_MOUS_list.append(temp_MOUS_list)
	except AttributeError:
		pass
		
	#Make plots showing which lines fall into which band, and how many projects have observed each spectral range. Added in atmospheric transmission overlays.
	if make_spec_plots:
		linecolor='k'
		atmcolor=colorscheme[4]
		datacolor=colorscheme[5]
		#Plot spectral coverage
		fig, axes = plt.subplots(2,4,figsize=(12,8))
		fig.suptitle(name+'(z='+str(zred)+')',fontsize=fs2,weight='bold')
		bigstep_x=np.arange(84.0,950.0,0.1)
		bigstep_y=np.zeros(len(bigstep_x))
		projlist=[]
		projvalues=[]
		for j in range(len(master_MOUS_list)):
			tempx1=master_MOUS_list[j]['min_freq_GHz']
			tempx2=master_MOUS_list[j]['max_freq_GHz']
			if master_MOUS_list[j]['project_code'] not in projlist:
				projlist.append(master_MOUS_list[j]['project_code'])
				tempy=np.zeros(len(bigstep_x))
				for jj in range(len(bigstep_x)):
					if bigstep_x[jj]>tempx1 and bigstep_x[jj]<tempx2:
						tempy[jj]=1
				projvalues.append(tempy)
			else:
				projindex=projlist.index(master_MOUS_list[j]['project_code'])
				tempy=projvalues[projindex]
				for jj in range(len(bigstep_x)):
					if bigstep_x[jj]>tempx1 and bigstep_x[jj]<tempx2:
						tempy[jj]=1
				projvalues[projindex]=tempy
		for j in range(len(projlist)):
			for jj in range(len(bigstep_x)):
				bigstep_y[jj]+=projvalues[j][jj]
		if projlist==[]:
			for plt_i in [0,1]:
				for plt_j in [0,1,2,3]:
					axes[plt_i,plt_j].set_ylim(-0.1,1.1)
		#Add lines
		for j in range(len(lines)):
			linefreq=lines[j][1]/(1+zred)
			if getsubplotnum(linefreq)[0]!=-1:
				axisnums=getsubplotnum(linefreq)
				axes[axisnums[0],axisnums[1]].axvline(linefreq,linestyle='dashed',c=linecolor)
				xmin, xmax, ymin, ymax = plt.axis()
				axes[axisnums[0],axisnums[1]].text(linefreq,(ymax-ymin)*0.5,lines[j][0],rotation=270,fontsize=8)
		#Add band names
		for j in range(8):
			axisnums=getsubplotnum2(j)
			axes[axisnums[0],axisnums[1]].plot(bigstep_x,bigstep_y,c=datacolor)
			axes[axisnums[0],axisnums[1]].set_ylim(0,max(bigstep_y)+1)
			axes[axisnums[0],axisnums[1]].set_xlim(almabands[j][1],almabands[j][2])
			xpostext=axes[axisnums[0],axisnums[1]].get_xlim()[0]+0.05*abs(axes[axisnums[0],axisnums[1]].get_xlim()[0]-axes[axisnums[0],axisnums[1]].get_xlim()[1])
			ypostext=axes[axisnums[0],axisnums[1]].get_ylim()[1]-0.09*abs(axes[axisnums[0],axisnums[1]].get_ylim()[0]-axes[axisnums[0],axisnums[1]].get_ylim()[1])
			axes[axisnums[0],axisnums[1]].text(xpostext,ypostext,'Band '+str(almabands[j][0]),bbox=dict(facecolor='white', alpha=0.8))
			axes[axisnums[0],axisnums[1]].yaxis.get_major_locator().set_params(integer=True)
			#
			ax2 = axes[axisnums[0],axisnums[1]].twinx()
			ax2.tick_params(axis='y', labelcolor=atmcolor)
			if j==0:
				ax2.plot(atm_data_3[0],atm_data_3[1],c=atmcolor,alpha=0.8)
			if j==1:
				ax2.plot(atm_data_4[0],atm_data_4[1],c=atmcolor,alpha=0.8)
			if j==2:
				ax2.plot(atm_data_5[0],atm_data_5[1],c=atmcolor,alpha=0.8)
			if j==3:
				ax2.plot(atm_data_6[0],atm_data_6[1],c=atmcolor,alpha=0.8)
			if j==4:
				ax2.plot(atm_data_7[0],atm_data_7[1],c=atmcolor,alpha=0.8)
			if j==5:
				ax2.plot(atm_data_8[0],atm_data_8[1],c=atmcolor,alpha=0.8)
			if j==6:
				ax2.plot(atm_data_9[0],atm_data_9[1],c=atmcolor,alpha=0.8)
			if j==7:
				ax2.plot(atm_data_10[0],atm_data_10[1],c=atmcolor,alpha=0.8)
			ax2.set_ylim(0.0,1.0)
		ax1=fig.add_subplot(111, frameon=False)
		plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
		plt.xlabel("Observed Frequency [GHz]",weight='bold')
		plt.ylabel("Number of Projects",weight='bold',color=datacolor)
		ax2 = fig.add_subplot(111, sharex=ax1, frameon=False)
		ax2.yaxis.set_label_position("right")
		plt.ylabel('Atmospheric Transmission', rotation=270,labelpad=50.0,color=atmcolor,weight='bold')
		plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
		plt.tight_layout(pad=3.0,w_pad=0.1,h_pad=0.1)
		#plt.show()
		plt.savefig('OUTPUTS/'+name+'/SpecCoverage_'+name+'.png',dpi=300,bbox_inches='tight')
		plt.close()

	#Make plots showing continuum sensitivity in each band.
	if make_cont_plots:
		try:
			atmcolor=colorscheme[4]
			datacolor=colorscheme[5]
			fig, axes = plt.subplots(2,4,figsize=(12,8))
			fig.suptitle(name+'(z='+str(zred)+')',fontsize=fs2,weight='bold')
			bigstep_all=[]
			for j in range(len(master_MOUS_list)):	
				Pp=master_MOUS_list[j]['project_code']
				Ss=master_MOUS_list[j]['cont_sens_bandwidth']
				Ff=master_MOUS_list[j]['central_freq_GHz']
				added=False
				for jj in range(len(bigstep_all)):
					if Pp==bigstep_all[jj]['P'] and Ss==bigstep_all[jj]['S']:
						bigstep_all[jj]['F'].append(Ff)
						added=True
				if not added:
					bigstep_all.append({'P':Pp,'S':Ss,'F':[Ff]})
			for j in range(8):
				axisnums=getsubplotnum2(j)
				max_S=-1
				for jj in range(len(bigstep_all)):
					axes[axisnums[0],axisnums[1]].plot(bigstep_all[jj]['F'],[bigstep_all[jj]['S'] for x in range(len(bigstep_all[jj]['F']))],c=datacolor,marker='v',linestyle='dashed')
					if bigstep_all[jj]['S']>max_S:
						max_S=bigstep_all[jj]['S']
				axes[axisnums[0],axisnums[1]].set_ylim(0,max_S*1.3)
				axes[axisnums[0],axisnums[1]].set_xlim(almabands[j][1],almabands[j][2])
				xpostext=axes[axisnums[0],axisnums[1]].get_xlim()[0]+0.05*abs(axes[axisnums[0],axisnums[1]].get_xlim()[0]-axes[axisnums[0],axisnums[1]].get_xlim()[1])
				ypostext=axes[axisnums[0],axisnums[1]].get_ylim()[1]-0.09*abs(axes[axisnums[0],axisnums[1]].get_ylim()[0]-axes[axisnums[0],axisnums[1]].get_ylim()[1])
				axes[axisnums[0],axisnums[1]].text(xpostext,ypostext,'Band '+str(almabands[j][0]),bbox=dict(facecolor='white', alpha=0.8))
				axes[axisnums[0],axisnums[1]].yaxis.get_major_locator().set_params(integer=True)
				#
				ax2 = axes[axisnums[0],axisnums[1]].twinx()
				ax2.tick_params(axis='y', labelcolor=atmcolor)
				if j==0:
					ax2.plot(atm_data_3[0],atm_data_3[1],c=atmcolor,alpha=0.8)
				if j==1:
					ax2.plot(atm_data_4[0],atm_data_4[1],c=atmcolor,alpha=0.8)
				if j==2:
					ax2.plot(atm_data_5[0],atm_data_5[1],c=atmcolor,alpha=0.8)
				if j==3:
					ax2.plot(atm_data_6[0],atm_data_6[1],c=atmcolor,alpha=0.8)
				if j==4:
					ax2.plot(atm_data_7[0],atm_data_7[1],c=atmcolor,alpha=0.8)
				if j==5:
					ax2.plot(atm_data_8[0],atm_data_8[1],c=atmcolor,alpha=0.8)
				if j==6:
					ax2.plot(atm_data_9[0],atm_data_9[1],c=atmcolor,alpha=0.8)
				if j==7:
					ax2.plot(atm_data_10[0],atm_data_10[1],c=atmcolor,alpha=0.8)
				ax2.set_ylim(0.0,1.0)
			ax1=fig.add_subplot(111, frameon=False)
			plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
			plt.xlabel("Observed Frequency [GHz]",weight='bold')
			plt.ylabel("Continuum Sensitivity [uJy/beam]",weight='bold',color=datacolor)
			ax2 = fig.add_subplot(111, sharex=ax1, frameon=False)
			ax2.yaxis.set_label_position("right")
			plt.ylabel('Atmospheric Transmission', rotation=270,labelpad=50.0,color=atmcolor,weight='bold')
			plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
			plt.tight_layout(pad=3.0,w_pad=0.1,h_pad=0.1)
			#plt.show()
			plt.savefig('OUTPUTS/'+name+'/ContCoverage_'+name+'.png',dpi=300,bbox_inches='tight')
			plt.close()
		except ValueError:
			print('No Continuum')

	if make_spat_plots:
		fig, axes = plt.subplots(1,1,figsize=(8,8))
		fig.suptitle(name+'(z='+str(zred)+')',fontsize=fs2,weight='bold')
		plt.axvline(0,linestyle='dashed',color='k',alpha=0.5)
		plt.axhline(0,linestyle='dashed',color='k',alpha=0.5)
		numcolor=0;projlist=[];colorlist=[];lslist=[]
		biggestr=40.
		for j in range(len(master_MOUS_list)):
			if master_MOUS_list[j]['project_code'] not in projlist:
				temp_ra=master_MOUS_list[j]['RAJ2000']
				temp_dec=master_MOUS_list[j]['DEJ2000']
				D_RA=(temp_ra-RA)*3600.
				D_DEC=(temp_dec-DEC)*3600.
				plotCircle(D_RA,D_DEC,master_MOUS_list[j]['ang_res_arcsec']/2.,True)
				colorlist.append(colorscheme[numcolor])
				if master_MOUS_list[j]['data_rights']=='Public':
					plotCircle(D_RA,D_DEC,master_MOUS_list[j]['FoV_arcsec']/2.,False)
					lslist.append('-')
				else:
					plotCircle(D_RA,D_DEC,master_MOUS_list[j]['FoV_arcsec']/2.,False,linestyle='dashed')
					lslist.append('dashed')
				numcolor+=1
				projlist.append(master_MOUS_list[j]['project_code'])
				if (0.5*master_MOUS_list[j]['FoV_arcsec'])>biggestr:
					biggestr=(0.5*master_MOUS_list[j]['FoV_arcsec'])+5.
		plt.xlabel("Relative R.A. [\"]",weight='bold')
		plt.ylabel("Relative Dec. [\"]",weight='bold')

		plt.xlim(biggestr,-1.*biggestr)
		plt.ylim(-1.*biggestr,biggestr)
		indexproj=sortpl(projlist)
		for j in indexproj:
			axes.plot(-1E+100,-1E+100,marker='.',color=colorlist[j],label=projlist[j],linestyle=lslist[j])
		plt.legend()
		plt.savefig('OUTPUTS/'+name+'/FoV_'+name+'.png',dpi=300,bbox_inches='tight')
		plt.close()

	if make_tables:
		line_table_txt=open('OUTPUTS/LineTable.txt','a')
		line_table_txt.write(name+'\n')
		line_table_txt.write('Line Project Ang_Res Int_Time Targ_Name\n')
		table_line_list=[]
		for j in range(len(lines)):
			linefreq=lines[j][1]/(1+zred)
			for k in range(len(master_MOUS_list)):
				tempx1=master_MOUS_list[k]['min_freq_GHz']
				tempx2=master_MOUS_list[k]['max_freq_GHz']
				if tempx1<linefreq and tempx2>linefreq:
					templineline=lines[j][0]+' '
					templineline+=str(master_MOUS_list[k]['project_code'])+' '
					templineline+=str(master_MOUS_list[k]['ang_res_arcsec'])+' '
					templineline+=str(master_MOUS_list[k]['t_exptime'])+' '
					templineline+=str(master_MOUS_list[k]['ALMA_source_name'])+'\n'
					if templineline not in table_line_list:
						line_table_txt.write(templineline)
						table_line_list.append(templineline)
		line_table_txt.write('-----\n')
		line_table_txt.close()
		#
		big_cont_list=[]
		for k in range(len(master_MOUS_list)):
			temp_cont={'freq':[float(master_MOUS_list[k]['central_freq_GHz'])], 'proj':str(master_MOUS_list[k]['project_code']), 'angres':str(master_MOUS_list[k]['ang_res_arcsec']), 'texp':str(master_MOUS_list[k]['t_exptime']), 'MOUS':str(master_MOUS_list[k]['member_ous_uid']), 'sourcename':str(master_MOUS_list[k]['ALMA_source_name'])}
			absorbed=False
			for ik in range(len(big_cont_list)):
				if big_cont_list[ik]['proj']==temp_cont['proj'] and big_cont_list[ik]['MOUS']==temp_cont['MOUS']:
					big_cont_list[ik]['freq'].append(temp_cont['freq'][0])
					absorbed=True
			if not absorbed:
				big_cont_list.append(temp_cont)
		cont_table_txt=open('OUTPUTS/ContTable.txt','a')
		cont_table_txt.write(name+'\n')
		cont_table_txt.write('Freq Project Ang_Res Int_Time Targ_Name Tunings MOUS\n')
		for k in range(len(big_cont_list)):
			tempcont=str(round(np.mean(big_cont_list[k]['freq']),3))+' '
			tempcont+=str(big_cont_list[k]['proj'])+' '
			tempcont+=str(big_cont_list[k]['angres'])+' '
			tempcont+=str(big_cont_list[k]['texp'])+' '
			tempcont+=str(big_cont_list[k]['sourcename'])+' '
			tempcont+=str(int(len(big_cont_list[k]['freq'])/4))+' '
			tempcont+=str(big_cont_list[k]['MOUS'])+'\n'
			cont_table_txt.write(tempcont)
		cont_table_txt.write('-----\n')
		cont_table_txt.close()
		#
		mous_table_txt=open('OUTPUTS/MOUSTable.txt','a')
		for k in range(len(master_MOUS_list)):
			mous_table_txt.write(str(master_MOUS_list[k]['member_ous_uid'])+'\n')
		mous_table_txt.write('-----\n')
		mous_table_txt.close()
		
	if make_line_plots:
		fig, axes = plt.subplots(3,4,figsize=(12,10))
		mlp=open('OUTPUTS/LineTable.txt','r')
		mlpf=mlp.readlines()
		line_lines=[-1,-1]
		for mlp_i in range(len(mlpf)):
			if mlpf[mlp_i].replace('\n','')==name:
				line_lines[0]=mlp_i
				for mlp_j in range(mlp_i,len(mlpf)):
					if '----' in mlpf[mlp_j]:
						line_lines[1]=mlp_j-1
						break
		linelist=[];numlines=0;L_INDEX=-1;bigvals=[]
		projlist=[];P_INDEX=-1
		if line_lines[1]-line_lines[0]!=1:
			for mlp_i in range(line_lines[0]+2,line_lines[1]+1):
				temp_mlp=mlpf[mlp_i].split(' ')
				LN=temp_mlp[0]+' '+temp_mlp[1]
				PR=str(temp_mlp[2]).replace('\n','')
				AR=float(temp_mlp[3])
				IT=float(temp_mlp[4])/3600.
				if LN not in linelist:
					linelist.append(LN)
					L_INDEX=numlines
					numlines+=1
					bigvals.append([IT,AR])
				else:
					L_INDEX=linelist.index(LN)
					if IT>bigvals[L_INDEX][0]:
						bigvals[L_INDEX][0]=IT
					if AR>bigvals[L_INDEX][1]:
						bigvals[L_INDEX][1]=AR
				if PR not in projlist:
					projlist.append(PR)
				PNUM=projlist.index(PR)
				plot_index=getsubplotnum3(L_INDEX)
				axes[plot_index[0],plot_index[1]].scatter([IT],[AR],marker='s',color=colorscheme[PNUM])
		try:
			for mlp_i in range(12):
				plot_index=getsubplotnum3(mlp_i)
				if mlp_i<len(linelist):
					axes[plot_index[0],plot_index[1]].set_title(linelist[mlp_i],weight='bold')
					axes[plot_index[0],plot_index[1]].set_xlim(0,1.1*bigvals[mlp_i][0])
					axes[plot_index[0],plot_index[1]].set_ylim(0,1.1*bigvals[mlp_i][1])
				else:
					axes[plot_index[0],plot_index[1]].set_visible(False)
		except TypeError:
			pass
		mlp.close()
		fig.add_subplot(111, frameon=False)
		fig.suptitle(name+'(z='+str(zred)+')',fontsize=fs2,weight='bold',x=0.56)
		plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
		plt.xlabel("Integration Time [hr]",weight='bold')
		plt.ylabel("Spatial Resolution [\"]",weight='bold')
		indexproj=sortpl(projlist)
		for PN_I in indexproj:
			plt.plot([], [], marker='s', label=projlist[PN_I], color=colorscheme[PN_I], linestyle=None)
		plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=4, fancybox=True, shadow=True, numpoints=1, handlelength=0)
		plt.tight_layout(pad=5.0,w_pad=0.1,h_pad=0.3)
		plt.savefig('OUTPUTS/'+name+'/LinePlot_'+name+'.png',dpi=300,bbox_inches='tight')
		plt.close()

	if download_products:
		#Make sure the data directory exists
		try:
			os.mkdir('OUTPUTS/'+name+'/data')
		except FileExistsError:
			pass
		#Get data
		alminer.download_data(myquery, fitsonly=False, dryrun=False, location='OUTPUTS/'+name+'/data', filename_must_include=['_of_'], print_urls=True)
		#Untar everything
		allfiles=os.listdir('OUTPUTS/'+name+'/data')
		for file in allfiles:
			if '.tar' in file:
				my_tar = tarfile.open('OUTPUTS/'+name+'/data/'+file)
				my_tar.extractall('OUTPUTS/'+name+'/data/'+file.replace('.tar','')) # specify which folder to extract to
				my_tar.close()
				'''
				#Get .fits
				res = []
				for (dir_path, dir_names, file_names) in os.walk('OUTPUTS/'+name+'/data/'+file.replace('.tar','')):
					res.extend(file_names)
				res2=[]
				for resi in res:
					if '.fits' in resi:
						res2.append(resi)
				#Check that our source is in each .fits
				for resi in range(len(res2)):
					hdul = fits.open(res2[resi])
					hdr = hdul[0].header
					wcs=WCS(hdr)
					print(wcs.all_pix2world())
				'''

f.close()

