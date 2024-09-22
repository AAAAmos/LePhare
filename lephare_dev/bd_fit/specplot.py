# Plot (file or window) the observed magnitudes
# of the objects, along with best-fit templates, 
# reading the info form .zsp and .pdz output 
# files of LePhare.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages 
import sys, getopt
#import subprocess
from matplotlib import rc
rc('text', usetex=True)
rc('xtick', labelsize=15)
rc('ytick', labelsize=15)


######## GET FILENAMES AND OPTIONS ########

nspec=len([s for s in sys.argv if s.find('.spec', -5) >=0 ]) # N of .spec files
odev='screen'  # print on screen as default 
sout='' ## use -o arg. to append a suffix to each file, just before .png/.eps/.ps

options, args=getopt.getopt(sys.argv[nspec+1:],"ho:d:c:",["help","output=","device=","context="])
for opt, arg in options:
    if opt in ('-o','--output'):
        sout=arg       
    if opt in ('-d','--device'): 
        if not (arg in ['pdf','png','eps','ps']):
            print("Please choose one of the following devices for output:")
            print("pdf, png, eps, ps")
            sys.exit()
        odev='.'+arg
        print( 'output device will be in '+odev+' format.')  #can be png,pdf (if none print on screen)
    if opt in ('-c','--context'): 
        print( 'context for models: ', arg)
        ctx=arg
    if opt in ('-h','--help'):
        print("""HELP still to be finished!!!

            SYNTAX:  python spec.py file[s].spec [OPTIONS]

            VARIABLES:
                        file[s] must be the output files from LePhare
                        (option SPEC_OUT = YES). Suffix .spec is 
                        compulsory.
            OPTIONS:
                    -d --device[=STR]:
                        select the output device; 
                        the string STR can be 'png','pdf','eps','ps'
                        (without quotes). If 'pdf' is chosen, all the plots
                        are collected in a single file. 
                        If the option is not set, then print on screen (with
                        a limit of 1 object/window).
                    -o --output[=STR]:
                        if --device = pdf, this option specifies the name of 
                        the pdf file. With any other --device values, a 
                        string STR is appended at the end of the filename,
                        just befor the extension (e.g. *STR.png). 
                        Do nothing when printing on screen.  
                    -c --context: 
                        TO BE IMPLEMENTED
                    -h --help:
                        print this help
            """)
        sys.exit()

if nspec==0:
    print('Please specify at least one .spec file')
    print('Try -h or --help options to get help.') 
    sys.exit()

if odev=='screen' and nspec > 1:
    print( """
Too many windows to be opened!
please reduce the number of .spec files
to one, or chose another output device.
""")
    sys.exit()

if odev=='.pdf': 
    print("All objects will be collected in a single pdf file named:")
    if len(sout)>0: 
        print( "--> ", sout+odev)
        pdp = PdfPages(sout+odev)
    else: 
        print( "--> MULTISPECv4.pdf")
        pdp = PdfPages('MULTISPECv4.pdf')
elif odev!='screen':
    if len(sout)>0:
        odev=sout+odev      
#when print on screen --output arg is not used

####### LOOP OVER .SPEC FILES ########

for k in range(nspec):

    ### Open .spec file[s] and read the parameters

    filename=sys.argv[1+k].replace('.spec','')
    fsp=open(filename+'.spec','r')

    bid=fsp.readline()    #header1
    line=fsp.readline()
    line=line.split()
    id=line[0]; zspec=line[1]; zphot=float(line[2])
    #z68low=float(line[3]); z68hig=float(line[4])

    bid=fsp.readline()    #header2
    line=fsp.readline()
    line=line.split()
    nfilt=int(line[1])
    
    bid=fsp.readline()    #header3
    line=fsp.readline()
    line=line.split()
    npdf=int(line[1])
    
    bid=fsp.readline()  
    #header4:  Type Nline Model Library Nband    Zphot Zinf Zsup Chi2  PDF     Extlaw EB-V Lir Age  Mass SFR SSFR
    models_info=[]
    for i in range(6):
        line=fsp.readline()
        model=line.split() 
        models_info.append(model)

    # Read observed mag, err_mag, filters' lambda_eff and width, models' mag
    mag=np.zeros(nfilt); em=np.zeros(nfilt); 
    lf=np.zeros(nfilt); dlf=np.zeros(nfilt); 
    mmod=np.zeros(nfilt); mfir=np.zeros(nfilt); mphys=np.zeros(nfilt)
    for i in range(nfilt):
        line=fsp.readline(); line=line.split() 
        mag[i]=float(line[0]); em[i]=float(line[1]); 
        lf[i]=float(line[2]); dlf[i]=float(line[3]); 
        mmod[i]=float(line[4]); mfir[i]=float(line[5]); mphys[i]=float(line[6])
        # print(mmod[i])
    
    #convert mag(AB syst.) in log(flux)
    ibad=np.where((mmod<=0) | (mmod>35))
    mmod=-0.4*(mmod-23.91)  # uJy
    mmod[ibad]=-10.
    print(mmod)
    ibad=np.where((mphys<=0) | (mphys>35))
    mphys=-0.4*(mphys-23.91)  # uJy
    mphys[ibad]=-10.
    ibad=np.where((mfir<=0) | (mfir>35))
    mfir=-0.4*(mfir-23.91)  # uJy
    mfir[ibad]=-10.
    
    zpdf=np.zeros([2,npdf])
    for i in range(npdf):
        line=fsp.readline()
        zpdf[:,i]=np.array(line.split())
        
    # Read spectra [lambda(A), Mag(AB)]
    # convert in log(F_nu) = -0.4*(mab-23.91) [uJy]
    # Loop over the 6 models (GAL-1, GAL-2, GAL-FIR, GAL-STOCH, QSO, STAR)
    lg=[]; mg=[]
    for m in range(6):
        nline=int(models_info[m][1])
        bid=np.zeros([2,nline])
        if nline>0:
            for i in range(nline): 
                line=fsp.readline()
                bid[:,i]=np.array(line.split())
                if (bid[1,i]>35):
                    bid[1,i]=-10.
                else:
                    bid[1,i]=-0.4*(bid[1,i]-23.91)
        lg.append(bid[0,:]/10000.)
        mg.append(bid[1,:])
    
    fsp.close() 
    
    ##############  PLOT  ###############
    
    ### Initialise the figure
    fig=plt.figure()
    
    ### Main panel
    ax1=fig.add_axes([.1,.1,.9,.6],
        xscale='log')
    
    ax1.set_xlabel('$\\lambda$ [$\\mu$m]', fontsize=18)
    ax1.set_ylabel('log(F$_{\\nu}$) [$\\mu$Jy]', fontsize=18)
    ax1.tick_params(axis='both', which='major', labelsize=16)
    ax1.set_xticks([0.5,1,2,3,4,5,6])
    ax1.set_xticklabels(['0.5','1','2','3','4','5','6'])
    
    # only the reliable obs mag will be plotted:
    em=em*2.
    dlf=dlf/2.
    mag1=mag[(mag>0.) & (mag<35) & (em>-3)]
    em1=em[(mag>0.) & (mag<35) & (em>-3)]
    lf1=lf[(mag>0.) & (mag<35) & (em>-3)]/10000.
    dlf1=dlf[(mag>0.) & (mag<35) & (em>-3)]/10000.
    
    ymin=max(mag1+2.); ymax=min(mag1-4.)
    if ymin>60:  ymin=30
    
    ic=[(em1>=0.) & (em1<2.)]
    # print(ic, lf1)
    lf2=lf1[ic[0]]
    mag2=-.4*(mag1[ic[0]]-23.91) 
    em2=0.4*em1[ic[0]]
    dlf2=dlf1[ic[0]]
    # low S/N bands:  
    ic2=[(em1>=2.) & (em1<8.)]
    lf2b=lf1[ic2[0]]
    mag2b=-.4*(mag1[ic2[0]]-23.91) 
    em2b=0.4*em1[ic2[0]]
    dlf2b=dlf1[ic2[0]]
    
    # normal fit
    
    # read star model file
    MODname = []
    with open('/home/yuan/LePhare/lephare_dev/sed/STAR/STAR_MOD_ALL2.list', 'r') as f:
        lines = f.readlines()
    for line in lines:
        MODname.append(line.split('\t')[0].split('/')[-1][:-4])
    
    # spec type fit
    
    # MODname = []
    # with open('/home/yuan/LePhare/lephare_dev/sed/STAR/bdtype_spec_ds.list', 'r') as f:
    #     lines = f.readlines()
    # for line in lines:
    #     MODname.append(line.split('\t')[0].split('/')[-1].split('.')[0])
    # Typename = []
    # with open('/home/yuan/LePhare/lephare_dev/sed/STAR/bdtype_specname.list', 'r') as f:
    #     lines = f.readlines()
    # for line in lines:
    #     Typename.append(line.split()[0])
    
    # set the plot aspect
    ax1.axis([min(lf1)*.85,max(lf1)*1.2,-0.4*(ymin-23.91),-0.4*(ymax-23.91)]) 
    ### plot SED and print info of best-fit models
    col_lst=['r','g','b','m','y','c']  #each one with a different color
    
    # plt.figtext(0.2,0.73,' Type:\n Chi$^2$:\n ModelName:', size=12)
    # plt.figtext(0.2,0.73,' Type:\n Modelnum:\n Chi$^2$:\n parameter:', size=10)
    # plt.figtext(0.2,0.73,' Type:\n Model:', size=10)
    
    parameter = MODname[int(models_info[5][2])-1]
    teff = parameter.split('_')[1].split('g')[0][1:]
    grav = parameter.split('_')[1].split('g')[1].split('nc')[0]
    if float(parameter.split('m')[1][:-1]) == 0.5:
        metal = '+0.5'
    elif float(parameter.split('m')[1][:-1]) == 0.0:
        metal = '0.0'
    else:
        metal = '-0.5'
    
    plt.figtext(0.1,0.73,
                '$T_{eff}$ = '+teff + 'K,  g = '+grav+'m/s$^2$,  [M/H] = '+metal, 
                size=18)
    
    #########################################plt.figtext(0.35,0.16,'MbsID=778',size=20)
    #  plt.figtext(0.35,0.11,'HSC ID=79671054930306067',size=20)
    
    model_name = ["Galaxy model", "Galaxy model 2", "Galaxy model (FIR)", "Galaxy model (Stochastic)", "QSO model", "Sonora-Bobcat model"]
    
    iml=0
    for im in range(6):
        if int(models_info[im][2])<0: continue   #print only models really used
        iml=iml+1  #counter of models used
        ax1.plot(lg[im],mg[im],color=col_lst[im],
                label=model_name[im]+', $\\chi^2$='+str(float(models_info[im][8])))  #plot the SED
        
    
    plt.legend(loc='upper left', fontsize=11)
    
    # plot absorption lines
    
    # H2O = [[1.18, 1.38], [1.7, 1.95], [2.4, 3.05], [4.4, 8]]
    # CH4 = [[1.5, 1.75], [2.2, 2.4], [3.1, 4.2]]
    # NH3 = [[0.85, 0.97], [1.32, 1.62], [2.25, 2.52]]
    
    # for i in range(len(H2O)):
    #     if i == 0:
    #         ax1.plot(H2O[i], [1.7, 1.7], color='k', linestyle='--', label='Strong H2O absorption region')
    #     else:
    #         ax1.plot(H2O[i], [1.7, 1.7], color='k', linestyle='--')
        
    # for i in range(len(CH4)):
    #     if i == 0:
    #         ax1.plot(CH4[i], [1.4, 1.4], color='k', linestyle='-.', label='Strong CH4 absorption region')
    #     else:
    #         ax1.plot(CH4[i], [1.4, 1.4], color='k', linestyle='-.')
            
    # for i in range(len(NH3)):
    #     if i == 0:
    #         ax1.plot(NH3[i], [1.1, 1.1], color='k', linestyle=':', label='Strong NH3 absorption region')
    #     else:
    #         ax1.plot(NH3[i], [1.1, 1.1], color='k', linestyle=':')
    
    # H2O = [1.15, 1.82, 2.7, 5.5]
    # CH4 = [1.6, 2.1, 3.7]
    # NH3 = [0.91, 1.45, 2.35]
    
    # for i in range(len(H2O)):
    #     if i == 0:
    #         ax1.arrow(H2O[i], 1-0.5, 0, -0.5, color='#00CCFF',
    #                 head_length=0.2, head_width=0.05*H2O[i], width=0.02*H2O[i], label='Strong H2O absorption region')
    #     else:
    #         ax1.arrow(H2O[i], 1-0.5, 0, -0.5, color='#00CCFF',
    #                 head_length=0.2, head_width=0.05*H2O[i], width=0.02*H2O[i])
            
    # for i in range(len(CH4)):
    #     if i == 0:
    #         ax1.arrow(CH4[i], 1-0.5, 0, -0.2, color='#FF99CC',
    #                 head_length=0.2, head_width=0.048*CH4[i], width=0.018*CH4[i], label='Strong CH4 absorption region')
    #     else:
    #         ax1.arrow(CH4[i], 1-0.5, 0, -0.2, color='#FF99CC',
    #                 head_length=0.2, head_width=0.048*CH4[i], width=0.018*CH4[i])
            
    # for i in range(len(NH3)):
    #     if i == 0:
    #         ax1.arrow([NH3[i], NH3[i]], [1.1, 1.1], color='#FFCC00', linestyle=':', label='Strong NH3 absorption region')
    #     else:
    #         ax1.arrow([NH3[i], NH3[i]], [1.1, 1.1], color='#FFCC00', linestyle=':')
    
    # plt.legend(loc='lower left', fontsize=10)
    
    # plot the obs mag...
    ax1.errorbar(lf2b,mag2b,yerr=em2b,xerr=dlf2b,fmt='o',color='0.6')
    ax1.errorbar(lf2,mag2,yerr=em2,xerr=dlf2,fmt='o',color='0.')
    #... and upper limits
    iu=np.where(em1<0)
    if len(iu[0])>0 :
        lf3=lf1[iu]
        mag3=-0.4*(mag1[iu]-23.91)
        ax1.quiver(lf3,mag3,0,-1,units='height',width=0.005,headwidth=5,color='k',pivot='tip', zorder=10)
        
    ### 2nd panel (inset) showing PDF(z)
    base=0.9-0.02*iml #starting position for the inset plot
    # if base>0.84: base=0.84
    # ax2=fig.add_axes([0.13,base-0.20,0.3,0.20],xlabel='z$\\_{phot}$') #,title='z\_spec='+zspec)
    # ax2.yaxis.set_label_position("right")
    # ax2.yaxis.set_ticks_position("right")
    # ax2.plot(zpdf[0,:],zpdf[1,:],color='r')
    
    # ax1.set_ylim([-2.5+0.5, 1+0.7])
    ax1.set_ylim([-2.5-0.3, 1+0.3])
    # ax1.set_ylim([-2.5-0.35, 1+0.05])
    # ax1.set_ylim([-2-1.5, 1.5+0.2])
    ax1.set_xlim([0.4, 5.5])
    #  ax2.set_xlim([5, 7])
    #  ax1.set_xticks([0.5,1,2])
    ax1.get_xaxis().get_major_formatter().labelOnlyBase = False
    #plot also z_phot with error bar
    #ax2.errorbar(zphot,0.5,fmt='ok',xerr=[[zphot-z68low],[z68hig-zphot]],mfc='none')
    
    # no chose  window/png for the moment TBI
    #plt.show()
    if  odev=='.pdf': 
        plt.savefig(pdp,format='pdf',dpi=300,bbox_inches='tight')
    elif odev=='screen':
        plt.show()
    else: 
        plt.savefig(filename+odev, bbox_inches='tight',dpi=300)

if odev=='.pdf': pdp.close()

print( 'end') 
